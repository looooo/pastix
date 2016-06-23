/**
 *
 * @file solver_compute_gpu_distrib.c
 *
 *  PaStiX analyse routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#if !defined(PASTIX_WITH_CUDA)
#error "This file should be compiled only with CUDA enabled"
#endif

#include <cuda.h>
#include <cuda_runtime_api.h>
#include "solver.h"
#include "flops.h"

#define PASTIX_GPU_MEMORY_WRITE  0.7

#ifndef PASTIX_GPU_MIN_NBPAGES
#  define PASTIX_GPU_MIN_NBPAGES 0
#endif
#ifndef PASTIX_GPU_MIN_UPDATES
#  define PASTIX_GPU_MIN_UPDATES 0
#endif
#ifndef PASTIX_GPU_MIN_FLOP
#  define PASTIX_GPU_MIN_FLOP 0
#endif

typedef struct _cblk_elem_ {
    SolverCblk *cblk;
    double      criterium;
    double      flops; /* in Mflop */
    int         cblk_size;
    int         updates;
} cblk_elem;

typedef struct _gpu_elem_ {
    int id;
    int avail_pages;
    int total_pages;
} gpu_elem;

int solverComputeGPUDistrib( SolverMatrix *solvmtx,
                             int           ngpus,
                             double        memory_percentage,
                             size_t        eltsize,
                             int           criterium,
                             int           factotype )
{
    SolverCblk    *cblktab = solvmtx->cblktab;
    SolverBlok    *bloktab = solvmtx->bloktab;
    SolverCblk    *cblk, *lcblk;
    SolverBlok    *blok;
    cblk_elem     *ecblktab, *ecblk;
    gpu_elem      *egputab, *egpu;
    pastix_queue_t cblkQueue;
    pastix_queue_t gpuQueue;
    pastix_int_t   cblkid, gpuid, i, j, m, n, k;
    int            pushgpu;
    int *devices_cblk, *devices_gemm;

    /**
     * Allocate arrays for statistics
     */
    devices_cblk = malloc((1+ngpus)*sizeof(int));
    devices_gemm = malloc((1+ngpus)*sizeof(int));
    memset(devices_cblk, 0, (1+ngpus)*sizeof(int));
    memset(devices_gemm, 0, (1+ngpus)*sizeof(int));

    /**
     * Create the list with all cblk that can potentially be computed on a GPU
     */
    MALLOC_INTERN(ecblktab, solvmtx->cblknbr, cblk_elem);
    pqueueInit(&cblkQueue, solvmtx->cblknbr);

    cblk  = cblktab;
    ecblk = ecblktab;
    for(cblkid = 0; cblkid < solvmtx->cblknbr; cblkid++, cblk++, ecblk++) {

        /* Let's start by setting the update on the CPU */
        cblk->gpuid = -2;
        ecblk->cblk      = cblk;
        ecblk->updates   = cblk[1].brownum - cblk[0].brownum;
        ecblk->cblk_size = cblk_colnbr( cblk ) * cblk->stride;
        ecblk->criterium = 0.;
        ecblk->flops     = 0.;

        /**
         * If it's a leaf or there is no GPUs available, we just make sure it
         * will go to the cpu
         */
        if ( ecblk->updates == 0 || ngpus <= 0 ) {
            devices_cblk[0]++;
        }
        else {

            /**
             * Let's compute the flops
             */
            double flops = 0.;
            double factor = 1.;

            if (factotype == PastixFactLU) factor = 2.;

            for(j=cblk[0].brownum; j<cblk[1].brownum-1; j++) {
                blok  = bloktab + solvmtx->browtab[j];
                lcblk = cblktab + blok->lcblknm;

                /* TODO: Incorrect in 2D !!!!! */
                m = lcblk->stride - blok->coefind;
                k = cblk_colnbr( lcblk );
                n = blok_rownbr( blok );
                flops += factor * FLOPS_DGEMM(m,n,k) / 1.e6;
            }
            /* Last block applied only on L for LU, and always for other */
            if (j < cblk[1].brownum)
            {
                blok  = bloktab + solvmtx->browtab[j];
                lcblk = cblktab + blok->lcblknm;

                /* TODO: Incorrect in 2D !!!!! */
                m = lcblk->stride - blok->coefind;
                k = cblk_colnbr( lcblk );
                n = blok_rownbr( blok );
                flops += FLOPS_DGEMM(m,n,k) / 1.e6;
            }
            ecblk->flops = flops;

            /* check for int overflow */
            /* assert(new_cblk->flops < 4294967296.); */

            /**
             * Compute the number of memory pages used by this cblk
             */
            ecblk->cblk_size = pastix_iceil( ecblk->cblk_size, eltsize );

            if (factotype == PastixFactLU){
                ecblk->cblk_size *= 2;
            }

            /* Sort criterium */
            switch(criterium){
            case API_GPU_CRITERION_UPDATES:
                ecblk->criterium = - (double)(ecblk->updates);
                break;
            case API_GPU_CRITERION_CBLKSIZE:
                ecblk->criterium = - (double)(ecblk->cblk_size);
                break;
            case API_GPU_CRITERION_FLOPS:
                ecblk->criterium = - (double)(ecblk->flops);
                break;
            default :
                ecblk->criterium = (double)cblkid;
            }
            pqueuePush1(&cblkQueue, cblkid, ecblk->criterium);
        }
    }

    /**
     * Create the list for the gpus
     */
    MALLOC_INTERN(egputab, ngpus, gpu_elem);
    pqueueInit(&gpuQueue, ngpus);
    egpu = egputab;
    for(i = 0; i < ngpus; i++, egpu++) {
        size_t total_mem, initial_free_mem;
        size_t how_much_we_allocate;

        cudaSetDevice(i);
        cudaMemGetInfo( &initial_free_mem, &total_mem );

        /**
         * PaRSEC allocated percent * initial_free_mem memory,
         * Unfortunatelly we can't get the information about how much memory has
         * been allocated, thus we base our computations on the total
         * memory. And for aour static mapping, we use only 70% of the allocated
         * memory.
         */
        how_much_we_allocate = ((double)memory_percentage * (double)total_mem) / 100.;
        how_much_we_allocate = how_much_we_allocate * PASTIX_GPU_MEMORY_WRITE;
        how_much_we_allocate = pastix_iceil( how_much_we_allocate, eltsize );

        egpu->id = i;
        egpu->total_pages = how_much_we_allocate;
        egpu->avail_pages = how_much_we_allocate;
        pqueuePush1(&gpuQueue, i, -(egpu->avail_pages));
    }

    /**
     * Dynamic algorithm to associate the highest priority cblk, to the less loaded GPU
     */
    while (pqueueSize(&cblkQueue) > 0) {
        pushgpu = 0;

        /* Get the cblk with the highest criterium */
        cblkid = pqueuePop(&cblkQueue);
        ecblk = ecblktab + cblkid;

        gpuid = ecblk->cblk->gpuid;
        egpu = NULL;

        /**
         * The cblk has already a prefered GPU
         * Let's try to use this GPU.
         */
        if (gpuid != -2) {
            assert(gpuid != -1);

            egpu = egputab + gpuid;
            if(egpu->avail_pages < ecblk->cblk_size) {
                egpu = NULL;
            }
        }

        /**
         * If no prefered GPU, let's take the one with the largest memory available
         */
        if (egpu == NULL) {
            pastix_int_t gpuid = pqueuePop(&gpuQueue);
            egpu = egputab + gpuid;
            pushgpu = 1;
        }

        gpuid = egpu->id;
        /**
         * If it doesn't match the minima, or if the GPU cannot handle it, we
         * force it to the CPU
         */
        if( (egpu->avail_pages < ecblk->cblk_size)      ||
            (ecblk->updates   < PASTIX_GPU_MIN_UPDATES) ||
            (ecblk->cblk_size < PASTIX_GPU_MIN_NBPAGES) ||
            (ecblk->flops     < PASTIX_GPU_MIN_FLOP) )
        {
            ecblk->cblk->gpuid = -2;
            devices_cblk[0] += 1;
            if (factotype == PastixFactLU) {
                devices_gemm[0] += 2 * ecblk->updates - 1;
            } else {
                devices_gemm[0] += ecblk->updates;
            }
        }
        else
        {
            /**
             * Let's link this cblk and gpu together, and add this GPU as the
             * prefered one to all the cblk that will contribute to this cblk.
             */
            egpu->avail_pages -= ecblk->cblk_size;
            ecblk->cblk->gpuid = gpuid;

            /* Add prediction to contributors */
            for(j=cblk[0].brownum; j<cblk[1].brownum; j++) {
                blok  = bloktab + solvmtx->browtab[j];
                lcblk = cblktab + blok->lcblknm;

                if (lcblk->gpuid != -2) {
                    lcblk->gpuid = gpuid;
                }
            }

            devices_cblk[gpuid+1] += 1;
            if (factotype == PastixFactLU) {
                devices_gemm[gpuid+1] += 2 * ecblk->updates - 1;
            } else {
                devices_gemm[gpuid+1] += ecblk->updates;
            }
        }

        if (pushgpu)
            pqueuePush1(&gpuQueue, gpuid, -(egpu->avail_pages));
    }

    fprintf(stdout,"cpu :  %d cblk, %d gemm\n", devices_cblk[0], devices_gemm[0]);
    egpu = egputab;
    for (i = 0 ; i < ngpus; i++, egpu++) {
        fprintf(stdout,"gpu %d:  %d cblk, %d gemm, memory : %d / %d \n",
                (int)i, devices_cblk[i+1], devices_gemm[i+1],
                egpu->avail_pages, egpu->total_pages );
    }

    pqueueExit(&cblkQueue);
    pqueueExit(&gpuQueue);
    return PASTIX_SUCCESS;
}
