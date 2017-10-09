/**
 *
 * @file kernels_trace.c
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * PaStiX trace and modelling routines
 *
 * @version 6.0.0
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2017-02-09
 *
 **/
#include "common.h"
#include "solver.h"
#include "kernels_trace.h"

volatile double kernels_flops[PastixKernelLvl1Nbr];

#if defined(PASTIX_WITH_EZTRACE)

int pastix_eztrace_level = 1;

#endif

#if defined(PASTIX_GENERATE_MODEL)

pastix_model_entry_t *model_entries     = NULL;
volatile int32_t      model_entries_nbr = -1;
int32_t               model_size        = 0;

#endif

/**
 *******************************************************************************
 *
 * @brief Start the trace module
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure of the problem to give input information
 *          to the different trace modes.
 *
 *******************************************************************************/
void
kernelsTraceStart( const pastix_data_t *pastix_data )
{
    const SolverMatrix *solvmtx = pastix_data->solvmatr;

#if defined(PASTIX_WITH_EZTRACE)
    {
        char *level = pastix_getenv("PASTIX_EZTRACE_LEVEL");
        if (level != NULL) {
            pastix_eztrace_level = atoi(level);
            pastix_cleanenv(level);
        }

        if ( pastix_data->dirtemp != NULL ) {
            pastix_setenv( "EZTRACE_TRACE_DIR", pastix_data->dirtemp, 1 );
        }
        eztrace_start ();
    }
#endif /* defined(PASTIX_WITH_EZTRACE) */

#if defined(PASTIX_GENERATE_MODEL)
    {
        pastix_int_t cblknbr   = solvmtx->cblknbr;
        pastix_int_t cblkmin2d = solvmtx->cblkmin2d;
        pastix_int_t total_number_of_tasks = 0;
        pastix_int_t nbfact, nbtrsm, nbgemm;
        pastix_int_t cblknum;
        SolverCblk  *cblk;

        /* Factorization kernels */
        nbfact = cblknbr;

        /* TRSM kernels */
        nbtrsm = cblkmin2d + (cblknbr - cblkmin2d) * solvmtx->cblkmaxblk;
        if ( solvmtx->factotype == PastixFactLU ) {
            nbtrsm *= 2;
        }

        /* GEMM kernels */
        nbgemm = solvmtx->bloknbr - cblknbr;
        if ( solvmtx->factotype == PastixFactLU ) {
            nbgemm *= 2;
        }

        cblk = solvmtx->cblktab+cblkmin2d;
        for(cblknum = cblkmin2d; cblknum < cblknbr; cblknum++, cblk++ ) {
            pastix_int_t nbodb = (cblk[1].fblokptr - cblk[0].fblokptr) - 1;

            if ( solvmtx->factotype == PastixFactLU ) {
                nbgemm += nbodb * nbodb;
            }
            else {
                nbgemm += (nbodb * (nbodb-1)) / 2;
            }
        }

        total_number_of_tasks = nbfact + nbtrsm + nbgemm;
        model_entries = malloc( total_number_of_tasks * sizeof(pastix_model_entry_t) );
        model_size = total_number_of_tasks;
    }
#endif

    memset( (void*)kernels_flops, 0, PastixKernelLvl1Nbr * sizeof(double) );

    (void)solvmtx;
    return;
}

/**
 *******************************************************************************
 *
 * @brief Stop the trace module
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure of the problem to get input information
 *          for the different trace modes, and store output statistics.
 *
 *******************************************************************************/
double
kernelsTraceStop( const pastix_data_t *pastix_data )
{
    double total_flops = 0.0;

#if defined(PASTIX_WITH_EZTRACE)
    eztrace_stop ();
#endif

#if defined(PASTIX_GENERATE_MODEL)
    {
        pastix_model_entry_t *entry = model_entries;
        pastix_int_t i, gpucase;
        FILE *f;

        f = fopen( "model.csv", "w" );
        if ( f == NULL ) {
            goto end_model;
        }

        gpucase = pastix_data->iparm[IPARM_GPU_NBR];
        if ( gpucase ) {
            fprintf(f, "# GPU Model data\n");
        }
        else {
            fprintf(f, "# CPU Model data\n");
        }

        for(i=0; i <= model_entries_nbr; i++, entry++ ) {
            switch( entry->ktype ) {
            case PastixKernelGETRF:
            case PastixKernelHETRF:
            case PastixKernelPOTRF:
            case PastixKernelPXTRF:
            case PastixKernelSYTRF:
            case PastixKernelSCALOCblk:
            case PastixKernelSCALOBlok:
            case PastixKernelTRSMCblk1d:
            case PastixKernelTRSMCblk2d:
            case PastixKernelTRSMCblkLR:
            case PastixKernelTRSMBlokLR:
            case PastixKernelGEMMCblk1d1d:
            case PastixKernelGEMMCblkFRLR:
            case PastixKernelGEMMCblkLRLR:
            case PastixKernelGEMMBlokLRLR:
                if ( gpucase ) {
                    continue;
                }

            default:
                fprintf( f, "%d;%d;%d;%d;%e\n",
                         entry->ktype, entry->m, entry->n, entry->k, entry->time );
            }
        }

        fclose( f );

        free( model_entries );

        /* Reinitialize values */
        model_entries     = NULL;
        model_entries_nbr = -1;
        model_size        = 0;
    }
  end_model:
#endif

    /* Disable as long as the counters are not correct */
    if (0)
    {
        int i;

        /* Compute the total number of flops */
        for( i=0; i<PastixKernelLvl1Nbr; i++ )
        {
            total_flops += kernels_flops[i];
        }

        fprintf(stderr, "The total number of flops excuted is %lf\n", total_flops);
    }

    (void)pastix_data;
    return total_flops;
}
