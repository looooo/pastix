/**
 *
 * @file sopalin_zgetrf.c
 *
 *  PaStiX factorization routines
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
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "isched.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_PARSEC)
#include <dague.h>
#include <dague/data.h>
#include <dague/data_distribution.h>
#include "parsec/sparse-matrix.h"
#endif

#if defined(PASTIX_WITH_HODLR)
static pastix_int_t compute_cblklevel( pastix_int_t cblknum )
{
    /* cblknum level has already been computed */
    pastix_int_t father = treetab[cblknum];
    if ( father == -1 ) {
        return 1;
    }
    else {
        return compute_cblklevel( father ) + 1;
    }
}
#endif /* defined(PASTIX_WITH_HODLR) */

void
sequential_zgetrf( pastix_data_t  *pastix_data,
                   sopalin_data_t *sopalin_data )
{
    SolverMatrix       *datacode = pastix_data->solvmatr;
    SolverCblk         *cblk;
    double              threshold = sopalin_data->diagthreshold;
    pastix_complex64_t *work;
    pastix_int_t  i;
    (void)pastix_data;

    MALLOC_INTERN( work, datacode->gemmmax, pastix_complex64_t );

    cblk = datacode->cblktab;
    current_cblk = 0;

    /* To apply contributions with a depth-first search */
    /* Warning: does not work with parallel implementations */
#if defined(PASTIX_WITH_HODLR)
    pastix_int_t  j;
    if (0){
        pastix_int_t max_level = 0;
        for (i=0; i<datacode->cblknbr; i++){
            pastix_int_t level = compute_cblklevel( current_cblk++ );
            if (level > max_level)
                max_level = level;
        }

        for (j=max_level; j>=1; j--){
            pastix_int_t id_b = 0;

            char *tol        = getenv("TOLERANCE");
            double tolerance = atof(tol);

            char env_char[256];
            double new_tol = 0.00001; //tolerance / (1 << j);
            sprintf(env_char, "%lf", new_tol);
            setenv("HODLR_TOLERANCE", &env_char, 1);

            current_cblk = 0;
            cblk         = datacode->cblktab;

            for (i=0; i<datacode->cblknbr; i++, cblk++){
                pastix_int_t level = compute_cblklevel( id_b++ );

                if (level == j){
                    printf("Supernode %ld level %ld Tolerance %lf\n", current_cblk, level, new_tol);
                    core_zgetrfsp1d( datacode, cblk, threshold, work );
                }
                else{
                    current_cblk++;
                }
            }
        }
    }

    else
#endif /* defined(PASTIX_WITH_HODLR) */
    {

        /* Compress information coming from A */
        for (i=0; i<datacode->cblknbr; i++, cblk++){
            pastix_int_t dima, dimb;

            SolverBlok *blok = cblk->fblokptr;
            pastix_complex64_t *U = cblk->ucoeftab;
            pastix_complex64_t *L = cblk->lcoeftab;

            /* Horizontal dimension */
            dima = cblk->lcolnum - cblk->fcolnum + 1;
            /* Vertical dimension */
            dimb = cblk->stride  - dima;

            /* Size of the diagonal blok */
            pastix_int_t totalsize = 0;
            pastix_int_t tot = 0;

            /* First pass to count the number of off-diagonal blocks */
            blok++;
            while (totalsize != dimb){
                totalsize += blok->lrownum - blok->frownum + 1;
                tot++;
                blok++;
            }

            /* printf("\nCBLK %ld got %ld off-diag (SIZE %ld %ld)\n", i, tot, dima, dimb); */

            /* Second pass to compress each off-diagonal block */
            blok = cblk->fblokptr + 1;
            totalsize = 0;
            tot = 0;

            while (totalsize != dimb)
            {
                pastix_int_t bloksize    = blok->lrownum - blok->frownum + 1;

                double *s;
                pastix_int_t dim_min = pastix_imin(dima, bloksize);
                s = malloc( dim_min * sizeof(double));
                memset(s, 0, dim_min * sizeof(double));

                /* printf("BLOK  starts at %ld with size %ld\n", dima + totalsize, bloksize); */

                pastix_complex64_t *uU, *vU;
                pastix_int_t rankU = -1;
                uU = malloc( bloksize * bloksize * sizeof(pastix_complex64_t));
                vU = malloc( dima     * dima     * sizeof(pastix_complex64_t));
                memset(uU, 0, bloksize * bloksize * sizeof(pastix_complex64_t));
                memset(vU, 0, dima     * dima     * sizeof(pastix_complex64_t));

                pastix_complex64_t *uL, *vL;
                pastix_int_t rankL = -1;
                uL = malloc( bloksize * bloksize * sizeof(pastix_complex64_t));
                vL = malloc( dima     * dima     * sizeof(pastix_complex64_t));
                memset(uL, 0, bloksize * bloksize * sizeof(pastix_complex64_t));
                memset(vL, 0, dima     * dima     * sizeof(pastix_complex64_t));

                if (bloksize > 1){
                    pastix_complex64_t *data = U + dima + totalsize;
                    rankU = z_compress_LR(data, dima, bloksize,
                                         s,
                                         uU, bloksize,
                                         vU, dima,
                                         cblk->stride);
                    /* if (rankU != -1){ */
                    /*     z_uncompress_LR(data, dima, bloksize, */
                    /*                     uU, bloksize, */
                    /*                     vU, dima, */
                    /*                     cblk->stride, rankU); */
                    /* } */

                    data = L + dima + totalsize;
                    rankL = z_compress_LR(data, dima, bloksize,
                                          s,
                                          uL, bloksize,
                                          vL, dima,
                                          cblk->stride);

                    /* if (rankL != -1){ */
                    /*     z_uncompress_LR(data, dima, bloksize, */
                    /*                     uL, bloksize, */
                    /*                     vL, dima, */
                    /*                     cblk->stride, rankL); */
                    /* } */

                }

                blok->coefU_u_LR = uU;
                blok->coefU_v_LR = vU;
                blok->rankU      = rankU;
                blok->coefL_u_LR = uL;
                blok->coefL_v_LR = vL;
                blok->rankL      = rankL;

                printf("Compress blok %ld %ld with ranks %ld %ld\n", dima, bloksize, rankU, rankL);

                totalsize += bloksize;
                tot++;
                blok++;
            }
        }

        cblk = datacode->cblktab;
        for (i=0; i<datacode->cblknbr; i++, cblk++){
            core_zgetrfsp1d( datacode, cblk, threshold, work );
        }
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( datacode, "getrf_L.txt" );
#endif

    memFree_null( work );
}

void
thread_pzgetrf( int rank, void *args )
{
    sopalin_data_t     *sopalin_data = (sopalin_data_t*)args;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    Task               *t;
    pastix_complex64_t *work;
    pastix_int_t  i, ii;
    pastix_int_t  tasknbr, *tasktab;

    MALLOC_INTERN( work, datacode->gemmmax, pastix_complex64_t );

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        /* Compute */
        core_zgetrfsp1d( datacode, cblk, sopalin_data->diagthreshold, work );
    }

#if defined(PASTIX_DEBUG_FACTO)
    isched_barrier_wait( &(((isched_t*)(sopalin_data->sched))->barrier) );
    if (rank == 0) {
        coeftab_zdump( datacode, "getrf_L.txt" );
    }
    isched_barrier_wait( &(((isched_t*)(sopalin_data->sched))->barrier) );
#endif

    memFree_null( work );
}

void
thread_zgetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    isched_parallel_call( pastix_data->isched, thread_pzgetrf, sopalin_data );
}

#if defined(PASTIX_WITH_PARSEC)
void
parsec_zgetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    sparse_matrix_desc_t desc;
    dague_context_t *ctx;
    int argc = 0;

    /* Start PaRSEC */
    if (pastix_data->parsec == NULL) {
        pastix_data->parsec = dague_init( -1, &argc, NULL );
    }
    ctx = pastix_data->parsec;

    /* Create the matrix descriptor */
    sparse_matrix_init( &desc, sopalin_data->solvmtx,
                        pastix_size_of( PastixComplex64 ), 1, 0 );

    /* Run the facto */
    dsparse_zgetrf_sp( ctx, &desc, sopalin_data );

    /* Destroy the decriptor */
    sparse_matrix_destroy( &desc );

    dague_fini( &(pastix_data->parsec) );


    /* TODO: add to the DAG */

   /* SolverMatrix       *datacode = sopalin_data->solvmtx; */
    /* SolverCblk         *cblk; */
    /* Task               *t; */
    /* pastix_int_t        tasknbr, *tasktab; */
    /* tasknbr = datacode->ttsknbr[0]; */
    /* tasktab = datacode->ttsktab[0]; */

    /* pastix_int_t i, ii; */
    /* for (ii=0; ii<tasknbr; ii++) { */
    /*     i = tasktab[ii]; */
    /*     t = datacode->tasktab + i; */
    /*     cblk = datacode->cblktab + t->cblknum; */

    /*     pastix_complex64_t *L = cblk->lcoeftab; */
    /*     pastix_complex64_t *U = cblk->ucoeftab; */

    /*     /\* Compute *\/ */
    /*     core_zgetrfsp1d_trsm2(cblk, L, U); */
    /* } */
}
#endif

static void (*zgetrf_table[4])(pastix_data_t *, sopalin_data_t *) = {
    sequential_zgetrf,
    thread_zgetrf,
#if defined(PASTIX_WITH_PARSEC)
    parsec_zgetrf,
#else
    NULL,
#endif
    NULL
};

void
sopalin_zgetrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zgetrf)(pastix_data_t *, sopalin_data_t *) = zgetrf_table[ sched ];

    if (zgetrf == NULL) {
        zgetrf = thread_zgetrf;
    }
    zgetrf( pastix_data, sopalin_data );
}
