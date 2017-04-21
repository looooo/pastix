/**
 *
 * @file sopalin_zsytrf.c
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
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_PARSEC)
#include <parsec.h>
#include <parsec/data.h>
#include <parsec/data_distribution.h>
#endif

void
sequential_zsytrf( pastix_data_t  *pastix_data,
                   sopalin_data_t *sopalin_data )
{
    SolverMatrix       *datacode = pastix_data->solvmatr;
    SolverCblk         *cblk;
    double              threshold = sopalin_data->diagthreshold;
    pastix_complex64_t *work1, *work2;
    pastix_int_t  i;
    (void)pastix_data;

    MALLOC_INTERN( work1, pastix_imax(datacode->gemmmax, datacode->diagmax),
                   pastix_complex64_t );
    MALLOC_INTERN( work2, datacode->gemmmax, pastix_complex64_t );

    cblk = datacode->cblktab;
    for (i=0; i<datacode->cblknbr; i++, cblk++){

        if ( cblk->cblktype & CBLK_IN_SCHUR )
            break;

        /* Compute */
        cpucblk_zsytrfsp1d( datacode, cblk, threshold,
                            work1, work2 );
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( datacode, "sytrf_L.txt" );
#endif

    memFree_null( work1 );
    memFree_null( work2 );
}

void
thread_pzsytrf( isched_thread_t *ctx, void *args )
{
    sopalin_data_t     *sopalin_data = (sopalin_data_t*)args;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    Task               *t;
    pastix_complex64_t *work1, *work2;
    pastix_int_t  i, ii;
    pastix_int_t  tasknbr, *tasktab;
    int rank = ctx->rank;

    MALLOC_INTERN( work1, pastix_imax(datacode->gemmmax, datacode->diagmax),
                   pastix_complex64_t );
    MALLOC_INTERN( work2, datacode->gemmmax, pastix_complex64_t );

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        if ( cblk->cblktype & CBLK_IN_SCHUR )
            continue;

        /* Wait */
        do {} while( cblk->ctrbcnt );

        /* Compute */
        cpucblk_zsytrfsp1d( datacode, cblk, sopalin_data->diagthreshold,
                            work1, work2 );
    }

#if defined(PASTIX_DEBUG_FACTO) && 0
    isched_barrier_wait( &(ctx->global_ctx->barrier) );
    if (rank == 0) {
        coeftab_zdump( datacode, "sytrf_L.txt" );
    }
    isched_barrier_wait( &(ctx->global_ctx->barrier) );
#endif

    memFree_null( work1 );
    memFree_null( work2 );
}

void
thread_zsytrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    isched_parallel_call( pastix_data->isched, thread_pzsytrf, sopalin_data );
}

#if defined(PASTIX_WITH_PARSEC) && 0
void
parsec_zsytrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    parsec_context_t *ctx;

    /* Start PaRSEC */
    if (pastix_data->parsec == NULL) {
        int argc = 0;
        pastix_parsec_init( pastix_data, &argc, NULL );
    }
    ctx = pastix_data->parsec;

    /* Run the facto */
    dsparse_zsytrf_sp( ctx, sopalin_data->solvmtx->parsec_desc, sopalin_data );
}
#endif

static void (*zsytrf_table[4])(pastix_data_t *, sopalin_data_t *) = {
    sequential_zsytrf,
    thread_zsytrf,
#if defined(PASTIX_WITH_PARSEC) && 0
    parsec_zsytrf,
#else
    NULL,
#endif
    NULL
};

void
sopalin_zsytrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zsytrf)(pastix_data_t *, sopalin_data_t *) = zsytrf_table[ sched ];

    if (zsytrf == NULL) {
        zsytrf = thread_zsytrf;
    }
    zsytrf( pastix_data, sopalin_data );
}
