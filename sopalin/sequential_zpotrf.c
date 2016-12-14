/**
 *
 * @file sopalin_zpotrf.c
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
#include <dague.h>
#include <dague/data.h>
#include <dague/data_distribution.h>

int dsparse_zpotrf_sp( dague_context_t *dague,
                       sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data );
#endif

void
sequential_zpotrf( pastix_data_t  *pastix_data,
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
    for (i=0; i<datacode->cblknbr; i++, cblk++){
        /* Compute */
        core_zpotrfsp1d( datacode, cblk, threshold, work );
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( datacode, "potrf_L.txt" );
#endif

    memFree_null( work );
}

void
thread_pzpotrf( isched_thread_t *ctx, void *args )
{
    sopalin_data_t     *sopalin_data = (sopalin_data_t*)args;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    Task               *t;
    pastix_complex64_t *work;
    pastix_int_t  i, ii;
    pastix_int_t  tasknbr, *tasktab;
    int rank = ctx->rank;

    MALLOC_INTERN( work, datacode->gemmmax, pastix_complex64_t );

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        /* Wait */
        do {} while( cblk->ctrbcnt );

        /* Compute */
        core_zpotrfsp1d( datacode, cblk, sopalin_data->diagthreshold, work );
    }

#if defined(PASTIX_DEBUG_FACTO) && 0
    isched_barrier_wait( &(ctx->global_ctx->barrier) );
    if (rank == 0) {
        coeftab_zdump( datacode, "potrf_L.txt" );
    }
    isched_barrier_wait( &(ctx->global_ctx->barrier) );
#endif

    memFree_null( work );
}

void
thread_zpotrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    isched_parallel_call( pastix_data->isched, thread_pzpotrf, sopalin_data );
}

#if defined(PASTIX_WITH_PARSEC)
void
parsec_zpotrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    dague_context_t *ctx;

    /* Start PaRSEC */
    if (pastix_data->parsec == NULL) {
        int argc = 0;
        pastix_parsec_init( pastix_data, &argc, NULL );
    }
    ctx = pastix_data->parsec;

    /* Run the facto */
    dsparse_zpotrf_sp( ctx, sopalin_data->solvmtx->parsec_desc, sopalin_data );
}
#endif

static void (*zpotrf_table[4])(pastix_data_t *, sopalin_data_t *) = {
    sequential_zpotrf,
    thread_zpotrf,
#if defined(PASTIX_WITH_PARSEC)
    parsec_zpotrf,
#else
    NULL,
#endif
    NULL
};

void
sopalin_zpotrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zpotrf)(pastix_data_t *, sopalin_data_t *) = zpotrf_table[ sched ];

    if (zpotrf == NULL) {
        zpotrf = thread_zpotrf;
    }
    zpotrf( pastix_data, sopalin_data );
}
