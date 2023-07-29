/**
 *
 * @file sequential_zdiag.c
 *
 * @copyright 2012-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Vincent Bridonneau
 * @date 2023-07-20
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_STARPU)
#include "starpu/pastix_zstarpu.h"
#endif

/**
 *******************************************************************************
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          TODO
 *
 * @param[in] sopalin_data
 *          TODO
 *
 * @param[in] nrhs
 *          TODO
 *
 * @param[in] b
 *          TODO
 *
 * @param[in] ldb
 *          TODO
 *
 *******************************************************************************/
void
sequential_zdiag( pastix_data_t      *pastix_data,
                  sopalin_data_t     *sopalin_data,
                  int                 nrhs,
                  pastix_complex64_t *b,
                  int                 ldb )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk   *cblk;
    pastix_int_t  i, cblknbr;
    pastix_solv_mode_t mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];

    cblk = datacode->cblktab;
    cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;
    for (i=0; i<cblknbr; i++, cblk++) {
        if ( cblk->ownerid != datacode->clustnum ) {
            continue;
        }
        solve_cblk_zdiag( cblk, nrhs,
                          b + cblk->lcolidx, ldb, NULL );
    }
}

/**
 * @brief TODO
 */
struct args_zdiag_t
{
    pastix_data_t  *pastix_data;
    sopalin_data_t *sopalin_data;
    int nrhs;
    pastix_complex64_t *b;
    int ldb;
    volatile int32_t taskcnt;
};

/**
 *******************************************************************************
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          TODO
 *
 * @param[in] args
 *          TODO
 *
 *******************************************************************************/
void
thread_zdiag_static( isched_thread_t *ctx,
                     void            *args )
{
    struct args_zdiag_t *arg = (struct args_zdiag_t*)args;
    pastix_data_t      *pastix_data  = arg->pastix_data;
    sopalin_data_t     *sopalin_data = arg->sopalin_data;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    pastix_complex64_t *b = arg->b;
    int nrhs  = arg->nrhs;
    int ldb   = arg->ldb;
    SolverCblk *cblk;
    Task       *t;
    pastix_int_t i, ii, cblknbr;
    pastix_int_t tasknbr, *tasktab;
    pastix_solv_mode_t mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];
    int rank = ctx->rank;

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];
    cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;

        if ( t->cblknum >= cblknbr ) {
            continue;
        }
        cblk = datacode->cblktab + t->cblknum;
        if ( cblk->ownerid != datacode->clustnum ) {
            continue;
        }
        solve_cblk_zdiag( cblk, nrhs,
                          b + cblk->lcolidx, ldb, NULL );
    }
}

/**
 *******************************************************************************
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          TODO
 *
 * @param[in] sopalin_data
 *          TODO
 *
 * @param[in] nrhs
 *          TODO
 *
 * @param[in] b
 *          TODO
 *
 * @param[in] ldb
 *          TODO
 *
 *******************************************************************************/
void
static_zdiag( pastix_data_t      *pastix_data,
              sopalin_data_t     *sopalin_data,
              int                 nrhs,
              pastix_complex64_t *b,
              int                 ldb )
{
    struct args_zdiag_t args_zdiag = {pastix_data, sopalin_data, nrhs, b, ldb, 0};
    isched_parallel_call( pastix_data->isched, thread_zdiag_static, &args_zdiag );
}

/**
 *******************************************************************************
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          TODO
 *
 * @param[in] args
 *          TODO
 *
 *******************************************************************************/
void
thread_zdiag_dynamic( isched_thread_t *ctx,
                      void            *args )
{
    struct args_zdiag_t *arg = (struct args_zdiag_t*)args;
    pastix_data_t       *pastix_data  = arg->pastix_data;
    sopalin_data_t      *sopalin_data = arg->sopalin_data;
    SolverMatrix        *datacode = sopalin_data->solvmtx;
    pastix_complex64_t  *b = arg->b;
    pastix_queue_t      *computeQueue;
    pastix_int_t         i, ii, cblknbr;//, lcblknbr;
    pastix_int_t         tasknbr, *tasktab, cblknum;
    pastix_solv_mode_t   mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];
    SolverCblk *cblk;
    Task       *t;
    int32_t     local_taskcnt = 0;
    int         nrhs = arg->nrhs;
    int         ldb  = arg->ldb;
    int         rank = ctx->rank;

    MALLOC_INTERN( datacode->computeQueue[rank], 1, pastix_queue_t );

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];
    computeQueue = datacode->computeQueue[rank];
    pqueueInit( computeQueue, tasknbr );
    cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;

        if ( t->cblknum >= cblknbr ) {
            continue;
        }
        cblk = datacode->cblktab + t->cblknum;
        pqueuePush1( computeQueue, t->cblknum, cblk->priority );
    }

    /* Make sure that all computeQueues are allocated */
    isched_barrier_wait( &(ctx->global_ctx->barrier) );

    /* When distributed -> only locals cblks, changes necessay */
    while( arg->taskcnt > 0 ){
        cblknum = pqueuePop(computeQueue);

        if( cblknum == -1 ){
            if ( local_taskcnt ) {
                pastix_atomic_sub_32b( &(arg->taskcnt), local_taskcnt );
                local_taskcnt = 0;
            }
            cblknum = stealQueue( datacode, rank,
                                  ctx->global_ctx->world_size );
        }
        if( cblknum != -1 ){
            cblk = datacode->cblktab + cblknum;
            solve_cblk_zdiag( cblk, nrhs,
                              b + cblk->lcolidx, ldb, NULL );
            local_taskcnt++;
        }
    }
    /* Make sure that everyone is done before freeing */
    isched_barrier_wait( &(ctx->global_ctx->barrier) );
    pqueueExit( computeQueue );
    memFree_null( computeQueue );
}

/**
 *******************************************************************************
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          TODO
 *
 * @param[in] sopalin_data
 *          TODO
 *
 * @param[in] nrhs
 *          TODO
 *
 * @param[in] b
 *          TODO
 *
 * @param[in] ldb
 *          TODO
 *
 *******************************************************************************/
void
dynamic_zdiag( pastix_data_t      *pastix_data,
               sopalin_data_t     *sopalin_data,
               int                 nrhs,
               pastix_complex64_t *b,
               int                 ldb )
{
    pastix_solv_mode_t  mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];
    pastix_int_t        tasknbr = (mode == PastixSolvModeSchur) ?
                                  sopalin_data->solvmtx->cblknbr : sopalin_data->solvmtx->cblkschur;
    struct args_zdiag_t args_zdiag = { pastix_data, sopalin_data, nrhs, b, ldb, tasknbr };
    /* Allocate the computeQueue */
    MALLOC_INTERN( sopalin_data->solvmtx->computeQueue,
                   pastix_data->isched->world_size, pastix_queue_t * );

    isched_parallel_call( pastix_data->isched, thread_zdiag_static, &args_zdiag );
    memFree_null( sopalin_data->solvmtx->computeQueue );
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static void (*zdiag_table[5])(pastix_data_t *, sopalin_data_t *,
                              int, pastix_complex64_t *, int) = {
    sequential_zdiag,
    static_zdiag,
#if defined(PASTIX_WITH_PARSEC)
    NULL, /* parsec_zdiag not yet implemented */
#else
    NULL,
#endif
#if defined(PASTIX_WITH_STARPU)
    starpu_zdiag,
#else
    NULL,
#endif
    dynamic_zdiag
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          TODO
 *
 * @param[in] sopalin_data
 *          TODO
 *
 * @param[in] nrhs
 *          TODO
 *
 * @param[in] b
 *          TODO
 *
 * @param[in] ldb
 *          TODO
 *
 *******************************************************************************/
void
sopalin_zdiag( pastix_data_t      *pastix_data,
               sopalin_data_t     *sopalin_data,
               int                 nrhs,
               pastix_complex64_t *b,
               int                 ldb )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zdiag)(pastix_data_t *, sopalin_data_t *, int, pastix_complex64_t *, int) = zdiag_table[ sched ];

    if (zdiag == NULL) {
        zdiag = static_zdiag;
    }
    zdiag( pastix_data, sopalin_data, nrhs, b, ldb );
}
