/**
 *
 * @file sequential_zsytrf.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Esragul Korkmaz
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2023-11-07
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "isched.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"
#include "sopalin/coeftab_z.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_PARSEC)
#include "parsec/pastix_zparsec.h"
#endif

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
 *******************************************************************************/
void
sequential_zsytrf( pastix_data_t  *pastix_data,
                   sopalin_data_t *sopalin_data )
{
    SolverMatrix       *datacode = pastix_data->solvmatr;
    SolverCblk         *cblk;
    pastix_complex64_t *work1, *work2;
    pastix_int_t  N, i, lwork1, lwork2;
    (void)sopalin_data;

    lwork1 = datacode->offdmax;
    lwork2 = pastix_imax( datacode->gemmmax, datacode->blokmax );
    if ( (datacode->lowrank.compress_when != PastixCompressNever) &&
         (datacode->lowrank.ilu_lvl < INT_MAX) )
    {
        lwork2 = pastix_imax( lwork2, 2 * datacode->blokmax );
    }
    MALLOC_INTERN( work1, lwork1, pastix_complex64_t );
    MALLOC_INTERN( work2, lwork2, pastix_complex64_t );

    cblk = datacode->cblktab;
    for (i=0; i<datacode->cblknbr; i++, cblk++){
        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            break;
        }

        /* Wait for incoming dependencies */
        if ( cpucblk_zincoming_deps( 0, PastixLCoef,
                                     datacode, cblk ) )
        {
            continue;
        }

        N = cblk_colnbr( cblk );

        /* Compute */
        cpucblk_zsytrfsp1d( datacode, cblk,
                            /*
                             * Workspace size has been computed without the
                             * diagonal block, thus in order to work with generic
                             * TRSM and GEMM kernels, we must shift the DLh workspace
                             * by the diagonal block size
                             */
                            work1 - (N*N), work2, lwork2 );
    }

    memFree_null( work1 );
    memFree_null( work2 );
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
thread_zsytrf_static( isched_thread_t *ctx,
                      void            *args )
{
    sopalin_data_t     *sopalin_data = (sopalin_data_t*)args;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    Task               *t;
    pastix_complex64_t *work1, *work2;
    pastix_int_t  N, i, ii, lwork1, lwork2;
    pastix_int_t  tasknbr, *tasktab;
    int rank = ctx->rank;

    lwork1 = datacode->offdmax;
    lwork2 = pastix_imax( datacode->gemmmax, datacode->blokmax );
    if ( (datacode->lowrank.compress_when != PastixCompressNever) &&
         (datacode->lowrank.ilu_lvl < INT_MAX) )
    {
        lwork2 = pastix_imax( lwork2, 2 * datacode->blokmax );
    }
    MALLOC_INTERN( work1, lwork1, pastix_complex64_t );
    MALLOC_INTERN( work2, lwork2, pastix_complex64_t );

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            continue;
        }

        /* Wait for incoming dependencies */
        if ( cpucblk_zincoming_deps( rank, PastixLCoef,
                                     datacode, cblk ) )
        {
            continue;
        }

        N = cblk_colnbr( cblk );

        /* Compute */
        cpucblk_zsytrfsp1d( datacode, cblk,
                            /*
                             * Workspace size has been computed without the
                             * diagonal block, thus in order to work with generic
                             * TRSM and GEMM kernels, we must shift the DLh workspace
                             * by the diagonal block size
                             */
                            work1 - (N*N), work2, lwork2 );
    }

    memFree_null( work1 );
    memFree_null( work2 );
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
 *******************************************************************************/
void
static_zsytrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    isched_parallel_call( pastix_data->isched, thread_zsytrf_static, sopalin_data );
}

/**
 * @brief TODO
 */
struct args_zsytrf_t
{
    sopalin_data_t     *sopalin_data;
    volatile int32_t    taskcnt;
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
thread_zsytrf_dynamic( isched_thread_t *ctx,
                       void            *args )
{
    struct args_zsytrf_t *arg = (struct args_zsytrf_t *)args;
    sopalin_data_t       *sopalin_data = arg->sopalin_data;
    SolverMatrix         *datacode = sopalin_data->solvmtx;
    SolverCblk           *cblk;
    SolverBlok           *blok;
    Task                 *t;
    pastix_queue_t       *computeQueue;
    pastix_complex64_t   *work1, *work2;
    pastix_int_t          N, i, ii, lwork1, lwork2;
    pastix_int_t          tasknbr, *tasktab, cblknum, bloknum;
    int32_t               local_taskcnt = 0;
    int                   rank = ctx->rank;

    lwork1 = datacode->offdmax;
    lwork2 = pastix_imax( datacode->gemmmax, datacode->blokmax );
    if ( (datacode->lowrank.compress_when != PastixCompressNever) &&
         (datacode->lowrank.ilu_lvl < INT_MAX) )
    {
        lwork2 = pastix_imax( lwork2, 2 * datacode->blokmax );
    }
    MALLOC_INTERN( work1, lwork1, pastix_complex64_t );
    MALLOC_INTERN( work2, lwork2, pastix_complex64_t );
    MALLOC_INTERN( datacode->computeQueue[rank], 1, pastix_queue_t );

    tasknbr      = datacode->ttsknbr[rank];
    tasktab      = datacode->ttsktab[rank];
    computeQueue = datacode->computeQueue[rank];
    pqueueInit( computeQueue, tasknbr );

    /* Initialize the local task queue with available cblks */
    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;

        if ( !(t->ctrbcnt) ) {
            cblk = datacode->cblktab + t->cblknum;
            pqueuePush1( computeQueue, t->cblknum, cblk->priority );
        }
    }

    /* Make sure that all computeQueues are allocated */
    isched_barrier_wait( &(ctx->global_ctx->barrier) );

    while( arg->taskcnt > 0 )
    {
        cblknum = pqueuePop(computeQueue);

#if defined(PASTIX_WITH_MPI)
        /* Nothing to do, let's make progress on communications */
        if( cblknum == -1 ) {
            cpucblk_zmpi_progress( PastixLCoef, datacode, rank );
            cblknum = pqueuePop( computeQueue );
        }
#endif

        /* No more local job, let's steal our neighbors */
        if( cblknum == -1 ) {
            if ( local_taskcnt ) {
                pastix_atomic_sub_32b( &(arg->taskcnt), local_taskcnt );
                local_taskcnt = 0;
            }
            cblknum = stealQueue( datacode, rank,
                                  ctx->global_ctx->world_size );
        }

        /* Still no job, let's loop again */
        if ( cblknum == -1 ) {
            continue;
        }

        if ( cblknum >= 0 ) {
            cblk = datacode->cblktab + cblknum;
            if ( cblk->cblktype & CBLK_IN_SCHUR ) {
                continue;
            }
            cblk->threadid = rank;

            N = cblk_colnbr( cblk );

            /* Compute */
            if ( cblk->cblktype & CBLK_TASKS_2D ) {
                cpucblk_zsytrfsp1dplus( datacode, cblk );
            }
            else {
                cpucblk_zsytrfsp1d( datacode, cblk,
                                    /*
                                    * Workspace size has been computed without the
                                    * diagonal block, thus in order to work with generic
                                    * TRSM and GEMM kernels, we must shift the DLh workspace
                                    * by the diagonal block size
                                    */
                                    work1 - (N*N), work2, lwork2 );
            }
        }
        else {
            bloknum = - cblknum - 1;
            blok    = datacode->bloktab + bloknum;
            cpucblk_zsytrfsp1dplus_update( datacode, blok, work2 );
        }
        local_taskcnt++;
    }
    memFree_null( work1 );
    memFree_null( work2 );

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
 *******************************************************************************/
void
dynamic_zsytrf( pastix_data_t  *pastix_data,
                 sopalin_data_t *sopalin_data )
{
    SolverMatrix         *datacode    = sopalin_data->solvmtx;
    int32_t               taskcnt     = datacode->tasknbr_1dp;
    struct args_zsytrf_t  args_zsytrf = { sopalin_data, taskcnt };

    /* Allocate the computeQueue */
    MALLOC_INTERN( datacode->computeQueue,
                   pastix_data->isched->world_size, pastix_queue_t * );

    isched_parallel_call( pastix_data->isched, thread_zsytrf_dynamic, &args_zsytrf );

    memFree_null( datacode->computeQueue );
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static void (*zsytrf_table[5])(pastix_data_t *, sopalin_data_t *) = {
    sequential_zsytrf,
    static_zsytrf,
/* #if defined(PASTIX_WITH_PARSEC) */
/*     parsec_zsytrf, */
/* #else */
    NULL,
/* #endif */
#if defined(PASTIX_WITH_STARPU)
    starpu_zsytrf,
#else
    NULL,
#endif
    dynamic_zsytrf
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
 *******************************************************************************/
void
sopalin_zsytrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zsytrf)(pastix_data_t *, sopalin_data_t *) = zsytrf_table[ sched ];

    if ( zsytrf == NULL ) {
        sched  = PastixSchedDynamic;
        zsytrf = dynamic_zsytrf;
    }

    if ( (sched == PastixSchedSequential) ||
         (sched == PastixSchedStatic)     ||
         (sched == PastixSchedDynamic) )
    {
        solverRequestInit( PastixFacto, sopalin_data->solvmtx );
        solverRecvInit( PastixLCoef, sopalin_data->solvmtx, PastixComplex64 );
    }

    zsytrf( pastix_data, sopalin_data );

    if ( (sched == PastixSchedSequential) ||
         (sched == PastixSchedStatic)     ||
         (sched == PastixSchedDynamic) )
    {
        cpucblk_zrequest_cleanup( PastixLCoef, sched, sopalin_data->solvmtx );
        solverRequestExit( sopalin_data->solvmtx );
        solverRecvExit( sopalin_data->solvmtx );
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( pastix_data, sopalin_data->solvmtx, "sytrf" );
#endif
}
