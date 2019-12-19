/**
 *
 * @file sequential_zsytrf.c
 *
 * @copyright 2012-2019 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "isched.h"
#include "solver.h"
#include "sopalin_data.h"
#include "sopalin/coeftab_z.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_PARSEC)
#include "parsec/pastix_zparsec.h"
#endif

#if defined(PASTIX_WITH_STARPU)
#include "starpu/pastix_zstarpu.h"
#endif

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
    if ( datacode->lowrank.compress_when == PastixCompressWhenBegin ) {
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

void
thread_zsytrf_static( isched_thread_t *ctx, void *args )
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
    if ( datacode->lowrank.compress_when == PastixCompressWhenBegin ) {
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
        if ( cpucblk_zincoming_deps( 1, PastixLCoef,
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

void
static_zsytrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    isched_parallel_call( pastix_data->isched, thread_zsytrf_static, sopalin_data );
}

struct args_zsytrf_t
{
    sopalin_data_t     *sopalin_data;
    volatile int32_t    taskcnt;
    pastix_complex64_t *recv;
};

void
thread_zsytrf_dynamic( isched_thread_t *ctx, void *args )
{
    struct args_zsytrf_t *arg = (struct args_zsytrf_t *)args;
    sopalin_data_t       *sopalin_data = arg->sopalin_data;
    SolverMatrix         *datacode = sopalin_data->solvmtx;
    SolverCblk           *cblk;
    Task                 *t;
    pastix_queue_t       *computeQueue;
    pastix_complex64_t   *work1, *work2;
    pastix_int_t          N, i, ii, lwork1, lwork2;
    pastix_int_t          tasknbr, *tasktab, cblknum;
    int32_t               local_taskcnt = 0;
    int                   rank = ctx->rank;
    int                   dest = (ctx->rank + 1)%ctx->global_ctx->world_size;
#if defined(PASTIX_WITH_MPI)
    pastix_complex64_t   *recv = arg->recv;
#endif

    lwork1 = datacode->offdmax;
    lwork2 = pastix_imax( datacode->gemmmax, datacode->blokmax );
    if ( datacode->lowrank.compress_when == PastixCompressWhenBegin ) {
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
            pqueuePush1( computeQueue, t->cblknum, t->prionum );
        }
    }

    /* Make sure that all computeQueues are allocated */
    isched_barrier_wait( &(ctx->global_ctx->barrier) );

    while( arg->taskcnt > 0 )
    {
        cblknum = pqueuePop(computeQueue);

#if defined(PASTIX_WITH_MPI)
        /* Nothing to do, let's make progress on comunications */
        if( cblknum == -1 ) {
            cpucblk_zmpi_progress( PastixLCoef, ctx, datacode, rank, recv );
            cblknum = pqueuePop(computeQueue);
        }
#endif

        /* No more local job, let's steal our neighbours */
        if( cblknum == -1 ) {
            if ( local_taskcnt ) {
                pastix_atomic_sub_32b( &(arg->taskcnt), local_taskcnt );
                local_taskcnt = 0;
            }
            cblknum = stealQueue( datacode, rank, &dest,
                                  ctx->global_ctx->world_size );
        }

        /* Still no job, let's loop again */
        if ( cblknum == -1 ) {
            continue;
        }

        cblk = datacode->cblktab + cblknum;
        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            continue;
        }
        cblk->threadid = rank;

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
        local_taskcnt++;
    }
    memFree_null( work1 );
    memFree_null( work2 );

    /* Make sure that everyone is done before freeing */
    isched_barrier_wait( &(ctx->global_ctx->barrier) );
    pqueueExit( computeQueue );
    memFree_null( computeQueue );
}

void
dynamic_zsytrf( pastix_data_t  *pastix_data,
                 sopalin_data_t *sopalin_data )
{
    SolverMatrix        *datacode = sopalin_data->solvmtx;
    int32_t              taskcnt = datacode->tasknbr;
    struct args_zsytrf_t args_zsytrf = { sopalin_data, taskcnt, NULL };
    /* Allocate the computeQueue */
    MALLOC_INTERN( datacode->computeQueue,
                   pastix_data->isched->world_size, pastix_queue_t * );

#if defined(PASTIX_WITH_MPI)
    MALLOC_INTERN( args_zsytrf.recv, datacode->maxrecv, pastix_complex64_t);

    /* Initialize the persistant communication */
    if( datacode->recvnbr ){
        MPI_Recv_init( args_zsytrf.recv, datacode->maxrecv,
                       PASTIX_MPI_COMPLEX64, MPI_ANY_SOURCE, MPI_ANY_TAG,
                       datacode->solv_comm, datacode->reqtab );
        MPI_Start( datacode->reqtab );
        datacode->reqnum++;
    }
#endif

    isched_parallel_call( pastix_data->isched, thread_zsytrf_dynamic, &args_zsytrf );

#if defined(PASTIX_WITH_MPI)
    memFree_null( args_zsytrf.recv );
    MPI_Barrier( pastix_data->inter_node_comm );
#endif
}

static void (*zsytrf_table[5])(pastix_data_t *, sopalin_data_t *) = {
    sequential_zsytrf,
    static_zsytrf,
#if defined(PASTIX_WITH_PARSEC)
    parsec_zsytrf,
#else
    NULL,
#endif
#if defined(PASTIX_WITH_STARPU)
    starpu_zsytrf,
#else
    NULL,
#endif
    dynamic_zsytrf
};

void
sopalin_zsytrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zsytrf)(pastix_data_t *, sopalin_data_t *) = zsytrf_table[ sched ];

    if (zsytrf == NULL) {
        sched = PastixSchedSequential;
        zsytrf = static_zsytrf;
    }

    if ( (sched == PastixSchedSequential) ||
         (sched == PastixSchedStatic)     ||
         (sched == PastixSchedDynamic) )
    {
        solverReqtabInit( sopalin_data->solvmtx, sched );
    }

    zsytrf( pastix_data, sopalin_data );

    if ( (sched == PastixSchedSequential) ||
         (sched == PastixSchedStatic)     ||
         (sched == PastixSchedDynamic) )
    {
        cpucblk_zrequest_cleanup( PastixLCoef, sched, sopalin_data->solvmtx );
        solverReqtabExit( sopalin_data->solvmtx );
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( pastix_data, sopalin_data->solvmtx, "sytrf.txt" );
#endif
}
