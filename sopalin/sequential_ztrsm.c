/**
 *
 * @file sequential_ztrsm.c
 *
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Vincent Bridonneau
 * @author Alycia Lisito
 * @date 2024-07-05
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "blend/solver.h"
#include "bcsc/bcsc.h"
#include "bcsc/bcsc_z.h"
#include "sopalin/sopalin_data.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_STARPU)
#include "starpu/pastix_zstarpu.h"
#endif

#if defined(PASTIX_WITH_MPI)
#include "sopalin/coeftab.h"
#endif

/**
 * @brief Arguments for the solve.
 */
struct args_ztrsm_t
{
    pastix_data_t      *pastix_data;
    const args_solve_t *enum_list;
    sopalin_data_t     *sopalin_data;
    pastix_rhs_t        rhsb;
    volatile int32_t    taskcnt;
};

/**
 *******************************************************************************
 *
 * @brief Applies the Sequential Forward or Backward solve.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] enums
 *          Enums needed for the solve.
 *
 * @param[in] sopalin_data
 *          The SolverMatrix structure from PaStiX.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
sequential_ztrsm( pastix_data_t  *pastix_data,
                  const args_solve_t   *enums,
                  sopalin_data_t *sopalin_data,
                  pastix_rhs_t    rhsb )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk   *cblk;
    pastix_int_t  i, cblknbr;

    /* Backward like */
    if ( enums->solve_step == PastixSolveBackward ) {
        cblknbr = (enums->mode == PastixSolvModeLocal) ? datacode->cblkschur : datacode->cblknbr;

        cblk = datacode->cblktab + cblknbr - 1;
        for (i=0; i<cblknbr; i++, cblk--){
            if( cblk->cblktype & CBLK_RECV ){
                cpucblk_zsend_rhs_backward( datacode, cblk, rhsb );
                continue;
            }

            if( cblk->cblktype & CBLK_FANIN ){
                cpucblk_zrecv_rhs_backward( datacode, cblk, rhsb );
            }

            solve_cblk_ztrsmsp_backward( enums, datacode, cblk, rhsb );
        }
    }
    /* Forward like */
    else {
        pastix_complex64_t *work;
        MALLOC_INTERN( work, datacode->colmax * rhsb->n, pastix_complex64_t );

        cblknbr = (enums->mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;
        cblk = datacode->cblktab;
        for (i=0; i<cblknbr; i++, cblk++){
            if( cblk->cblktype & CBLK_FANIN ){
                cpucblk_zsend_rhs_forward( datacode, cblk, rhsb );
                continue;
            }

            if( cblk->cblktype & CBLK_RECV ) {
                cpucblk_zrecv_rhs_forward( datacode, cblk, work, rhsb );
                continue;
            }

            solve_cblk_ztrsmsp_forward( enums, datacode, cblk, rhsb );
        }

        memFree_null(work);
    }

#if !defined(NDEBUG)
    {
        pastix_int_t nbbuffers = datacode->faninnbr + datacode->recvnbr;
        int i;
        for( i=0; i<nbbuffers; i++ ) {
            assert( rhsb->cblkb[i] == NULL );
        }
    }
#endif
    (void)pastix_data;
}

/**
 *******************************************************************************
 *
 * @brief Applies the Static Forward or Backward solve.
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          Thread structure of the execution context of one instance of the
 *          scheduler.
 *
 * @param[in] args
 *          Arguments for the Static solve.
 *
 *******************************************************************************/
void
thread_ztrsm_static( isched_thread_t *ctx,
                     void            *args )
{
    struct args_ztrsm_t *arg          = (struct args_ztrsm_t*)args;
    sopalin_data_t      *sopalin_data = arg->sopalin_data;
    SolverMatrix        *datacode     = sopalin_data->solvmtx;
    pastix_rhs_t         rhsb         = arg->rhsb;
    const args_solve_t  *enums        = arg->enum_list;
    pastix_int_t         thrd_size    = (pastix_int_t)ctx->global_ctx->world_size;
    pastix_int_t         thrd_rank    = (pastix_int_t)ctx->rank;
    SolverCblk          *cblk;
    Task                *t;
    pastix_int_t         i, ii;
    pastix_int_t         tasknbr, *tasktab;
    pastix_int_t         cblkfirst, cblklast;

    /* Computes range to update the ctrbnbr */
    cblkfirst = (datacode->cblknbr / thrd_size ) * thrd_rank;
    cblklast  = (datacode->cblknbr / thrd_size ) * (thrd_rank + 1);
    if ( thrd_rank == (thrd_size-1) ) {
        cblklast = datacode->cblknbr;
    }

    tasknbr = datacode->ttsknbr[thrd_rank];
    tasktab = datacode->ttsktab[thrd_rank];

    /* Backward like */
    if ( enums->solve_step == PastixSolveBackward ) {
        /* Init ctrbcnt in parallel */
        cblk = datacode->cblktab + cblkfirst;
        for (ii=cblkfirst; ii<cblklast; ii++, cblk++) {
            if ( (cblk->cblktype & CBLK_IN_SCHUR) && (enums->mode != PastixSolvModeSchur) ) {
                cblk->ctrbcnt = 0;
            }
            else {
                cblk->ctrbcnt = cblk[1].fblokptr - cblk[0].fblokptr - 1;
            }
        }
        isched_barrier_wait( &(ctx->global_ctx->barrier) );

        for (ii=tasknbr-1; ii>=0; ii--) {
            i = tasktab[ii];
            t = datacode->tasktab + i;
            cblk = datacode->cblktab + t->cblknum;

            /* Wait for incoming dependencies */
            if ( cpucblk_zincoming_rhs_bwd_deps( thrd_rank, enums, datacode, cblk, rhsb ) ) {
                continue;
            }

            /* Computes */
            solve_cblk_ztrsmsp_backward( enums, datacode, cblk, rhsb );
        }
    }
    /* Forward like */
    else {
        /* Init ctrbcnt in parallel */
        cblk = datacode->cblktab + cblkfirst;
        for (ii=cblkfirst; ii<cblklast; ii++, cblk++) {
            cblk->ctrbcnt = cblk[1].brownum - cblk[0].brownum;
        }
        isched_barrier_wait( &(ctx->global_ctx->barrier) );

        for (ii=0; ii<tasknbr; ii++) {
            i = tasktab[ii];
            t = datacode->tasktab + i;
            cblk = datacode->cblktab + t->cblknum;

            if ( (cblk->cblktype & CBLK_IN_SCHUR) &&
                 (enums->mode != PastixSolvModeSchur) ) {
                continue;
            }

            /* Wait for incoming dependencies */
            if ( cpucblk_zincoming_rhs_fwd_deps( thrd_rank, enums,
                                                 datacode, cblk, rhsb ) ) {
                continue;
            }
            /* Computes */
            solve_cblk_ztrsmsp_forward( enums, datacode, cblk, rhsb );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Applies the Static Forward or Backward solve.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] enums
 *          Enums needed for the solve.
 *
 * @param[in] sopalin_data
 *          The SolverMatrix structure from PaStiX.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
static_ztrsm( pastix_data_t      *pastix_data,
              const args_solve_t *enums,
              sopalin_data_t     *sopalin_data,
              pastix_rhs_t        rhsb  )
{
    struct args_ztrsm_t args_ztrsm = { pastix_data, enums, sopalin_data, rhsb, 0 };
    isched_parallel_call( pastix_data->isched, thread_ztrsm_static, &args_ztrsm );
}

/**
 *******************************************************************************
 *
 * @brief Applies the Dynamic Forward or Backward solve.
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          Thread structure of the execution context of one instance of the
 *          scheduler.
 *
 * @param[in] args
 *          Arguments for the Static solve.
 *
 *******************************************************************************/
void
thread_ztrsm_dynamic( isched_thread_t *ctx,
                      void            *args )
{
    struct args_ztrsm_t *arg           = (struct args_ztrsm_t*)args;
    pastix_data_t       *pastix_data   = arg->pastix_data;
    sopalin_data_t      *sopalin_data  = arg->sopalin_data;
    SolverMatrix        *datacode      = sopalin_data->solvmtx;
    const args_solve_t  *enums         = arg->enum_list;
    pastix_rhs_t         rhsb          = arg->rhsb;
    pastix_int_t         thrd_size     = (pastix_int_t)ctx->global_ctx->world_size;
    pastix_int_t         thrd_rank     = (pastix_int_t)ctx->rank;
    int32_t              local_taskcnt = 0;
    SolverCblk          *cblk;
    pastix_queue_t      *computeQueue;
    pastix_int_t         ii;
    pastix_int_t         tasknbr;
    pastix_int_t         cblkfirst, cblklast, cblknum;

    /* Computes range to update the ctrbnbr */
    cblkfirst = (datacode->cblknbr / thrd_size ) * thrd_rank;
    cblklast  = (datacode->cblknbr / thrd_size ) * (thrd_rank + 1);
    if ( thrd_rank == (thrd_size-1) ) {
        cblklast = datacode->cblknbr;
    }

    MALLOC_INTERN( datacode->computeQueue[thrd_rank], 1, pastix_queue_t );

    tasknbr      = datacode->ttsknbr[thrd_rank];
    computeQueue = datacode->computeQueue[thrd_rank];
    pqueueInit( computeQueue, tasknbr );

    /* Backward like */
    if ( enums->solve_step == PastixSolveBackward ) {
        /* Init ctrbcnt in parallel */
        cblk = datacode->cblktab + cblkfirst;
        for (ii=cblkfirst; ii<cblklast; ii++, cblk++) {
            if ( (cblk->cblktype & CBLK_IN_SCHUR) && (enums->mode != PastixSolvModeSchur) ) {
                cblk->ctrbcnt = 0;
            }
            else {
                cblk->ctrbcnt = cblk[1].fblokptr - cblk[0].fblokptr - 1;
            }
            if ( !(cblk->ctrbcnt) && !(cblk->cblktype & (CBLK_FANIN | CBLK_RECV)) ) {
                pqueuePush1( computeQueue, ii, - cblk->priority );
            }
        }
        isched_barrier_wait( &(ctx->global_ctx->barrier) );

        while( arg->taskcnt > 0 ) {
            cblknum = pqueuePop(computeQueue);

#if defined(PASTIX_WITH_MPI)
            /* Nothing to do, let's make progress on communications */
            if ( ( pastix_data->inter_node_procnbr > 1 ) && ( cblknum == -1 ) ) {
                cpucblk_zmpi_rhs_bwd_progress( enums, datacode, rhsb, thrd_rank );
                cblknum = pqueuePop(computeQueue);
            }
#endif

            /* No more local job, let's steal our neighbors */
            if ( cblknum == -1 ) {
                if ( local_taskcnt ) {
                    pastix_atomic_sub_32b( &(arg->taskcnt), local_taskcnt );
                    local_taskcnt = 0;
                }
                cblknum = stealQueue( datacode, thrd_rank, thrd_size );
            }

            /* Still no job, let's loop again */
            if ( cblknum == -1 ) {
                continue;
            }

            cblk           = datacode->cblktab + cblknum;
            cblk->threadid = thrd_rank;

            /* Computes */
            solve_cblk_ztrsmsp_backward( enums, datacode, cblk, rhsb );
            local_taskcnt++;
        }
    }
    /* Forward like */
    else {
        /* Init ctrbcnt in parallel */
        cblk = datacode->cblktab + cblkfirst;
        for (ii=cblkfirst; ii<cblklast; ii++, cblk++) {
            cblk->ctrbcnt = cblk[1].brownum - cblk[0].brownum;
            if ( !(cblk->ctrbcnt) ) {
                if  (!(cblk->cblktype & (CBLK_FANIN|CBLK_RECV)) ) {
                    pqueuePush1( computeQueue, ii, cblk->priority );
                }
            }
        }
        isched_barrier_wait( &(ctx->global_ctx->barrier) );

        while( arg->taskcnt > 0 ) {
            cblknum = pqueuePop(computeQueue);

#if defined(PASTIX_WITH_MPI)
            /* Nothing to do, let's make progress on communications */
            if ( ( pastix_data->inter_node_procnbr > 1 ) && ( cblknum == -1 ) ) {
                cpucblk_zmpi_rhs_fwd_progress( enums, datacode, rhsb, thrd_rank );
                cblknum = pqueuePop(computeQueue);
            }
#endif

            /* No more local job, let's steal our neighbors */
            if ( cblknum == -1 ) {
                if ( local_taskcnt ) {
                    pastix_atomic_sub_32b( &(arg->taskcnt), local_taskcnt );
                    local_taskcnt = 0;
                }
                cblknum = stealQueue( datacode, thrd_rank,
                                      thrd_size );
            }

            /* Still no job, let's loop again */
            if ( cblknum == -1 ) {
                continue;
            }

            cblk           = datacode->cblktab + cblknum;
            cblk->threadid = thrd_rank;

            if ( (cblk->cblktype & CBLK_IN_SCHUR) &&
                 (enums->mode != PastixSolvModeSchur) ) {
                continue;
            }

            /* Computes */
            solve_cblk_ztrsmsp_forward( enums, datacode, cblk, rhsb );
            local_taskcnt++;
        }
    }
    /* Make sure that everyone is done before freeing */
    isched_barrier_wait( &(ctx->global_ctx->barrier) );
    assert( computeQueue->used == 0 );
    pqueueExit( computeQueue );
    memFree_null( computeQueue );

    (void)pastix_data;
}

/**
 *******************************************************************************
 *
 * @brief Applies the Dynamic Forward or Backward solve.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] enums
 *          Enums needed for the solve.
 *
 * @param[in] sopalin_data
 *          The SolverMatrix structure from PaStiX.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
dynamic_ztrsm( pastix_data_t      *pastix_data,
               const args_solve_t *enums,
               sopalin_data_t     *sopalin_data,
               pastix_rhs_t        rhsb  )
{
    SolverMatrix        *datacode   = sopalin_data->solvmtx;
    int32_t              taskcnt    = datacode->tasknbr - (datacode->cblknbr - datacode->cblkschur);
    struct args_ztrsm_t  args_ztrsm = { pastix_data, enums, sopalin_data, rhsb, taskcnt };

    /* Reintroduce Schur tasks in the counter for backward */
    if ( enums->solve_step == PastixSolveBackward ) {
        args_ztrsm.taskcnt = datacode->cblknbr - datacode->recvnbr;
    }

    /* Allocates the computeQueue */
    MALLOC_INTERN( datacode->computeQueue,
                   pastix_data->isched->world_size, pastix_queue_t * );

    isched_parallel_call( pastix_data->isched, thread_ztrsm_dynamic, &args_ztrsm );

    memFree_null( datacode->computeQueue );
}

#if defined(PASTIX_WITH_MPI)
/**
 *******************************************************************************
 *
 * @brief Applies the Reuntime Forward or Backward solve.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] enums
 *          Enums needed for the solve.
 *
 * @param[in] sopalin_data
 *          The SolverMatrix structure from PaStiX.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
runtime_ztrsm( pastix_data_t      *pastix_data,
               const args_solve_t *enums,
               sopalin_data_t     *sopalin_data,
               pastix_rhs_t        rhsb )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk   *cblk;
    pastix_int_t  i, cblknbr;

    /* Collect the matrix on node 0 */
    coeftab_gather( datacode, datacode->solv_comm, 0, PastixComplex64 );

    if ( sopalin_data->solvmtx->clustnum == 0 ) {

        /* Backward like */
        if ( enums->solve_step == PastixSolveBackward ) {
            cblknbr = (enums->mode == PastixSolvModeLocal) ? datacode->cblkschur : datacode->cblknbr;

            cblk = datacode->cblktab + cblknbr - 1;
            for ( i=0; i<cblknbr; i++, cblk-- ) {
                assert( !(cblk->cblktype & (CBLK_FANIN | CBLK_RECV)) );
                solve_cblk_ztrsmsp_backward( enums, datacode, cblk, rhsb );
            }
        }
        /* Forward like */
        else {
            cblknbr = (enums->mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;
            cblk = datacode->cblktab;
            for (i=0; i<cblknbr; i++, cblk++){
                solve_cblk_ztrsmsp_forward( enums, datacode, cblk, rhsb );
            }
        }

        /* Free the gathered coefficients of the matrix */
        coeftab_nullify( datacode );
    }
    else {
        memset( rhsb->b, 0, rhsb->ld * rhsb->n * sizeof(pastix_complex64_t) );
    }

    bvec_zallreduce( pastix_data, rhsb->b );
}
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static void (*ztrsm_table[5])(pastix_data_t *, const args_solve_t *,
                              sopalin_data_t *, pastix_rhs_t) =
{
    sequential_ztrsm,
    static_ztrsm,
#if defined(PASTIX_WITH_PARSEC)
    NULL, /* parsec_ztrsm not yet implemented */
#else
    NULL,
#endif
#if defined(PASTIX_WITH_STARPU)
    starpu_ztrsm,
#else
    NULL,
#endif
    dynamic_ztrsm
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief Calls the sequential, static, dynamic or runtime solve according to
 * scheduler.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] side
 *          Specify whether the off-diagonal blocks appear on the left or right
 *          in the equation. It has to be either PastixLeft or PastixRight.
 *
 * @param[in] uplo
 *          Specify whether the off-diagonal blocks are upper or lower
 *          triangular. It has to be either PastixUpper or PastixLower.
 *
 * @param[in] trans
 *          Specify the transposition used for the off-diagonal blocks. It has
 *          to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 * @param[in] sopalin_data
 *          The SolverMatrix structure from PaStiX.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
sopalin_ztrsm( pastix_data_t  *pastix_data,
               pastix_side_t   side,
               pastix_uplo_t   uplo,
               pastix_trans_t  trans,
               pastix_diag_t   diag,
               sopalin_data_t *sopalin_data,
               pastix_rhs_t    rhsb  )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*ztrsm)( pastix_data_t *, const args_solve_t *,
                   sopalin_data_t *, pastix_rhs_t ) = ztrsm_table[ sched ];
    solve_step_t  solve_step = compute_solve_step( side, uplo, trans );
    args_solve_t *enum_list = malloc( sizeof( args_solve_t ) );

    enum_list->solve_step = solve_step;
    enum_list->mode       = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];
    enum_list->side       = side;
    enum_list->uplo       = uplo;
    enum_list->trans      = trans;
    enum_list->diag       = diag;

    if (ztrsm == NULL) {
        ztrsm = static_ztrsm;
    }

    /* parsec_ztrsm and starpu_ztrsm not implemented yet, runtime_ztrsm works only for starpu and
       parsec with mpi in distributed and replicated cases */
#if defined ( PASTIX_WITH_MPI )
    if( pastix_data->inter_node_procnbr > 1 ) {
        if( sched == PastixSchedParsec ) {
            ztrsm = runtime_ztrsm;
        }
    }
#endif

    if ( (sched == PastixSchedStatic)     ||
         (sched == PastixSchedDynamic) )
    {
        solverRequestInit( solve_step, sopalin_data->solvmtx );
        solverRhsRecvInit( solve_step, sopalin_data->solvmtx, PastixComplex64, rhsb );
    }

    enum_list->sched = sched;
    ztrsm( pastix_data, enum_list, sopalin_data, rhsb );

    if ( (sched == PastixSchedStatic)     ||
         (sched == PastixSchedDynamic) )
    {
        if ( solve_step == PastixSolveForward ) {
            cpucblk_zrequest_rhs_fwd_cleanup( enum_list, sched, sopalin_data->solvmtx, rhsb );
        }
        else {
            cpucblk_zrequest_rhs_bwd_cleanup( enum_list, sched, sopalin_data->solvmtx, rhsb );
        }
        solverRequestExit( sopalin_data->solvmtx );
        solverRhsRecvExit( sopalin_data->solvmtx );
    }

#if defined(PASTIX_WITH_MPI)
   MPI_Barrier( pastix_data->inter_node_comm );
#endif
    free(enum_list);
}
