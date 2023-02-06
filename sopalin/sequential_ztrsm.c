/**
 *
 * @file sequential_ztrsm.c
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
 * @date 2022-12-03
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
 *******************************************************************************
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          TODO
 *
 * @param[in] side
 *          TODO
 *
 * @param[in] uplo
 *          TODO
 *
 * @param[in] trans
 *          TODO
 *
 * @param[in] diag
 *          TODO
 *
 * @param[in] sopalin_data
 *          TODO
 *
 * @param[in] rhsb
 *          TODO
 *
 *******************************************************************************/
void
sequential_ztrsm( pastix_data_t  *pastix_data,
                  int             side,
                  int             uplo,
                  int             trans,
                  int             diag,
                  sopalin_data_t *sopalin_data,
                  pastix_rhs_t    rhsb )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk   *cblk;
    pastix_int_t i, cblknbr;

    pastix_solv_mode_t mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];

    /* Backward like */
    if ( ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
         ( (side == PastixLeft)  && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
         ( (side == PastixRight) && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
         ( (side == PastixRight) && (uplo == PastixLower) && (trans == PastixNoTrans) ) )
    {
        cblknbr = (mode == PastixSolvModeLocal) ? datacode->cblkschur : datacode->cblknbr;

        cblk = datacode->cblktab + cblknbr - 1;
        for (i=0; i<cblknbr; i++, cblk--){
            if( cblk->cblktype & CBLK_RECV ){
                cpucblk_zsend_rhs_backward( datacode, cblk, rhsb );
                continue;
            }

            if( cblk->cblktype & CBLK_FANIN ){
                cpucblk_zrecv_rhs_backward( datacode, cblk, rhsb );
            }

            solve_cblk_ztrsmsp_backward( mode, side, uplo, trans, diag,
                                         datacode, cblk, rhsb );
        }
    }
    /* Forward like */
    else
        /**
         * ( (side == PastixRight) && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
         * ( (side == PastixRight) && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
         * ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
         * ( (side == PastixLeft)  && (uplo == PastixLower) && (trans == PastixNoTrans) )
         */
    {
        pastix_complex64_t *work;
        MALLOC_INTERN( work, datacode->colmax * rhsb->n, pastix_complex64_t );

        cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;
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

            solve_cblk_ztrsmsp_forward( mode, side, uplo, trans, diag,
                                        datacode, cblk, rhsb );
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
}

/**
 * @brief TODO
 */
struct args_ztrsm_t
{
    pastix_data_t  *pastix_data;
    int side, uplo, trans, diag;
    sopalin_data_t *sopalin_data;
    pastix_rhs_t    rhsb;
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
thread_ztrsm_static( isched_thread_t *ctx,
                     void            *args )
{
    struct args_ztrsm_t *arg          = (struct args_ztrsm_t*)args;
    pastix_data_t       *pastix_data  = arg->pastix_data;
    sopalin_data_t      *sopalin_data = arg->sopalin_data;
    SolverMatrix        *datacode     = sopalin_data->solvmtx;
    pastix_rhs_t         rhsb         = arg->rhsb;
    int                  side         = arg->side;
    int                  uplo         = arg->uplo;
    int                  trans        = arg->trans;
    int                  diag         = arg->diag;
    pastix_solv_mode_t   mode         = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];
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
    if ( ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
         ( (side == PastixLeft)  && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
         ( (side == PastixRight) && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
         ( (side == PastixRight) && (uplo == PastixLower) && (trans == PastixNoTrans) ) )
    {
        /* Init ctrbcnt in parallel */
        cblk = datacode->cblktab + cblkfirst;
        for (ii=cblkfirst; ii<cblklast; ii++, cblk++) {
            if ( (cblk->cblktype & CBLK_IN_SCHUR) && (mode != PastixSolvModeSchur) ) {
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
            if ( cpucblk_zincoming_rhs_deps( thrd_rank, PastixSolveBackward,
                                             datacode, cblk, rhsb ) ) {
                continue;
            }

            /* Computes */
            solve_cblk_ztrsmsp_backward( mode, side, uplo, trans, diag,
                                         datacode, cblk, rhsb );
        }
    }
    /* Forward like */
    else
        /**
         * ( (side == PastixRight) && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
         * ( (side == PastixRight) && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
         * ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
         * ( (side == PastixLeft)  && (uplo == PastixLower) && (trans == PastixNoTrans) )
         */
    {
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
                 (mode != PastixSolvModeSchur) ) {
                continue;
            }

            /* Wait for incoming dependencies */
            if ( cpucblk_zincoming_rhs_deps( thrd_rank, PastixSolveForward,
                                             datacode, cblk, rhsb ) ) {
                continue;
            }
            /* Computes */
            solve_cblk_ztrsmsp_forward( mode, side, uplo, trans, diag,
                                        datacode, cblk, rhsb );
        }
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
 * @param[in] side
 *          TODO
 *
 * @param[in] uplo
 *          TODO
 *
 * @param[in] trans
 *          TODO
 *
 * @param[in] diag
 *          TODO
 *
 * @param[in] sopalin_data
 *          TODO
 *
 * @param[in] rhsb
 *          TODO
 *
 *******************************************************************************/
void
static_ztrsm( pastix_data_t  *pastix_data,
              int             side,
              int             uplo,
              int             trans,
              int             diag,
              sopalin_data_t *sopalin_data,
              pastix_rhs_t    rhsb  )
{
    struct args_ztrsm_t args_ztrsm = { pastix_data, side, uplo, trans, diag, sopalin_data, rhsb };
    isched_parallel_call( pastix_data->isched, thread_ztrsm_static, &args_ztrsm );
}

#if defined(PASTIX_WITH_MPI)
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
 * @param[in] side
 *          TODO
 *
 * @param[in] uplo
 *          TODO
 *
 * @param[in] trans
 *          TODO
 *
 * @param[in] diag
 *          TODO
 *
 * @param[in] sopalin_data
 *          TODO
 *
 * @param[in] rhsb
 *          TODO
 *
 *******************************************************************************/
void
runtime_ztrsm( pastix_data_t  *pastix_data,
               int             side,
               int             uplo,
               int             trans,
               int             diag,
               sopalin_data_t *sopalin_data,
               pastix_rhs_t    rhsb  )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk *cblk;
    pastix_int_t i, cblknbr;

    /* Collect the matrix on node 0 */
    coeftab_gather( datacode, datacode->solv_comm, 0, PastixComplex64 );

    if ( sopalin_data->solvmtx->clustnum == 0 ) {
        pastix_solv_mode_t mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];

        /* Backward like */
        if ( ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
             ( (side == PastixLeft)  && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
             ( (side == PastixRight) && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
             ( (side == PastixRight) && (uplo == PastixLower) && (trans == PastixNoTrans) ) )
        {
            cblknbr = (mode == PastixSolvModeLocal) ? datacode->cblkschur : datacode->cblknbr;

            cblk = datacode->cblktab + cblknbr - 1;
            for ( i=0; i<cblknbr; i++, cblk-- ) {
                assert( !(cblk->cblktype & (CBLK_FANIN | CBLK_RECV)) );
                solve_cblk_ztrsmsp_backward( mode, side, uplo, trans, diag,
                                             datacode, cblk, rhsb );
            }
        }
        /* Forward like */
        else
            /**
             * ( (side == PastixRight) && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
             * ( (side == PastixRight) && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
             * ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
             * ( (side == PastixLeft)  && (uplo == PastixLower) && (trans == PastixNoTrans) )
             */
        {
            cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;
            cblk = datacode->cblktab;
            for (i=0; i<cblknbr; i++, cblk++){
                solve_cblk_ztrsmsp_forward( mode, side, uplo, trans, diag,
                                            datacode, cblk, rhsb );
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
static void (*ztrsm_table[5])(pastix_data_t *, int, int, int, int,
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
    static_ztrsm
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
 * @param[in] side
 *          TODO
 *
 * @param[in] uplo
 *          TODO
 *
 * @param[in] trans
 *          TODO
 *
 * @param[in] diag
 *          TODO
 *
 * @param[in] sopalin_data
 *          TODO
 *
 * @param[in] rhsb
 *          TODO
 *
 *******************************************************************************/
void
sopalin_ztrsm( pastix_data_t  *pastix_data,
               int             side,
               int             uplo,
               int             trans,
               int             diag,
               sopalin_data_t *sopalin_data,
               pastix_rhs_t    rhsb  )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    int solve_step = ( ( (side == PastixLeft)  && (uplo == PastixUpper) && (trans == PastixNoTrans) ) ||
                       ( (side == PastixLeft)  && (uplo == PastixLower) && (trans != PastixNoTrans) ) ||
                       ( (side == PastixRight) && (uplo == PastixUpper) && (trans != PastixNoTrans) ) ||
                       ( (side == PastixRight) && (uplo == PastixLower) && (trans == PastixNoTrans) ) ) ?
                      PastixSolveBackward : PastixSolveForward;
    void (*ztrsm)( pastix_data_t *, int, int, int, int,
                   sopalin_data_t *, pastix_rhs_t ) = ztrsm_table[ sched ];

    if (ztrsm == NULL) {
        ztrsm = static_ztrsm;
    }

#if defined(PASTIX_WITH_MPI)
    if( pastix_data->inter_node_procnbr > 1 ) {
        if( (sched == PastixSchedStarPU) || (sched == PastixSchedParsec) ) {
            ztrsm = runtime_ztrsm;
        }
        else {
            if ( sched == PastixSchedStatic ) {
                if ( solve_step == PastixSolveBackward ) {
                    ztrsm = sequential_ztrsm;
                }
                else {
                    ztrsm = static_ztrsm;
                }
            }
            else {
                ztrsm = sequential_ztrsm;
            }
        }
    }
#endif


#if defined(PASTIX_WITH_MPI)
    if ( ( pastix_data->inter_node_procnbr > 1 ) &&
         ( sched      == PastixSchedStatic )     &&
         ( solve_step == PastixSolveForward ) ) {
        solverRequestInit( sopalin_data->solvmtx );
        solverRhsRecvInit( sopalin_data->solvmtx, PastixComplex64 );
    }
#endif

    ztrsm( pastix_data, side, uplo, trans, diag, sopalin_data, rhsb );

#if defined(PASTIX_WITH_MPI)
    if ( ( pastix_data->inter_node_procnbr > 1 ) &&
         ( sched      == PastixSchedStatic )     &&
         ( solve_step == PastixSolveForward ) ) {
        cpucblk_zrequest_rhs_cleanup( PastixSolveForward, sched, sopalin_data->solvmtx, rhsb );
        solverRequestExit( sopalin_data->solvmtx );
        solverRhsRecvExit( sopalin_data->solvmtx );
    }
#endif

#if defined(PASTIX_WITH_MPI)
   MPI_Barrier( pastix_data->inter_node_comm );
#endif

    (void)solve_step;
}
