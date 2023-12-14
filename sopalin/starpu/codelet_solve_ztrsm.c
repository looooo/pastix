/**
 *
 * @file codelet_solve_ztrsm.c
 *
 * StarPU codelet for TRSM function
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Vincent Bridonneau
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2023-07-21
 *
 * @precisions normal z -> z c d s
 *
 * @addtogroup starpu_trsm_solve
 * @{
 *
 **/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define _GNU_SOURCE
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#include "common.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_starpu.h"
#include "pastix_zstarpu.h"
#include "codelets.h"

#if defined( PASTIX_STARPU_PROFILING )
/**
 * @brief Block version
 */
starpu_profile_t solve_blok_ztrsm_profile = {
    .next = NULL,
    .name = "solve_blok_ztrsm"
};

/**
 * @brief Profiling registration function
 */
void solve_blok_ztrsm_profile_register( void ) __attribute__( ( constructor ) );
void
solve_blok_ztrsm_profile_register( void )
{
    profiling_register_cl( &solve_blok_ztrsm_profile );
}
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct cl_solve_blok_ztrsm_args_s {
    profile_data_t    profile_data;
    pastix_side_t     side;
    pastix_uplo_t     uplo;
    pastix_trans_t    trans;
    pastix_diag_t     diag;
    const SolverCblk *cblk;
};

static struct starpu_perfmodel starpu_solve_blok_ztrsm_model = {
    .type   = STARPU_HISTORY_BASED,
    .symbol = "solve_blok_ztrsm",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void
fct_solve_blok_ztrsm_cpu( void *descr[], void *cl_arg )
{
    struct cl_solve_blok_ztrsm_args_s *args = (struct cl_solve_blok_ztrsm_args_s *) cl_arg;
    void                              *A;
    pastix_complex64_t                *B;
    pastix_int_t                       nrhs, ldb;

    A    = pastix_starpu_cblk_get_ptr( descr[0] );
    B    = (pastix_complex64_t *)STARPU_MATRIX_GET_PTR( descr[1] );
    ldb  = (pastix_int_t)        STARPU_MATRIX_GET_LD( descr[1] );
    nrhs = (pastix_int_t)        STARPU_MATRIX_GET_NY( descr[1] );

    solve_blok_ztrsm( args->side, args->uplo,
                      args->trans, args->diag, args->cblk,
                      nrhs, A, B, ldb );
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( solve_blok_ztrsm, 2 );
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief Submit a task to do a trsm related to a diagonal
 * block of the matrix A.
 *
 *******************************************************************************
 *
 * @param[in] sopalin_data
 *          The data that provide the SolverMatrix structure from PaStiX, and
 *          descriptor of b (providing nrhs, b and ldb).
 *
 * @param[inout] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 * @param[in] coef
 *          Specify whether the computation are made with the L part, or the U
 *          part of A. It has to be either PastixLCoef, or PastixUCoef.
 *
 * @param[in] side
 *          Specify the side parameter of the TRSM.
 *
 * @param[in] uplo
 *          Specify the uplo parameter of the TRSM.
 *
 * @param[in] trans
 *          Specify the transposition used for the matrix A in the
 *          computation. It has to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] diag
 *          Specify if the off-diagonal blocks are unit triangular. It has to be
 *          either PastixUnit or PastixNonUnit.
 *
 * @param[in] cblk
 *          The cblk structure that corresponds to the A and B matrix.

 * @param[in] prio
 *          The priority of the task in th DAG.
 *
 *******************************************************************************/
void
starpu_stask_blok_ztrsm( sopalin_data_t   *sopalin_data,
                         pastix_rhs_t      rhsb,
                         pastix_coefside_t coef,
                         pastix_side_t     side,
                         pastix_uplo_t     uplo,
                         pastix_trans_t    trans,
                         pastix_diag_t     diag,
                         const SolverCblk *cblk,
                         pastix_int_t      prio )
{
    struct cl_solve_blok_ztrsm_args_s *cl_arg;
    starpu_data_handle_t               handle  = cblk->handler[coef];
    SolverMatrix                      *solvmtx = sopalin_data->solvmtx;
    pastix_int_t                       cblknum = cblk - solvmtx->cblktab;
#if defined(PASTIX_DEBUG_STARPU)
    char                              *task_name;
#endif

    /*
     * Check if it needs to be submitted
     */
#if defined(PASTIX_WITH_MPI)
    {
        int need_submit = 0;
        if ( cblk->ownerid == sopalin_data->solvmtx->clustnum ) {
            need_submit = 1;
        }
        if ( starpu_mpi_cached_receive( rhsb->starpu_desc->handletab[cblknum] ) ) {
            need_submit = 1;
        }
        if ( !need_submit ) {
            return;
        }
    }
#endif

    /*
     * Create the arguments array
     */
    cl_arg                        = malloc( sizeof(struct cl_solve_blok_ztrsm_args_s) );
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->profile_data.measures = solve_blok_ztrsm_profile.measures;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->side                  = side;
    cl_arg->trans                 = trans;
    cl_arg->uplo                  = uplo;
    cl_arg->diag                  = diag;
    cl_arg->cblk                  = cblk;

#if defined(PASTIX_DEBUG_STARPU)
    asprintf( &task_name, "%s( %ld )",
              cl_solve_blok_ztrsm_cpu.name,
              (long)(cblknum) );
#endif

    pastix_starpu_insert_task(
        &cl_solve_blok_ztrsm_cpu,
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_solve_blok_ztrsm_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
        STARPU_R,                       handle,
        STARPU_RW,                      rhsb->starpu_desc->handletab[cblknum],
#if defined(PASTIX_DEBUG_STARPU)
        STARPU_NAME,                    task_name,
#endif
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketSolveTRSM,
#endif
        0);
    (void)prio;
}

/**
 * @}
 */
