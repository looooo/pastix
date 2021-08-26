/**
 *
 * @file codelet_solve_zgemm.c
 *
 * StarPU codelet for gemm function
 *
 * @copyright 2016-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Vincent Bridonneau
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2021-06-21
 *
 * @precisions normal z -> z c d s
 *
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#define _GNU_SOURCE
#include "common.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_starpu.h"
#include "pastix_zstarpu.h"
#include "codelets.h"

/**
 * Block version
 */
#if defined( PASTIX_STARPU_PROFILING )
starpu_profile_t solve_blok_zgemm_profile = {
    .next = NULL,
    .name = "solve_blok_zgemm"
};

/**
 * @brief Profiling registration function
 */
void solve_blok_zgemm_profile_register( void ) __attribute__( ( constructor ) );
void
solve_blok_zgemm_profile_register( void )
{
    profiling_register_cl( &solve_blok_zgemm_profile );
}
#endif

struct cl_solve_blok_zgemm_args_s {
    profile_data_t    profile_data;
    pastix_side_t     side;
    pastix_trans_t    trans;
    const SolverCblk *cblk;
    const SolverBlok *blok;
    SolverCblk       *fcbk;
};

static struct starpu_perfmodel starpu_solve_blok_zgemm_model =
{
    .type = STARPU_HISTORY_BASED,
    .symbol = "solve_blok_zgemm",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_solve_blok_zgemm_cpu(void *descr[], void *cl_arg)
{
    const void                        *dataA = NULL;
    const pastix_lrblock_t            *lrA;
    const pastix_complex64_t          *A;
    pastix_complex64_t                *B, *C;
    pastix_int_t                       nrhs, ldb, ldc;
    struct cl_solve_blok_zgemm_args_s *args = (struct cl_solve_blok_zgemm_args_s *) cl_arg;

    dataA = (const void *)STARPU_VECTOR_GET_PTR(descr[0]);
    B     = (pastix_complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    ldb   = (pastix_int_t)        STARPU_MATRIX_GET_LD (descr[1]);
    nrhs  = (pastix_int_t)        STARPU_MATRIX_GET_NY (descr[1]);
    C     = (pastix_complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    ldc   = (pastix_int_t)        STARPU_MATRIX_GET_LD (descr[2]);

    /*
     * Make sure we get the correct pointer to the lrA, or to the right position in [lu]coeftab
     */
    if ( (args->side == PastixLeft) && (args->cblk->cblktype & CBLK_COMPRESSED) ) {
        lrA = dataA;
        lrA += (args->blok - args->cblk->fblokptr);
        dataA = lrA;
    }
    else if ( (args->side == PastixRight) && (args->fcbk->cblktype & CBLK_COMPRESSED) ) {
        lrA = dataA;
        lrA += (args->blok - args->fcbk->fblokptr);
        dataA = lrA;
    }
    else {
        A = dataA;
        A += args->blok->coefind;
        dataA = A;
    }

    solve_blok_zgemm( args->side, args->trans, nrhs,
                      args->cblk, args->blok, args->fcbk, dataA, B, ldb, C, ldc );
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( solve_blok_zgemm, 3 );

/**
 *******************************************************************************
 *
 * @brief Submit a task to perform a gemm.
 *
 *******************************************************************************
 *
 * @param[in] coef
 *          Specify whether the computation are made with the L part, or the U
 *          part of A. It has to be either PastixLCoef, or PastixUCoef.
 *
 * @param[in] side
 *          Specify the side parameter of the TRSM.
 *
 * @param[in] trans
 *          Specify the transposition used for the matrix A in the
 *          computation. It has to be either PastixTrans or PastixConjTrans.
 *
 * @param[in] cblk
 *          The cblk structure that corresponds to the A and B matrix.
 *
 * @param[in] blok
 *          The blok structure that corresponds to the A matrix, and that
 *          belongs either to cblk or fcbk depending on the side parameter.
 *
 * @param[inout] fcbk
 *          The cblk structure that corresponds to the C matrix.

 *
 * @param[in] sopalin_data
 *          The data that provide the SolverMatrix structure from PaStiX, and
 *          descriptor of b (providing nrhs, b and ldb).

 * @param[in] prio
 *          The priority of the task in th DAG.
 *
 *******************************************************************************/
void
starpu_stask_blok_zgemm( sopalin_data_t   *sopalin_data,
                         pastix_coefside_t coef,
                         pastix_side_t     side,
                         pastix_trans_t    trans,
                         const SolverCblk *cblk,
                         const SolverBlok *blok,
                         SolverCblk       *fcbk,
                         pastix_int_t      prio )
{
    struct cl_solve_blok_zgemm_args_s *cl_arg;
    SolverMatrix                      *solvmtx = sopalin_data->solvmtx;
    pastix_int_t                       cblknum = cblk - solvmtx->cblktab;
    pastix_int_t                       fcbknum = fcbk - solvmtx->cblktab;
    starpu_data_handle_t               handle;
#if defined(PASTIX_DEBUG_STARPU)
    char                             *task_name;
    asprintf( &task_name, "%s( %ld, %ld, %ld )",
              cl_solve_blok_zgemm_cpu.name,
              (long) (( side == PastixRight ) ? fcbknum : cblknum), (long) cblknum, (long) fcbknum );

#endif
    
    /*
     * Create the arguments array
     */
    cl_arg                        = malloc( sizeof(struct cl_solve_blok_zgemm_args_s) );
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->profile_data.measures = solve_blok_zgemm_profile.measures;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->side                  = side;
    cl_arg->trans                 = trans;
    cl_arg->cblk                  = cblk;
    cl_arg->blok                  = blok;
    cl_arg->fcbk                  = fcbk;

    if ( side == PastixRight ) {
        handle = fcbk->handler[coef];
    }
    else {
        handle = cblk->handler[coef];
    }

    starpu_insert_task(
        pastix_codelet(&cl_solve_blok_zgemm_cpu),
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_solve_blok_zgemm_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
        STARPU_R,                       handle,
        STARPU_R,                       solvmtx->starpu_desc_rhs->handletab[cblknum],
        STARPU_RW,                      solvmtx->starpu_desc_rhs->handletab[fcbknum],
#if defined(PASTIX_DEBUG_STARPU)
        STARPU_NAME,                    task_name,
#endif
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketSolveGEMM,
#endif
        0);
    (void)prio;
}

/**
 * @}
 */
