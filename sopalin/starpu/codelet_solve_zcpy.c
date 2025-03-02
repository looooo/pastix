/**
 *
 * @file codelet_solve_zcpy.c
 *
 * StarPU codelet for add function
 *
 * @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Vincent Bridonneau
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2024-07-05
 *
 * @precisions normal z -> z c d s
 *
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define _GNU_SOURCE
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#include "common.h"
#include <lapacke.h>
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
starpu_profile_t solve_blok_zcpy_profile = {
    .next = NULL,
    .name = "solve_blok_zcpy"
};

/**
 * @brief Profiling registration function
 */
void solve_blok_zcpy_profile_register( void ) __attribute__( ( constructor ) );
void
solve_blok_zcpy_profile_register( void )
{
    profiling_register_cl( &solve_blok_zcpy_profile );
}
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct cl_solve_blok_zcpy_args_s {
    profile_data_t    profile_data;
    SolverCblk       *cblk;
    const SolverCblk *fcblk;
};

static struct starpu_perfmodel starpu_solve_blok_zcpy_model = {
    .type   = STARPU_HISTORY_BASED,
    .symbol = "solve_blok_zcpy",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void
fct_solve_blok_zcpy_cpu( void *descr[], void *cl_arg )
{
    const pastix_complex64_t         *A;
    pastix_complex64_t               *B;
    pastix_int_t                      m, n, lda, ldb;
    struct cl_solve_blok_zcpy_args_s *args = (struct cl_solve_blok_zcpy_args_s *)cl_arg;
    int                               rc;

#if defined(PASTIX_DEBUG_STARPU)
    fprintf( stderr, "[pastix][%s] Add copy cblk into recv = %d ok\n", __func__, args->cblk->gcblknum );
#endif

    A    = (pastix_complex64_t *)STARPU_MATRIX_GET_PTR( descr[0] );
    lda  = (pastix_int_t)        STARPU_MATRIX_GET_LD( descr[0] );

    m    = (pastix_int_t)        STARPU_MATRIX_GET_NX( descr[1] );
    n    = (pastix_int_t)        STARPU_MATRIX_GET_NY( descr[1] );
    ldb  = (pastix_int_t)        STARPU_MATRIX_GET_LD( descr[1] );
    B    = (pastix_complex64_t *)STARPU_MATRIX_GET_PTR( descr[1] );

    assert( m <= ((pastix_int_t)STARPU_MATRIX_GET_NX( descr[0] )) );

    A += (args->cblk->lcolidx - args->fcblk->lcolidx);

    rc = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, B, ldb );
    assert( rc == 0 );

    (void)rc;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( solve_blok_zcpy, 2 );
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief Insert the task to add a fanin cblk on the receiver side (The fanin is
 * seen on this side as the RECV cblk).  Note that the caller always execute the
 * task.
 *
 *******************************************************************************
 *
 * @param[in] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 * @param[in] cblk
 *          The column block of the matrix.
 *
 * @param[in] fcblk
 *          The facing column block of the matrix.
 *
 * @param[in] prio
 *          The task priority.
 *
 *******************************************************************************/
void
starpu_stask_blok_zcpy_bwd_recv( sopalin_data_t   *sopalin_data,
                                 pastix_rhs_t      rhsb,
                                 SolverCblk       *cblk,
                                 const SolverCblk *fcblk,
                                 int               prio )
{
    struct cl_solve_blok_zcpy_args_s *cl_arg   = NULL;
    SolverMatrix                     *solvmtx  = sopalin_data->solvmtx;
    pastix_int_t                      cblknum  = cblk - solvmtx->cblktab;
    pastix_int_t                      fcblknum = fcblk - solvmtx->cblktab;
#if defined(PASTIX_DEBUG_STARPU)
    char                       *task_name;
#endif

#if defined(PASTIX_DEBUG_STARPU)
    fprintf( stderr, "[%2d][%s] cblk = %d, ownerid = %d, handler = %p, size = %ld\n",
             solvmtx->clustnum, __func__, cblk->gcblknum, cblk->ownerid, rhb->starpu_desc->handletab[cblknum],
             cblk_colnbr( cblk ) * sizeof(pastix_complex64_t) );
#endif

    /*
     * Create the arguments array
     */
    cl_arg                        = malloc( sizeof( struct cl_solve_blok_zcpy_args_s) );
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->profile_data.measures = solve_blok_zcpy_profile.measures;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->cblk                  = cblk;
    cl_arg->fcblk                 = fcblk;

#if defined(PASTIX_DEBUG_STARPU)
    /* This actually generates a memory leak */
    asprintf( &task_name, "%s( %ld )",
              cl_solve_blok_zcpy_cpu.name,
              (long)(cblk - sopalin_data->solvmtx->cblktab) );
#endif

    assert( cblk->cblktype & CBLK_RECV );
    assert( !(fcblk->cblktype & (CBLK_RECV|CBLK_FANIN)) );

    pastix_starpu_insert_task(
        &cl_solve_blok_zcpy_cpu,
        STARPU_CL_ARGS,                 cl_arg,                 sizeof( struct cl_solve_blok_zcpy_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, solve_blok_zcpy_callback, cl_arg,
#endif
        STARPU_R,                       rhsb->starpu_desc->handletab[fcblknum],
        STARPU_W,                       rhsb->starpu_desc->handletab[cblknum],
#if defined(PASTIX_DEBUG_STARPU)
        STARPU_NAME,                    task_name,
#endif
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketSolveGEMM,
#else
        STARPU_PRIORITY,                prio,
#endif
        0);

    (void)prio;
}

/**
 * @}
 */
