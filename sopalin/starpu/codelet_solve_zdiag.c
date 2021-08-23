/**
 *
 * @file codelet_solve_zdiag.c
 *
 * StarPU codelets for diag functions.
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
 * @addtogroup starpu_diag_solve
 * @{
 *
 **/
#include "common.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_starpu.h"
#include "pastix_zstarpu.h"
#include "codelets.h"

/**
 * Cblk version
 */
#if defined( PASTIX_STARPU_PROFILING )
measure_t solve_cblk_zdiag_perf[STARPU_NMAXWORKERS];
#endif

struct cl_solve_cblk_zdiag_args_s {
    profile_data_t  profile_data;
    SolverCblk     *cblk;
};

static struct starpu_perfmodel starpu_solve_cblk_zdiag_model =
{
    .type = STARPU_HISTORY_BASED,
    .symbol = "solve_cblk_zdiag",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_solve_cblk_zdiag_cpu(void *descr[], void *cl_arg)
{
    int                                nrhs;
    pastix_complex64_t                *b;
    int                                ldb;
    struct cl_solve_cblk_zdiag_args_s *args = (struct cl_solve_cblk_zdiag_args_s *) cl_arg;

    b    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(descr[1]);
    ldb  = (int)STARPU_MATRIX_GET_LD(descr[1]);
    nrhs = (int)STARPU_MATRIX_GET_NY(descr[1]);

    solve_cblk_zdiag( args->cblk, nrhs,
                      b, ldb, NULL );
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( solve_cblk_zdiag, 2 );

/**
 *******************************************************************************
 *
 * @brief Submit a task to perform a diagonal solve related to one cblk
 * to all the right hand side.
 *
 *******************************************************************************
 *
 * @param[inout] sopalin_data
 *          The structure providing descriptor about b and the A matrix.
 *
 * @param[in] cblk
 *          The cblk structure to which diagonal block belongs to.
 *
 * @param[in] prio
 *          The priority of the task in the DAG.
 *
 *******************************************************************************/
void
starpu_stask_cblk_zdiag( sopalin_data_t *sopalin_data,
                         SolverCblk     *cblk,
                         int             prio )
{
    struct cl_solve_cblk_zdiag_args_s *cl_arg;
    starpu_data_handle_t               handle;
    SolverMatrix                      *solvmtx = sopalin_data->solvmtx;
    pastix_int_t                       cblknum = cblk - solvmtx->cblktab;

    /*
     * Create the arguments array
     */
    cl_arg                         = malloc( sizeof(struct cl_solve_cblk_zdiag_args_s) );
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->profile_data.measures  = solve_cblk_zdiag_perf;
    cl_arg->profile_data.flops     = NAN;
#endif
    cl_arg->cblk                   = cblk;

    /* if ( cblk->cblktype & CBLK_TASKS_2D ) { */
    /*     handle = cblk->fblokptr->handler[0]; */
    /* } */
    /* else { */
        handle = cblk->handler[0];
    /* } */

    starpu_insert_task(
        pastix_codelet(&cl_solve_cblk_zdiag_cpu),
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_solve_cblk_zdiag_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
        STARPU_R,                       handle,
        STARPU_RW,                      solvmtx->starpu_desc_rhs->handletab[cblknum],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME,                    "solve_cblk_zdiag",
#endif
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketSolveDiag,
#endif
        0);
    (void) prio;
}
/**
 * @}
 */
