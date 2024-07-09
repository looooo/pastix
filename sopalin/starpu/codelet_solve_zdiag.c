/**
 *
 * @file codelet_solve_zdiag.c
 *
 * StarPU codelets for diag functions.
 *
 * @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Vincent Bridonneau
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Alycia Lisito
 * @author Nolan Bredel
 * @author Tom Moenne-Loccoz
 * @date 2023-11-07
 *
 * @precisions normal z -> z c d s
 *
 * @addtogroup starpu_diag_solve
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
 * @brief Cblk version
 */
starpu_profile_t solve_cblk_zdiag_profile = {
    .next = NULL,
    .name = "solve_cblk_zdiag"
};

/**
 * @brief Profiling registration function
 */
void solve_cblk_zdiag_profile_register( void ) __attribute__( ( constructor ) );
void
solve_cblk_zdiag_profile_register( void )
{
    profiling_register_cl( &solve_cblk_zdiag_profile );
}
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct cl_solve_cblk_zdiag_args_s {
    profile_data_t  profile_data;
    SolverCblk     *cblk;
};

static struct starpu_perfmodel starpu_solve_cblk_zdiag_model = {
    .type   = STARPU_HISTORY_BASED,
    .symbol = "solve_cblk_zdiag",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void
fct_solve_cblk_zdiag_cpu( void *descr[], void *cl_arg )
{
    struct cl_solve_cblk_zdiag_args_s *args = (struct cl_solve_cblk_zdiag_args_s *) cl_arg;
    void                              *A;
    int                                nrhs;
    pastix_complex64_t                *b;
    int                                ldb;

    A    = pastix_starpu_cblk_get_ptr( descr[0] );
    b    = (pastix_complex64_t *)STARPU_MATRIX_GET_PTR( descr[1] );
    ldb  = (int)STARPU_MATRIX_GET_LD( descr[1] );
    nrhs = (int)STARPU_MATRIX_GET_NY( descr[1] );

    solve_cblk_zdiag( args->cblk, A, nrhs, b, ldb, NULL );
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( solve_cblk_zdiag, 2 );
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

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
 * @param[inout] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
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
                         pastix_rhs_t    rhsb,
                         SolverCblk     *cblk,
                         int             prio )
{
    struct cl_solve_cblk_zdiag_args_s *cl_arg;
    starpu_data_handle_t               handle;
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
    cl_arg                        = malloc( sizeof(struct cl_solve_cblk_zdiag_args_s) );
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->profile_data.measures = solve_cblk_zdiag_profile.measures;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->cblk                  = cblk;

#if defined(PASTIX_DEBUG_STARPU)
    asprintf( &task_name, "%s( %ld )",
              cl_solve_cblk_zdiag_cpu.name,
              (long)cblknum );
#endif

    handle = cblk->handler[0];

    pastix_starpu_insert_task(
        &cl_solve_cblk_zdiag_cpu,
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_solve_cblk_zdiag_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
        STARPU_R,                       handle,
        STARPU_RW,                      rhsb->starpu_desc->handletab[cblknum],
#if defined(PASTIX_DEBUG_STARPU)
        STARPU_NAME,                    task_name,
#endif
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketSolveDiag,
#endif
        0);
    (void)prio;
}
/**
 * @}
 */
