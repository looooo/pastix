/**
 *
 * @file codelet_zsytrfsp.c
 *
 * StarPU codelets for LDL^t functions
 *
 * @copyright 2016-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
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
#include "common.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_starpu.h"
#include "pastix_zstarpu.h"
#include "codelets.h"
#include "pastix_starpu_model.h"

/**
 * Cblk version
 */
#if defined( PASTIX_STARPU_PROFILING )
measure_t cblk_zsytrfsp_perf[STARPU_NMAXWORKERS];
#endif

struct cl_cblk_zsytrfsp_args_s {
    profile_data_t  profile_data;
    SolverMatrix   *solvmtx;
    SolverCblk     *cblk;
};

static struct starpu_perfmodel starpu_cblk_zsytrfsp1d_panel_model =
{
#if defined(PASTIX_STARPU_COST_PER_ARCH)
    .type = STARPU_PER_ARCH,
    .arch_cost_function = cblk_sytrf_cost,
#else
    .type = STARPU_HISTORY_BASED,
#endif
    .symbol = "cblk_zsytrf",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_cblk_zsytrfsp1d_panel_cpu(void *descr[], void *cl_arg)
{
    SolverMatrix                 *solvmtx;
    pastix_complex64_t           *L;
    pastix_complex64_t           *DL;
    int                           nbpivot;
    struct cl_cblk_zsytrfsp_args_s *args = (struct cl_cblk_zsytrfsp_args_s *)cl_arg;

    L  = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[0]);
    DL = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[1]);

    solvmtx = args->solvmtx;
    nbpivot = cpucblk_zsytrfsp1d_panel( solvmtx, args->cblk, L, DL );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( cblk_zsytrfsp1d_panel, 2 );

void
starpu_task_cblk_zsytrfsp1d_panel( sopalin_data_t *sopalin_data,
                                   SolverCblk     *cblk,
                                   int             prio )
{
    struct cl_cblk_zsytrfsp_args_s *cl_arg;
    starpu_data_handle_t           *handler = (starpu_data_handle_t *)( cblk->handler );
    pastix_int_t                    N       = cblk_colnbr( cblk );
    pastix_int_t                    M       = cblk->stride;

    if ( M-N > 0 ) {
        starpu_vector_data_register( handler + 1, -1, (uintptr_t)NULL, M * N,
                                     sopalin_data->solvmtx->starpu_desc->typesze );
    }
    else {
        starpu_vector_data_register( handler + 1, -1, (uintptr_t)NULL, 0,
                                     sopalin_data->solvmtx->starpu_desc->typesze );
    }

#if defined(PASTIX_WITH_MPI)
    {
        int64_t tag_desc = sopalin_data->solvmtx->starpu_desc->mpitag;
        int64_t tag_cblk = 2 * cblk->gcblknum + 1;

        starpu_mpi_data_register( *(handler+1),
                                  tag_desc | tag_cblk,
                                  cblk->ownerid );
    }
#endif /* PASTIX_WITH_MPI */

    cl_arg                        = malloc( sizeof(struct cl_cblk_zsytrfsp_args_s) );
    cl_arg->solvmtx               = sopalin_data->solvmtx;
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->profile_data.measures = cblk_zsytrfsp_perf;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->cblk                  = cblk;

    starpu_insert_task(
        pastix_codelet(&cl_cblk_zsytrfsp1d_panel_cpu),
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_cblk_zsytrfsp_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
        STARPU_RW,                      cblk->handler[0],
        STARPU_W,                       cblk->handler[1],
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketFacto1D,
#else
        STARPU_PRIORITY,                prio,
#endif
        0);
    (void) prio;
}

/**
 * Blok version
 */
#if defined( PASTIX_STARPU_PROFILING )
measure_t blok_zsytrfsp_perf[STARPU_NMAXWORKERS];
#endif

struct cl_blok_zsytrfsp_args_s {
    profile_data_t  profile_data;
    SolverMatrix   *solvmtx;
    SolverCblk     *cblk;
};

static struct starpu_perfmodel starpu_blok_zsytrfsp_model =
{
#if defined(PASTIX_STARPU_COST_PER_ARCH)
    .type = STARPU_PER_ARCH,
    .arch_cost_function = blok_sytrf_cost,
#else
    .type = STARPU_HISTORY_BASED,
#endif
    .symbol = "blok_zsytrfsp",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_blok_zsytrfsp_cpu(void *descr[], void *cl_arg)
{
    SolverMatrix                   *solvmtx;
    pastix_complex64_t             *L;
    int                             nbpivot;
    struct cl_blok_zsytrfsp_args_s *args = (struct cl_blok_zsytrfsp_args_s *)cl_arg;

    L = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[0]);

    assert(args->cblk->cblktype & CBLK_TASKS_2D);

    solvmtx = args->solvmtx;
    nbpivot = cpucblk_zsytrfsp1d_sytrf( solvmtx, args->cblk, L );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( blok_zsytrfsp, 1 );

void
starpu_task_blok_zsytrf( sopalin_data_t *sopalin_data,
                         SolverCblk     *cblk,
                         int             prio )
{
    struct cl_blok_zsytrfsp_args_s *cl_arg;

    /*
     * Create the arguments array
     */
    cl_arg                        = malloc( sizeof(struct cl_blok_zsytrfsp_args_s) );
    cl_arg->solvmtx               = sopalin_data->solvmtx;
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->profile_data.measures = blok_zsytrfsp_perf;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->cblk                  = cblk;

    starpu_insert_task(
        pastix_codelet(&cl_blok_zsytrfsp_cpu),
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_blok_zsytrfsp_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
        STARPU_RW,                      cblk->fblokptr->handler[0],
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketFacto2D,
#else
        STARPU_PRIORITY,                prio,
#endif
        0);
    (void) prio;
}

/**
 * @}
 */
