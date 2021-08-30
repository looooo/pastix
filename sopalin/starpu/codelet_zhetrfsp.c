/**
 *
 * @file codelet_zhetrfsp.c
 *
 * StarPU codelets for LDL^h functions
 *
 * @copyright 2016-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Ian Masliah
 * @date 2021-06-21
 *
 * @precisions normal z -> z c
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
#include "pastix_starpu_model.h"

/**
 * Cblk version
 */
#if defined( PASTIX_STARPU_PROFILING )
starpu_profile_t cblk_zhetrfsp_profile = {
    .next = NULL,
    .name = "cblk_zhetrfsp"
};

/**
 * @brief Profiling registration function
 */
void cblk_zhetrfsp_profile_register( void ) __attribute__( ( constructor ) );
void
cblk_zhetrfsp_profile_register( void )
{
    profiling_register_cl( &cblk_zhetrfsp_profile );
}
#endif

struct cl_cblk_zhetrfsp_args_s {
    profile_data_t  profile_data;
    SolverMatrix   *solvmtx;
    SolverCblk     *cblk;
};

static struct starpu_perfmodel starpu_cblk_zhetrfsp1d_panel_model = {
#if defined(PASTIX_STARPU_COST_PER_ARCH)
    .type               = STARPU_PER_ARCH,
    .arch_cost_function = cblk_hetrf_cost,
#else
    .type               = STARPU_HISTORY_BASED,
#endif
    .symbol             = "cblk_zhetrf",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void
fct_cblk_zhetrfsp1d_panel_cpu( void *descr[], void *cl_arg )
{
    SolverMatrix                   *solvmtx;
    void                           *L;
    void                           *DL;
    int                             nbpivot;
    struct cl_cblk_zhetrfsp_args_s *args = (struct cl_cblk_zhetrfsp_args_s *)cl_arg;

    L  = (void *)STARPU_VECTOR_GET_PTR( descr[0] );
    DL = (void *)STARPU_VECTOR_GET_PTR( descr[1] );

    solvmtx = args->solvmtx;
    nbpivot = cpucblk_zhetrfsp1d_panel( solvmtx, args->cblk, L, DL );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( cblk_zhetrfsp1d_panel, 2 );

void
starpu_task_cblk_zhetrfsp1d_panel( sopalin_data_t *sopalin_data,
                                   SolverCblk     *cblk,
                                   int             prio )
{
    struct cl_cblk_zhetrfsp_args_s *cl_arg;
    starpu_data_handle_t           *handler = (starpu_data_handle_t *)( cblk->handler );
    pastix_int_t                    N       = cblk_colnbr( cblk );
    pastix_int_t                    M       = cblk->stride;
#if defined(PASTIX_DEBUG_STARPU)
    char                           *task_name;
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
        if ( starpu_mpi_cached_receive( cblk->handler[0] ) ) {
            need_submit = 1;
        }
        if ( !need_submit ) {
            return;
        }
    }
#endif

    if ( (M - N) > 0 ) {
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

        starpu_mpi_data_register( *(handler + 1),
                                  tag_desc | tag_cblk,
                                  cblk->ownerid );
    }
#endif /* PASTIX_WITH_MPI */

    /*
     * Create the arguments array
     */
    cl_arg                        = malloc( sizeof(struct cl_cblk_zhetrfsp_args_s) );
    cl_arg->solvmtx               = sopalin_data->solvmtx;
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->profile_data.measures = cblk_zhetrfsp_profile.measures;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->cblk                  = cblk;

#if defined(PASTIX_DEBUG_STARPU)
    asprintf( &task_name, "%s( %ld )",
              cl_cblk_zhetrfsp1d_panel_cpu.name,
              (long)( cblk - sopalin_data->solvmtx->cblktab ) );
#endif

    starpu_insert_task(
        pastix_codelet(&cl_cblk_zhetrfsp1d_panel_cpu),
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_cblk_zhetrfsp_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
        STARPU_RW,                      cblk->handler[0],
        STARPU_W,                       cblk->handler[1],
#if defined(PASTIX_DEBUG_STARPU)
        STARPU_NAME,                    task_name,
#endif
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketFacto1D,
#else
        STARPU_PRIORITY,                prio,
#endif
        0);
    (void)prio;
}

/**
 * Blok version
 */
#if defined( PASTIX_STARPU_PROFILING )
starpu_profile_t blok_zhetrfsp_profile = {
    .next = NULL,
    .name = "blok_zhetrfsp"
};

/**
 * @brief Profiling registration function
 */
void blok_zhetrfsp_profile_register( void ) __attribute__( ( constructor ) );
void
blok_zhetrfsp_profile_register( void )
{
    profiling_register_cl( &blok_zhetrfsp_profile );
}
#endif

struct cl_blok_zhetrfsp_args_s {
    profile_data_t  profile_data;
    SolverMatrix   *solvmtx;
    SolverCblk     *cblk;
};

static struct starpu_perfmodel starpu_blok_zhetrfsp_model = {
#if defined(PASTIX_STARPU_COST_PER_ARCH)
    .type               = STARPU_PER_ARCH,
    .arch_cost_function = blok_hetrf_cost,
#else
    .type               = STARPU_HISTORY_BASED,
#endif
    .symbol             = "blok_zhetrfsp",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void
fct_blok_zhetrfsp_cpu( void *descr[], void *cl_arg )
{
    SolverMatrix                   *solvmtx;
    void                           *L;
    int                             nbpivot;
    struct cl_blok_zhetrfsp_args_s *args = (struct cl_blok_zhetrfsp_args_s *)cl_arg;

    L = (void *)STARPU_VECTOR_GET_PTR( descr[0] );

    assert( args->cblk->cblktype & CBLK_TASKS_2D );

    solvmtx = args->solvmtx;
    nbpivot = cpucblk_zhetrfsp1d_hetrf( solvmtx, args->cblk, L );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( blok_zhetrfsp, 1 );

void
starpu_task_blok_zhetrf( sopalin_data_t *sopalin_data,
                         SolverCblk     *cblk,
                         int             prio )
{
    struct cl_blok_zhetrfsp_args_s *cl_arg;
#if defined(PASTIX_DEBUG_STARPU)
    char                           *task_name;
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
        if ( starpu_mpi_cached_receive( cblk->fblokptr->handler[0] ) ) {
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
    cl_arg                        = malloc( sizeof(struct cl_blok_zhetrfsp_args_s) );
    cl_arg->solvmtx               = sopalin_data->solvmtx;
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->profile_data.measures = blok_zhetrfsp_profile.measures;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->cblk                  = cblk;

#if defined(PASTIX_DEBUG_STARPU)
    asprintf( &task_name, "%s( %ld )",
              cl_blok_zhetrfsp_cpu.name,
              (long)(cblk - sopalin_data->solvmtx->cblktab) );
#endif

    starpu_insert_task(
        pastix_codelet(&cl_blok_zhetrfsp_cpu),
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_blok_zhetrfsp_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
        STARPU_RW,                      cblk->fblokptr->handler[0],
#if defined(PASTIX_DEBUG_STARPU)
        STARPU_NAME,                    task_name,
#endif
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketFacto2D,
#else
        STARPU_PRIORITY,                prio,
#endif
        0);
    (void)prio;
}

/**
 * @}
 */
