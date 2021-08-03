/**
 *
 * @file codelet_zgetrfsp.c
 *
 * StarPU codelets for LU functions
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
#include "pastix_starpu_model.h"

/**
 * Cblk version
 */
#if defined( PASTIX_STARPU_PROFILING )
starpu_profile_t cblk_zgetrfsp_profile = {
    .next = NULL,
    .name = "cblk_zgetrfsp"
};

/**
 * @brief Profiling registration function
 */
void cblk_zgetrfsp_profile_register( void ) __attribute__( ( constructor ) );
void
cblk_zgetrfsp_profile_register( void )
{
    profiling_register_cl( &cblk_zgetrfsp_profile );
}
#endif

struct cl_cblk_zgetrfsp_args_s {
    profile_data_t  profile_data;
    SolverMatrix   *solvmtx;
    SolverCblk     *cblk;
};

static struct starpu_perfmodel starpu_cblk_zgetrfsp1d_panel_model = {
#if defined( PASTIX_STARPU_COST_PER_ARCH )
    .type               = STARPU_PER_ARCH,
    .arch_cost_function = cblk_getrf_cost,
#else
    .type               = STARPU_HISTORY_BASED,
#endif
    .symbol             = "cblk_zgetrfsp1d_panel",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void
fct_cblk_zgetrfsp1d_panel_cpu( void *descr[], void *cl_arg )
{
    struct cl_cblk_zgetrfsp_args_s *args = (struct cl_cblk_zgetrfsp_args_s *)cl_arg;
    void                           *L;
    void                           *U;
    int                             nbpivot;

    L = pastix_starpu_cblk_get_ptr( descr[0] );
    U = pastix_starpu_cblk_get_ptr( descr[1] );

    nbpivot = cpucblk_zgetrfsp1d_panel( args->solvmtx, args->cblk, L, U );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( cblk_zgetrfsp1d_panel, 2 );

void
starpu_task_cblk_zgetrfsp1d_panel( sopalin_data_t *sopalin_data,
                                   SolverCblk     *cblk,
                                   int             prio )
{
    struct cl_cblk_zgetrfsp_args_s *cl_arg;
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
        if ( starpu_mpi_cached_receive( cblk->handler[1] ) ) {
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
    cl_arg                        = malloc( sizeof( struct cl_cblk_zgetrfsp_args_s) );
    cl_arg->solvmtx               = sopalin_data->solvmtx;
#if defined( PASTIX_STARPU_PROFILING )
    cl_arg->profile_data.measures = cblk_zgetrfsp_profile.measures;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->cblk                  = cblk;

#if defined(PASTIX_DEBUG_STARPU)
    asprintf( &task_name, "%s( %ld )",
              cl_cblk_zgetrfsp1d_panel_cpu.name,
              (long)(cblk - sopalin_data->solvmtx->cblktab) );
#endif

    starpu_insert_task(
        pastix_codelet(&cl_cblk_zgetrfsp1d_panel_cpu),
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_cblk_zgetrfsp_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
        STARPU_RW,                      cblk->handler[0],
        STARPU_RW,                      cblk->handler[1],
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
starpu_profile_t blok_zgetrfsp_profile = {
    .next = NULL,
    .name = "blok_zgetrfsp"
};

/**
 * @brief Profiling registration function
 */
void blok_zgetrfsp_profile_register( void ) __attribute__( ( constructor ) );
void
blok_zgetrfsp_profile_register( void )
{
    profiling_register_cl( &blok_zgetrfsp_profile );
}
#endif

struct cl_blok_zgetrfsp_args_s {
    profile_data_t  profile_data;
    SolverMatrix   *solvmtx;
    SolverCblk     *cblk;
};

static struct starpu_perfmodel starpu_blok_zgetrfsp_model = {
#if defined(PASTIX_STARPU_COST_PER_ARCH)
    .type               = STARPU_PER_ARCH,
    .arch_cost_function = blok_getrf_cost,
#else
    .type               = STARPU_HISTORY_BASED,
#endif
    .symbol             = "blok_zgetrfsp",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void
fct_blok_zgetrfsp_cpu( void *descr[], void *cl_arg )
{
    struct cl_blok_zgetrfsp_args_s *args = (struct cl_blok_zgetrfsp_args_s *)cl_arg;
    void                           *L, *U;
    int                             nbpivot;

    L = pastix_starpu_blok_get_ptr( descr[0] );
    U = pastix_starpu_blok_get_ptr( descr[1] );

    assert( args->cblk->cblktype & CBLK_TASKS_2D );

    nbpivot = cpucblk_zgetrfsp1d_getrf( args->solvmtx, args->cblk, L, U );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( blok_zgetrfsp, 2 );

void
starpu_task_blok_zgetrf( sopalin_data_t *sopalin_data,
                         SolverCblk     *cblk,
                         int             prio )
{
    struct cl_blok_zgetrfsp_args_s *cl_arg;
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
        if ( starpu_mpi_cached_receive( cblk->fblokptr->handler[1] ) ) {
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
    cl_arg                        = malloc( sizeof( struct cl_blok_zgetrfsp_args_s) );
    cl_arg->solvmtx               = sopalin_data->solvmtx;
#if defined( PASTIX_STARPU_PROFILING )
    cl_arg->profile_data.measures = blok_zgetrfsp_profile.measures;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->cblk                  = cblk;

#if defined(PASTIX_DEBUG_STARPU)
    asprintf( &task_name, "%s( %ld )",
              cl_blok_zgetrfsp_cpu.name,
              (long)(cblk - sopalin_data->solvmtx->cblktab) );
#endif

    starpu_insert_task(
        pastix_codelet(&cl_blok_zgetrfsp_cpu),
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_blok_zgetrfsp_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
        STARPU_RW,                      cblk->fblokptr->handler[0],
        STARPU_RW,                      cblk->fblokptr->handler[1],
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
