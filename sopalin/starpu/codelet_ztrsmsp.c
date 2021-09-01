/**
 *
 * @file codelet_ztrsmsp.c
 *
 * StarPU codelets for blas-like functions
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
#if defined(PASTIX_WITH_CUDA)
#include "pastix_zcuda.h"
#endif
#include "pastix_starpu.h"
#include "pastix_zstarpu.h"
#include "codelets.h"
#include "pastix_starpu_model.h"

/**
 * Block version
 */
#if defined( PASTIX_STARPU_PROFILING )
measure_t blok_ztrsmsp_perf[STARPU_NMAXWORKERS];

starpu_profile_t blok_ztrsmsp_profile = {
    .next = NULL,
    .name = "cblk_ztrsmsp"
};

/**
 * @brief Profiling registration function
 */
void blok_ztrsmsp_profile_register( void ) __attribute__( ( constructor ) );
void
blok_ztrsmsp_profile_register( void )
{
    profiling_register_cl( &blok_ztrsmsp_profile );
}
#endif

struct cl_blok_ztrsmsp_args_s {
    profile_data_t    profile_data;
    sopalin_data_t   *sopalin_data;
    pastix_side_t     side;
    pastix_uplo_t     uplo;
    pastix_trans_t    trans;
    pastix_diag_t     diag;
    const SolverCblk *cblk;
    pastix_int_t      blok_m;
};

#if defined(PASTIX_STARPU_PROFILING_LOG)
void
cl_profiling_cb_blok_ztrsmsp( void *callback_arg )
{
    cl_profiling_callback( callback_arg );

    struct starpu_task                *task     = starpu_task_get_current();
    struct starpu_profiling_task_info *info     = task->profiling_info;

    if ( info == NULL ) {
        return;
    }
    struct cl_blok_ztrsmsp_args_s     *args     = (struct cl_blok_ztrsmsp_args_s *) callback_arg;
    double                             flops    = args->profile_data.flops;
    double                             duration = starpu_timing_timespec_delay_us( &info->start_time, &info->end_time );
    double                             speed    = flops / ( 1000.0 * duration );
    pastix_int_t                       N        = cblk_colnbr( args->cblk );
    pastix_int_t                       full_m   = 0;
    const SolverBlok                  *lblok    = args->cblk[1].fblokptr;
    const SolverBlok                  *blok     = args->cblk[0].fblokptr + args->blok_m;
    pastix_int_t                       cblk_m   = blok->fcblknm;
    for (; (blok < lblok) && (blok->fcblknm == cblk_m); blok++) {
        full_m += blok_rownbr(blok);
    }
    cl_profiling_log_register(task->name, "blok_ztrsmsp", full_m, N, 0, flops, speed);
}
#endif

static struct starpu_perfmodel starpu_blok_ztrsmsp_model =
{
#if defined(PASTIX_STARPU_COST_PER_ARCH)
    .type               = STARPU_PER_ARCH,
    .arch_cost_function = cblk_gemmsp_cost,
#else
    .type               = STARPU_HISTORY_BASED,
#endif
    .symbol             = "blok_ztrsmsp",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void
fct_blok_ztrsmsp_cpu( void *descr[], void *cl_arg )
{
    struct cl_blok_ztrsmsp_args_s *args = (struct cl_blok_ztrsmsp_args_s *)cl_arg;
    const void                    *A;
    void                          *C;

    A = pastix_starpu_blok_get_ptr( descr[0] );
    C = pastix_starpu_blok_get_ptr( descr[1] );

    assert( args->cblk->cblktype & CBLK_TASKS_2D );

    args->profile_data.flops = cpublok_ztrsmsp( args->side, args->uplo,
                                                args->trans, args->diag,
                                                args->cblk, args->blok_m, A, C,
                                                &(args->sopalin_data->solvmtx->lowrank) );
}

#if defined(PASTIX_WITH_CUDA)
static void
fct_blok_ztrsmsp_gpu( void *descr[], void *cl_arg )
{
    struct cl_blok_ztrsmsp_args_s *args = (struct cl_blok_ztrsmsp_args_s *)cl_arg;
    const void                    *A;
    void                          *C;

    A = pastix_starpu_blok_get_ptr( descr[0] );
    C = pastix_starpu_blok_get_ptr( descr[1] );

    assert( args->cblk->cblktype & CBLK_TASKS_2D );

    args->profile_data.flops = gpublok_ztrsmsp( args->side, args->uplo,
                                                args->trans, args->diag,
                                                args->cblk, args->blok_m, A, C,
                                                &(args->sopalin_data->solvmtx->lowrank),
                                                starpu_cuda_get_local_stream() );
}
#endif /* defined(PASTIX_WITH_CUDA) */
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_ANY( blok_ztrsmsp, 2, STARPU_CUDA_ASYNC );

void
starpu_task_blok_ztrsmsp( sopalin_data_t   *sopalin_data,
                          pastix_coefside_t coef,
                          pastix_side_t     side,
                          pastix_uplo_t     uplo,
                          pastix_trans_t    trans,
                          pastix_diag_t     diag,
                          const SolverCblk *cblk,
                          SolverBlok       *blok,
                          int               prio )
{
    struct cl_blok_ztrsmsp_args_s *cl_arg;
    long long                      execute_where = cl_blok_ztrsmsp_any.where;
    pastix_int_t                   blok_m        = blok - cblk->fblokptr;
#if defined(PASTIX_DEBUG_STARPU) || defined(PASTIX_STARPU_PROFILING_LOG)
    char                          *task_name;
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
        if ( starpu_mpi_cached_receive( blok->handler[coef] ) ) {
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
    cl_arg                        = malloc( sizeof(struct cl_blok_ztrsmsp_args_s) );
    cl_arg->sopalin_data          = sopalin_data;
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->profile_data.measures = blok_ztrsmsp_profile.measures;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->side                  = side;
    cl_arg->uplo                  = uplo;
    cl_arg->trans                 = trans;
    cl_arg->diag                  = diag;
    cl_arg->cblk                  = cblk;
    cl_arg->blok_m                = blok_m;

#if defined(PASTIX_WITH_CUDA)
    if ( (cblk->cblktype & CBLK_COMPRESSED) ) {
        execute_where &= (~STARPU_CUDA);
    }
#endif

#if defined(PASTIX_DEBUG_STARPU) || defined(PASTIX_STARPU_PROFILING_LOG)
    asprintf( &task_name, "%s( %ld, %ld, %ld )",
              cl_blok_ztrsmsp_any.name,
              (long)(cblk - sopalin_data->solvmtx->cblktab),
              (long)(blok - sopalin_data->solvmtx->bloktab),
              (long)coef );
#endif

    starpu_insert_task(
        pastix_codelet(&cl_blok_ztrsmsp_any),
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_blok_ztrsmsp_args_s ),
        STARPU_EXECUTE_WHERE,           execute_where,
#if defined(PASTIX_STARPU_PROFILING)
#if defined(PASTIX_STARPU_PROFILING_LOG)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_cb_blok_ztrsmsp, cl_arg,
#else
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
#endif
        STARPU_R,                       cblk->fblokptr->handler[coef],
        STARPU_RW,                      blok->handler[coef],
#if defined(PASTIX_DEBUG_STARPU) || defined(PASTIX_STARPU_PROFILING_LOG)
        STARPU_NAME,                    task_name,
#endif
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketTRSM2D,
#else
        STARPU_PRIORITY,                prio,
#endif
        0);
    (void)prio;
}

/**
 * @}
 */
