/**
 *
 * @file codelet_zgemmsp.c
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
 * Cblk version
 */
#if defined( PASTIX_STARPU_PROFILING )
starpu_profile_t cblk_zgemmsp_profile = { .next = NULL, .name = "cblk_zgemmsp" };

/**
 * @brief Profiling registration function
 */
void cblk_zgemmsp_profile_register( void ) __attribute__( ( constructor ) );
void
cblk_zgemmsp_profile_register( void )
{
    profiling_register_cl( &cblk_zgemmsp_profile );
}
#endif

struct cl_cblk_zgemmsp_args_s {
    profile_data_t    profile_data;
    sopalin_data_t   *sopalin_data;
    pastix_coefside_t sideA;
    pastix_trans_t    trans;
    const SolverCblk *cblk;
    const SolverBlok *blok;
    SolverCblk       *fcblk;
};

#if defined(PASTIX_STARPU_PROFILING_LOG)
void
cl_profiling_cb_cblk_zgemmsp( void *callback_arg )
{
    cl_profiling_callback( callback_arg );

    struct starpu_task                *task     = starpu_task_get_current();
    struct starpu_profiling_task_info *info     = task->profiling_info;

    if ( info == NULL ) {
        return;
    }
    struct cl_cblk_zgemmsp_args_s     *args     = (struct cl_cblk_zgemmsp_args_s *) callback_arg;
    double                             flops    = args->profile_data.flops;
    double                             duration = starpu_timing_timespec_delay_us( &info->start_time, &info->end_time );
    double                             speed    = flops / ( 1000.0 * duration );

    pastix_int_t                       k        = cblk_colnbr( args->cblk );
    pastix_int_t                       n        = blok_rownbr( args->blok );
    pastix_int_t                       m        = args->cblk->stride;

    m -= (args->cblk->cblktype & CBLK_LAYOUT_2D) ? args->blok->coefind / k : args->blok->coefind;
    m -= (args->sideA == PastixUCoef) ? blok_rownbr( args->blok ) : 0;

    cl_profiling_log_register(task->name, "cblk_zgemmsp", m, n, k, flops, speed);
}
#endif

static struct starpu_perfmodel starpu_cblk_zgemmsp_model =
{
#if defined(PASTIX_STARPU_COST_PER_ARCH)
    .type               = STARPU_PER_ARCH,
    .arch_cost_function = cblk_gemmsp_cost,
#else
    .type               = STARPU_HISTORY_BASED,
#endif
    .symbol             = "cblk_zgemmsp",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void
fct_cblk_zgemmsp_cpu( void *descr[], void *cl_arg )
{
    const void                    *A;
    const void                    *B;
    void                          *C;
    struct cl_cblk_zgemmsp_args_s *args = (struct cl_cblk_zgemmsp_args_s *)cl_arg;

    A = (const void *)STARPU_VECTOR_GET_PTR( descr[0] );
    B = (const void *)STARPU_VECTOR_GET_PTR( descr[1] );
    C = (void *)      STARPU_VECTOR_GET_PTR( descr[2] );

    /* Check layout due to NULL workspace for now */
    assert( args->cblk->cblktype  & CBLK_LAYOUT_2D );
    assert( args->fcblk->cblktype & CBLK_LAYOUT_2D );

    args->profile_data.flops = cpucblk_zgemmsp( args->sideA, args->trans,
                                                args->cblk, args->blok, args->fcblk,
                                                A, B, C, NULL, 0,
                                                &( args->sopalin_data->solvmtx->lowrank ) );
}

#if defined(PASTIX_WITH_CUDA)
static void
fct_cblk_zgemmsp_gpu( void *descr[], void *cl_arg )
{
    const void                    *A;
    const void                    *B;
    void                          *C;
    struct cl_cblk_zgemmsp_args_s *args  = (struct cl_cblk_zgemmsp_args_s *)cl_arg;

    A = (const void *)STARPU_VECTOR_GET_PTR( descr[0] );
    B = (const void *)STARPU_VECTOR_GET_PTR( descr[1] );
    C = (void *)      STARPU_VECTOR_GET_PTR( descr[2] );

    args->profile_data.flops = gpucblk_zgemmsp( args->sideA, args->trans,
                                                args->cblk, args->blok, args->fcblk,
                                                A, B, C,
                                                &( args->sopalin_data->solvmtx->lowrank ),
                                                starpu_cuda_get_local_stream() );
}
#endif /* defined(PASTIX_WITH_CUDA) */
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_ANY( cblk_zgemmsp, 3, STARPU_CUDA_ASYNC );

void
starpu_task_cblk_zgemmsp( sopalin_data_t   *sopalin_data,
                          pastix_coefside_t sideA,
                          pastix_coefside_t sideB,
                          pastix_trans_t    trans,
                          const SolverCblk *cblk,
                          const SolverBlok *blok,
                          SolverCblk       *fcblk,
                          int               prio )
{
    struct cl_cblk_zgemmsp_args_s *cl_arg = NULL;
    long long                      execute_where;
    int                            need_exec = 1;
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
        if ( fcblk->ownerid == sopalin_data->solvmtx->clustnum ) {
            need_submit = 1;
        }
        else {
            need_exec = 0;
        }
        if ( starpu_mpi_cached_receive( fcblk->handler[sideA] ) ) {
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
    if ( need_exec ) {
        cl_arg                        = malloc( sizeof( struct cl_cblk_zgemmsp_args_s ) );
        cl_arg->sopalin_data          = sopalin_data;
#if defined( PASTIX_STARPU_PROFILING )
        cl_arg->profile_data.measures = cblk_zgemmsp_profile.measures;
        cl_arg->profile_data.flops    = NAN;
#endif
        cl_arg->sideA                 = sideA;
        cl_arg->trans                 = trans;
        cl_arg->cblk                  = cblk;
        cl_arg->blok                  = blok;
        cl_arg->fcblk                 = fcblk;

        execute_where = cl_cblk_zgemmsp_any.where;
#if defined(PASTIX_WITH_CUDA)
        if ( (cblk->cblktype  & CBLK_COMPRESSED) ||
             (fcblk->cblktype & CBLK_COMPRESSED) )
        {
            /* Disable CUDA */
            execute_where &= (~STARPU_CUDA);
        }
#endif
    }

#if defined(PASTIX_DEBUG_STARPU) || defined(PASTIX_STARPU_PROFILING_LOG)
    asprintf( &task_name, "%s( %ld )",
              cl_cblk_zgemmsp_any.name,
              (long)(cblk - sopalin_data->solvmtx->cblktab) );
#endif

    starpu_insert_task(
        pastix_codelet(&cl_cblk_zgemmsp_any),
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_cblk_zgemmsp_args_s ),
        STARPU_EXECUTE_WHERE,           execute_where,
#if defined(PASTIX_STARPU_PROFILING)
#if defined(PASTIX_STARPU_PROFILING_LOG)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_cb_cblk_zgemmsp, cl_arg,
#else
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
#endif
        STARPU_R,                       cblk->handler[sideA],
        STARPU_R,                       cblk->handler[sideB],
        STARPU_RW,                      fcblk->handler[sideA],
#if defined(PASTIX_STARPU_DEBUG) || defined(PASTIX_STARPU_PROFILING_LOG)
        STARPU_NAME,                    task_name,
#endif
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketGEMM1D,
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
starpu_profile_t blok_zgemmsp_profile = {
    .next = NULL,
    .name = "blok_zgemmsp"
};

/**
 * @brief Profiling registration function
 */
void blok_zgemmsp_profile_register( void ) __attribute__( ( constructor ) );
void
blok_zgemmsp_profile_register( void )
{
    profiling_register_cl( &blok_zgemmsp_profile );
}
#endif

struct cl_blok_zgemmsp_args_s {
    profile_data_t    profile_data;
    sopalin_data_t   *sopalin_data;
    pastix_trans_t    trans;
    const SolverCblk *cblk;
    SolverCblk       *fcblk;
    pastix_int_t      blok_mk;
    pastix_int_t      blok_nk;
    pastix_int_t      blok_mn;
};

#if defined(PASTIX_STARPU_PROFILING_LOG)
void
cl_profiling_cb_blok_zgemmsp( void *callback_arg )
{
    cl_profiling_callback( callback_arg );

    struct starpu_task                *task = starpu_task_get_current();
    struct starpu_profiling_task_info *info = task->profiling_info;

    if ( info == NULL ) {
        return;
    }
    struct cl_blok_zgemmsp_args_s *args     = (struct cl_blok_zgemmsp_args_s *) callback_arg;
    double          flops                   = args->profile_data.flops;
    double          duration                = starpu_timing_timespec_delay_us( &info->start_time, &info->end_time );
    double          speed                   = flops / ( 1000.0 * duration );


    pastix_int_t      K      = cblk_colnbr( args->cblk );
    const SolverBlok *blokA  = args->cblk[0].fblokptr + args->blok_mk;
    const SolverBlok *blokB  = args->cblk[0].fblokptr + args->blok_nk;
    pastix_int_t      cblk_m = blokA->fcblknm;
    pastix_int_t      cblk_n = blokB->fcblknm;
    const SolverBlok *lblokK = args->cblk[1].fblokptr;
    pastix_int_t      full_m = 0;
    pastix_int_t      full_n = 0;
    for (; (blokA < lblokK) && (blokA->fcblknm == cblk_m); blokA++) {
        full_m += blok_rownbr(blokA);
    }
    for (; (blokB < lblokK) && (blokB->fcblknm == cblk_n); blokB++) {
        full_n += blok_rownbr( blokB );
    }
    cl_profiling_log_register(task->name, "blok_zgemmsp", full_m, full_n, K, flops, speed);
}
#endif

static struct starpu_perfmodel starpu_blok_zgemmsp_model =
{
#if defined(PASTIX_STARPU_COST_PER_ARCH)
    .type               = STARPU_PER_ARCH,
    .arch_cost_function = blok_gemmsp_cost,
#else
    .type               = STARPU_HISTORY_BASED,
#endif
    .symbol             = "blok_zgemmsp",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void
fct_blok_zgemmsp_cpu( void *descr[], void *cl_arg )
{
    struct cl_blok_zgemmsp_args_s *args = (struct cl_blok_zgemmsp_args_s *)cl_arg;
    const void *                   A;
    const void *                   B;
    void *                         C;

    A = (const void *)STARPU_VECTOR_GET_PTR( descr[0] );
    B = (const void *)STARPU_VECTOR_GET_PTR( descr[1] );
    C = (void *)      STARPU_VECTOR_GET_PTR( descr[2] );

    assert( args->cblk->cblktype  & CBLK_TASKS_2D );
    assert( args->fcblk->cblktype & CBLK_TASKS_2D );

    args->profile_data.flops = cpublok_zgemmsp( args->trans,
                                                args->cblk, args->fcblk,
                                                args->blok_mk, args->blok_nk, args->blok_mn,
                                                A, B, C,
                                                &(args->sopalin_data->solvmtx->lowrank) );
}

#if defined(PASTIX_WITH_CUDA)
static void
fct_blok_zgemmsp_gpu( void *descr[], void *cl_arg )
{
    const void                    *A;
    const void                    *B;
    void                          *C;
    struct cl_blok_zgemmsp_args_s *args = (struct cl_blok_zgemmsp_args_s *)cl_arg;

    A = (const void *)STARPU_VECTOR_GET_PTR( descr[0] );
    B = (const void *)STARPU_VECTOR_GET_PTR( descr[1] );
    C = (void *)      STARPU_VECTOR_GET_PTR( descr[2] );

    assert( args->cblk->cblktype  & CBLK_TASKS_2D );
    assert( args->fcblk->cblktype & CBLK_TASKS_2D );

    args->profile_data.flops = gpublok_zgemmsp( args->trans,
                                                args->cblk, args->fcblk,
                                                args->blok_mk, args->blok_nk, args->blok_mn,
                                                A, B, C,
                                                &(args->sopalin_data->solvmtx->lowrank),
                                                starpu_cuda_get_local_stream() );
}
#endif /* defined(PASTIX_WITH_CUDA) */
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_ANY( blok_zgemmsp, 3, STARPU_CUDA_ASYNC );

void
starpu_task_blok_zgemmsp( sopalin_data_t   *sopalin_data,
                          pastix_coefside_t sideA,
                          pastix_coefside_t sideB,
                          pastix_trans_t    trans,
                          const SolverCblk *cblk,
                          SolverCblk       *fcblk,
                          const SolverBlok *blokA,
                          const SolverBlok *blokB,
                          int               prio )
{
    pastix_int_t frownum;
    pastix_int_t lrownum;
    pastix_int_t blok_mn = 0, j = 0;
    pastix_int_t blok_mk = blokA - cblk->fblokptr;
    pastix_int_t blok_nk = blokB - cblk->fblokptr;
    SolverBlok  *blokC   = fcblk->fblokptr;

    struct cl_blok_zgemmsp_args_s *cl_arg = NULL;
    long long                      execute_where;
    int                            need_exec = 1;
#if defined(PASTIX_DEBUG_STARPU) || defined(PASTIX_STARPU_PROFILING_LOG)
    char                          *task_name;
#endif

    assert( blok_nk <= blok_mk );

    do {
        frownum = blokC->frownum;
        lrownum = blokC->lrownum;
        blok_mn += j;
        j = 1;

        /* Increase lrownum as long as blocks are facing the same cblk */
        while( (blokC < fcblk[1].fblokptr-1) &&
               (blokC[0].fcblknm == blokC[1].fcblknm) &&
               (blokC[0].lcblknm == blokC[1].lcblknm) )
        {
            blokC++; j++;
            lrownum = blokC->lrownum;
        }
        blokC++;
    }
    while( !((blokA->frownum >= frownum) &&
             (blokA->lrownum <= lrownum)) );

    blokC = fcblk->fblokptr + blok_mn;

    assert( blokA->lcblknm == blokB->lcblknm );
    assert( blokB->fcblknm == blokC->lcblknm );
    assert( blokC->frownum <= blokA->frownum );
    assert( blokA[-1].fcblknm != blokA[0].fcblknm );
    assert( blokB[-1].fcblknm != blokB[0].fcblknm );
    assert( (blok_mn == 0) || (blokC[-1].fcblknm != blokC[0].fcblknm) );

    /*
     * Check if it needs to be submitted
     */
#if defined(PASTIX_WITH_MPI)
    {
        int need_submit = 0;
        if ( cblk->ownerid == sopalin_data->solvmtx->clustnum ) {
            need_submit = 1;
        }
        if ( fcblk->ownerid == sopalin_data->solvmtx->clustnum ) {
            need_submit = 1;
        }
        else {
            need_exec = 0;
        }
        if ( starpu_mpi_cached_receive( blokC->handler[sideA] ) ) {
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
    if ( need_exec ) {
        cl_arg                        = malloc( sizeof(struct cl_blok_zgemmsp_args_s) );
        cl_arg->sopalin_data          = sopalin_data;
#if defined(PASTIX_STARPU_PROFILING)
        cl_arg->profile_data.measures = blok_zgemmsp_profile.measures;
        cl_arg->profile_data.flops    = NAN;
#endif
        cl_arg->trans                 = trans;
        cl_arg->cblk                  = cblk;
        cl_arg->fcblk                 = fcblk;
        cl_arg->blok_mk               = blok_mk;
        cl_arg->blok_nk               = blok_nk;
        cl_arg->blok_mn               = blok_mn;

        execute_where = cl_blok_zgemmsp_any.where;
#if defined(PASTIX_WITH_CUDA)
        if ( (cblk->cblktype  & CBLK_COMPRESSED) ||
             (fcblk->cblktype & CBLK_COMPRESSED) )
        {
            /* Disable CUDA */
            execute_where &= (~STARPU_CUDA);
        }
#endif
    }

#if defined(PASTIX_DEBUG_STARPU) || defined(PASTIX_STARPU_PROFILING_LOG)
    asprintf( &task_name, "%s( %ld, %ld, %ld )",
              cl_blok_zgemmsp_any.name,
              (long)(blokA - sopalin_data->solvmtx->bloktab),
              (long)(blokB - sopalin_data->solvmtx->bloktab),
              (long)(blokC - sopalin_data->solvmtx->bloktab) );
#endif

    starpu_insert_task(
        pastix_codelet(&cl_blok_zgemmsp_any),
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_blok_zgemmsp_args_s ),
        STARPU_EXECUTE_WHERE,           execute_where,
#if defined(PASTIX_STARPU_PROFILING)
#if defined(PASTIX_STARPU_PROFILING_LOG)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_cb_blok_zgemmsp, cl_arg,
#else
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
#endif
        STARPU_R,                       blokA->handler[sideA],
        STARPU_R,                       blokB->handler[sideB],
        STARPU_RW,                      blokC->handler[sideA],
#if defined(PASTIX_DEBUG_STARPU)
        STARPU_NAME,                    task_name,
#endif
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketGEMM2D,
#else
        STARPU_PRIORITY,                prio,
#endif
        0);
    (void)prio;
}

/**
 * @}
 */
