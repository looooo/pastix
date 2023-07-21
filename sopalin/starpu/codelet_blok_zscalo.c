/**
 *
 * @file codelet_blok_zscalo.c
 *
 * StarPU codelets for blas-like functions
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Ian Masliah
 * @date 2023-06-07
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
starpu_profile_t blok_zscalo_profile = {
    .next = NULL,
    .name = "blok_zscalo"
};

/**
 * @brief Profiling registration function
 */
void blok_zscalo_profile_register( void ) __attribute__( ( constructor ) );
void
blok_zscalo_profile_register( void )
{
    profiling_register_cl( &blok_zscalo_profile );
}
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct cl_blok_zscalo_args_s {
    profile_data_t  profile_data;
    pastix_trans_t  trans;
    SolverCblk     *cblk;
    pastix_int_t    blok_m;
};

static struct starpu_perfmodel starpu_blok_zscalo_model = {
    .type   = STARPU_HISTORY_BASED,
    .symbol = "blok_zscalo",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void
fct_blok_zscalo_cpu( void *descr[], void *cl_arg )
{
    struct cl_blok_zscalo_args_s *args = (struct cl_blok_zscalo_args_s *)cl_arg;
    const void                   *A;
    const void                   *D;
    void                         *B;

    A = pastix_starpu_blok_get_ptr( descr[0] );
    D = pastix_starpu_blok_get_ptr( descr[1] );
    B = pastix_starpu_blok_get_ptr( descr[2] );

    assert( args->cblk->cblktype & CBLK_TASKS_2D );

    cpublok_zscalo( args->trans, args->cblk, args->blok_m, A, D, B );
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( blok_zscalo, 3 );
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] sopalin_data
 *          TODO
 *
 * @param[in] trans
 *          TODO
 *
 * @param[in] cblk
 *          TODO
 *
 * @param[in] blok
 *          TODO
 *
 * @param[in] prio
 *          TODO
 *
 *******************************************************************************/
void
starpu_task_blok_zscalo( sopalin_data_t *sopalin_data,
                         pastix_trans_t  trans,
                         SolverCblk     *cblk,
                         SolverBlok     *blok,
                         int             prio )
{
    struct cl_blok_zscalo_args_s *cl_arg;
    starpu_data_handle_t         *handler = (starpu_data_handle_t *)( blok->handler );
    SolverBlok                   *blokA   = blok;
    pastix_int_t                  blok_m  = blok - cblk->fblokptr;
    pastix_int_t                  M       = blok_rownbr_ext( blokA );
#if defined(PASTIX_DEBUG_STARPU)
    char                         *task_name;
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
        if ( !need_submit ) {
            return;
        }
    }
#endif

    starpu_vector_data_register( handler + 1, -1, (uintptr_t)NULL, M * cblk_colnbr( cblk ),
                                 sopalin_data->solvmtx->starpu_desc->typesze );

#if defined(PASTIX_WITH_MPI)
    {
        int64_t tag_desc = sopalin_data->solvmtx->starpu_desc->mpitag;
        int64_t bloknum  = blok - sopalin_data->solvmtx->bloktab;
        int64_t tag_blok = 2 * (sopalin_data->solvmtx->cblknbr + bloknum) + 1;

        starpu_mpi_data_register( *(handler + 1),
                                  tag_desc | tag_blok,
                                  cblk->ownerid );
    }
#endif /* PASTIX_WITH_MPI */

    /*
     * Create the arguments array
     */
    cl_arg                        = malloc( sizeof(struct cl_blok_zscalo_args_s) );
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->profile_data.measures = blok_zscalo_profile.measures;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->trans                 = trans;
    cl_arg->cblk                  = cblk;
    cl_arg->blok_m                = blok_m;

#if defined(PASTIX_DEBUG_STARPU)
    asprintf( &task_name, "%s( %ld, %ld )",
              cl_blok_zscalo_cpu.name,
              (long)(cblk - sopalin_data->solvmtx->cblktab),
              (long)(blok - sopalin_data->solvmtx->bloktab) );
#endif

    pastix_starpu_insert_task(
        &cl_blok_zscalo_cpu,
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_blok_zscalo_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
        STARPU_R,                       blok->handler[0],
        STARPU_R,                       cblk->fblokptr->handler[0],
        STARPU_W,                       blok->handler[1],
#if defined(PASTIX_DEBUG_STARPU)
        STARPU_NAME,                    task_name,
#endif
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketScalo,
#else
        STARPU_PRIORITY,                prio,
#endif
        0);
    (void)prio;
}
/**
 * @}
 */
