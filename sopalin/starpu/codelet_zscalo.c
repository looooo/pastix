/**
 *
 * @file codelet_zscalo.c
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
measure_t blok_zscalo_perf[STARPU_NMAXWORKERS];
#endif

struct cl_blok_zscalo_args_s {
    profile_data_t  profile_data;
    pastix_trans_t  trans;
    SolverCblk     *cblk;
    pastix_int_t    blok_m;
};

static struct starpu_perfmodel starpu_blok_zscalo_model =
{
    .type = STARPU_HISTORY_BASED,
    .symbol = "blok_zscalo",
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void fct_blok_zscalo_cpu(void *descr[], void *cl_arg)
{
    const void                   *A;
    const void                   *D;
    void                         *B;
    struct cl_blok_zscalo_args_s *args = (struct cl_blok_zscalo_args_s *) cl_arg;

    A = (const void *)STARPU_VECTOR_GET_PTR(descr[0]);
    D = (const void *)STARPU_VECTOR_GET_PTR(descr[1]);
    B = (void *)STARPU_VECTOR_GET_PTR(descr[2]);

    assert( args->cblk->cblktype & CBLK_TASKS_2D );

    cpublok_zscalo( args->trans, args->cblk, args->blok_m, A, D, B );
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( blok_zscalo, 3 );

void
starpu_task_blok_zscalo( sopalin_data_t   *sopalin_data,
                         pastix_trans_t    trans,
                         SolverCblk       *cblk,
                         SolverBlok       *blok,
                         int               prio )
{
    struct cl_blok_zscalo_args_s *cl_arg;
    starpu_data_handle_t         *handler = (starpu_data_handle_t*)(blok->handler);
    SolverBlok                   *blokA = blok;
    pastix_int_t                  blok_m = blok - cblk->fblokptr;
    pastix_int_t                  M = blok_rownbr_ext( blokA );

    starpu_vector_data_register( handler + 1, -1, (uintptr_t)NULL, M * cblk_colnbr( cblk ),
                                 sopalin_data->solvmtx->starpu_desc->typesze );

#if defined(PASTIX_WITH_MPI)
    {
        int64_t tag_desc = sopalin_data->solvmtx->starpu_desc->mpitag;
        int64_t bloknum  = blok - sopalin_data->solvmtx->bloktab;
        int64_t tag_blok = 2 * (sopalin_data->solvmtx->cblknbr + bloknum) + 1;

        starpu_mpi_data_register( *(handler+1),
                                  tag_desc | tag_blok,
                                  cblk->ownerid );
    }
#endif /* PASTIX_WITH_MPI */

    /*
     * Create the arguments array
     */
    cl_arg                        = malloc( sizeof(struct cl_blok_zscalo_args_s) );
#if defined(PASTIX_STARPU_PROFILING)
    cl_arg->profile_data.measures = blok_zscalo_perf;
    cl_arg->profile_data.flops    = NAN;
#endif
    cl_arg->trans                 = trans;
    cl_arg->cblk                  = cblk;
    cl_arg->blok_m                = blok_m;

    starpu_insert_task(
        pastix_codelet(&cl_blok_zscalo_cpu),
        STARPU_CL_ARGS,                 cl_arg,                sizeof( struct cl_blok_zscalo_args_s ),
#if defined(PASTIX_STARPU_PROFILING)
        STARPU_CALLBACK_WITH_ARG_NFREE, cl_profiling_callback, cl_arg,
#endif
        STARPU_R,                       blok->handler[0],
        STARPU_R,                       cblk->fblokptr->handler[0],
        STARPU_W,                       blok->handler[1],
#if defined(PASTIX_STARPU_HETEROPRIO)
        STARPU_PRIORITY,                BucketScalo,
#else
        STARPU_PRIORITY,                prio,
#endif
        0);
    (void) prio;
}
/**
 * @}
 */
