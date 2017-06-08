/**
 *
 * @file codelet_zgemmsp.c
 *
 * StarPU codelets for Cholesky functions
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2013-06-24
 *
 * @precisions normal z -> z c d s
 *
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"
#include "sopalin/starpu/codelets.h"

#if !defined(PASTIX_STARPU_SIMULATION)
static void cl_cblk_zgemmsp_cpu(void *descr[], void *cl_arg)
{
    pastix_coefside_t sideA;
    pastix_coefside_t sideB;
    pastix_trans_t    trans;
    SolverCblk       *cblk;
    SolverBlok       *blok;
    SolverCblk       *fcblk;
    sopalin_data_t   *sopalin_data;
    const pastix_complex64_t *A;
    const pastix_complex64_t *B;
    pastix_complex64_t *C;

    A = (const pastix_complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (const pastix_complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (pastix_complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);

    /* Check layout due to workspace */
    starpu_codelet_unpack_args(cl_arg, &sideA, &sideB, &trans, &cblk, &blok, &fcblk, &sopalin_data);

    assert( cblk->cblktype & CBLK_LAYOUT_2D );
    assert( fcblk->cblktype & CBLK_LAYOUT_2D );

    cpucblk_zgemmsp( sideA, sideB, trans,
                     cblk, blok, fcblk,
                     A, B, C, NULL,
                     &(sopalin_data->solvmtx->lowrank) );
}

#if defined(PASTIX_WITH_CUDA)
static void cl_cblk_zgemmsp_gpu(void *descr[], void *cl_arg)
{
    pastix_coefside_t sideA;
    pastix_coefside_t sideB;
    pastix_trans_t    trans;
    SolverCblk       *cblk;
    SolverBlok       *blok;
    SolverCblk       *fcblk;
    sopalin_data_t   *sopalin_data;
    const cuDoubleComplex *A;
    const cuDoubleComplex *B;
    cuDoubleComplex *C;

    A = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    B = (const cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]);
    C = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[2]);

    /* Check layout due to workspace */
    assert( cblk->cblktype & CBLK_LAYOUT_2D );
    assert( fcblk->cblktype & CBLK_LAYOUT_2D );

    starpu_codelet_unpack_args(cl_arg, &sideA, &sideB, &trans, &cblk, &blok, &fcblk, &sopalin_data);

    gpucblk_zgemmsp( sideA, sideB, trans,
                     cblk, blok, fcblk,
                     A, B, C,
                     &(sopalin_data->solvmtx->lowrank),
                     starpu_cuda_get_local_stream() );
}
#endif /* defined(PASTIX_WITH_CUDA) */
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_GPU( cblk_zgemmsp, 3, STARPU_CUDA_ASYNC )

void
starpu_task_cblk_zgemmsp( pastix_coefside_t sideA,
                          pastix_coefside_t sideB,
                          pastix_trans_t    trans,
                          const SolverCblk *cblk,
                          const SolverBlok *blok,
                          SolverCblk       *fcblk,
                          sopalin_data_t   *sopalin_data )
{
    starpu_insert_task(
        pastix_codelet(&cl_cblk_zgemmsp),
        STARPU_VALUE, &sideA,             sizeof(pastix_coefside_t),
        STARPU_VALUE, &sideB,             sizeof(pastix_coefside_t),
        STARPU_VALUE, &trans,             sizeof(pastix_trans_t),
        STARPU_VALUE, &cblk,              sizeof(SolverCblk*),
        STARPU_VALUE, &blok,              sizeof(SolverBlok*),
        STARPU_VALUE, &fcblk,             sizeof(SolverCblk*),
        STARPU_R,      cblk->handler[sideA],
        STARPU_R,      cblk->handler[sideB],
        STARPU_RW,     fcblk->handler[sideA],
        STARPU_VALUE, &sopalin_data,     sizeof(sopalin_data_t*),
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "cblk_zgemmsp",
#endif
        0);
}

/**
 * @}
 */
