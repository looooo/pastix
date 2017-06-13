/**
 *
 * @file codelet_zpotrfsp.c
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
static void cl_cblk_zpotrfsp1d_panel_cpu(void *descr[], void *cl_arg)
{
    sopalin_data_t *sopalin_data;
    SolverCblk *cblk;
    pastix_complex64_t *L;
    int nbpivot;

    L = (pastix_complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &cblk, &sopalin_data);

    assert( !(cblk->cblktype & CBLK_TASKS_2D) );

    nbpivot = cpucblk_zpotrfsp1d_panel( cblk, L, sopalin_data->diagthreshold,
                                        &(sopalin_data->solvmtx->lowrank) );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( cblk_zpotrfsp1d_panel, 1 )

void
starpu_task_cblk_zpotrfsp1d_panel( sopalin_data_t *sopalin_data,
                                   SolverCblk     *cblk )
{
    starpu_insert_task(
        pastix_codelet(&cl_cblk_zpotrfsp1d_panel),
        STARPU_VALUE, &cblk,             sizeof(SolverCblk*),
        STARPU_VALUE, &sopalin_data,     sizeof(sopalin_data_t*),
        STARPU_RW,     cblk->handler[0],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "cblk_zpotrfsp1d_panel",
#endif
        0);
}

#if !defined(PASTIX_STARPU_SIMULATION)
static void cl_blok_zpotrfsp_cpu(void *descr[], void *cl_arg)
{
    sopalin_data_t *sopalin_data;
    SolverCblk *cblk;
    pastix_complex64_t *L;
    int nbpivot;

    L = (pastix_complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &cblk, &sopalin_data);

    assert(cblk->cblktype & CBLK_TASKS_2D);

    nbpivot = cpucblk_zpotrfsp1d_potrf( cblk, L, sopalin_data->diagthreshold );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( blok_zpotrfsp, 1 )

void
starpu_task_blok_zpotrf( sopalin_data_t *sopalin_data,
                         SolverCblk     *cblk )
{
    starpu_insert_task(
        pastix_codelet(&cl_blok_zpotrfsp),
        STARPU_VALUE, &cblk,             sizeof(SolverCblk*),
        STARPU_VALUE, &sopalin_data,     sizeof(sopalin_data_t*),
        STARPU_RW,     cblk->fblokptr->handler[0],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "blok_zpotrfsp",
#endif
        0);
}

/**
 * @}
 */
