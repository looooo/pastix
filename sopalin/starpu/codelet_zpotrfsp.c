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
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#include "common.h"
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
    (void)nb;

    starpu_insert_task(
        pastix_codelet(&cl_cblk_zpotrfsp1d_panel),
        STARPU_VALUE, &cblk,             sizeof(SolverCblk*),
        STARPU_VALUE, &sopalin_data,     sizeof(sopalin_data*),
        STARPU_RW,     cblk->handler[0],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "cblk_zpotrfsp1d_panel",
#endif
        0);
}

/**
 * @}
 */
