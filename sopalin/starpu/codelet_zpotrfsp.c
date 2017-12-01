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
 * @date 2017-06-24
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
#include "pastix_starpu.h"
#include "codelets.h"
#include "pastix_starpu_model.h"

/**
 * Cblk version
 */
static struct starpu_perfmodel starpu_cblk_zpotrfsp1d_panel_model =
{
    .type = STARPU_PER_ARCH,
    .symbol = "cblk_zpotrfsp1d_panel",
    .arch_cost_function = cblk_potrf_cost,
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void cl_cblk_zpotrfsp1d_panel_cpu(void *descr[], void *cl_arg)
{
    sopalin_data_t *sopalin_data;
    SolverCblk *cblk;
    pastix_complex64_t *L;
    int nbpivot;

    L = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[0]);

    starpu_codelet_unpack_args( cl_arg, &cblk, &sopalin_data );

    assert( !(cblk->cblktype & CBLK_TASKS_2D) );

    nbpivot = cpucblk_zpotrfsp1d_panel( cblk, L, sopalin_data->diagthreshold,
                                        &(sopalin_data->solvmtx->lowrank) );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( cblk_zpotrfsp1d_panel, 1 )

void
starpu_task_cblk_zpotrfsp1d_panel( sopalin_data_t *sopalin_data,
                                   SolverCblk     *cblk,
                                   int             prio )
{
    starpu_insert_task(
        pastix_codelet(&cl_cblk_zpotrfsp1d_panel),
        STARPU_VALUE, &cblk,         sizeof(SolverCblk*),
        STARPU_VALUE, &sopalin_data, sizeof(sopalin_data_t*),
        STARPU_RW,     cblk->handler[0],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "cblk_zpotrfsp1d_panel",
#endif
        STARPU_PRIORITY, prio,
        0);
}

/**
 * Blok version
 */
static struct starpu_perfmodel starpu_blok_zpotrfsp_model =
{
    .type = STARPU_PER_ARCH,
    .symbol = "blok_zpotrfsp",
    .arch_cost_function = blok_potrf_cost,
};

#if !defined(PASTIX_STARPU_SIMULATION)
static void cl_blok_zpotrfsp_cpu(void *descr[], void *cl_arg)
{
    sopalin_data_t *sopalin_data;
    SolverCblk *cblk;
    pastix_complex64_t *L;
    int nbpivot;

    L = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[0]);

    starpu_codelet_unpack_args( cl_arg, &cblk, &sopalin_data );

    assert(cblk->cblktype & CBLK_TASKS_2D);

    nbpivot = cpucblk_zpotrfsp1d_potrf( cblk, L, sopalin_data->diagthreshold );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU( blok_zpotrfsp, 1 )

void
starpu_task_blok_zpotrf( sopalin_data_t *sopalin_data,
                         SolverCblk     *cblk,
                         int             prio )
{
    starpu_insert_task(
        pastix_codelet(&cl_blok_zpotrfsp),
        STARPU_VALUE, &cblk,         sizeof(SolverCblk*),
        STARPU_VALUE, &sopalin_data, sizeof(sopalin_data_t*),
        STARPU_RW,     cblk->fblokptr->handler[0],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "blok_zpotrfsp",
#endif
        STARPU_PRIORITY, prio,
        0);
}

/**
 * @}
 */
