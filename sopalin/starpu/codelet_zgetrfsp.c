/**
 *
 * @file codelet_zgetrfsp.c
 *
 * StarPU codelets for LU functions
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
 * StarPU models
 */

static struct starpu_perfmodel starpu_cblk_zgetrfsp1d_panel_model =
{
  .type = STARPU_PER_ARCH,
  .arch_cost_function = cblk_getrfsp1d_cost,
  .symbol = "cblk_zgetrfsp1d_panel",
};

static struct starpu_perfmodel starpu_blok_zgetrfsp1d_panel_model =
{
  .type = STARPU_PER_ARCH,
  .arch_cost_function = blok_getrfsp1d_cost,
  .symbol = "blok_zgetrfsp1d_panel",
};

/**
 * Cblk version
 */

#if !defined(PASTIX_STARPU_SIMULATION)
static void cl_cblk_zgetrfsp1d_panel_cpu(void *descr[], void *cl_arg)
{
    sopalin_data_t *sopalin_data;
    SolverCblk *cblk;
    pastix_complex64_t *L;
    pastix_complex64_t *U;
    int nbpivot;

    L = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[0]);
    U = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &cblk, &sopalin_data);

    nbpivot = cpucblk_zgetrfsp1d_panel( cblk, L, U, sopalin_data->diagthreshold,
                                        &(sopalin_data->solvmtx->lowrank) );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU_MODEL( cblk_zgetrfsp1d_panel, 2,starpu_cblk_zgetrfsp1d_panel_model)

void
starpu_task_cblk_zgetrfsp1d_panel( sopalin_data_t *sopalin_data,
                                   SolverCblk     *cblk,
                                   int             prio )
{
    starpu_insert_task(
        pastix_codelet(&cl_cblk_zgetrfsp1d_panel),
        STARPU_VALUE, &cblk,             sizeof(SolverCblk*),
        STARPU_VALUE, &sopalin_data,     sizeof(sopalin_data_t*),
        STARPU_RW,     cblk->handler[0],
        STARPU_RW,     cblk->handler[1],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "cblk_zgetrfsp1d_panel",
#endif
        STARPU_PRIORITY, prio,
        0);
}

/**
 * Blok version
 */

#if !defined(PASTIX_STARPU_SIMULATION)
static void cl_blok_zgetrfsp_cpu(void *descr[], void *cl_arg)
{
    sopalin_data_t *sopalin_data;
    SolverCblk *cblk;
    pastix_complex64_t *L, *U;
    int nbpivot;

    L = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[0]);
    U = (pastix_complex64_t *)STARPU_VECTOR_GET_PTR(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &cblk, &sopalin_data);

    assert(cblk->cblktype & CBLK_TASKS_2D);

    nbpivot = cpucblk_zgetrfsp1d_getrf( cblk, L, U, sopalin_data->diagthreshold );

    (void)nbpivot;
}
#endif /* !defined(PASTIX_STARPU_SIMULATION) */

CODELETS_CPU_MODEL( blok_zgetrfsp, 2, starpu_blok_zgetrfsp1d_panel_model )

void
starpu_task_blok_zgetrf( sopalin_data_t *sopalin_data,
                         SolverCblk     *cblk,
                         int             prio)
{
    starpu_insert_task(
        pastix_codelet(&cl_blok_zgetrfsp),
        STARPU_VALUE, &cblk,             sizeof(SolverCblk*),
        STARPU_VALUE, &sopalin_data,     sizeof(sopalin_data_t*),
        STARPU_RW,     cblk->fblokptr->handler[0],
        STARPU_RW,     cblk->fblokptr->handler[1],
#if defined(PASTIX_STARPU_CODELETS_HAVE_NAME)
        STARPU_NAME, "blok_zgetrfsp",
#endif
        STARPU_PRIORITY, prio,
        0);
}

/**
 * @}
 */
