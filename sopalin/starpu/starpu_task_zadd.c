/**
 *
 * @file starpu_task_zadd.c
 *
 * PaStiX zadd StarPU wrapper.
 *
 * @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.1
 * @author Mathieu Faverge
 * @author Alycia Lisito
 * @date 2023-12-01
 * @precisions normal z -> s d c
 *
 **/

#include "common.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_starpu.h"
#include "pastix_zstarpu.h"

/**
 *******************************************************************************
 *
 * @brief Submits starpu zadd task to send the fanin cblk.
 *
 *******************************************************************************
 *
 * @param[in] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 * @param[in] side
 *          The side of the update.
 *
 * @param[in] cblk
 *          The fanin column block.
 *
 * @param[in] prio
 *          The task priority.
 *
 *******************************************************************************/
void
starpu_task_zadd_1dp_fanin( sopalin_data_t    *sopalin_data,
                            pastix_coefside_t  side,
                            const SolverCblk  *cblk,
                            int                prio )
{
    /* If this is a fanin, let's submit the send */
    assert( cblk->cblktype & CBLK_FANIN );
    starpu_task_cblk_zadd_fanin( sopalin_data, side, cblk, prio );
}

/**
 *******************************************************************************
 *
 * @brief Submits starpu zadd task to receive and add the recv cblk.
 *
 *******************************************************************************
 *
 * @param[in] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 * @param[in] side
 *          The side of the update.
 *
 * @param[in] cblk
 *          The recv column block.
 *
 * @param[in] prio
 *          The task priority.
 *
 *******************************************************************************/
void
starpu_task_zadd_1dp_recv( sopalin_data_t    *sopalin_data,
                           pastix_coefside_t  side,
                           const SolverCblk  *cblk,
                           int                prio )
{
    SolverCblk *fcblk = sopalin_data->solvmtx->cblktab + cblk->fblokptr->fcblknm;

    assert( cblk->cblktype & CBLK_RECV );
    starpu_task_cblk_zadd_recv( sopalin_data, side, cblk, fcblk, prio );
}

/**
 *******************************************************************************
 *
 * @brief Submits starpu zadd task to send the fanin block.
 *
 *******************************************************************************
 *
 * @param[in] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 * @param[in] side
 *          The side of the update.
 *
 * @param[in] cblk
 *          The fanin column block.
 *
 * @param[in] prio
 *          The task priority.
 *
 *******************************************************************************/
void
starpu_task_zadd_2d_fanin( sopalin_data_t    *sopalin_data,
                           pastix_coefside_t  side,
                           const SolverCblk  *cblk,
                           int                prio )
{
    const SolverBlok *blok = cblk[0].fblokptr;
    const SolverBlok *lblk = cblk[1].fblokptr;

    assert( cblk->cblktype & CBLK_FANIN );

    for ( ; blok < lblk; blok++ ) {

        starpu_task_blok_zadd_fanin( sopalin_data, side, cblk, blok, prio );

        /* Skip blocks facing the same cblk */
        while ( ( blok < lblk ) &&
                ( blok[0].fcblknm == blok[1].fcblknm ) &&
                ( blok[0].lcblknm == blok[1].lcblknm ) )
        {
            blok++;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Submits starpu zadd task to receive and add the recv block.
 *
 *******************************************************************************
 *
 * @param[in] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 * @param[in] side
 *          The side of the update.
 *
 * @param[in] cblk
 *          The recv column block.
 *
 * @param[in] prio
 *          The task priority.
 *
 *******************************************************************************/
void
starpu_task_zadd_2d_recv( sopalin_data_t    *sopalin_data,
                          pastix_coefside_t  side,
                          const SolverCblk  *cblk,
                          int                prio )
{
    const SolverBlok *blok = cblk[0].fblokptr;
    const SolverBlok *lblk = cblk[1].fblokptr;
    SolverCblk *fcblk  = sopalin_data->solvmtx->cblktab + cblk->fblokptr->fcblknm;
    SolverBlok *fblok  = fcblk[0].fblokptr;
    SolverBlok *lfblok = fcblk[1].fblokptr;

    assert( cblk->cblktype & CBLK_RECV );

    for ( ; blok < lblk; blok++ ) {
        while ( (blok->fcblknm != fblok->fcblknm) && (fblok < lfblok) ) {
            fblok++;
        }
        assert( fblok->lcblknm == cblk->fblokptr->fcblknm );
        starpu_task_blok_zadd_recv( sopalin_data, side, cblk, blok, fcblk, fblok, prio );

        /* Skip blocks facing the same cblk */
        while ( ( blok < lblk ) &&
                ( blok[0].fcblknm == blok[1].fcblknm ) &&
                ( blok[0].lcblknm == blok[1].lcblknm ) )
        {
            blok++;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Submits starpu zadd task to send the recv cblk.
 *
 *******************************************************************************
 *
 * @param[in] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 * @param[in] side
 *          The side of the update.
 *
 * @param[in] cblk
 *          The recv column block.
 *
 * @param[in] prio
 *          The task priority.
 *
 *******************************************************************************/
void
starpu_task_zadd_recv( sopalin_data_t    *sopalin_data,
                       pastix_coefside_t  side,
                       const SolverCblk  *cblk,
                       int                prio )
{
    if ( cblk->cblktype & CBLK_TASKS_2D ) {
        starpu_task_zadd_2d_recv( sopalin_data, side, cblk, prio );
        return;
    }
    else {
        starpu_task_zadd_1dp_recv( sopalin_data, side, cblk, prio );
        return;
    }
}

/**
 *******************************************************************************
 *
 * @brief Submits starpu zadd task to send the fanin cblk.
 *
 *******************************************************************************
 *
 * @param[in] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 * @param[in] side
 *          The side of the update.
 *
 * @param[in] cblk
 *          The fanin column block.
 *
 * @param[in] prio
 *          The task priority.
 *
 *******************************************************************************/
void
starpu_task_zadd_fanin( sopalin_data_t    *sopalin_data,
                        pastix_coefside_t  side,
                        const SolverCblk  *cblk,
                        int                prio )
{
    if ( cblk->cblktype & CBLK_TASKS_2D ) {
        starpu_task_zadd_2d_fanin( sopalin_data, side, cblk, prio );
        return;
    }
    else {
        starpu_task_zadd_1dp_fanin( sopalin_data, side, cblk, prio );
        return;
    }
}
