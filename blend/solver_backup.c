/**
 *
 * @file solver_backup.c
 *
 * PaStiX solver structure routines.
 *
 * @copyright 2004-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2024-07-05
 *
 **/
#include "common.h"
#include "blend/solver.h"

/**
 * @ingroup blend_dev_solver
 * @brief Structure to store backup of counter modified during numerical factorization and solve steps.
 */
struct SolverBackup_s {
    pastix_int_t *task_ctrbcnt;   /**< Number of contribution: the counter is decreased to 0 during factorization  */
    pastix_int_t *fanin_ctrbnbr;  /**< Number of contribution to FanIn: decreased during facto and solve           */
    pastix_int_t *fanin_prionum;  /**< Replaced by the number of msg packed during factorization sends             */
    pastix_int_t *symbol_cblknum; /**< Replaced by the negative FanIn index during facto and solve                 */
    pastix_int_t  symbol_nodenbr; /**< ???                                                                         */
    pastix_int_t  recvcnt;
    pastix_int_t  fanincnt;
};

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Initialize the backup structure.
 *
 * This function saves the initial values of the counters before
 * modifications. See structure description for more information about what is
 * saved and why.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix structure holding information for factorization
 *          and solve steps.
 *
 *******************************************************************************
 *
 * @retval Pointer to the allocated backup structure with the copy of
 *         information that might be destroyed.
 *
 *******************************************************************************/
SolverBackup_t *
solverBackupInit( const SolverMatrix *solvmtx )
{
    SolverBackup_t *b;
    pastix_int_t i;

    MALLOC_INTERN( b, 1, SolverBackup_t );
    memset( b, 0, sizeof(SolverBackup_t) );

    if (solvmtx->tasknbr)
    {
        Task *task = solvmtx->tasktab;

        MALLOC_INTERN(b->task_ctrbcnt, solvmtx->tasknbr, pastix_int_t);

        for (i=0; i<solvmtx->tasknbr; i++, task++)
        {
            b->task_ctrbcnt[i] = task->ctrbcnt;
        }
    }

    if (solvmtx->bloknbr) {
        SolverBlok *blok = solvmtx->bloktab;

        MALLOC_INTERN(b->symbol_cblknum, solvmtx->bloknbr, pastix_int_t);

        for (i=0; i<solvmtx->bloknbr; i++, blok++) {
            b->symbol_cblknum[i] = blok->fcblknm;
        }
    }

    b->symbol_nodenbr = solvmtx->nodenbr;

    {
        SolverCblk *cblk = solvmtx->cblktab;
        for (i=0; i<solvmtx->cblknbr; i++, cblk++)
        {
            cblk->ctrbcnt = cblk[1].brownum - cblk[0].brownum;
            cblk->partitioned = 0;
        }
    }

    b->recvcnt  = solvmtx->recvcnt;
    b->fanincnt = solvmtx->fanincnt;

    return b;
}

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Restore initial values.
 *
 * Restore counter values to be able to call a second factorization or solve
 * step. The amount of information restored depends on the value of
 * solvmtx->restore. If it is equal to:
 *     - 0: Nothing is restored
 *     - 1: A solve step has been performed and partial information is restored.
 *     - 2: A factorization step has been performed and full information is restored.
 * The value of solvmtx->restore is noramally initialized to 0 during the
 * structure creation, and then set to the correct value by the routine
 * modifying the solvmtx structure.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix structure holding information for factorization
 *          and solve steps.
 *          On exit, the counters have been restored to their original value
 *          stored in the backup structure.
 *
 * @param[in] b
 *          The backup structure pointer returned by the call to solverBackupInit().
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the data has been restored successfuly.
 * @retval PASTIX_ERR_BADPARAMETER if one of the parameter is incorrect.
 *
 *******************************************************************************/
int
solverBackupRestore( SolverMatrix         *solvmtx,
                     const SolverBackup_t *b       )
{
    pastix_int_t i;

    if ( solvmtx == NULL || b == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    if ( solvmtx->restore == 0 ) {
        return PASTIX_SUCCESS;
    }

    /* After factorization */
    if ( solvmtx->restore == 2 ) {
        if (solvmtx->tasknbr)
        {
            Task *task = solvmtx->tasktab;

            for (i=0; i<solvmtx->tasknbr; i++, task++)
            {
                task->ctrbcnt = b->task_ctrbcnt[i];
            }
        }
    }

    if (solvmtx->bloknbr) {
        SolverBlok *blok = solvmtx->bloktab;

        for (i=0; i<solvmtx->bloknbr; i++, blok++) {
            blok->fcblknm = b->symbol_cblknum[i];
        }
    }

    solvmtx->nodenbr  = b->symbol_nodenbr;
    solvmtx->recvcnt  = b->recvcnt;
    solvmtx->fanincnt = b->fanincnt;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver
 *
 * @brief Free the solver backup data structure.
 *
 * Clean the data structure holding the information backup and free the given
 * pointer because it has necessarily been allocated in solverBackupInit().
 *
 *******************************************************************************
 *
 * @param[inout] b
 *          The backup structure to destroy. On exit, b cannot be used anymore.
 *
 *******************************************************************************/
void
solverBackupExit( SolverBackup_t *b )
{
    if (b->task_ctrbcnt)
    {
        memFree_null(b->task_ctrbcnt);
    }

    if (b->fanin_ctrbnbr)
    {
        memFree_null(b->fanin_ctrbnbr);
        memFree_null(b->fanin_prionum);
    }

    if (b->symbol_cblknum) {
        memFree_null(b->symbol_cblknum);
    }
    memFree_null(b);
}
