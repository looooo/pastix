/**
 *
 * @file starpu_zdiag.c
 *
 * PaStiX zdiag StarPU wrapper.
 *
 * @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Vincent Bridonneau
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2023-07-21
 * @precisions normal z -> s d c
 *
 * @addtogroup starpu_diag_solve
 * @{
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
 * @brief Apply the diagonal solve on a problem of 1 dimension.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          Provide informations about starpu and the schur solving mode.
 *
 * @param[in] sopalin_data
 *          Provide descriptor of b and the sparse matrix and the SolverMatrix
 *          structure.
 *
 * @param[inout] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
starpu_zdiag_sp1dplus( pastix_data_t  *pastix_data,
                       sopalin_data_t *sopalin_data,
                       pastix_rhs_t    rhsb )
{
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    pastix_int_t        k, cblknbr;
    pastix_solv_mode_t  mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];

    cblk = datacode->cblktab;
    cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;
    for (k=0; k<cblknbr; k++, cblk++) {
        if ( cblk->ownerid != datacode->clustnum ) {
            continue;
        }
        starpu_stask_cblk_zdiag( sopalin_data, rhsb, cblk, cblknbr - k );
    }
}

/**
 *******************************************************************************
 *
 * @brief Apply the diagonal solve (StarPU version)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          Provide informations about starpu and the schur solving mode.
 *
 * @param[in] sopalin_data
 *          Provide descriptor of b and the sparse matrix and the SolverMatrix
 *          structure.
 *
 * @param[inout] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
starpu_zdiag( pastix_data_t  *pastix_data,
              sopalin_data_t *sopalin_data,
              pastix_rhs_t    rhsb )
{
    starpu_sparse_matrix_desc_t *sdesc = sopalin_data->solvmtx->starpu_desc;
    starpu_rhs_desc_t           *ddesc = rhsb->starpu_desc;

    /*
     * Start StarPU if not already started
     */
    if (pastix_data->starpu == NULL) {
        int argc = 0;
        pastix_starpu_init( pastix_data, &argc, NULL, NULL );
    }

    if ( sdesc == NULL ) {
        /* Create the sparse matrix descriptor */
        starpu_sparse_matrix_init( sopalin_data->solvmtx,
                                   PastixHermitian,
                                   pastix_data->inter_node_procnbr,
                                   pastix_data->inter_node_procnum,
                                   PastixComplex64 );
        sdesc = sopalin_data->solvmtx->starpu_desc;
    }

    if ( ddesc == NULL ) {
        /* Create the dense matrix descriptor */
        starpu_rhs_init( pastix_data->solvmatr, rhsb,
                         PastixComplex64,
                         pastix_data->inter_node_procnbr,
                         pastix_data->inter_node_procnum );
        ddesc = rhsb->starpu_desc;
    }

#if defined(STARPU_USE_FXT)
    if (pastix_data->iparm[IPARM_TRACE] & PastixTraceSolve) {
        starpu_fxt_start_profiling();
    }
#endif
    starpu_resume();
    starpu_zdiag_sp1dplus( pastix_data, sopalin_data, rhsb );

    starpu_sparse_matrix_getoncpu( sdesc );
    starpu_rhs_getoncpu( ddesc );
    starpu_task_wait_for_all();
#if defined(PASTIX_WITH_MPI)
    starpu_mpi_wait_for_all( pastix_data->pastix_comm );
    starpu_mpi_barrier(pastix_data->inter_node_comm);
#endif
    starpu_pause();
#if defined(STARPU_USE_FXT)
    if (pastix_data->iparm[IPARM_TRACE] & PastixTraceSolve) {
        starpu_fxt_stop_profiling();
    }
#endif

    return;
}

/**
 *@}
 */
