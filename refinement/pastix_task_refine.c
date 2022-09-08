/**
 *
 * @file pastix_task_refine.c
 *
 * PaStiX refinement functions implementations.
 *
 * @copyright 2015-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @author Tony Delarue
 * @date 2022-10-11
 *
 **/
#include "common.h"
#include "bcsc/bcsc.h"
#include "refinement/z_refine_functions.h"
#include "refinement/c_refine_functions.h"
#include "refinement/d_refine_functions.h"
#include "refinement/s_refine_functions.h"
#include "pastix/order.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_dev_refine
 *
 * @brief Select the refinement function to call depending on the matrix type
 * and the precision
 *
 *******************************************************************************/
typedef pastix_int_t (*refine_fct_t)( pastix_data_t *, pastix_rhs_t, pastix_rhs_t );
static refine_fct_t sopalinRefine[4][4] =
{
    /* PastixRefineGMRES */
    {
        s_gmres_smp,
        d_gmres_smp,
        c_gmres_smp,
        z_gmres_smp
    },
    /* PastixRefineCG */
    {
        s_grad_smp,
        d_grad_smp,
        c_grad_smp,
        z_grad_smp
    },
    /* PastixRefineSR */
    {
        s_pivot_smp,
        d_pivot_smp,
        c_pivot_smp,
        z_pivot_smp
    },
    /* PastixRefineBiCGSTAB */
    {
        s_bicgstab_smp,
        d_bicgstab_smp,
        c_bicgstab_smp,
        z_bicgstab_smp
    }
};

/**
 *******************************************************************************
 *
 * @ingroup pastix_refine
 *
 * @brief Perform the iterative refinement without apply the permutations.
 *
 * This routine is affected by the following parameters:
 *   IPARM_REFINEMENT, DPARM_EPSILON_REFINEMENT
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 * @param[in] n
 *          The size of system to solve, and the number of rows of both
 *          matrices b and x.
 *
 * @param[in] nrhs
 *          The number of right hand side members, and the number of columns of
 *          b and x.
 *
 * @param[inout] b
 *          The right hand side matrix of size ldb-by-nrhs.
 *          B is noted as inout, as permutation might be performed on the
 *          matrix. On exit, the matrix is restored as it was on entry.
 *
 * @param[in] ldb
 *          The leading dimension of the matrix b. ldb >= n.
 *
 * @param[inout] x
 *          The matrix x of size ldx-by-nrhs.
 *          On entry, the initial guess x0 for the refinement step, that may be
 *          the solution returned by the solve step or any other initial guess.
 *          On exit, contains the final solution after the iterative refinement.
 *
 * @param[in] ldx
 *          The leading dimension of the matrix x. ldx >= n.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 *
 *******************************************************************************/
int
pastix_subtask_refine( pastix_data_t *pastix_data,
                       pastix_rhs_t   Bp,
                       pastix_rhs_t   Xp )
{
    pastix_int_t  *iparm = pastix_data->iparm;
    pastix_bcsc_t *bcsc  = pastix_data->bcsc;
    double         timer;
    pastix_int_t   n    = Bp->m;
    pastix_int_t   nrhs = Bp->n;
    pastix_int_t   ldb  = Bp->ld;
    pastix_int_t   ldx  = Xp->ld;
    const void    *b    = Bp->b;
    void          *x    = Xp->b;

    if (nrhs > 1)
    {
        pastix_print_warning( "Refinement works only with 1 rhs, We will iterate on each RHS one by one\n" );
    }

    if ( (pastix_data->schur_n > 0) && (iparm[IPARM_SCHUR_SOLV_MODE] != PastixSolvModeLocal))
    {
        fprintf(stderr, "Refinement is not available with Schur complement when non local solve is required\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Prepare the refinement threshold, if not set by the user */
    if ( pastix_data->dparm[DPARM_EPSILON_REFINEMENT] < 0. ) {
        int isDouble = (bcsc->flttype == PastixDouble) || (bcsc->flttype == PastixComplex64);
        if ( (!isDouble) ) {
            pastix_data->dparm[DPARM_EPSILON_REFINEMENT] = 1e-6;
        }
        else {
            pastix_data->dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
        }
    }

    clockSyncStart( timer, pastix_data->inter_node_comm );
    {
        refine_fct_t refinefct = sopalinRefine[iparm[IPARM_REFINEMENT]][pastix_data->bcsc->flttype -2];
        char *xptr = (char *)x;
        char *bptr = (char *)b;
        size_t shiftx, shiftb;
        int i;

        shiftx = ldx * pastix_size_of( Xp->flttype );
        shiftb = ldb * pastix_size_of( Bp->flttype );
        Bp->n  = 1;
        Xp->n  = 1;

        for(i=0; i<nrhs; i++, xptr += shiftx, bptr += shiftb ) {
            pastix_int_t it;
            Bp->b = bptr;
            Xp->b = xptr;
            it = refinefct( pastix_data, Xp, Bp );
            pastix_data->iparm[IPARM_NBITER] = pastix_imax( it, pastix_data->iparm[IPARM_NBITER] );
        }

        Bp->n = nrhs;
        Bp->b = (void*)b;
        Xp->n = nrhs;
        Xp->b = (void*)x;
    }
    clockSyncStop( timer, pastix_data->inter_node_comm );

    pastix_data->dparm[DPARM_REFINE_TIME] = clockVal(timer);
    if ( iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
        pastix_print( pastix_data->inter_node_procnum,
                      0, OUT_TIME_REFINE,
                      pastix_data->dparm[DPARM_REFINE_TIME] );
    }

    (void)n;
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_users
 *
 * @brief Perform iterative refinement.
 *
 * This routine performs the permutation of x, and b before and after the
 * iterative refinement solution. To prevent extra permuation to happen, see
 * pastix_subtask_refine().
 * This routine is affected by the following parameters:
 *   IPARM_REFINEMENT, DPARM_EPSILON_REFINEMENT
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The PaStiX data structure that describes the solver instance.
 *
 * @param[in] m
 *          The size of system to solve, and the number of rows of both
 *          matrices b and x.
 *
 * @param[in] nrhs
 *          The number of right hand side members, and the number of columns of
 *          b and x.
 *
 * @param[inout] b
 *          The right hand side matrix of size ldb-by-nrhs.
 *          B is noted as inout, as permutation might be performed on the
 *          matrix. On exit, the matrix is restored as it was on entry.
 *
 * @param[in] ldb
 *          The leading dimension of the matrix b. ldb >= n.
 *
 * @param[inout] x
 *          The matrix x of size ldx-by-nrhs.
 *          On entry, the initial guess x0 for the refinement step, that may be
 *          the solution returned by the solve step or any other initial guess.
 *          On exit, contains the final solution after the iterative refinement.
 *
 * @param[in] ldx
 *          The leading dimension of the matrix x. ldx >= n.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 *
 *******************************************************************************/
int
pastix_task_refine( pastix_data_t *pastix_data,
                    pastix_int_t   m,
                    pastix_int_t   nrhs,
                    void          *b,
                    pastix_int_t   ldb,
                    void          *x,
                    pastix_int_t   ldx )
{
    pastix_int_t  *iparm = pastix_data->iparm;
    pastix_bcsc_t *bcsc  = pastix_data->bcsc;
    pastix_rhs_t   Bp, Xp;
    void          *bglob, *btmp = NULL;
    void          *xglob, *xtmp = NULL;
    int            rc;

    if ( (pastix_data->schur_n > 0) && (iparm[IPARM_SCHUR_SOLV_MODE] != PastixSolvModeLocal))
    {
        fprintf(stderr, "Refinement is not available with Schur complement when non local solve is required\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* The spm is distributed, we have to gather it for the moment */
    if ( pastix_data->csc->loc2glob != NULL ) {
        if( iparm[IPARM_VERBOSE] > PastixVerboseNo ) {
            pastix_print( pastix_data->procnum, 0, "pastix_task_refine: the RHS has to be centralized for the moment\n" );
        }
        btmp = b;
        xtmp = x;

        m = pastix_data->csc->gNexp;

        bglob = malloc( m * nrhs * pastix_size_of( bcsc->flttype ) );
        xglob = malloc( m * nrhs * pastix_size_of( bcsc->flttype ) );

        spmGatherRHS( nrhs, pastix_data->csc, b, ldb, -1, bglob, m );
        spmGatherRHS( nrhs, pastix_data->csc, x, ldx, -1, xglob, m );

        b = bglob;
        x = xglob;
        ldb = m;
        ldx = m;
    }

    /* Compute P * b */
    rc = pastixRhsInit( pastix_data, &Bp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    rc = pastix_subtask_applyorder( pastix_data, PastixDirForward, m, nrhs, b, ldb, Bp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    /* Compute P * x */
    rc = pastixRhsInit( pastix_data, &Xp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    rc = pastix_subtask_applyorder( pastix_data, PastixDirForward, m, nrhs, x, ldx, Xp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    /* Performe the iterative refinement */
    rc = pastix_subtask_refine( pastix_data, Bp, Xp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    /* Compute P * b */
    rc = pastix_subtask_applyorder( pastix_data, PastixDirBackward, m, nrhs, b, ldb, Bp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }
    rc = pastixRhsFinalize( pastix_data, Bp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    /* Compute P * x */
    rc = pastix_subtask_applyorder( pastix_data, PastixDirBackward, m, nrhs, x, ldx, Xp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    rc = pastixRhsFinalize( pastix_data, Xp );
    if( rc != PASTIX_SUCCESS ) {
        return rc;
    }

    if ( btmp != NULL) {
        pastix_int_t ldbglob = ldb;
        ldb = pastix_data->csc->nexp;
        b   = btmp;

        spmExtractLocalRHS( nrhs, pastix_data->csc, bglob, ldbglob, b, ldb );
        memFree_null( bglob );
    }
    if ( xtmp != NULL ) {
        pastix_int_t ldxglob = ldx;
        ldx = pastix_data->csc->nexp;
        x   = xtmp;

        spmExtractLocalRHS( nrhs, pastix_data->csc, xglob, ldxglob, x, ldx );
        memFree_null( xglob );
    }

    (void)m;
    return PASTIX_SUCCESS;
}
