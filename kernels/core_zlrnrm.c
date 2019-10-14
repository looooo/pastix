/**
 *
 * @file core_zlrnrm.c
 *
 * PaStiX low-rank kernel to compute the norms of a low-rank block.
 *
 * @copyright 2016-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @date 2019-10-10
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <lapacke.h>
#include "pastix_zlrcores.h"

/**
 *******************************************************************************
 *
 * @brief Compute the norm of a low-rank matrix.
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          The matrix norm to compute.
 *
 * @param[in] A
 *          The low-rank matrix
 *
 *******************************************************************************
 *
 * @return The norm of the matrix A
 *
 *******************************************************************************/
double
core_zlrnrm( pastix_normtype_t ntype, int transV,
             pastix_int_t M, pastix_int_t N,
             const pastix_lrblock_t *A )
{
    if ( ntype != PastixFrobeniusNorm ) {
        fprintf( stderr, "core_zlrnrm: Only the Frobenieu norm is available for now\n");
        ntype = PastixFrobeniusNorm;
    }

    if ( A->rk == -1 ) {
        assert( transV == PastixNoTrans );
        return LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f',
                                    M, N, A->u, M, NULL );
    }
    else if ( A->rk == 0 ) {
        return 0.;
    }
    else {
        double normU, normV;

        normU = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f',
                                     M, A->rk, A->u, M, NULL );
        if ( transV == PastixNoTrans ) {
            normV = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f',
                                         A->rkmax, N, A->v, A->rkmax, NULL );
        }
        else {
            normV = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f',
                                         N, A->rk, A->v, N, NULL );
        }
        return normU * normV;
    }
}
