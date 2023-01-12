/**
 *
 * @file cpucblk_zmpi_rhs.c
 *
 * Precision dependent routines to manag communications for the solve part.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2023-01-10
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common/common.h"
#include <lapacke.h>
#include "blend/solver.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

/**
 *******************************************************************************
 *
 * @brief Send the rhs associated to a cblk->lcolidx to the remote node.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix holding the communicator.
 *
 * @param[in] cblk
 *          The cblk which defines the part to sent.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
cpucblk_zsend_rhs_forward( const SolverMatrix *solvmtx,
                           SolverCblk         *cblk,
                           pastix_rhs_t        rhsb )
{
#if defined(PASTIX_WITH_MPI)
    pastix_complex64_t *b;
    pastix_int_t        colnbr = cblk_colnbr(cblk);
    pastix_int_t        size   = colnbr * rhsb->n;
    pastix_int_t        idx    = - cblk->bcscnum - 1;
    int                 rc;

    assert( colnbr <= solvmtx->colmax );
    assert( cblk->cblktype & CBLK_FANIN );
    assert( rhsb->cblkb[ idx ] != NULL );

#if defined (PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Fwd: Send cblk %ld to %2d at index %ld of size %ld\n",
             solvmtx->clustnum, (long)cblk->gcblknum, cblk->ownerid,
             (long)cblk->lcolidx, (long)colnbr );
#endif

    b = (pastix_complex64_t*)(rhsb->cblkb[ idx ]);
    assert( b != NULL );

    rc = MPI_Send( b, size, PASTIX_MPI_COMPLEX64,
                   cblk->ownerid, cblk->gcblknum, solvmtx->solv_comm );
    assert( rc == MPI_SUCCESS );

    memFree_null( rhsb->cblkb[ idx ] );

    (void)rc;
#else
    (void)solvmtx;
    (void)cblk;
    (void)rhsb;
#endif
}

/**
 *******************************************************************************
 *
 * @brief Send the rhs associated to a cblk->lcolidx to the remote node.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix holding the communicator.
 *
 * @param[in] cblk
 *          The cblk which defines the part to sent.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
cpucblk_zsend_rhs_backward( const SolverMatrix *solvmtx,
                            SolverCblk         *cblk,
                            pastix_rhs_t        rhsb )
{
#if defined(PASTIX_WITH_MPI)
    pastix_complex64_t *b   = rhsb->b;
    pastix_int_t        colnbr = cblk_colnbr(cblk);
    pastix_int_t        idx  = - cblk->bcscnum - 1;
    int rc;

    assert( colnbr <= solvmtx->colmax );
    assert( cblk->cblktype & CBLK_RECV );
    assert( rhsb->cblkb[ idx ] == NULL );

#if defined (PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Bwd: Send cblk %ld to %2d at index %ld of size %ld\n",
             solvmtx->clustnum, (long)cblk->gcblknum, cblk->ownerid,
             (long)cblk->lcolidx, (long)colnbr );
#endif

    b += cblk->lcolidx;
    if ( rhsb->n > 1 ) {
        rhsb->cblkb[ idx ] = malloc( colnbr * rhsb->n * sizeof(pastix_complex64_t) );
        rc = LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', colnbr, rhsb->n, b, rhsb->ld, rhsb->cblkb[ idx ], colnbr );
        assert( rc == 0 );
        b = rhsb->cblkb[ idx ];
    }

    rc = MPI_Send( b, colnbr * rhsb->n, PASTIX_MPI_COMPLEX64,
                   cblk->ownerid, cblk->gcblknum, solvmtx->solv_comm );
    assert( rc == MPI_SUCCESS );

    if ( rhsb->n > 1 ) {
        memFree_null( rhsb->cblkb[ idx ] );
    }

    (void)rc;
#else
    (void)solvmtx;
    (void)cblk;
    (void)rhsb;
#endif
}

/**
 *******************************************************************************
 *
 * @brief Receive the rhs associated to a cblk->lcolidx to the remote node.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix holding the communicator.
 *
 * @param[in] cblk
 *          The cblk which may define the part to sent.
 *
 * @param[inout] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
cpucblk_zrecv_rhs_backward( const SolverMatrix *solvmtx,
                            SolverCblk         *cblk,
                            pastix_rhs_t        rhsb )
{
#if defined(PASTIX_WITH_MPI)
    MPI_Status   status;
    pastix_int_t colnbr = cblk_colnbr(cblk);
    pastix_int_t idx  = - cblk->bcscnum - 1;
    int rc;

    assert( colnbr <= solvmtx->colmax );
    assert( cblk->cblktype & CBLK_FANIN );

#if defined (PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Bwd: Recv cblk %ld from %ld at index %ld of size %ld\n",
             solvmtx->clustnum, (long)cblk->gcblknum, (long)cblk->ownerid,
             (long)cblk->lcolidx, (long)colnbr );
#endif

    assert( rhsb->cblkb[ idx ] == NULL );
    rhsb->cblkb[ idx ] = malloc( colnbr * rhsb->n * sizeof(pastix_complex64_t) );

    rc = MPI_Recv( rhsb->cblkb[ idx ], colnbr * rhsb->n, PASTIX_MPI_COMPLEX64,
                   cblk->ownerid, cblk->gcblknum, solvmtx->solv_comm, &status );
    assert( rc == MPI_SUCCESS );

#if defined (PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Bwd: Received cblk %ld from %2d\n",
             solvmtx->clustnum, (long)cblk->gcblknum, status.MPI_SOURCE );
#endif

    (void)rc;
#else
    (void)solvmtx;
    (void)cblk;
    (void)rhsb;
#endif
}

/**
 *******************************************************************************
 *
 * @brief Receive the rhs associated to a cblk->lcolidx to the remote node.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix holding the communicator.
 *
 * @param[in] cblk
 *          The cblk which may define the part to sent.
 *
 * @param[inout] work
 *          The temporary buffer to receive the remote data
 *
 * @param[inout] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
cpucblk_zrecv_rhs_forward( const SolverMatrix *solvmtx,
                           SolverCblk         *cblk,
                           pastix_complex64_t *work,
                           pastix_rhs_t        rhsb )
{
#if defined(PASTIX_WITH_MPI)
    pastix_complex64_t *b      = rhsb->b;
    pastix_int_t        colnbr = cblk_colnbr(cblk);
    MPI_Status          status;
    int rc;

    assert( colnbr <= solvmtx->colmax );
    assert( cblk->cblktype & CBLK_RECV );

#if defined (PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Fwd: Recv cblk %ld from %ld at index %ld of size %ld\n",
             solvmtx->clustnum, (long)cblk->gcblknum, (long)cblk->ownerid,
             (long)cblk->lcolidx, (long)colnbr );
#endif

    rc = MPI_Recv( work, colnbr * rhsb->n, PASTIX_MPI_COMPLEX64,
                   cblk->ownerid, cblk->gcblknum, solvmtx->solv_comm, &status );
    assert( rc == MPI_SUCCESS );

#if defined (PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Fwd: Received cblk %ld from %2d\n",
                     solvmtx->clustnum, (long)cblk->gcblknum, status.MPI_SOURCE );
#endif

    b += cblk->lcolidx;
    core_zgeadd( PastixNoTrans, colnbr, rhsb->n,
                 1., work, colnbr,
                 1., b,    rhsb->ld );

    (void)rc;
#else
    (void)solvmtx;
    (void)cblk;
    (void)work;
    (void)rhsb;
#endif
}
