/**
 *
 * @file cpucblk_zmpi_rhs_fwd.c
 *
 * Precision dependent routines to manag communications for the solve part.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.1
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2023-09-20
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common/common.h"
#include "common/pastixdata.h"
#include <lapacke.h>
#include "blend/solver.h"
#include "blend/solver_comm_matrix.h"
#include "kernels.h"
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

#if defined( PASTIX_WITH_MPI )
/**
 *******************************************************************************
 *
 * @brief Asynchronously sends the rhs associated to cblk->ownerid
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 * @param[in] cblk
 *          The column block that will be sent.
 *
 *******************************************************************************/
static inline void
cpucblk_zisend_rhs_fwd( SolverMatrix *solvmtx,
                        pastix_rhs_t  rhsb,
                        SolverCblk   *cblk )
{
    MPI_Request         request;
    pastix_complex64_t *b;
    pastix_int_t        colnbr = cblk_colnbr(cblk);
    pastix_int_t        size   = colnbr * rhsb->n;
    pastix_int_t        idx    = - cblk->bcscnum - 1;
    int                 rc;

    assert( colnbr <= solvmtx->colmax );
    assert( cblk->cblktype & CBLK_FANIN );
    assert( rhsb->cblkb[ idx ] != NULL );

#if defined(PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] RHS Fwd: Post Isend cblk %ld to %2d at index %ld of size %ld\n",
             solvmtx->clustnum, (long)cblk->gcblknum, cblk->ownerid,
             (long)cblk->lcolidx, (long)(size * sizeof(pastix_complex64_t) ) );
#endif

    b = (pastix_complex64_t*)(rhsb->cblkb[ idx ]);
    assert( b != NULL );

    rc = MPI_Isend( b, size, PASTIX_MPI_COMPLEX64, cblk->ownerid, cblk->gcblknum,
                    solvmtx->solv_comm, &request );
    assert( rc == MPI_SUCCESS );

    solverCommMatrixAdd( solvmtx, cblk->ownerid, size * sizeof(pastix_complex64_t) );

    /* Register the request to make it progress */
    pastix_atomic_lock( &(solvmtx->reqlock) );

    assert( solvmtx->reqidx[ solvmtx->reqnum ] == -1 );
    assert( solvmtx->reqnum >= 0 );
    assert( solvmtx->reqnum < solvmtx->reqnbr );

    solvmtx->reqtab[ solvmtx->reqnum ] = request;
    solvmtx->reqidx[ solvmtx->reqnum ] = cblk - solvmtx->cblktab;
    solvmtx->reqnum++;

    pastix_atomic_unlock( &(solvmtx->reqlock) );

    (void)rc;
}

/**
 *******************************************************************************
 *
 * @brief Handle a finished request on a fanin
 *
 *******************************************************************************
 *
 * @param[in] enums
 *          Enums needed for the solve.
 *
 * @param[in] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 * @param[inout] cblk
 *          The cblk concerned by the computation.
 *
 *******************************************************************************/
void
cpucblk_zrequest_rhs_fwd_handle_send( const args_solve_t *enums,
                                      SolverMatrix       *solvmtx,
                                      pastix_rhs_t        rhsb,
                                      const SolverCblk   *cblk )
{
    pastix_int_t idx = - cblk->bcscnum - 1;

    assert( cblk->cblktype & CBLK_FANIN );
    assert( enums->solve_step == PastixSolveForward );

#if defined(PASTIX_DEBUG_MPI)
    {
        size_t cblksize = cblk_colnbr( cblk ) * rhsb->n * sizeof(pastix_complex64_t);

        fprintf( stderr, "[%2d] RHS Fwd: Isend for cblk %ld toward %2d ( %ld Bytes ) (DONE)\n",
                 solvmtx->clustnum, (long)cblk->gcblknum, cblk->ownerid, (long)cblksize );
    }
#endif

    memFree_null( rhsb->cblkb[ idx ] );
    (void)solvmtx;
}

/**
 *******************************************************************************
 *
 * @brief Handle a finished request on a recv cblk.
 *
 *******************************************************************************
 *
 * @param[in] enums
 *          Enums needed for the solve.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 * @param[inout] threadid
 *          Id of the thread calling this method.
 *
 * @param[inout] status
 *          Statuses of the completed request.
 *
 * @param[inout] recvbuf
 *          Reception buffer of the completed request.
 *
 *******************************************************************************/
static inline void
cpucblk_zrequest_rhs_fwd_handle_recv( const args_solve_t *enums,
                                      SolverMatrix       *solvmtx,
                                      pastix_rhs_t        rhsb,
                                      int                 threadid,
                                      const MPI_Status   *status,
                                      pastix_complex64_t *recvbuf )
{
    SolverCblk   *cblk, *fcbk;
    pastix_int_t  colnbr;
    int           src  = status->MPI_SOURCE;
    int           tag  = status->MPI_TAG;
    int           nrhs = rhsb->n;

    assert( ( 0 <= src ) && ( src < solvmtx->clustnbr ) );
    assert( ( 0 <= tag ) && ( tag < solvmtx->gcblknbr ) );

    /*
     * Let's look for the local cblk
     */
    fcbk = solvmtx->cblktab + solvmtx->gcbl2loc[ tag ];
    cblk = fcbk-1;

    /* Get through source */
    while( cblk->ownerid != src ) {
        cblk--;
        assert( cblk >= solvmtx->cblktab );
        assert( cblk->gcblknum == tag );
        assert( cblk->cblktype & CBLK_RECV );
    }
    assert( fcbk == (solvmtx->cblktab + cblk->fblokptr->fcblknm) );

    colnbr = cblk_colnbr( cblk );
#if defined(PASTIX_DEBUG_MPI)
    {
        int          rc;
        int          count = 0;
        pastix_int_t size  = colnbr * nrhs * sizeof(pastix_complex64_t);

        rc = MPI_Get_count( status, MPI_CHAR, &count );
        assert( rc == MPI_SUCCESS );
        assert( count == size );

        /* We can't know the sender easily, so we don't print it */
        fprintf( stderr, "[%2d] RHS Fwd  : recv of size %d/%ld for cblk %ld (DONE)\n",
                 solvmtx->clustnum, count, (long)size, (long)cblk->gcblknum );
    }
#endif

    /* Initialize the cblk with the reception buffer */
    cblk->threadid = (fcbk->threadid == -1) ? threadid : fcbk->threadid;

    {
        pastix_complex64_t *b = rhsb->b;
        b += cblk->lcolidx;
        pastix_cblk_lock( fcbk );
        core_zgeadd( PastixNoTrans, colnbr, nrhs,
                     1., recvbuf, colnbr,
                     1., b,       rhsb->ld );
        pastix_cblk_unlock( fcbk );
    }

    /* Receptions cblks contribute to themselves */
    cpucblk_zrelease_rhs_fwd_deps( enums, solvmtx, rhsb, cblk, fcbk );

    /* Free the CBLK_RECV */
    memFree_null( recvbuf );
}

/**
 *******************************************************************************
 *
 * @brief Handle a finished request.
 *
 * If cblktype & CBLK_FANIN : Will deallocate the coeftab
 * If cblktype & CBLK_RECV  : Will add cblk and deallocate the coeftab
 *
 *******************************************************************************
 *
 * @param[in] enums
 *          Enums needed for the solve.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 * @param[in] threadid
 *          Id of the thread calling this method.
 *
 * @param[in] outcount
 *          Amount of finshed requests.
 *
 * @param[in] indexes
 *          Array of completed requests.
 *
 * @param[in] statuses
 *          Array of statuses for the completed requests.
 *
 *******************************************************************************
 *
 * @retval The amount of finished requests updated.
 *
 *******************************************************************************/
static inline int
cpucblk_zrequest_rhs_fwd_handle( const args_solve_t *enums,
                                 SolverMatrix       *solvmtx,
                                 pastix_rhs_t        rhsb,
                                 int                 threadid,
                                 int                 outcount,
                                 const int          *indexes,
                                 const MPI_Status   *statuses )
{
    pastix_int_t i, reqid;
    int          nbrequest = outcount;

    for ( i = 0; i < outcount; i++ ) {
        reqid = indexes[i];

        /*
         * Handle the reception
         */
        if ( solvmtx->reqidx[reqid] == -1 ) {
            /* We're on a cblk recv, copy datas and restart communications */
            pastix_complex64_t *recvbuf;
            MPI_Status status;
            int        size;

            memcpy( &status, statuses + i, sizeof(MPI_Status) );
            MPI_Get_count( &status, PASTIX_MPI_COMPLEX64, &size );

            MALLOC_INTERN( recvbuf, size, pastix_complex64_t );
            memcpy( recvbuf, solvmtx->rcoeftab, size * sizeof(pastix_complex64_t) );

            solvmtx->recvcnt--;

            /* Let's restart the communication */
            assert( solvmtx->recvcnt >= 0 );
            if ( solvmtx->recvcnt > 0 ) {
                MPI_Start( solvmtx->reqtab + reqid );
                nbrequest--;
            }
            else {
                MPI_Request_free( solvmtx->reqtab + reqid );
                solvmtx->reqtab[reqid] = MPI_REQUEST_NULL;
            }

            cpucblk_zrequest_rhs_fwd_handle_recv( enums, solvmtx, rhsb,
                                                  threadid, &status, recvbuf );
        }
        /*
         * Handle the emission
         */
        else {
            SolverCblk *cblk = solvmtx->cblktab + solvmtx->reqidx[ reqid ];
            assert( cblk->cblktype & CBLK_FANIN );

            cpucblk_zrequest_rhs_fwd_handle_send( enums, solvmtx, rhsb, cblk );

#if !defined(NDEBUG)
            solvmtx->reqidx[ reqid ] = -1;
#endif
            solvmtx->fanincnt--;
        }
    }

    return nbrequest;
}

/**
 *******************************************************************************
 *
 * @ingroup kernel_fact
 * @brief Progress communications for one process
 *
 * If a communication is completed, it will be treated.
 * If cblktype & CBLK_FANIN : Will deallocate coeftab
 * If cblktype & CBLK_RECV  : Will add cblk to fcblk
 *
 *******************************************************************************
 *
 * @param[in] enums
 *          Enums needed for the solve.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 * @param[in] threadid
 *          Id of the thread calling this method.
 *
 *******************************************************************************/
void
cpucblk_zmpi_rhs_fwd_progress( const args_solve_t *enums,
                               SolverMatrix       *solvmtx,
                               pastix_rhs_t        rhsb,
                               int                 threadid )
{
    pthread_t  tid = pthread_self();
    int        outcount = 1;
    int        nbrequest, nbfree;
    int        indexes[ solvmtx->reqnbr ];
    MPI_Status statuses[ solvmtx->reqnbr ];

    /* Check if someone is already communicating or not */
    pthread_mutex_lock( &pastix_comm_lock );
    if ( pastix_comm_tid == (pthread_t)-1 ) {
        pastix_comm_tid = tid;
    }
    pthread_mutex_unlock( &pastix_comm_lock );

    if ( tid != pastix_comm_tid ) {
        pastix_yield();
        return;
    }

    /*
     * Let's register the number of active requests.
     * We now suppose that the current thread is working on the first nbrequest
     * active in the reqtab array. Additional requests can be posted during this
     * progression, but it will be with a larger index. Thus, we do not need to
     * protect every changes in these requests.
     * When this is done, the requests arrays is locked to be packed, and the
     * number of requests is updated for the next round.
     */
    pastix_atomic_lock( &(solvmtx->reqlock) );
    nbrequest = solvmtx->reqnum;
    pastix_atomic_unlock( &(solvmtx->reqlock) );

    while( (outcount > 0) && (nbrequest > 0) )
    {
        MPI_Testsome( nbrequest, solvmtx->reqtab, &outcount, indexes, statuses );
        nbfree = 0;

        /* Handle all the completed requests */
        if ( outcount > 0 ) {
            nbfree = cpucblk_zrequest_rhs_fwd_handle( enums, solvmtx, rhsb, threadid,
                                                      outcount, indexes, statuses );
        }

        /*
         * Pack the request arrays, and update the number of active requests by
         * removing the completed ones
         */
        pastix_atomic_lock( &(solvmtx->reqlock) );
        if ( nbfree > 0 ) {
            cpucblk_zupdate_reqtab( solvmtx );
        }
        nbrequest = solvmtx->reqnum;
        pastix_atomic_unlock( &(solvmtx->reqlock) );
    }

    pastix_comm_tid = (pthread_t)-1;
    pastix_yield();
}
#endif /* defined(PASTIX_WITH_MPI) */

/**
 *******************************************************************************
 *
 * @brief Wait for incoming dependencies, and return when cblk->ctrbcnt has
 * reached 0.
 *
 *******************************************************************************
 *
 * @param[in] rank
 *          The rank of the current thread.
 *
 * @param[in] enums
 *          Enums needed for the solve.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[inout] cblk
 *          The column block that contribute to fcblk.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************
 *
 * @return 1 if the cblk is a fanin, 0 otherwise
 *
 *******************************************************************************/
int
cpucblk_zincoming_rhs_fwd_deps( int                 rank,
                                const args_solve_t *enums,
                                SolverMatrix       *solvmtx,
                                SolverCblk         *cblk,
                                pastix_rhs_t        rhsb )
{
#if defined(PASTIX_WITH_MPI)
    if ( cblk->cblktype & CBLK_FANIN ) {
        /*
         * We are in the sequential case, we progress on communications and
         * return if nothing.
         */
        //cpucblk_ztestsome( side, solvmtx );
        return 1;
    }

    if ( cblk->cblktype & CBLK_RECV ) {
        return 1;
    }

    /* Make sure we receive every contribution */
    while( cblk->ctrbcnt > 0 ) {
        cpucblk_zmpi_rhs_fwd_progress( enums, solvmtx, rhsb, rank );
    }
#else
    assert( !(cblk->cblktype & (CBLK_FANIN | CBLK_RECV)) );
    do { pastix_yield(); } while( cblk->ctrbcnt > 0 );
#endif

    (void)rank;
    (void)enums;
    (void)solvmtx;
    (void)rhsb;

    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Release the dependencies of the given cblk after an update.
 *
 *******************************************************************************
 *
 * @param[in] enums
 *          Enums needed for the solve.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 * @param[in] cblk
 *          The column block that contribute to fcblk.
 *
 * @param[inout] fcbk
 *          The facing column block that is updated by cblk.
 *
 *******************************************************************************/
void
cpucblk_zrelease_rhs_fwd_deps( const args_solve_t *enums,
                               SolverMatrix       *solvmtx,
                               pastix_rhs_t        rhsb,
                               const SolverCblk   *cblk,
                               SolverCblk         *fcbk )
{
    int32_t ctrbcnt;
    ctrbcnt = pastix_atomic_dec_32b( &(fcbk->ctrbcnt) );
    if ( !ctrbcnt ) {
#if defined(PASTIX_WITH_MPI)
        if ( ( fcbk->cblktype & CBLK_FANIN ) &&
             ( enums->solve_step == PastixSolveForward ) ) {
                cpucblk_zisend_rhs_fwd( solvmtx, rhsb, fcbk );
                return;
        }
#else
        (void)enums;
        (void)rhsb;
#endif
        if ( solvmtx->computeQueue ) {
            pastix_queue_t *queue = solvmtx->computeQueue[ cblk->threadid ];
            assert( fcbk->priority != -1 );
            pqueuePush1( queue, fcbk - solvmtx->cblktab, fcbk->priority );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief  Waitall routine for current cblk request
 *
 * It may be possible that some cblk will not be deallocated with the static
 * scheduler. So a cleanup may be necessary.
 *
 *******************************************************************************
 *
 * @param[in] enums
 *          Enums needed for the solve.
 *
 * @param[in] sched
 *          Define which sched is used
 *          @arg PastixSchedSequential if sequential
 *          @arg PastixSchedStatic if multi-threaded static scheduler
 *          @arg PastixSchedDynamic if multi-threaded dynamic scheduler
 *          No other scheduler is supported.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
cpucblk_zrequest_rhs_fwd_cleanup( const args_solve_t *enums,
                                  pastix_int_t        sched,
                                  SolverMatrix       *solvmtx,
                                  pastix_rhs_t        rhsb )
{
    if ( (sched != PastixSchedSequential) &&
         (sched != PastixSchedStatic)     &&
         (sched != PastixSchedDynamic) )
    {
        return;
    }
#if defined(PASTIX_WITH_MPI)
    pastix_int_t i;
    int          rc;
    SolverCblk  *cblk;
    int          reqnbr = solvmtx->reqnum;
    MPI_Status   status;

#if defined(PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] Wait for all pending communications\n",
             solvmtx->clustnum );
#endif

    for ( i=0; i<reqnbr; i++ ) {
        if ( solvmtx->reqtab[i] == MPI_REQUEST_NULL ) {
            assert( 0 /* MPI_REQUEST_NULL should have been pushed to the end */ );
            solvmtx->reqnum--;
            continue;
        }

        /* Make sure that we don't have an already cleaned request in dynamic */
        assert( solvmtx->reqidx[i] != -1 );

        rc = MPI_Wait( solvmtx->reqtab + i, &status );
        assert( rc == MPI_SUCCESS );

        cblk = solvmtx->cblktab + solvmtx->reqidx[i];

        /* We should wait only for fanin */
        assert( cblk->cblktype & CBLK_FANIN );

        cpucblk_zrequest_rhs_fwd_handle_send( enums, solvmtx, rhsb, cblk );

        solvmtx->reqnum--;
    }
    assert( solvmtx->reqnum == 0 );
    (void)rc;
#else
    (void)enums;
    (void)solvmtx;
    (void)rhsb;
#endif
}
