/**
 *
 * @file cpucblk_zmpi_coeftab.c
 *
 * Precision dependent routines to send and receive cblks coeftab.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.1
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Nolan Bredel
 * @date 2023-09-20
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common/common.h"
#include "blend/solver.h"
#include "blend/solver_comm_matrix.h"
#include "kernels.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"
#include "cpucblk_zpack.h"
#include <lapacke.h>

#if defined( PASTIX_WITH_MPI )
/**
 *******************************************************************************
 *
 * @brief Asynchronously send a cblk to cblk->ownerid
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] cblk
 *          The column block that will be sent.
 *
 *******************************************************************************/
void
cpucblk_zisend( pastix_coefside_t side,
                SolverMatrix     *solvmtx,
                SolverCblk       *cblk )
{
    MPI_Request request;
    size_t      bufsize;
    void       *buffer;
    int         rc;

    assert( cblk->cblktype & CBLK_FANIN );

    bufsize = cpucblk_zcompute_size( side, cblk );
    buffer  = cpucblk_zpack( side, cblk, bufsize );

    /* Backup buffer pointer when packed */
    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        cpucblk_zfree( side, cblk );
        if ( side != PastixUCoef ) {
            cblk->lcoeftab = buffer;
        }
        else {
            cblk->ucoeftab = buffer;
        }
    }

#if defined(PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] Post Isend for cblk %ld toward %2d ( %ld Bytes )\n",
             solvmtx->clustnum, (long)cblk->gcblknum, cblk->ownerid, bufsize );
#endif

    rc = MPI_Isend( buffer, bufsize, MPI_CHAR,
                    cblk->ownerid, cblk->gcblknum, solvmtx->solv_comm, &request );
    assert( rc == MPI_SUCCESS );

    solverCommMatrixAdd( solvmtx, cblk->ownerid, bufsize );

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
 * @param[in] side
 *          Define which side of the cblk are concerned.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] solvmtx
 *          The solver matrix structure.
 *
 * @param[inout] cblk
 *          The cblk concerned by the computation.
 *
 *******************************************************************************/
static inline void
cpucblk_zrequest_handle_fanin( pastix_coefside_t   side,
                               const SolverMatrix *solvmtx,
                               SolverCblk         *cblk )
{
    assert( cblk->cblktype & CBLK_FANIN );

#if defined(PASTIX_DEBUG_MPI)
    {
        size_t cblksize = ( side == PastixLUCoef ) ? 2 : 1;
        cblksize = cblksize * cblk_colnbr( cblk ) * cblk->stride * sizeof(pastix_complex64_t);

        fprintf( stderr, "[%2d] Isend for cblk %ld toward %2d ( %ld Bytes ) (DONE)\n",
                 solvmtx->clustnum, (long)cblk->gcblknum, cblk->ownerid, (long)cblksize );
    }
#endif

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        /* Free the packed buffer */
        if ( side != PastixUCoef ) {
            memFree_null( cblk->lcoeftab );
        }
        else {
            memFree_null( cblk->ucoeftab );
        }
    }
    else {
        /* Free the cblk */
        cpucblk_zfree( side, cblk );
    }

    (void)solvmtx;
}

/**
 *******************************************************************************
 *
 * @brief Handle a finished request on a recv cblk.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk are concerned.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[inout] threadid
 *          TODO
 *
 * @param[inout] status
 *          TODO
 *
 * @param[inout] recvbuf
 *          TODO
 *
 *******************************************************************************/
static inline void
cpucblk_zrequest_handle_recv( pastix_coefside_t   side,
                              SolverMatrix       *solvmtx,
                              int                 threadid,
                              const MPI_Status   *status,
                              char               *recvbuf )
{
    SolverCblk *cblk, *fcbk;
    int src = status->MPI_SOURCE;
    int tag = status->MPI_TAG;

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

#if defined(PASTIX_DEBUG_MPI)
    {
        int          rc;
        int          count = 0;
        pastix_int_t size  = (cblk_colnbr(cblk) * cblk->stride) * sizeof(pastix_complex64_t);

        if ( side != PastixLCoef ) {
            size *= 2;
        }

        rc = MPI_Get_count( status, MPI_CHAR, &count );
        assert( rc == MPI_SUCCESS );
        assert( (cblk->cblktype & CBLK_COMPRESSED) ||
                (!(cblk->cblktype & CBLK_COMPRESSED) && (count == size)) );

        /* We can't know the sender easily, so we don't print it */
        fprintf( stderr, "[%2d] Irecv of size %d/%ld for cblk %ld (DONE)\n",
                 solvmtx->clustnum, count, (long)size, (long)cblk->gcblknum );
        (void)rc;
    }
#endif

    /* Initialize the cblk with the reception buffer */
    cblk->threadid = (fcbk->threadid == -1) ? threadid : fcbk->threadid;
    cpucblk_zunpack( side, cblk, recvbuf );

    fcbk = solvmtx->cblktab + cblk->fblokptr->fcblknm;
    cpucblk_zadd( PastixLCoef, 1., cblk, fcbk, &solvmtx->lowrank );

    /* If side is LU, let's add the U part too */
    if ( side != PastixLCoef ) {
        cpucblk_zadd( PastixUCoef, 1., cblk, fcbk, &solvmtx->lowrank );
    }

    /* Receptions cblks contribute to themselves */
    cpucblk_zrelease_deps( side, solvmtx, cblk, fcbk );

    /* Free the CBLK_RECV */
    cpucblk_zfree( side, cblk );
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
 * @param[in] side
 *          Define which side of the cblk are concerned.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] threadid
 *          Id of the thread calling this method.
 *
 * @param[in] outcount
 *          Amount of finshed requests
 *
 * @param[in] indexes
 *          Array of completed requests
 *
 * @param[in] statuses
 *          Array of statuses for the completed requests
 *
 *******************************************************************************
 *
 * @retval TODO
 *
 *******************************************************************************/
static inline int
cpucblk_zrequest_handle( pastix_coefside_t  side,
                         SolverMatrix      *solvmtx,
                         int                threadid,
                         int                outcount,
                         const int         *indexes,
                         const MPI_Status  *statuses )
{
    pastix_int_t i, reqid;
    int          nbrequest = outcount;

    for( i = 0; i < outcount; i++ ){
        reqid = indexes[i];

        /*
         * Handle the reception
         */
        if ( solvmtx->reqidx[reqid] == -1 ) {
            /* We're on a cblk recv, copy datas and restart communications */
            char       *recvbuf;
            MPI_Status status;
            int        size;

            memcpy( &status, statuses + i, sizeof(MPI_Status) );
            MPI_Get_count( &status, MPI_CHAR, &size );

            MALLOC_INTERN( recvbuf, size, char );
            memcpy( recvbuf, solvmtx->rcoeftab, size );

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

            cpucblk_zrequest_handle_recv( side, solvmtx, threadid,
                                          &status, recvbuf );
        }
        /*
         * Handle the emission
         */
        else {
            SolverCblk *cblk = solvmtx->cblktab + solvmtx->reqidx[ reqid ];
            assert( cblk->cblktype & CBLK_FANIN );

            cpucblk_zrequest_handle_fanin( side, solvmtx, cblk );

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
 * @brief Update Request array ands Request indexes in a contiguous way.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix structure with the updated arrays.
 *
 *******************************************************************************/
void
cpucblk_zupdate_reqtab( SolverMatrix *solvmtx )
{
    /* Pointer to the compressed array of request */
    MPI_Request  *outrequest = solvmtx->reqtab;
    pastix_int_t *outreqloc  = solvmtx->reqidx;
    int           outreqnbr  = 0;

    /* Pointer to the input array of request */
    MPI_Request  *inrequest = solvmtx->reqtab;
    pastix_int_t *inreqloc  = solvmtx->reqidx;
    int           inreqnbr  = 0;

    /* Look for the first completed request */
    while( (outreqnbr < solvmtx->reqnum) &&
           (*outrequest != MPI_REQUEST_NULL) )
    {
        outrequest++;
        outreqnbr++;
        outreqloc++;
    }

    inrequest = outrequest;
    inreqloc  = outreqloc;
    inreqnbr  = outreqnbr;
    for( ; inreqnbr < solvmtx->reqnum;
         inreqnbr++, inrequest++, inreqloc++ )
    {
        if ( *inrequest == MPI_REQUEST_NULL )
        {
            continue;
        }

        /* Pack the uncompleted request */
        *outrequest = *inrequest;
        *outreqloc  = *inreqloc;

        /* Move to the next one */
        outrequest++;
        outreqloc++;
        outreqnbr++;
    }

#if !defined(NDEBUG)
    /* Set to -1 remaining of the array */
    memset( outreqloc, 0xff, (solvmtx->reqnbr - outreqnbr) * sizeof(pastix_int_t) );
#endif

#if defined(PASTIX_DEBUG_MPI)
    int  i;
    for( i = outreqnbr; i < solvmtx->reqnbr; i++ )
    {
        solvmtx->reqtab[i] = MPI_REQUEST_NULL;
    }
#endif
    assert( outreqnbr < solvmtx->reqnum );
    solvmtx->reqnum = outreqnbr;
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
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] threadid
 *          Id of the thread calling this method.
 *
 *******************************************************************************/
void
cpucblk_zmpi_progress( pastix_coefside_t   side,
                       SolverMatrix       *solvmtx,
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
            nbfree = cpucblk_zrequest_handle( side, solvmtx, threadid,
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
 * @param[in] side
 *          Define which side of the cblk must be released.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[inout] cblk
 *          The column block that contribute to fcblk.
 *
 *******************************************************************************
 *
 * @return 1 if the cblk is a fanin, 0 otherwise
 *
 *******************************************************************************/
int
cpucblk_zincoming_deps( int                rank,
                        pastix_coefside_t  side,
                        SolverMatrix      *solvmtx,
                        SolverCblk        *cblk )
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
        cpucblk_zmpi_progress( side, solvmtx, rank );
    }
#else
    assert( !(cblk->cblktype & (CBLK_FANIN | CBLK_RECV)) );
    do { pastix_yield(); } while( cblk->ctrbcnt > 0 );
#endif

    (void)rank;
    (void)side;
    (void)solvmtx;

    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Release the dependencies of the given cblk after an update.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be released.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] cblk
 *          The column block that contribute to fcblk.
 *
 * @param[inout] fcbk
 *          The facing column block that is updated by cblk.
 *
 *******************************************************************************/
void
cpucblk_zrelease_deps( pastix_coefside_t  side,
                       SolverMatrix      *solvmtx,
                       const SolverCblk  *cblk,
                       SolverCblk        *fcbk )
{
    int32_t ctrbcnt;
    ctrbcnt = pastix_atomic_dec_32b( &(fcbk->ctrbcnt) );
    if ( !ctrbcnt ) {
#if defined(PASTIX_WITH_MPI)
        if ( fcbk->cblktype & CBLK_FANIN ) {
            cpucblk_zisend( side, solvmtx, fcbk );
            return;
        }
#else
        (void)side;
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
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
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
 *******************************************************************************/
void
cpucblk_zrequest_cleanup( pastix_coefside_t side,
                          pastix_int_t      sched,
                          SolverMatrix     *solvmtx )
{
    if ( (sched != PastixSchedSequential) &&
         (sched != PastixSchedStatic)     &&
         (sched != PastixSchedDynamic) )
    {
        return;
    }
#if defined(PASTIX_WITH_MPI)
    pastix_int_t i;
    int rc;
    SolverCblk  *cblk;
    int          reqnbr = solvmtx->reqnum;
    MPI_Status   status;

#if defined(PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] Wait for all pending communications\n",
             solvmtx->clustnum );
#endif

    for( i=0; i<reqnbr; i++ )
    {
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

        cpucblk_zrequest_handle_fanin( side, solvmtx, cblk );

        solvmtx->reqnum--;
    }
    assert( solvmtx->reqnum == 0 );
    (void)rc;
#else
    (void)side;
    (void)solvmtx;
#endif
}
