/**
 *
 * @file cpucblk_zmpi_coeftab.c
 *
 * Precision dependent routines to send and receive cblks coeftab.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2019-04-11
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "kernels.h"
#include "solver.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

#if defined(PASTIX_WITH_MPI)
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
 * @param[inout] cblk
 *          The column block that will be sent.
 *
 *******************************************************************************/
void
cpucblk_zisend( pastix_coefside_t side,
                SolverMatrix     *solvmtx,
                const SolverCblk *cblk )
{
    pastix_int_t cblksize = cblk->stride * cblk_colnbr(cblk);
    MPI_Request  request;
    int rc;

    assert( !(cblk->cblktype & CBLK_COMPRESSED) );
    assert(   cblk->cblktype & CBLK_FANIN       );

    if ( side == PastixLUCoef ) {
        cblksize *= 2;
    }

#if defined(PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] Post Isend for cblk %ld toward %2d\n",
             solvmtx->clustnum, (long)cblk->gcblknum, cblk->ownerid );
#endif

    if ( side == PastixUCoef ) {
        rc = MPI_Isend( cblk->ucoeftab, cblksize, PASTIX_MPI_COMPLEX64,
                        cblk->ownerid, cblk->gcblknum, solvmtx->solv_comm, &request );
    }
    else {
        rc = MPI_Isend( cblk->lcoeftab, cblksize, PASTIX_MPI_COMPLEX64,
                        cblk->ownerid, cblk->gcblknum, solvmtx->solv_comm, &request );
    }
    assert( rc == MPI_SUCCESS );

    /* Register the request to make it progress */
    pastix_atomic_lock( &pastix_mpi_lock );
    if ( solvmtx->computeQueue ){
        solvmtx->reqtab[ solvmtx->reqnum ] = request;
        solvmtx->reqlocal[ solvmtx->reqnum ] = cblk - solvmtx->cblktab;
        solvmtx->reqnum++;
    }
    else{
        solvmtx->reqtab[ cblk->reqindex ] = request;
    }
    pastix_atomic_unlock( &pastix_mpi_lock );

    (void)rc;
}

/**
 *******************************************************************************
 *
 * @brief  Asynchronously receive a cblk from cblk->ownerid
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] solvmtx
 *          The solver matrix structure.
 *
 * @param[inout] cblk
 *          The column block that will be sent.
 *
 *******************************************************************************/
void
cpucblk_zirecv( pastix_coefside_t   side,
                const SolverMatrix *solvmtx,
                SolverCblk         *cblk )
{
    pastix_int_t cblksize = cblk->stride * cblk_colnbr(cblk);
    MPI_Request *request  = solvmtx->reqtab + cblk->reqindex;
    int rc;

    /* Compression is not handled yet */
    assert( !(cblk->cblktype & CBLK_COMPRESSED) );
    assert(   cblk->cblktype & CBLK_RECV        );

    if ( side == PastixLUCoef ) {
        cblksize *= 2;
    }

    /* Init rcoeftab only when it's called */
    if ( cblk->rcoeftab == NULL ) {
        MALLOC_INTERN( cblk->rcoeftab, cblksize, pastix_complex64_t );
    }

    assert( cblk->rcoeftab );
#if defined(PASTIX_DEBUG_MPI)
    memset( cblk->rcoeftab, 0, cblksize * sizeof(pastix_complex64_t) );

    fprintf( stderr, "[%2d] Post Irecv for cblk %ld from any source\n",
             solvmtx->clustnum, (long)cblk->gcblknum );
#endif

    rc = MPI_Irecv( cblk->rcoeftab, cblksize, PASTIX_MPI_COMPLEX64,
                    MPI_ANY_SOURCE, cblk->gcblknum, solvmtx->solv_comm, request );
    assert( rc == MPI_SUCCESS );

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
    fprintf( stderr, "[%2d] Isend for cblk %ld toward %2d (DONE)\n",
             solvmtx->clustnum, (long)cblk->gcblknum, cblk->ownerid );
#endif
    cpucblk_zfree( side, cblk );

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
 * @param[inout] cblk
 *          The cblk concerned by the computation.
 *
 *******************************************************************************/
static inline void
cpucblk_zrequest_handle_recv( pastix_coefside_t  side,
                              SolverMatrix      *solvmtx,
                              SolverCblk        *cblk )
{
    assert( cblk->cblktype & CBLK_RECV );

#if defined(PASTIX_DEBUG_MPI)
    /* We can't know the sender easily, so we don't print it */
    fprintf( stderr, "[%2d] Irecv for cblk %ld (DONE)\n",
             solvmtx->clustnum, (long)cblk->gcblknum );
#endif

    cpucblk_zadd_recv( PastixLCoef, 1., cblk );

    /* If side is LU, let's add the U part too */
    if ( side != PastixLCoef ) {
        cpucblk_zadd_recv( PastixUCoef, 1., cblk );
    }

    /* Receptions cblks contribute to themselves */
    cpucblk_zrelease_deps( side, solvmtx, cblk, cblk );

    /* We are done with this update, we can now receive another one */
    cblk->rcoeftab = NULL;
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
 * @param[in] outcount
 *          Amount of finshed requests
 *
 * @param[inout] indexes
 *          Array of completed requests
 *
 *******************************************************************************/
static inline void
cpucblk_zrequest_handle( pastix_coefside_t  side,
                         SolverMatrix      *solvmtx,
                         int                outcount,
                         int               *indexes  )
{
    pastix_int_t i, index;
    SolverCblk *cblk;

    for( i = 0; i < outcount; i++ ){
        index = indexes[i];
        assert( solvmtx->reqtab[index] == MPI_REQUEST_NULL );

        cblk = solvmtx->cblktab + solvmtx->reqlocal[index];

        if ( cblk->cblktype & CBLK_FANIN ) {
            cpucblk_zrequest_handle_fanin( side, solvmtx, cblk );
            continue;
        }

        if ( cblk->cblktype & CBLK_RECV ) {
            cpucblk_zrequest_handle_recv( side, solvmtx, cblk );

            cblk->recvcnt--;
            /* Relaunch the reception if there is still more work to do */
            if ( cblk->recvcnt ) {
                cpucblk_zirecv( side, solvmtx, cblk );
            }
            else {
                memFree_null( cblk->rcoeftab );
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Wait some active MPI requests and process the finished ones.
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
 *******************************************************************************/
void
cpucblk_zwaitsome( pastix_coefside_t  side,
                   SolverMatrix      *solvmtx )
{
    int        outcount, rc;
    int        reqnbr = solvmtx->reqnbr;
    int        indexes[reqnbr];
    MPI_Status statuses[reqnbr];

    pastix_atomic_lock( &pastix_mpi_lock );
    rc = MPI_Waitsome( reqnbr, solvmtx->reqtab, &outcount, indexes, statuses );
    pastix_atomic_unlock( &pastix_mpi_lock );
    assert( rc == MPI_SUCCESS );

    if ( outcount != MPI_UNDEFINED ) {
        cpucblk_zrequest_handle( side, solvmtx, outcount, indexes  );
    }

    (void)rc;
}

/**
 *******************************************************************************
 *
 * @brief Test some active MPI requests and process the finished ones.
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
 *******************************************************************************/
void
cpucblk_ztestsome( pastix_coefside_t  side,
                   SolverMatrix      *solvmtx )
{
    int reqnbr = solvmtx->reqnbr;
    int indexes[reqnbr];
    int outcount = 0;
    int rc;
    MPI_Status statuses[reqnbr];

    pastix_atomic_lock( &pastix_mpi_lock );
    rc = MPI_Testsome( reqnbr, solvmtx->reqtab, &outcount, indexes, statuses );
    pastix_atomic_unlock( &pastix_mpi_lock );
    assert( rc == MPI_SUCCESS );

    if ( outcount > 0 ) {
        cpucblk_zrequest_handle( side, solvmtx, outcount, indexes );
    }

    (void)rc;
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
static inline void
cpucblk_zupdate_reqtab( SolverMatrix *solvmtx )
{
    /* If there is only one reception request */
    if ( solvmtx->reqnbr == 1 ) {
        return;
    }

    pastix_atomic_lock( &pastix_mpi_lock );
    {
        /* Pointer to the compressed array of request */
        MPI_Request  *outrequest = solvmtx->reqtab;
        pastix_int_t *outreqloc  = solvmtx->reqlocal;
        int           outreqnbr  = 0;

        /* Pointer to the input array of request */
        MPI_Request  *inrequest = solvmtx->reqtab;
        pastix_int_t *inreqloc  = solvmtx->reqlocal;
        int           inreqnbr  = solvmtx->reqnum;

        /* Look for the first completed request */
        while( (outreqnbr < solvmtx->reqnbr) &&
               (*outrequest != MPI_REQUEST_NULL) )
        {
            outrequest++;
            outreqnbr++;
            outreqloc++;
        }

        inrequest = outrequest;
        inreqloc  = outreqloc;
        while( outreqnbr < inreqnbr )
        {
            /*
             * Skip all completed requests
             * until the next non completed one
             */
            while( *inrequest == MPI_REQUEST_NULL )
            {
                inrequest++;
                inreqloc++;
            }
            /* Pack the uncompleted request */
            *outrequest = *inrequest;
            *outreqloc  = *inreqloc;

            /* Move to the next one */
            outrequest++;
            outreqloc++;
            outreqnbr++;

            inrequest++;
            inreqloc++;
        }

        /* Set to -1 remaining of the array */
        memset( outreqloc, 0xff, (solvmtx->reqnbr - outreqnbr) * sizeof(pastix_int_t) );

#if defined(PASTIX_DEBUG_MPI)
        int  i;
        for( i = outreqnbr; i < solvmtx->reqnbr; i++ )
        {
            solvmtx->reqtab[i] = MPI_REQUEST_NULL;
        }
#endif
    }
    pastix_atomic_unlock( &pastix_mpi_lock );
}

/**
 *******************************************************************************
 *
 * @brief Progress communications with the dynamic scheduler.
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
 * @param[in] ctx
 *          Context to assure thread-safety.
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] threadid
 *          Id of the thread calling this method.
 *
 * @param[inout] recv
 *          Reception buffer of the persistant communication.
 *
 *******************************************************************************/
void
cpucblk_zmpi_progress( pastix_coefside_t   side,
                       isched_thread_t    *ctx,
                       SolverMatrix       *solvmtx,
                       pastix_int_t        threadid,
                       pastix_complex64_t *recv )
{
    int        ii, outcount = 0;
    int        indexes[ solvmtx->reqnbr ];
    MPI_Status statuses[ solvmtx->reqnbr ];

    /* Check if someone is already communicating or not */
    pthread_mutex_lock( &(ctx->global_ctx->commlock) );
    if ( ctx->global_ctx->commid == -1 ) {
        ctx->global_ctx->commid = threadid;
    }
    pthread_mutex_unlock( &(ctx->global_ctx->commlock) );

    if ( threadid != ctx->global_ctx->commid ) {
        return;
    }

    /*
     * We suppose that we have at least one reception if needed, or some send
     */
    assert( solvmtx->reqnum >= 0 );
    MPI_Testsome( solvmtx->reqnum, solvmtx->reqtab, &outcount, indexes, statuses );

    for( ii=0; ii<outcount; ii++ )
    {
        SolverCblk *cblk;
        int         reqid  = indexes[ii];

        if ( solvmtx->reqlocal[ reqid ] == -1 ) {
            int source = statuses[ii].MPI_SOURCE;
            int tag    = statuses[ii].MPI_TAG;

            assert( ( 0 <= tag    ) && ( tag    < solvmtx->gcblknbr ) );
            assert( ( 0 <= source ) && ( source < solvmtx->clustnbr ) );

            /* Get the associated recv cblk */
            cblk = solvmtx->cblktab + solvmtx->gcbl2loc[ tag ];
            assert( cblk != NULL );
            (void)source;
        }
        else {
            cblk = solvmtx->cblktab + solvmtx->reqlocal[ reqid ];
        }

        if ( cblk->cblktype & CBLK_FANIN ) {
            cpucblk_zrequest_handle_fanin( side, solvmtx, cblk );

            solvmtx->reqnum--;
            solvmtx->fanincnt--;
            continue;
        }

        if ( cblk->cblktype & CBLK_RECV ) {
            cblk->threadid = threadid;
            assert( cblk->rcoeftab == NULL );
            cblk->rcoeftab = recv;

            cpucblk_zrequest_handle_recv( side, solvmtx, cblk );
            solvmtx->recvcnt--;

            /* Let's restart the communication */
            if ( solvmtx->recvcnt > 0 ) {
                MPI_Start( solvmtx->reqtab + reqid );
            }
            else {
                MPI_Request_free( solvmtx->reqtab + reqid );
                solvmtx->reqtab[reqid] = MPI_REQUEST_NULL;
                solvmtx->reqnum--;
            }
            continue;
        }
        assert(0);
    }

    if ( outcount > 0 ) {
        cpucblk_zupdate_reqtab( solvmtx );
    }
    ctx->global_ctx->commid = -1;
}
#endif /* defined(PASTIX_WITH_MPI) */

/**
 *******************************************************************************
 *
 * @brief Wait for incoming dependencies, and return when cblk->ctrbcnt has reached 0.
 *
 *******************************************************************************
 *
 * @param[in] mt_flag
 *          @arg 0, the function is called in a sequential environment, and we
 *                  can wait on each communication.
 *          @arg 1, the function is called in a multi-threaded environment, and
 *                  we need to test the communication to avoid dead locks.
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
 * @return 1 if the cblk is a fanin, 0 otherwise
 *
 *******************************************************************************/
int
cpucblk_zincoming_deps( int                mt_flag,
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
        cpucblk_ztestsome( side, solvmtx );
        return 1;
    }

    if ( cblk->cblktype & CBLK_RECV ) {
        cpucblk_zirecv( side, solvmtx, cblk );
    }

    /* Make sure we receive every contribution */
    while( cblk->ctrbcnt > 0 ) {
        if ( mt_flag ) {
            /*
             * Change for test to make sure we do not dead lock in a wait on
             * something received by another thread
             */
            cpucblk_ztestsome( side, solvmtx );
        }
        else {
            cpucblk_zwaitsome( side, solvmtx );
        }
    }
#else
    assert( !(cblk->cblktype & (CBLK_FANIN | CBLK_RECV)) );
    do { } while( cblk->ctrbcnt > 0 );
#endif

    (void)mt_flag;
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
#endif
        if ( solvmtx->computeQueue ) {
            pastix_queue_t *queue = solvmtx->computeQueue[ cblk->threadid ];
            pqueuePush1( queue, fcbk - solvmtx->cblktab, queue->size );
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
    int          reqnbr = (sched == PastixSchedDynamic) ?
                           solvmtx->reqnum : solvmtx->reqnbr;
    MPI_Status   status;

#if defined(PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] Wait for all pending communications\n",
             solvmtx->clustnum );
#endif

    for( i=0; i<reqnbr; i++ )
    {
        if ( solvmtx->reqtab[i] == MPI_REQUEST_NULL ) {
            continue;
        }

        /* Make sure that we don't have an already cleaned request in dynamic */
        assert( solvmtx->reqlocal[i] != -1 );

        rc = MPI_Wait( solvmtx->reqtab + i, &status );
        assert( rc == MPI_SUCCESS );

        cblk = solvmtx->cblktab + solvmtx->reqlocal[i];

        /* We should wait only for fanin */
        assert( cblk->cblktype & CBLK_FANIN );

        cpucblk_zrequest_handle_fanin( side, solvmtx, cblk );

        solvmtx->reqnum = (sched == PastixSchedDynamic) ?
                           solvmtx->reqnum - 1 : solvmtx->reqnum;
    }
#else
    (void)side;
    (void)solvmtx;
#endif
}
