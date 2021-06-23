/**
 *
 * @file cpucblk_zmpi_coeftab.c
 *
 * Precision dependent routines to send and receive cblks coeftab.
 *
 * @copyright 2015-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.0
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2021-04-07
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common/common.h"
#include "blend/solver.h"
#include "kernels.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"
#include <lapacke.h>
#if defined( PASTIX_WITH_MPI )

/**
 *******************************************************************************
 *
 * @brief Compute the size of a block to send in LR
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] cblk
 *          The column block that will be sent.
 *
 * @param[in] blok
 *          A block that will be sent.
 *
 *******************************************************************************
 *
 * @return Size of a blok to send in LR
 *
 *******************************************************************************/
pastix_uint_t
cpublok_zcompute_size_lr( pastix_coefside_t side, const SolverCblk *cblk, const SolverBlok *blok )
{
    assert( cblk->cblktype & CBLK_COMPRESSED );

    pastix_uint_t M = blok_rownbr( blok );
    pastix_uint_t N = cblk_colnbr( cblk );

    pastix_uint_t suv  = 0;
    pastix_uint_t coef = 0;
    /* Add lower part size */
    if ( side != PastixUCoef ) {
        if ( blok->LRblock[0].rk != -1 ) {
            suv += blok->LRblock[0].rk * ( M + N );
        }
        else {
            suv += M * N;
        }
        coef++;
    }

    /* Add upper part size */
    if ( side != PastixLCoef ) {
        if ( blok->LRblock[1].rk != -1 ) {
            suv += blok->LRblock[1].rk * ( M + N );
        }
        else {
            suv += M * N;
        }
        coef++;
    }

    /* size of rk(int) + sizeof(u+v) + size of u + size of v */
    return ( sizeof( int ) + sizeof( pastix_uint_t ) ) * coef + suv * sizeof( pastix_complex64_t );
}

/**
 *******************************************************************************
 *
 * @brief Compute the size of a cblk to send in LR
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] cblk
 *          The column block that will be sent.
 *
 *******************************************************************************
 *
 * @return Size of a cblk to send in LR
 *
 *******************************************************************************/

pastix_uint_t
cpucblk_zcompute_size_lr( pastix_coefside_t side, const SolverCblk *cblk )
{
    assert( cblk->cblktype & CBLK_COMPRESSED );

    pastix_uint_t     size  = 0;
    const SolverBlok *blok  = cblk->fblokptr;
    const SolverBlok *lblok = cblk[1].fblokptr;
    for ( ; blok < lblok; blok++ ) {
        size += cpublok_zcompute_size_lr( side, cblk, blok );
    }

    /* size of all bloks */
    return size;
}

/**
 *******************************************************************************
 *
 * @brief Compute the size of the buffer to send.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] cblk
 *          The column block that will be sent.
 *
 *******************************************************************************
 *
 * @return Size of the buffer to send.
 *
 *******************************************************************************/
pastix_uint_t
cpucblk_zcompute_size( pastix_coefside_t side, const SolverCblk *cblk )
{
    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        return cpucblk_zcompute_size_lr( side, cblk );
    }
    else {
        pastix_int_t cblksize = cblk->stride * cblk_colnbr( cblk );
        if ( side == PastixLUCoef ) {
            cblksize *= 2;
        }
        return cblksize;
    }
}

/**
 *******************************************************************************
 *
 * @brief Pack low-rank data by side
 *
 *******************************************************************************
 *
 * @param[in] shift
 *          Define which side of the blok must be tested.
 *          @arg 0 if the side is PastixLCoef
 *          @arg 1 if the side is PastixUCoef
 *
 * @param[in] blok
 *          Block that will be sent.
 *
 * @param[in] M
 *          Number of rows of the matrix A.
 *
 * @param[in] N
 *          Number of columns of the matrix A.
 *
 * @param[inout] buffer
 *          Pointer on packed data
 *
 *******************************************************************************
 *
 * @return Pointer to the data buffer.
 *
 *******************************************************************************/
void *
cpublok_zpack_lr_by_side( pastix_int_t      shift,
                          const SolverBlok *blok,
                          pastix_uint_t     M,
                          pastix_uint_t     N,
                          void *            buffer )
{
    int   rk    = blok->LRblock[shift].rk;
    int   rkmax = blok->LRblock[shift].rkmax;
    void *u     = blok->LRblock[shift].u;
    void *v     = blok->LRblock[shift].v;

    buffer = (char *)buffer;

    memcpy( buffer, &rk, sizeof( int ) );
    buffer += sizeof( int );

    pastix_uint_t suv;
    if ( rk != -1 ) {
        suv = rk * ( M + N );
    }
    else {
        suv = M + N;
    }
    memcpy( buffer, &suv, sizeof( pastix_uint_t ) );
    buffer += sizeof( pastix_uint_t );

    if ( rk != -1 ) {
        memcpy( buffer, u, rk * M * sizeof( pastix_complex64_t ) );
        buffer += rk * M * sizeof( pastix_complex64_t );
        if ( rk == rkmax ) {
            memcpy( buffer, v, rk * N * sizeof( pastix_complex64_t ) );
            buffer += rk * N * sizeof( pastix_complex64_t );
        }
        else {
            LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', rk, N, v, rkmax, buffer, rk );
            buffer += rk * N * sizeof( pastix_complex64_t );
        }
    }
    else {
        memcpy( buffer, u, N * M * sizeof( pastix_complex64_t ) );
        buffer += M * N * sizeof( pastix_complex64_t );
    }
    return (void *)buffer;
}

/**
 *******************************************************************************
 *
 * @brief Pack low-rank data for a block
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[in] blok
 *          Block that will be sent.
 *
 * @param[in] M
 *          Number of rows of the matrix A.
 *
 * @param[in] N
 *          Number of columns of the matrix A.
 *
 * @param[inout] buffer
 *          Pointer on packed data
 *
 *******************************************************************************
 *
 * @return Pointer to the data buffer.
 *
 *******************************************************************************/
void *
cpublok_zpack_lr( pastix_coefside_t side,
                  const SolverBlok *blok,
                  pastix_uint_t     M,
                  pastix_uint_t     N,
                  void *            buffer )
{
    if ( side == PastixLCoef ) {
        buffer = cpublok_zpack_lr_by_side( 0, blok, M, N, buffer );
    }
    else if ( side == PastixUCoef ) {
        buffer = cpublok_zpack_lr_by_side( 1, blok, M, N, buffer );
    }
    else {
        buffer = cpublok_zpack_lr_by_side( 0, blok, M, N, buffer );
        buffer = cpublok_zpack_lr_by_side( 1, blok, M, N, buffer );
    }

    return buffer;
}

/**
 *******************************************************************************
 *
 * @brief Pack low-rank data
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] cblk
 *          The column block that will be sent.
 *
 * @param[in] size
 *          Size of packed data
 *
 *******************************************************************************
 *
 * @return Pointer to the data buffer.
 *
 *******************************************************************************/
void *
cpucblk_zpack_lr( pastix_coefside_t side, const SolverCblk *cblk, pastix_uint_t size )
{
    assert( cblk->cblktype & CBLK_COMPRESSED );

    char *buffer = malloc( size );
    char *tmp    = buffer;

    const SolverBlok *blok  = cblk->fblokptr;
    const SolverBlok *lblok = cblk[1].fblokptr;

    pastix_int_t M = 0;
    pastix_int_t N = 0;

    for ( ; blok < lblok; blok++ ) {
        M   = blok_rownbr( blok );
        N   = cblk_colnbr( cblk );
        tmp = cpublok_zpack_lr( side, blok, M, N, tmp );
    }

    return buffer;
}

/**
 *******************************************************************************
 *
 * @brief Unpack low rank data and fill the cblk concerned by the computation
 *
 *******************************************************************************
 *
 * @param[in] shift
 *          Define which side of the blok must be tested.
 *          @arg 0 if the side is PastixLCoef
 *          @arg 1 if the side is PastixUCoef
 *
 * @param[inout] blok
 *          The blok concerned by the computation.
 *
 * @param[in] M
 *          Number of rows of the matrix.
 *
 * @param[in] N
 *          Number of columns of the matrix.
 *
 * @param[inout] buffer
 *          Pointer on packed data
 *
 *******************************************************************************
 *
 * @return Pointer to the data buffer.
 *
 *******************************************************************************/
void *
cpublok_zunpack_lr_by_side( pastix_int_t  shift,
                            SolverBlok *  blok,
                            pastix_uint_t M,
                            pastix_uint_t N,
                            void *        buffer )
{
    int rk;
    memcpy( &rk, buffer, sizeof( int ) );
    buffer += sizeof( int );
    pastix_uint_t suv;
    memcpy( &suv, buffer, sizeof( pastix_uint_t ) );
    buffer += sizeof( pastix_int_t );

    blok->LRblock[shift].rk    = rk;
    blok->LRblock[shift].rkmax = rk;

    if ( rk != -1 ) {
        memcpy( blok->LRblock[shift].u, buffer, M * rk * sizeof( pastix_complex64_t ) );
        buffer += M * rk * sizeof( pastix_complex64_t );
        memcpy( blok->LRblock[shift].v, buffer, suv - M * rk * sizeof( pastix_complex64_t ) );
        buffer += suv - M * rk * sizeof( pastix_complex64_t );
    }
    else {
        memcpy( blok->LRblock[shift].u, buffer, M * N * sizeof( pastix_complex64_t ) );
        buffer += M * N * sizeof( pastix_complex64_t );
    }
    return buffer;
}

/**
 *******************************************************************************
 *
 * @brief Unpack low rank data and fill the cblk concerned by the computation
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix structure.
 *
 * @param[inout] blok
 *          The blok concerned by the computation.
 *
 * @param[in] M
 *          Number of rows of the matrix.
 *
 * @param[in] N
 *          Number of columns of the matrix.
 *
 * @param[inout] buffer
 *          Pointer on packed data
 *
 *******************************************************************************
 *
 * @return Pointer to the data buffer.
 *
 *******************************************************************************/
void *
cpublok_zunpack_lr( pastix_coefside_t side,
                    SolverBlok *      blok,
                    pastix_uint_t     M,
                    pastix_uint_t     N,
                    void *            buffer )
{
    if ( side == PastixLCoef ) {
        buffer = cpublok_zunpack_lr_by_side( 0, blok, M, N, buffer );
    }
    else if ( side == PastixUCoef ) {
        buffer = cpublok_zunpack_lr_by_side( 1, blok, M, N, buffer );
    }
    else {
        buffer = cpublok_zunpack_lr_by_side( 0, blok, M, N, buffer );
        buffer = cpublok_zunpack_lr_by_side( 1, blok, M, N, buffer );
    }
    return buffer;
}

/**
 *******************************************************************************
 *
 * @brief Unpack low rank data and fill the cblk concerned by the computation
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] cblk
 *          The cblk concerned by the computation.
 *
 * @param[inout] buffer
 *          Pointer on packed data
 *
 *******************************************************************************
 *
 * @return Pointer to the data buffer.
 *
 *******************************************************************************/
void
cpucblk_zunpack_lr( pastix_coefside_t side, SolverCblk *cblk, void *buffer )
{
    assert( cblk->cblktype & CBLK_COMPRESSED );

    SolverBlok *blok  = cblk->fblokptr;
    SolverBlok *lblok = cblk[1].fblokptr;

    for ( ; blok < lblok; blok++ ) {
        pastix_uint_t M = blok_rownbr( blok );
        pastix_uint_t N = cblk_colnbr( cblk );
        buffer          = cpublok_zunpack_lr( side, blok, M, N, buffer );
    }
}

/**
 *******************************************************************************
 *
 * @brief Pack data in full rank
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
 *******************************************************************************
 *
 * @return Pointer to the data buffer.
 *
 *******************************************************************************/
void *
cpucblk_zpack_fr( pastix_coefside_t side, const SolverCblk *cblk )
{
    assert( !( cblk->cblktype & CBLK_COMPRESSED ) );

    return side == PastixUCoef ? cblk->ucoeftab : cblk->lcoeftab;
}

/**
 *******************************************************************************
 *
 * @brief Unpack data in full rank
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
cpucblk_zunpack_fr( pastix_coefside_t side, SolverMatrix *solvmtx, SolverCblk *cblk )
{
    assert( !( cblk->cblktype & CBLK_COMPRESSED ) );

    cblk->lcoeftab = solvmtx->rcoeftab;
    if ( side != PastixLCoef ) {
        pastix_complex64_t *recv = cblk->lcoeftab;
        cblk->ucoeftab           = recv + ( cblk_colnbr( cblk ) * cblk->stride );
    }
}

/**
 *******************************************************************************
 *
 * @brief Pack data (Full rank or low rank)
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] cblk
 *          The column block that will be sent.it a
 *
 * @param[in] size
 *          Size of the packed data
 *
 *******************************************************************************
 *
 * @return Pointer to the data buffer
 *
 *******************************************************************************/
void *
cpucblk_zpack( pastix_coefside_t side, const SolverCblk *cblk, const pastix_uint_t size )
{
    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        return cpucblk_zpack_lr( side, cblk, size );
    }
    else {
        return cpucblk_zpack_fr( side, cblk );
    }
}

/**
 *******************************************************************************
 *
 * @brief Unpack data and fill the cblk concerned by the computation
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
 *          The cblk concerned by the computation.
 *
 *******************************************************************************/
void
cpucblk_zunpack( pastix_coefside_t side, SolverMatrix *solvmtx, SolverCblk *cblk )
{
    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        cpucblk_zunpack_lr( side, cblk, solvmtx->rcoeftab );
    }
    else {
        cpucblk_zunpack_fr( side, solvmtx, cblk );
    }
}

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
                const SolverCblk *cblk )
{
    MPI_Request  request;
    int rc;

    assert( cblk->cblktype & CBLK_FANIN );


#if defined(PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] Post Isend for cblk %ld toward %2d ( %ld Bytes )\n",
             solvmtx->clustnum, (long)cblk->gcblknum, cblk->ownerid,
             (long)(cblksize * sizeof(pastix_complex64_t)) );
#endif

    pastix_uint_t bufsize = cpucblk_zcompute_size( side, cblk );
    void *buffer = cpucblk_zpack( side, cblk, bufsize );

    rc = MPI_Isend( buffer, bufsize, PASTIX_MPI_COMPLEX64,
                        cblk->ownerid, cblk->gcblknum, solvmtx->solv_comm, &request );

    if (cblk->cblktype & CBLK_COMPRESSED) {
        free(buffer);
    }

    assert( rc == MPI_SUCCESS );

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
                              int threadid, const MPI_Status *status )
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
    cblk = fcbk--;

    /* Get through source */
    while( cblk->ownerid != src ) {
        cblk--;
        assert( cblk >= solvmtx->cblktab );
        assert( cblk->gcblknum == tag );
        assert( cblk->cblktype & CBLK_RECV );
    }

#if defined(PASTIX_DEBUG_MPI)
    {
        pastix_int_t size = (cblk_colnbr(cblk) * cblk->stride);
        int count = 0;

        if ( side != PastixLCoef ) {
            size *= 2;
        }

        MPI_Get_count( status, PASTIX_MPI_COMPLEX64, &count );
        assert( count == size );

        /* We can't know the sender easily, so we don't print it */
        fprintf( stderr, "[%2d] Irecv of size %d/%ld for cblk %ld (DONE)\n",
                 solvmtx->clustnum, count, (long)size, (long)cblk->gcblknum );
    }
#endif

    /* Initialize the cblk with the reception buffer */
    cblk->threadid = (fcbk->threadid == -1) ? threadid : fcbk->threadid;

    cpucblk_zunpack( side, solvmtx, cblk );

    fcbk = solvmtx->cblktab + cblk->fblokptr->fcblknm;

    cpucblk_zadd( PastixLCoef, 1., cblk, fcbk, NULL );

    /* If side is LU, let's add the U part too */
    if ( side != PastixLCoef ) {
        cpucblk_zadd( PastixUCoef, 1., cblk, fcbk, NULL );
    }

    /* Receptions cblks contribute to themselves */
    cpucblk_zrelease_deps( side, solvmtx, cblk, fcbk );
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
            cpucblk_zrequest_handle_recv( side, solvmtx,
                                          threadid, statuses + i );
            solvmtx->recvcnt--;

            /* Let's restart the communication */
            if ( solvmtx->recvcnt > 0 ) {
                MPI_Start( solvmtx->reqtab + reqid );
                nbrequest--;
            }
            else {
                MPI_Request_free( solvmtx->reqtab + reqid );
                solvmtx->reqtab[reqid] = MPI_REQUEST_NULL;
            }
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
static inline void
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

    pastix_comm_tid = -1;
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
    do { } while( cblk->ctrbcnt > 0 );
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
    int          reqnbr =  solvmtx->reqnum;
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
