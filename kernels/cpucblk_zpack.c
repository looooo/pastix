/**
 *
 * @file cpucblk_zpack.c
 *
 * Precision dependent routines to pack and unpack cblks.
 *
 * @copyright 2021-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Nolan Bredel
 * @date 2023-07-21
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common/common.h"
#include "blend/solver.h"
#include "cpucblk_zpack.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

/**
 *******************************************************************************
 *
 * @brief Compute the size of a block to send in LR.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] N
 *          The number of columns of the block.
 *
 * @param[in] blok
 *          The block for which the size is computed.
 *
 *******************************************************************************
 *
 * @return size of the LR block to send in bytes
 *
 *******************************************************************************/
size_t
cpublok_zcompute_size_lr( pastix_coefside_t  side,
                          pastix_int_t       N,
                          const SolverBlok  *blok )
{
    pastix_int_t M    = blok_rownbr( blok );
    pastix_int_t suv  = 0;
    pastix_int_t coef = 0;

    /* Add lower part size */
    if ( side != PastixUCoef ) {
        suv += core_zlrgetsize( M, N, blok->LRblock[0] );
        coef++;
    }

    /* Add upper part size */
    if ( side != PastixLCoef ) {
        suv += core_zlrgetsize( M, N, blok->LRblock[1] );
        coef++;
    }

    /* size of rk(int) + size of u + size of v */
    return coef * sizeof( int ) + suv * sizeof( pastix_complex64_t );
}

/**
 *******************************************************************************
 *
 * @brief Compute the size of a column block to send in LR.
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
 *          The low-rank column block for which the size is computed.
 *
 *******************************************************************************
 *
 * @return size of the LR column block to send in bytes.
 *
 *******************************************************************************/
pastix_uint_t
cpucblk_zcompute_size_lr( pastix_coefside_t  side,
                          const SolverCblk  *cblk )
{
    assert( cblk->cblktype & CBLK_COMPRESSED );

    pastix_uint_t     N     = cblk_colnbr( cblk );
    pastix_uint_t     size  = 0;
    const SolverBlok *blok  = cblk->fblokptr;
    const SolverBlok *lblok = cblk[1].fblokptr;

    for ( ; blok < lblok; blok++ ) {
        size += cpublok_zcompute_size_lr( side, N, blok );
    }

    /* size of all blocks */
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
 *          The column block for which the size is computed.
 *
 *******************************************************************************
 *
 * @return Size of the buffer to send.
 *
 *******************************************************************************/
size_t
cpucblk_zcompute_size( pastix_coefside_t  side,
                       const SolverCblk  *cblk )
{
    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        return cpucblk_zcompute_size_lr( side, cblk );
    }
    else {
        size_t cblksize = cblk->stride * cblk_colnbr( cblk );
        if ( side == PastixLUCoef ) {
            cblksize *= 2;
        }
        return cblksize * sizeof( pastix_complex64_t );
    }
}

/**
 *******************************************************************************
 *
 * @brief Pack low-rank data for a block.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] N
 *          Number of columns of the block.
 *
 * @param[in] blok
 *          The solver block to pack.
 *
 * @param[inout] buffer
 *          Pointer to the buffer where to pack the data.
 *
 *******************************************************************************
 *
 * @return Pointer to the end of the packed data.
 *
 *******************************************************************************/
char *
cpublok_zpack_lr( pastix_coefside_t  side,
                  pastix_uint_t      N,
                  const SolverBlok  *blok,
                  char              *buffer )
{
    pastix_int_t M = blok_rownbr( blok );

    if ( side != PastixUCoef ) {
        buffer = core_zlrpack( M, N, blok->LRblock[0], buffer );
    }

    if ( side != PastixLCoef ) {
        buffer = core_zlrpack( M, N, blok->LRblock[1], buffer );
    }

    return buffer;
}

/**
 *******************************************************************************
 *
 * @brief Pack low-rank data for column block.
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
 *          The column block to pack.
 *
 * @param[in] size
 *          Size to allocate the buffer.
 *
 *******************************************************************************
 *
 * @return Pointer to the packed data.
 *
 *******************************************************************************/
void *
cpucblk_zpack_lr( pastix_coefside_t side,
                  SolverCblk       *cblk,
                  size_t            size )
{
    assert( cblk->cblktype & CBLK_COMPRESSED );

    char *buffer = malloc(size);
    char *tmp    = buffer;

    const SolverBlok *blok  = cblk->fblokptr;
    const SolverBlok *lblok = cblk[1].fblokptr;

    pastix_int_t N = cblk_colnbr( cblk );

    for ( ; blok < lblok; blok++ ) {
        tmp = cpublok_zpack_lr( side, N, blok, tmp );
    }

    assert( (side == PastixUCoef) || (cblk->lcoeftab == (void*)-1) );
    assert( (side == PastixLCoef) || (cblk->ucoeftab == (void*)-1) );

    return buffer;
}

/**
 *******************************************************************************
 *
 * @brief Unpack low rank data and fill the block concerned by the computation.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] N
 *          Number of columns of the block.
 *
 * @param[inout] blok
 *          The block concerned by the computation.
 *
 * @param[inout] buffer
 *           Pointer on the packed data.
 *
 *******************************************************************************
 *
 * @return Pointer to the next block to unpack.
 *
 *******************************************************************************/
char *
cpublok_zunpack_lr( pastix_coefside_t side,
                    pastix_int_t      N,
                    SolverBlok       *blok,
                    char             *buffer )
{
    pastix_int_t M = blok_rownbr( blok );

    if ( side != PastixUCoef ) {
        buffer = core_zlrunpack( M, N, blok->LRblock[0], buffer );
    }

    if ( side != PastixLCoef ) {
        buffer = core_zlrunpack( M, N, blok->LRblock[1], buffer );
    }

    return buffer;
}

/**
 *******************************************************************************
 *
 * @brief Unpack low rank data and fill the column block concerned by the computation.
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
 *          The column block to fill.
 *
 * @param[inout] buffer
 *          Pointer on packed data.
 *
 *******************************************************************************/
void
cpucblk_zunpack_lr( pastix_coefside_t  side,
                    SolverCblk        *cblk,
                    void              *buffer )
{
    assert( cblk->cblktype & CBLK_COMPRESSED );

    SolverBlok *blok  = cblk->fblokptr;
    SolverBlok *lblok = cblk[1].fblokptr;
    pastix_int_t N    = cblk_colnbr( cblk );

    cpucblk_zalloc_lr( side, cblk, 0 );

    for ( ; blok < lblok; blok++ ) {
        buffer = cpublok_zunpack_lr( side, N, blok, buffer );
    }
}

/**
 *******************************************************************************
 *
 * @brief Pack data in full rank.
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
 * @return Pointer to the data buffer. Note that in full-rank, the data pointer is directly returned
 * without data copy.
 *
 *******************************************************************************/
void *
cpucblk_zpack_fr( pastix_coefside_t  side,
                  const SolverCblk  *cblk )
{
    assert( !( cblk->cblktype & CBLK_COMPRESSED ) );

    return side == PastixUCoef ? cblk->ucoeftab : cblk->lcoeftab;
}

/**
 *******************************************************************************
 *
 * @brief Unpack data in full rank and fill the column block concerned by the computation.
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
 *          The column block to unpack.
 *
 * @param[in] buffer
 *          Pointer on packed data.
 *
 *******************************************************************************/
void
cpucblk_zunpack_fr( pastix_coefside_t   side,
                    SolverCblk         *cblk,
                    pastix_complex64_t *buffer )
{
    assert( !( cblk->cblktype & CBLK_COMPRESSED ) );

    cblk->lcoeftab = buffer;
    if ( side != PastixLCoef ) {
        cblk->ucoeftab = buffer + ( cblk_colnbr( cblk ) * cblk->stride );
    }
}

/**
 *******************************************************************************
 *
 * @brief Pack a column block (Full rank or low rank).
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
 *          TODO
 *
 *******************************************************************************
 *
 * @return Pointer to the data buffer.
 *
 *******************************************************************************/
void *
cpucblk_zpack( pastix_coefside_t  side,
               SolverCblk        *cblk,
               size_t             size )
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
 * @brief Unpack data and fill the column block concerned by the computation.
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
 *          The cblk to unpack.
 *
 * @param[in] buffer
 *          Pointer on packed data.
 *
 *******************************************************************************/
void
cpucblk_zunpack( pastix_coefside_t  side,
                 SolverCblk        *cblk,
                 void              *buffer )
{
    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        cpucblk_zunpack_lr( side, cblk, buffer );

        /*
         * In low-rank, we created a copy so we need to directly free the
         * reception buffer
         */
        free( buffer );
    }
    else {
        cpucblk_zunpack_fr( side, cblk, buffer );
    }
}
