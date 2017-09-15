/**
 *
 * @file core_zscalo.c
 *
 * PaStiX kernel routines
 *
 * @copyright 2010-2015 Univ. of Tennessee, Univ. of California Berkeley and
 *                      Univ. of Colorado Denver. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "blend/solver.h"
#include "pastix_zcores.h"
#include "cblas.h"

/**
 ******************************************************************************
 *
 * @brief Scale a matrix by a diagonal out of place
 *
 * Perform the operation: B <- op(A) * D, where A is a general matrix, and D a
 * diagonal matrix.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *         @arg PastixNoTrans:   No transpose, op( A ) = A;
 *         @arg PastixTrans:     Transpose, op( A ) = A;
 *         @arg PastixConjTrans: Conjugate Transpose, op( A ) = conj(A).
 *
 * @param[in] M
 *          Number of rows of the matrix B.
 *          Number of rows of the matrix A.
 *
 * @param[in] N
 *          Number of columns of the matrix B.
 *          Number of columns of the matrix A.
 *
 * @param[in] A
 *          Matrix of size lda-by-N.
 *
 * @param[in] lda
 *          Leading dimension of the array A. lda >= max(1,M).
 *
 * @param[in] D
 *          Diagonal matrix of size ldd-by-N.
 *
 * @param[in] ldd
 *          Leading dimension of the array D. ldd >= 1.
 *
 * @param[inout] B
 *          Matrix of size LDB-by-N.
 *
 * @param[in] ldb
 *          Leading dimension of the array B. ldb >= max(1,M)
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS successful exit
 * @retval <0 if -i, the i-th argument had an illegal value
 * @retval 1, not yet implemented
 *
 ******************************************************************************/
int
core_zscalo( pastix_trans_t            trans,
             pastix_int_t              M,
             pastix_int_t              N,
             const pastix_complex64_t *A,
             pastix_int_t              lda,
             const pastix_complex64_t *D,
             pastix_int_t              ldd,
             pastix_complex64_t       *B,
             pastix_int_t              ldb )
{
    pastix_complex64_t alpha;
    pastix_int_t i, j;

#if !defined(NDEBUG)
    if ((trans < PastixNoTrans)   ||
        (trans > PastixConjTrans))
    {
        return -1;
    }

    if (M < 0) {
        return -2;
    }
    if (N < 0) {
        return -3;
    }
    if ( lda < pastix_imax(1,M) )
    {
        return -5;
    }
    if ( ldd < 1 )
    {
        return -7;
    }
    if ( ldb < pastix_imax(1,M) ) {
        return -9;
    }
#endif

    if (trans == PastixConjTrans) {
        for( j=0; j<N; j++, D += ldd ) {
            alpha = *D;
            for( i=0; i<M; i++, B++, A++ ) {
                *B = conj(*A) * alpha;
            }
            A += lda - M;
            B += ldb - M;
        }
    }
    else {
        for( j=0; j<N; j++, D += ldd ) {
            alpha = *D;
            for( i=0; i<M; i++, B++, A++ ) {
                *B = (*A) * alpha;
            }
            A += lda - M;
            B += ldb - M;
        }
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Copy the L term with scaling for the two-terms algorithm
 *
 * Performs LD = op(L) * D
 *
 *******************************************************************************
 *
 * @param[in] trans
 *         @arg PastixNoTrans:   No transpose, op( L ) = L;
 *         @arg PastixTrans:     Transpose, op( L ) = L;
 *         @arg PastixConjTrans: Conjugate Transpose, op( L ) = conj(L).
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[inout] LD
 *          The pointer to workspace of size cblk->stride - by -
 *          cblk_colnbr(cblk) that will store L * D on exit, where L is the
 *          lower coeftab array and D the diagonal matrix taken from the
 *          diagonal of the diagonal block.
 *
 *******************************************************************************/
void
cpucblk_zscalo( pastix_trans_t      trans,
                SolverCblk         *cblk,
                pastix_complex64_t *LD )
{
    const SolverBlok *blok, *lblk;

    blok = cblk->fblokptr + 1; /* Firt off-diagonal block */
    lblk = cblk[1].fblokptr;   /* Next diagonal block     */

    /* if there are off-diagonal supernodes in the column */
    if ( blok < lblk )
    {
        const pastix_complex64_t *D;
        const pastix_complex64_t *L;
        pastix_complex64_t *B;
        pastix_int_t M, N, ldl, ldd, ldb;

        N = cblk_colnbr( cblk );

        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            D   = cblk->fblokptr->LRblock[0].u;
            ldd = N+1;

            for(; blok < lblk; blok++) {
                M = blok_rownbr( blok );

                memcpy( blok->LRblock + 1, blok->LRblock, sizeof(pastix_lrblock_t) );

                if ( blok->LRblock[1].rk == -1 ) {
                    assert( M == blok->LRblock[1].rkmax );

                    blok->LRblock[1].u = LD + blok->coefind;

                    L = blok->LRblock[0].u;
                    B = blok->LRblock[1].u;
                    ldl = M;
                    ldb = M;
                }
                else {
                    blok->LRblock[1].v = LD + blok->coefind;
                    L = blok->LRblock[0].v;
                    B = blok->LRblock[1].v;
                    M = blok->LRblock[0].rkmax;
                }

                ldl = M;
                ldb = M;

                /* Compute B = LD */
                core_zscalo( trans, M, N,
                             L, ldl, D, ldd,
                             B, ldb );
            }
        }
        else if ( cblk->cblktype & CBLK_LAYOUT_2D ) {
            L = D = cblk->lcoeftab;
            ldd = N+1;

            for(; blok < lblk; blok++) {
                M = blok_rownbr( blok );

                /* Compute B = LD */
                core_zscalo( trans, M, N,
                             L  + blok->coefind, M, D, ldd,
                             LD + blok->coefind, M );
            }
        }
        else {
            L = D = cblk->lcoeftab;
            ldl = cblk->stride;
            ldd = cblk->stride+1;

            M   = cblk->stride - N;
            B   = LD + blok->coefind;
            ldb = cblk->stride;

            core_zscalo( trans, M, N, L + blok->coefind, ldl, D, ldd, B, ldb );
        }
    }
}
