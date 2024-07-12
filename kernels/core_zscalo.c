/**
 *
 * @file core_zscalo.c
 *
 * PaStiX kernel routines
 *
 * @copyright 2010-2015 Univ. of Tennessee, Univ. of California Berkeley and
 *                      Univ. of Colorado Denver. All rights reserved.
 * @copyright 2012-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Mathieu Faverge
 * @author Nolan Bredel
 * @date 2024-07-05
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "blend/solver.h"
#include "pastix_zcores.h"
#include "cblas.h"
#include "kernels_trace.h"

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

#if defined(PRECISION_z) || defined(PRECISION_c)
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
    else
#endif
    {
        for( j=0; j<N; j++, D += ldd ) {
            alpha = *D;
            for( i=0; i<M; i++, B++, A++ ) {
                *B = (*A) * alpha;
            }
            A += lda - M;
            B += ldb - M;
        }
    }

    (void)trans;
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
 * @param[inout] dataL
 *          The pointer to the correct representation of lower part of the data.
 *          - coeftab if the block is in full rank. Must be of size cblk.stride -by- cblk.width.
 *          - pastix_lr_block if the block is compressed.
 *
 * @param[inout] dataLD
 *          The pointer to the correct representation of LD.
 *          - coeftab if the block is in full rank. Must be of size cblk.stride -by- cblk.width.
 *          - pastix_lr_block if the block is compressed.
 *
 *******************************************************************************/
void
cpucblk_zscalo( pastix_trans_t    trans,
                const SolverCblk *cblk,
                void             *dataL,
                void             *dataLD )
{
    const SolverBlok *blok, *lblk;
    pastix_int_t M, N;
    pastix_lrblock_t *lrL, *lrLD;
    pastix_fixdbl_t time;
    pastix_complex64_t *LD;

    time = kernel_trace_start( PastixKernelSCALOCblk );

    N = cblk_colnbr( cblk );

    blok = cblk->fblokptr + 1; /* Firt off-diagonal block */
    lblk = cblk[1].fblokptr;   /* Next diagonal block     */

    /* if there are off-diagonal supernodes in the column */
    if ( blok < lblk )
    {
        const pastix_complex64_t *L;
        const pastix_complex64_t *D;
        pastix_int_t ldl, ldd, ldld;

        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            lrL  = (pastix_lrblock_t *)dataL;
            lrLD = (pastix_lrblock_t *)dataLD;
            D   = lrL->u;
            ldd = N+1;

            lrL++; lrLD++;
            for(; blok < lblk; blok++, lrL++, lrLD++) {
                M = blok_rownbr( blok );

                assert( lrLD->rk == -1 );

                /* Copy L in LD */
                lrLD->rk    = lrL->rk;
                lrLD->rkmax = lrL->rkmax;

                if ( lrL->rk == -1 ) {
                    assert( M == lrL->rkmax );

                    /* Initialize the workspace */
                    memcpy( lrLD->u, lrL->u, lrL->rkmax * N * sizeof(pastix_complex64_t) );
                    lrLD->v = NULL;

                    L  = lrL->u;
                    LD = lrLD->u;
                }
                else {
                    /*
                     * Initialize the workspace
                     */
                    memcpy( lrLD->u, lrL->u, M * lrL->rk    * sizeof(pastix_complex64_t) );
                    lrLD->v = ((pastix_complex64_t *)lrLD->u) + M * lrL->rk;
                    memcpy( lrLD->v, lrL->v, N * lrL->rkmax * sizeof(pastix_complex64_t) );

                    L  = lrL->v;
                    LD = lrLD->v;
                    M  = lrLD->rkmax;
                }

                ldl  = M;
                ldld = M;

                /* Compute LD = L * D */
                core_zscalo( trans, M, N,
                             L, ldl, D, ldd,
                             LD, ldld );
            }
        }
        else if ( cblk->cblktype & CBLK_LAYOUT_2D ) {
            L = D = (pastix_complex64_t *)dataL;
            LD = (pastix_complex64_t *)dataLD;
            ldd = N+1;

            for(; blok < lblk; blok++) {
                M = blok_rownbr( blok );

                /* Compute LD = L * D */
                core_zscalo( trans, M, N,
                             L  + blok->coefind, M, D, ldd,
                             LD + blok->coefind, M );
            }
        }
        else {
            L = D = (pastix_complex64_t *)dataL;
            LD = (pastix_complex64_t *)dataLD;
            ldl = cblk->stride;
            ldd = cblk->stride+1;

            M    = cblk->stride - N;
            LD   = LD + blok->coefind;
            ldld = cblk->stride;

            core_zscalo( trans, M, N, L + blok->coefind, ldl, D, ldd, LD, ldld );
        }
    }

    M = cblk->stride - N;
    kernel_trace_stop( cblk->fblokptr->inlast, PastixKernelSCALOCblk, M, N, 0, (pastix_fixdbl_t)(M*N), time );
}

/**
 *******************************************************************************
 *
 * @brief Copy the lower terms of the block with scaling for the two-terms
 * algorithm.
 *
 * Performs B = op(A) * D
 *
 *******************************************************************************
 *
 * @param[in] trans
 *         @arg PastixNoTrans:   No transpose, op( A ) = A;
 *         @arg PastixTrans:     Transpose, op( A ) = A;
 *         @arg PastixConjTrans: Conjugate Transpose, op( A ) = conj(A).
 *
 * @param[in] cblk
 *          Pointer to the structure representing the panel to factorize in the
 *          cblktab array.  Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok_m
 *          Index of the off-diagonal block to be solved in the cblk. All blocks
 *          facing the same cblk, in the current column block will be solved.
 *
 * @param[in] dataA
 *          The pointer to the correct representation of data of A.
 *          - coeftab if the block is in full rank. Must be of size cblk.stride -by- cblk.width.
 *          - pastix_lr_block if the block is compressed.
 *
 * @param[in] dataD
 *          The pointer to the correct representation of data of D.
 *          - coeftab if the block is in full rank. Must be of size cblk.stride -by- cblk.width.
 *          - pastix_lr_block if the block is compressed.
 *
 * @param[inout] dataB
 *          The pointer to the correct representation of data of B.
 *          - coeftab if the block is in full rank. Must be of size cblk.stride -by- cblk.width.
 *          - pastix_lr_block if the block is compressed.
 *
 *******************************************************************************/
void
cpublok_zscalo( pastix_trans_t    trans,
                const SolverCblk *cblk,
                pastix_int_t      blok_m,
                const void       *dataA,
                const void       *dataD,
                void             *dataB )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, ldd, offset, cblk_m;
    const pastix_complex64_t *lA;
    pastix_lrblock_t *lrD, *lrB, *lrA;
    pastix_complex64_t *D, *B, *A;
    pastix_complex64_t *lB;

    N     = cblk_colnbr( cblk );
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */
    ldd   = blok_rownbr( fblok ) + 1;

    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_LAYOUT_2D );

    blok   = fblok + blok_m;
    offset = blok->coefind;
    cblk_m = blok->fcblknm;

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        lrA = (pastix_lrblock_t *)dataA;
        lrD = (pastix_lrblock_t *)dataD;
        lrB = (pastix_lrblock_t *)dataB;
        D = lrD->u;
        for (; (blok < lblok) && (blok->fcblknm == cblk_m); blok++, lrA++, lrB++) {
            M = blok_rownbr( blok );

            /* Copy A in B */
            lrB->rk    = lrA->rk;
            lrB->rkmax = lrA->rkmax;

            if ( lrB->rk == -1 ) {
                assert( M == lrA->rkmax );
                assert( NULL == lrA->v );

                /* Initialize the workspace */
                memcpy( lrB->u, lrA->u, lrA->rkmax * N * sizeof(pastix_complex64_t) );
                lrB->v = NULL;

                lA = lrA->u;
                lB = lrB->u;
            }
            else {
                /*
                 * Initialize the workspace
                 */
                memcpy( lrB->u, lrA->u, M * lrA->rk    * sizeof(pastix_complex64_t) );
                lrB->v = ((pastix_complex64_t *)lrB->u) + M * lrA->rk;
                memcpy( lrB->v, lrA->v, N * lrA->rkmax * sizeof(pastix_complex64_t) );

                lA = lrA->v;
                lB = lrB->v;
                M  = lrA->rkmax;
            }

            /* Compute B = op(A) * D */
            core_zscalo( trans, M, N,
                         lA, M, D, ldd, lB, M );
        }
    }
    else {
        A = (pastix_complex64_t *)dataA;
        D = (pastix_complex64_t *)dataD;
        B = (pastix_complex64_t *)dataB;

        for (; (blok < lblok) && (blok->fcblknm == cblk_m); blok++) {
            lA = A + blok->coefind - offset;
            lB = B + blok->coefind - offset;
            M  = blok_rownbr(blok);

            /* Compute B = op(A) * D */
            core_zscalo( trans, M, N,
                         lA, M, D, ldd, lB, M );
        }
    }
}
