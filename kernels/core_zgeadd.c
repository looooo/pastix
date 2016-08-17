/**
 *
 * @file core_zgeadd.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.6.0
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
 * @ingroup pastix_kernel
 *
 *  core_zgeadd adds to matrices together.
 *
 *       B <- alpha * op(A)  + B
 *
 *******************************************************************************
 *
 * @param[in] trans
 *         @arg CblasNoTrans   :  No transpose, op( A ) = A;
 *         @arg CblasTrans     :  Transpose, op( A ) = A';
 *         @arg CblasConjTrans :  Conjugate Transpose, op( A ) = conj(A').
 *
 * @param[in] M
 *          Number of rows of the matrix B.
 *          Number of rows of the matrix A, if trans == CblasNoTrans, number of
 *          columns of A otherwise.
 *
 * @param[in] N
 *          Number of columns of the matrix B.
 *          Number of columns of the matrix A, if trans == CblasNoTrans, number
 *          of rows of A otherwise.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N, if trans == CblasNoTrans, LDA-by-M,
 *          otherwise.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,K).
 *          K = M if trans == CblasNoTrans, K = N otherwise.
 *
 * @param[in,out] B
 *          Matrix of size LDB-by-N.
 *
 * @param[in] LDB
 *          Leading dimension of the array B. LDB >= max(1,M)
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *          \retval 1, not yet implemented
 *
 ******************************************************************************/
int core_zgeadd( pastix_int_t trans, pastix_int_t M, pastix_int_t N,
                       pastix_complex64_t  alpha,
                 const pastix_complex64_t *A, pastix_int_t LDA,
                       pastix_complex64_t  beta,
                       pastix_complex64_t *B, pastix_int_t LDB)
{
    int i, j;

#if defined(PASTIX_DEBUG)
    if ((trans != PastixNoTrans) &&
        (trans != PastixTrans)   &&
        (trans != PastixConjTrans))
    {
        //coreblas_error(1, "illegal value of trans");
        return -1;
    }

    if (M < 0) {
        //coreblas_error(2, "Illegal value of M");
        return -2;
    }
    if (N < 0) {
        coreblas_error(3, "Illegal value of N");
        return -3;
    }
    if ( ((trans == PastixNoTrans) && (LDA < max(1,M)) && (M > 0)) ||
         ((trans != PastixNoTrans) && (LDA < max(1,N)) && (N > 0)) )
    {
        //coreblas_error(6, "Illegal value of LDA");
        return -6;
    }
    if ( (LDB < max(1,M)) && (M > 0) ) {
        //coreblas_error(8, "Illegal value of LDB");
        return -8;
    }
#endif

    switch( trans ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PastixConjTrans:
        if ( alpha == 0. ) {
            for (j=0; j<N; j++) {
                for(i=0; i<M; i++, B++) {
                    *B = beta * (*B);
                }
                B += LDB-M;
            }
        }
        else if ( beta == 0. ) {
            for (j=0; j<N; j++, A++) {
                for(i=0; i<M; i++, B++) {
                    *B = alpha * conj(A[LDA*i]);
                }
                B += LDB-M;
            }
        }
        else {
            for (j=0; j<N; j++, A++) {
                for(i=0; i<M; i++, B++) {
                    *B = beta * (*B) + alpha * conj(A[LDA*i]);
                }
                B += LDB-M;
            }
        }
        break;
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    case PastixTrans:
        if ( alpha == 0. ) {
            for (j=0; j<N; j++) {
                for(i=0; i<M; i++, B++) {
                    *B = beta * (*B);
                }
                B += LDB-M;
            }
        }
        else if ( beta == 0. ) {
            for (j=0; j<N; j++, A++) {
                for(i=0; i<M; i++, B++) {
                    *B = alpha * A[LDA*i];
                }
                B += LDB-M;
            }
        }
        else {
            for (j=0; j<N; j++, A++) {
                for(i=0; i<M; i++, B++) {
                    *B = beta * (*B) + alpha * A[LDA*i];
                }
                B += LDB-M;
            }
        }
        break;

    case PastixNoTrans:
    default:
        if ( alpha == 0. ) {
            for (j=0; j<N; j++) {
                for(i=0; i<M; i++, B++) {
                    *B = beta * (*B);
                }
                B += LDB-M;
            }
        }
        else if ( beta == 0. ) {
            for (j=0; j<N; j++) {
                for(i=0; i<M; i++, B++, A++) {
                    *B = alpha * (*A);
                }
                A += LDA-M;
                B += LDB-M;
            }
        }
        else {
            for (j=0; j<N; j++) {
                for(i=0; i<M; i++, B++, A++) {
                    *B = beta * (*B) + alpha * (*A);
                }
                A += LDA-M;
                B += LDB-M;
            }
        }
    }
    return PASTIX_SUCCESS;
}

/**
 ******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 *  core_zgeaddsp1d adds two column blocks together.
 *
 *       cblk2 <- cblk1  + cblk2
 *
 *******************************************************************************
 *
 * @param[in] cblk1
 *          The pointer to the data structure that describes the panel to add.
 *          Next column blok must be accessible through cblk1[1].
 *
 * @param[in] cblk2
 *          The pointer to the data structure that describes the panel in which
 *          we add.
 *          Next column blok must be accessible through cblk2[1].
 *
 * @param[in] L
 *          The pointer to the lower matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] Cl
 *          The pointer to the lower matrix storing the coefficients of the
 *          updated panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] U
 *          The pointer to the upper matrix storing the coefficients of the
 *          panel. Must be of size cblk.stride -by- cblk.width. Ignored if
 *          NULL.
 *
 *
 * @param[in,out] Cu
 *          The pointer to the upper matrix storing the coefficients of the
 *          updated panel. Must be of size cblk.stride -by- cblk.width
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS successful exit
 *
 ******************************************************************************/
int
core_zgeaddsp1d(SolverCblk * cblk1,
                SolverCblk * cblk2,
                pastix_complex64_t * L,
                pastix_complex64_t * Cl,
                pastix_complex64_t * U,
                pastix_complex64_t * Cu)
{
    SolverBlok *iterblok;
    SolverBlok *firstblok;
    SolverBlok *lastblok;
    SolverBlok *fblok;
    pastix_int_t ncol1 = cblk1->lcolnum - cblk1->fcolnum + 1;
    pastix_complex64_t *ga, *gb;

    assert(0 /* Outdated */);

    firstblok = cblk1->fblokptr;
    lastblok  = cblk1[1].fblokptr;
    fblok = cblk2->fblokptr;

    assert(cblk1->fcolnum >= cblk2->fcolnum);

    for (iterblok = firstblok; iterblok < lastblok; iterblok++) {
        pastix_int_t nrow;
        /* Find facing bloknum */
        while (!is_block_inside_fblock( iterblok, fblok )) {
            fblok++;
            assert( fblok < cblk2[1].fblokptr );
        }
        ga = L + iterblok->coefind;
        gb = Cl + cblk2->stride*(cblk1->fcolnum-cblk2->fcolnum) +
            fblok->coefind +
            iterblok->frownum - fblok->frownum;
        nrow = iterblok->lrownum - iterblok->frownum + 1;
        core_zgeadd( CblasNoTrans,
                     nrow, ncol1,
                     1.0, ga, cblk1->stride,
		     1.0, gb, cblk2->stride );
        if (U != NULL) {
            ga = U + iterblok->coefind;
            gb = Cu + cblk2->stride*(cblk1->fcolnum-cblk2->fcolnum) +
                fblok->coefind +
                iterblok->frownum - fblok->frownum;
            core_zgeadd( CblasNoTrans,
                         nrow, ncol1,
                         1.0, ga, cblk1->stride,
                         1.0, gb, cblk2->stride );
        }
    }
    return PASTIX_SUCCESS;
}
