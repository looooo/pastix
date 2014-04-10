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
#include <cblas.h>

/**
 ******************************************************************************
 *
 * @ingroup CORE_pastix_complex64_t
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
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int core_zgeadd(int trans, int M, int N, pastix_complex64_t alpha,
                const pastix_complex64_t *A, int LDA,
                      pastix_complex64_t *B, int LDB)
{
    int j;

#if defined(PASTIX_DEBUG)
    if (M < 0) {
        //coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        //coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if ( (LDA < pastix_imax(1,M)) && (M > 0) ) {
        //coreblas_error(5, "Illegal value of LDA");
        return -5;
    }
    if ( (LDB < pastix_imax(1,M)) && (M > 0) ) {
        //coreblas_error(7, "Illegal value of LDB");
        return -7;
    }
#endif

    if (trans == CblasNoTrans) {
        if (M == LDA && M == LDB)
            cblas_zaxpy(M*N, CBLAS_SADDR(alpha), A, 1, B, 1);
        else {
            for (j = 0; j < N; j++)
                cblas_zaxpy(M, CBLAS_SADDR(alpha), A + j*LDA, 1, B + j*LDB, 1);
        }
    }
    else if (trans == CblasTrans ) {

    }
    /* trans == CblasConjTrans */
    else {

    }

    return 0;
}
