/**
 *
 * @file core_zgemdm.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.5
 * @author Dulceneia Becker
 * @author Mathieu Faverge
 * @date 2011-1-18
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <cblas.h>

/***************************************************************************//**
 *
 * @ingroup CORE_pastix_complex64_t
 *
 * CORE_zgemdm performs one of the matrix-matrix operations
 *
 *       C := alpha*op( A )*D*op( B ) + beta*C,
 *
 * where op( X ) is one of
 *
 *       op( X ) = X   or   op( X ) = X',
 *
 * alpha and beta are scalars, and A, B, C and D are matrices, with
 *
 *       op( A ) an m by k matrix,
 *       op( B ) an k by n matrix,
 *       C an m by n matrix and
 *       D an k by k matrix.
 *
 *******************************************************************************
 *
 * @param[in] TRANSA
 *         INTEGER
 *         @arg PlasmaNoTrans   :  No transpose, op( A ) = A;
 *         @arg PlasmaConjTrans :  Transpose, op( A ) = A'.
 *
 * @param[in] TRANSB
 *         INTEGER
 *         @arg PlasmaNoTrans   :  No transpose, op( B ) = B;
 *         @arg PlasmaConjTrans :  Transpose, op( B ) = B'.
 *
 * @param[in] M
 *         INTEGER
 *         The number of rows  of the  matrix op( A ) and of the
 *         matrix C.  M  must  be at least  zero.
 *
 * @param[in] N
 *         INTEGER
 *         The number of columns of the matrix  op( B ) and the
 *         number of columns of the matrix C. N must be at least zero.
 *
 * @param[in] K
 *         INTEGER
 *         The number of columns of the matrix op( A ), the number of
 *         rows of the matrix op( B ), and the number of rows and columns
 *         of matrix D. K must be at least  zero.
 *
 * @param[in] ALPHA
 *         pastix_complex64_t.
 *         On entry, ALPHA specifies the scalar alpha.
 *         Unchanged on exit.
 *
 * @param[in] A
 *         pastix_complex64_t array of DIMENSION ( LDA, ka ), where ka is
 *         k  when  TRANSA = PlasmaTrans, and is  m  otherwise.
 *         Before entry with  TRANSA = PlasmaTrans,  the leading  m by k
 *         part of the array  A  must contain the matrix  A,  otherwise
 *         the leading  k by m  part of the array  A  must contain  the
 *         matrix A.
 *         Unchanged on exit.
 *
 * @param[in] LDA
 *        INTEGER.
 *        On entry, LDA specifies the first dimension of A as declared
 *        in the calling (sub) program. When TRANSA = PlasmaTrans then
 *        LDA must be at least  max( 1, m ), otherwise  LDA must be at
 *        least  max( 1, k ).
 *        Unchanged on exit.
 *
 * @param[in] B
 *        pastix_complex64_t array of DIMENSION ( LDB, kb ), where kb is
 *        n  when TRANSB = PlasmaTrans, and is k otherwise.
 *        Before entry with TRANSB = PlasmaTrans, the leading  k by n
 *        part of the array  B  must contain the matrix B, otherwise
 *        the leading n by k part of the array B must contain  the
 *        matrix B.
 *        Unchanged on exit.
 *
 * @param[in] LDB
 *       INTEGER.
 *       On entry, LDB specifies the first dimension of B as declared
 *       in the calling (sub) program. When  TRANSB = PlasmaTrans then
 *       LDB must be at least  max( 1, k ), otherwise  LDB must be at
 *       least  max( 1, n ).
 *       Unchanged on exit.
 *
 * @param[in] BETA
 *       pastix_complex64_t.
 *       On entry,  BETA  specifies the scalar  beta.  When  BETA  is
 *       supplied as zero then C need not be set on input.
 *       Unchanged on exit.
 *
 * @param[in] C
 *       pastix_complex64_t array of DIMENSION ( LDC, n ).
 *       Before entry, the leading  m by n  part of the array  C must
 *       contain the matrix  C,  except when  beta  is zero, in which
 *       case C need not be set on entry.
 *       On exit, the array  C  is overwritten by the  m by n  matrix
 *       ( alpha*op( A )*D*op( B ) + beta*C ).
 *
 * @param[in] LDC
 *       INTEGER
 *       On entry, LDC specifies the first dimension of C as declared
 *       in  the  calling  (sub)  program.   LDC  must  be  at  least
 *       max( 1, m ).
 *       Unchanged on exit.
 *
 * @param[in] D
 *        pastix_complex64_t array of DIMENSION ( LDD, k ).
 *        Before entry, the leading  k by k part of the array  D
 *        must contain the matrix D.
 *        Unchanged on exit.
 *
 * @param[in] LDD
 *       INTEGER.
 *       On entry, LDD specifies the first dimension of D as declared
 *       in  the  calling  (sub)  program.   LDD  must  be  at  least
 *       max( 1, k ).
 *       Unchanged on exit.
 *
 * @param[workspace] WORK
 *       pastix_complex64_t array, dimension (MAX(1,LWORK))
 *
 * @param[in] LWORK
 *       INTEGER
 *       The length of WORK.
 *       On entry, if TRANSA = PlasmaTrans and TRANSB = PlasmaTrans then
 *       LWORK >= max(1, K*N). Otherwise LWORK >= max(1, M*K).
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 ******************************************************************************/
int CORE_zgemdm(int transA, int transB,
                int M, int N, int K,
                pastix_complex64_t alpha, pastix_complex64_t *A, int LDA,
                pastix_complex64_t *B, int LDB,
                pastix_complex64_t beta, pastix_complex64_t *C, int LDC,
                pastix_complex64_t *D, int incD,
                pastix_complex64_t *WORK, int LWORK)
{
    int j, Am, Bm;
    pastix_complex64_t delta;
    pastix_complex64_t *wD, *w;

    Am = (transA == CblasNoTrans ) ? M : K;
    Bm = (transB == CblasNoTrans ) ? K : N;

    /* Check input arguments */
    if ((transA != CblasNoTrans) && (transA != CblasTrans) && (transA != CblasConjTrans)) {
        //coreblas_error(1, "Illegal value of transA");
        return -1;
    }
    if ((transB != CblasNoTrans) && (transB != CblasTrans) && (transB != CblasConjTrans)) {
        //coreblas_error(2, "Illegal value of transB");
        return -2;
    }
    if (M < 0) {
        //coreblas_error(3, "Illegal value of M");
        return -3;
    }
    if (N < 0) {
        //coreblas_error(4, "Illegal value of N");
        return -4;
    }
    if (K < 0) {
        //coreblas_error(5, "Illegal value of K");
        return -5;
    }
    if ((LDA < pastix_imax(1,Am)) && (Am > 0)) {
        //coreblas_error(8, "Illegal value of LDA");
        return -8;
    }
    if ((LDB < pastix_imax(1,Bm)) && (Bm > 0)) {
        //coreblas_error(10, "Illegal value of LDB");
        return -10;
    }
    if ((LDC < pastix_imax(1,M)) && (M > 0)) {
        //coreblas_error(13, "Illegal value of LDC");
        return -13;
    }
    if ( incD < 0 ) {
        //coreblas_error(15, "Illegal value of incD");
        return -15;
    }
    if ( ( ( transA == CblasNoTrans ) && ( LWORK < (M+1)*K) ) ||
         ( ( transA != CblasNoTrans ) && ( LWORK < (N+1)*K) ) ){
        //coreblas_error(17, "Illegal value of LWORK");
        return -17;
    }

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == 0.0 || K == 0) && beta == 1.0) ) {
        return PASTIX_SUCCESS;
    }

    if ( incD == 1 ) {
        wD = D;
    } else {
        wD = WORK;
        cblas_zcopy(K, D, incD, wD, 1);
    }
    w = WORK + K;

    /*
     * transA == CblasNoTrans
     */
    if ( transA == CblasNoTrans )
    {
        /* WORK = A * D */
      for (j=0; j<K; j++, wD++) {
            delta = *wD;
            cblas_zcopy(M, &A[LDA*j], 1,       &w[M*j], 1);
            cblas_zscal(M, CBLAS_SADDR(delta), &w[M*j], 1);
        }

        /* C = alpha * WORK * op(B) + beta * C */
        cblas_zgemm(CblasColMajor, CblasNoTrans, transB,
                    M, N, K,
                    CBLAS_SADDR(alpha), w, M,
                                        B, LDB,
                    CBLAS_SADDR(beta),  C, LDC);
    }
    else
    {
        if ( transB == CblasNoTrans ) /* Worst case*/
        {
            /* WORK = (D * B)' */
          for (j=0; j<K; j++, wD++) {
                delta = *wD;
                cblas_zcopy(N, &B[j],     LDB,     &w[N*j], 1);
                cblas_zscal(N, CBLAS_SADDR(delta), &w[N*j], 1);
            }

            /* C = alpha * op(A) * WORK' + beta * C */
            cblas_zgemm(CblasColMajor, transA, CblasTrans,
                        M, N, K,
                        CBLAS_SADDR(alpha), A, LDA,
                                            w, N,
                        CBLAS_SADDR(beta),  C, LDC);
        }
        else
        {
#ifdef COMPLEX
            if ( transB == CblasConjTrans )
            {
                /* WORK = D * B' */
              for (j=0; j<K; j++, wD++) {
                    delta = *wD;
                    cblas_zcopy(N, &B[LDB*j], 1,       &w[N*j], 1);
                    LAPACKE_zlacgv_work(N,             &w[N*j], 1);
                    cblas_zscal(N, CBLAS_SADDR(delta), &w[N*j], 1);
                }
            }
            else
#endif
            {
                /* WORK = D * B' */
              for (j=0; j<K; j++, wD++) {
                    delta = *wD;
                    cblas_zcopy(N, &B[LDB*j], 1,       &w[N*j], 1);
                    cblas_zscal(N, CBLAS_SADDR(delta), &w[N*j], 1);
                }
            }

            /* C = alpha * op(A) * WORK + beta * C */
            cblas_zgemm(CblasColMajor, transA, CblasNoTrans,
                        M, N, K,
                        CBLAS_SADDR(alpha), A, LDA,
                                            w, N,
                        CBLAS_SADDR(beta),  C, LDC);
        }
    }
    return PASTIX_SUCCESS;
}
