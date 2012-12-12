/**
 *
 * @file core_zgeadd.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zgeadd adds to matrices together.
 *
 *       B <- alpha * A  + B
 *
 *******************************************************************************
 *
 * @param[in] M
 *          Number of rows of the matrices A and B.
 *
 * @param[in] N
 *          Number of columns of the matrices A and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Matrix of size LDA-by-N.
 *
 * @param[in] LDA
 *          Leading dimension of the array A. LDA >= max(1,M)
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
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeadd = PCORE_zgeadd
#define CORE_zgeadd PCORE_zgeadd
#endif
int CORE_zgeadd(int M, int N, PLASMA_Complex64_t alpha,
                const PLASMA_Complex64_t *A, int LDA,
                      PLASMA_Complex64_t *B, int LDB)
{
    int j;

    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (N < 0) {
        coreblas_error(2, "Illegal value of N");
        return -2;
    }
    if ( (LDA < max(1,M)) && (M > 0) ) {
        coreblas_error(5, "Illegal value of LDA");
        return -5;
    }
    if ( (LDB < max(1,M)) && (M > 0) ) {
        coreblas_error(7, "Illegal value of LDB");
        return -7;
    }

    if (M == LDA && M == LDB)
        cblas_zaxpy(M*N, CBLAS_SADDR(alpha), A, 1, B, 1);
    else {
        for (j = 0; j < N; j++)
            cblas_zaxpy(M, CBLAS_SADDR(alpha), A + j*LDA, 1, B + j*LDB, 1);
    }

    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgeadd(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb, PLASMA_Complex64_t alpha,
                       const PLASMA_Complex64_t *A, int lda,
                             PLASMA_Complex64_t *B, int ldb)
{
    DAG_CORE_GEADD;
    QUARK_Insert_Task(quark, CORE_zgeadd_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeadd_quark = PCORE_zgeadd_quark
#define CORE_zgeadd_quark PCORE_zgeadd_quark
#endif
void CORE_zgeadd_quark(Quark *quark)
{
    int M;
    int N;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;

    quark_unpack_args_7(quark, M, N, alpha, A, LDA, B, LDB);
    CORE_zgeadd(M, N, alpha, A, LDA, B, LDB);
    return;
}

