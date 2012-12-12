/**
 *
 * @file core_zsyr2k.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Jakub Kurzak
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zsyr2k = PCORE_zsyr2k
#define CORE_zsyr2k PCORE_zsyr2k
#endif
void CORE_zsyr2k(PLASMA_enum uplo, PLASMA_enum trans,
                 int N, int K,
                 PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                 const PLASMA_Complex64_t *B, int LDB,
                 PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC)
{
    cblas_zsyr2k(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        CBLAS_SADDR(alpha), A, LDA, B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zsyr2k(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum trans,
                       int n, int k, int nb,
                       PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                       const PLASMA_Complex64_t *B, int ldb,
                       PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int ldc)
{
    DAG_CORE_SYR2K;
    QUARK_Insert_Task(quark, CORE_zsyr2k_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(PLASMA_Complex64_t),         &beta,      VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zsyr2k_quark = PCORE_zsyr2k_quark
#define CORE_zsyr2k_quark PCORE_zsyr2k_quark
#endif
void CORE_zsyr2k_quark(Quark *quark)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    int n;
    int k;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_Complex64_t *B;
    int ldb;
    PLASMA_Complex64_t beta;
    PLASMA_Complex64_t *C;
    int ldc;

    quark_unpack_args_12(quark, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    cblas_zsyr2k(CblasColMajor, (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
                 n, k, CBLAS_SADDR(alpha), A, lda, B, ldb, CBLAS_SADDR(beta), C, ldc);
}
