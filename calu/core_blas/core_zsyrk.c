/**
 *
 * @file core_zsyrk.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
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
#pragma weak CORE_zsyrk = PCORE_zsyrk
#define CORE_zsyrk PCORE_zsyrk
#endif
void CORE_zsyrk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC)
{
    cblas_zsyrk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        CBLAS_SADDR(alpha), A, LDA,
        CBLAS_SADDR(beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zsyrk(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum uplo, PLASMA_enum trans,
                      int n, int k, int nb,
                      PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                      PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int ldc)
{
    DAG_CORE_SYRK;
    QUARK_Insert_Task(quark, CORE_zsyrk_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex64_t),         &beta,      VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zsyrk_quark = PCORE_zsyrk_quark
#define CORE_zsyrk_quark PCORE_zsyrk_quark
#endif
void CORE_zsyrk_quark(Quark *quark)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    int n;
    int k;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_Complex64_t beta;
    PLASMA_Complex64_t *C;
    int ldc;

    quark_unpack_args_10(quark, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
    cblas_zsyrk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        n, k,
        CBLAS_SADDR(alpha), A, lda,
        CBLAS_SADDR(beta), C, ldc);
}
