/**
 *
 * @file core_zherk.c
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
 * @precisions normal z -> c
 *
 **/
#include "common.h"

#undef REAL
#define COMPLEX
#ifdef COMPLEX
/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zherk = PCORE_zherk
#define CORE_zherk PCORE_zherk
#endif
void CORE_zherk(PLASMA_enum uplo, PLASMA_enum trans,
                int N, int K,
                double alpha, const PLASMA_Complex64_t *A, int LDA,
                double beta, PLASMA_Complex64_t *C, int LDC)
{
    cblas_zherk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        N, K,
        alpha, A, LDA,
        beta, C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zherk(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum uplo, PLASMA_enum trans,
                      int n, int k, int nb,
                      double alpha, const PLASMA_Complex64_t *A, int lda,
                      double beta, PLASMA_Complex64_t *C, int ldc)
{
    DAG_CORE_HERK;
    QUARK_Insert_Task(quark, CORE_zherk_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &trans,     VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(double),                     &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double),                     &beta,      VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    C,                 INOUT,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zherk_quark = PCORE_zherk_quark
#define CORE_zherk_quark PCORE_zherk_quark
#endif
void CORE_zherk_quark(Quark *quark)
{
    PLASMA_enum uplo;
    PLASMA_enum trans;
    int n;
    int k;
    double alpha;
    PLASMA_Complex64_t *A;
    int lda;
    double beta;
    PLASMA_Complex64_t *C;
    int ldc;

    quark_unpack_args_10(quark, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
    cblas_zherk(
        CblasColMajor,
        (CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,
        n, k,
        alpha, A, lda,
        beta, C, ldc);
}
#endif
