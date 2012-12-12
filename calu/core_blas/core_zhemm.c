/**
 *
 * @file core_zhemm.c
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
#pragma weak CORE_zhemm = PCORE_zhemm
#define CORE_zhemm PCORE_zhemm
#endif
void CORE_zhemm(PLASMA_enum side, PLASMA_enum uplo,
                int M, int N,
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                const PLASMA_Complex64_t *B, int LDB,
                PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC)
{
    cblas_zhemm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zhemm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo,
                      int m, int n, int nb,
                      PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                      const PLASMA_Complex64_t *B, int ldb,
                      PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int ldc)
{
    DAG_CORE_HEMM;
    QUARK_Insert_Task(quark, CORE_zhemm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,    VALUE,
        sizeof(PLASMA_enum),                &uplo,    VALUE,
        sizeof(int),                        &m,       VALUE,
        sizeof(int),                        &n,       VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,               INPUT,
        sizeof(int),                        &lda,     VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    B,               INPUT,
        sizeof(int),                        &ldb,     VALUE,
        sizeof(PLASMA_Complex64_t),         &beta,    VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    C,               INOUT,
        sizeof(int),                        &ldc,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zhemm_quark = PCORE_zhemm_quark
#define CORE_zhemm_quark PCORE_zhemm_quark
#endif
void CORE_zhemm_quark(Quark *quark)
{
    PLASMA_enum side;
    PLASMA_enum uplo;
    int M;
    int N;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;
    PLASMA_Complex64_t beta;
    PLASMA_Complex64_t *C;
    int LDC;

    quark_unpack_args_12(quark, side, uplo, M, N, alpha, A, LDA, B, LDB, beta, C, LDC);
    cblas_zhemm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}
#endif
