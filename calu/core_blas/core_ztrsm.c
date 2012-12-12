/**
 *
 * @file core_ztrsm.c
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
#pragma weak CORE_ztrsm = PCORE_ztrsm
#define CORE_ztrsm PCORE_ztrsm
#endif
void CORE_ztrsm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                PLASMA_Complex64_t *B, int LDB)
{
    cblas_ztrsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ztrsm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                      PLASMA_Complex64_t *B, int ldb)
{
    DAG_CORE_TRSM;
    QUARK_Insert_Task(quark, CORE_ztrsm_quark, task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    B,                 INOUT,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztrsm_quark = PCORE_ztrsm_quark
#define CORE_ztrsm_quark PCORE_ztrsm_quark
#endif
void CORE_ztrsm_quark(Quark *quark)
{
    PLASMA_enum side;
    PLASMA_enum uplo;
    PLASMA_enum transA;
    PLASMA_enum diag;
    int m;
    int n;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_Complex64_t *B;
    int ldb;

    quark_unpack_args_11(quark, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
    cblas_ztrsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        m, n,
        CBLAS_SADDR(alpha), A, lda,
        B, ldb);
}
