/**
 *
 * @file core_ztrmm.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Julien Langou
 * @author Henricus Bouwmeester
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
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztrmm = PCORE_ztrmm
#define CORE_ztrmm PCORE_ztrmm
#endif
void CORE_ztrmm(PLASMA_enum side, PLASMA_enum uplo,
                PLASMA_enum transA, PLASMA_enum diag,
                int M, int N,
                PLASMA_Complex64_t alpha,
                const PLASMA_Complex64_t *A, int LDA,
                PLASMA_Complex64_t *B, int LDB)
{
    cblas_ztrmm(
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
void QUARK_CORE_ztrmm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                      int m, int n, int nb,
                      PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                      PLASMA_Complex64_t *B, int ldb)
{
    DAG_CORE_TRMM;
    QUARK_Insert_Task(quark, CORE_ztrmm_quark, task_flags,
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
#pragma weak CORE_ztrmm_quark = PCORE_ztrmm_quark
#define CORE_ztrmm_quark PCORE_ztrmm_quark
#endif
void CORE_ztrmm_quark(Quark *quark)
{
    PLASMA_enum side;
    PLASMA_enum uplo;
    PLASMA_enum transA;
    PLASMA_enum diag;
    int M;
    int N;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;

    quark_unpack_args_11(quark, side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB);
    cblas_ztrmm(
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
void QUARK_CORE_ztrmm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum transA, PLASMA_enum diag,
                         int m, int n, int nb,
                         PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                         PLASMA_Complex64_t **B, int ldb)
{
    DAG_CORE_TRMM;
    QUARK_Insert_Task(quark, CORE_ztrmm_p2_quark, task_flags,
        sizeof(PLASMA_enum),                &side,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex64_t*),         B,                 INOUT,
        sizeof(int),                        &ldb,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztrmm_p2_quark = PCORE_ztrmm_p2_quark
#define CORE_ztrmm_p2_quark PCORE_ztrmm_p2_quark
#endif
void CORE_ztrmm_p2_quark(Quark *quark)
{
    PLASMA_enum side;
    PLASMA_enum uplo;
    PLASMA_enum transA;
    PLASMA_enum diag;
    int M;
    int N;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t **B;
    int LDB;

    quark_unpack_args_11(quark, side, uplo, transA, diag, M, N, alpha, A, LDA, B, LDB);
    cblas_ztrmm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        M, N,
        CBLAS_SADDR(alpha), A, LDA,
        *B, LDB);
}
