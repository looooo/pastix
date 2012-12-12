/**
 *
 * @file core_zgemm.c
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
#pragma weak CORE_zgemm = PCORE_zgemm
#define CORE_zgemm PCORE_zgemm
#endif
void CORE_zgemm(PLASMA_enum transA, int transB,
                int M, int N, int K,
                PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int LDA,
                                          const PLASMA_Complex64_t *B, int LDB,
                PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int LDC)
{
    cblas_zgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgemm(Quark *quark, Quark_Task_Flags *task_flags,
                      PLASMA_enum transA, int transB,
                      int m, int n, int k, int nb,
                      PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                                                const PLASMA_Complex64_t *B, int ldb,
                      PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int ldc)
{
    DAG_CORE_GEMM;
    QUARK_Insert_Task(quark, CORE_zgemm_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
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
void QUARK_CORE_zgemm2( Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum transA, int transB,
                        int m, int n, int k, int nb,
                        PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                        const PLASMA_Complex64_t *B, int ldb,
                        PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int ldc)
{
    DAG_CORE_GEMM;
    QUARK_Insert_Task(quark, CORE_zgemm_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(PLASMA_Complex64_t),         &beta,      VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    C,                 INOUT | LOCALITY | GATHERV,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgemm_quark = PCORE_zgemm_quark
#define CORE_zgemm_quark PCORE_zgemm_quark
#endif
void CORE_zgemm_quark(Quark *quark)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    int m;
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

    quark_unpack_args_13(quark, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
    cblas_zgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        m, n, k,
        CBLAS_SADDR(alpha), A, lda,
        B, ldb,
        CBLAS_SADDR(beta), C, ldc);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgemm_f2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, int transB,
                         int m, int n, int k, int nb,
                         PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                                                   const PLASMA_Complex64_t *B, int ldb,
                         PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int ldc,
                         PLASMA_Complex64_t *fake1, int szefake1, int flag1,
                         PLASMA_Complex64_t *fake2, int szefake2, int flag2)
{
    DAG_CORE_GEMM;
    QUARK_Insert_Task(quark, CORE_zgemm_f2_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(PLASMA_Complex64_t),         &beta,      VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    C,                 INOUT | LOCALITY,
        sizeof(int),                        &ldc,       VALUE,
        sizeof(PLASMA_Complex64_t)*szefake1, fake1,             flag1,
        sizeof(PLASMA_Complex64_t)*szefake2, fake2,             flag2,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgemm_f2_quark = PCORE_zgemm_f2_quark
#define CORE_zgemm_f2_quark PCORE_zgemm_f2_quark
#endif
void CORE_zgemm_f2_quark(Quark* quark)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    int M;
    int N;
    int K;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;
    PLASMA_Complex64_t beta;
    PLASMA_Complex64_t *C;
    int LDC;
    void *fake1, *fake2;

    quark_unpack_args_15(quark, transA, transB, M, N, K, alpha,
                         A, LDA, B, LDB, beta, C, LDC, fake1, fake2);
    cblas_zgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgemm_p2(Quark *quark, Quark_Task_Flags *task_flags,
                         PLASMA_enum transA, int transB,
                         int m, int n, int k, int nb,
                         PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                         const PLASMA_Complex64_t **B, int ldb,
                         PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int ldc)
{
    DAG_CORE_GEMM;
    QUARK_Insert_Task(quark, CORE_zgemm_p2_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex64_t*),         B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(PLASMA_Complex64_t),         &beta,      VALUE,
        sizeof(PLASMA_Complex64_t)*ldc*nb,    C,                 INOUT | LOCALITY,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgemm_p2_quark = PCORE_zgemm_p2_quark
#define CORE_zgemm_p2_quark PCORE_zgemm_p2_quark
#endif
void CORE_zgemm_p2_quark(Quark* quark)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    int M;
    int N;
    int K;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t **B;
    int LDB;
    PLASMA_Complex64_t beta;
    PLASMA_Complex64_t *C;
    int LDC;

    quark_unpack_args_13(quark, transA, transB, M, N, K, alpha,
                         A, LDA, B, LDB, beta, C, LDC);
    cblas_zgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        CBLAS_SADDR(alpha), A, LDA,
        *B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgemm_p3(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, int transB,
                           int m, int n, int k, int nb,
                           PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                           const PLASMA_Complex64_t *B, int ldb,
                           PLASMA_Complex64_t beta, PLASMA_Complex64_t **C, int ldc)
{
    DAG_CORE_GEMM;
    QUARK_Insert_Task(quark, CORE_zgemm_p3_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex64_t)*ldb*nb,   B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(PLASMA_Complex64_t),         &beta,      VALUE,
        sizeof(PLASMA_Complex64_t*),         C,                 INOUT | LOCALITY,
        sizeof(int),                        &ldc,       VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgemm_p3_quark = PCORE_zgemm_p3_quark
#define CORE_zgemm_p3_quark PCORE_zgemm_p3_quark
#endif
void CORE_zgemm_p3_quark(Quark* quark)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    int M;
    int N;
    int K;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t *B;
    int LDB;
    PLASMA_Complex64_t beta;
    PLASMA_Complex64_t **C;
    int LDC;

    quark_unpack_args_13(quark, transA, transB, M, N, K, alpha,
                         A, LDA, B, LDB, beta, C, LDC);
    cblas_zgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        CBLAS_SADDR(alpha), A, LDA,
        B, LDB,
        CBLAS_SADDR(beta), *C, LDC);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgemm_p2f1(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum transA, int transB,
                           int m, int n, int k, int nb,
                           PLASMA_Complex64_t alpha, const PLASMA_Complex64_t *A, int lda,
                           const PLASMA_Complex64_t **B, int ldb,
                           PLASMA_Complex64_t beta, PLASMA_Complex64_t *C, int ldc,
                           PLASMA_Complex64_t *fake1, int szefake1, int flag1)
{
    DAG_CORE_GEMM;
    QUARK_Insert_Task(quark, CORE_zgemm_p2f1_quark, task_flags,
        sizeof(PLASMA_enum),                &transA,    VALUE,
        sizeof(PLASMA_enum),                &transB,    VALUE,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(int),                        &k,         VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha,     VALUE,
        sizeof(PLASMA_Complex64_t)*lda*nb,   A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex64_t*),         B,                 INPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(PLASMA_Complex64_t),         &beta,      VALUE,
        sizeof(PLASMA_Complex64_t)*ldc*nb,    C,                 INOUT | LOCALITY,
        sizeof(int),                        &ldc,       VALUE,
        sizeof(PLASMA_Complex64_t)*szefake1, fake1,             flag1,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgemm_p2f1_quark = PCORE_zgemm_p2f1_quark
#define CORE_zgemm_p2f1_quark PCORE_zgemm_p2f1_quark
#endif
void CORE_zgemm_p2f1_quark(Quark* quark)
{
    PLASMA_enum transA;
    PLASMA_enum transB;
    int M;
    int N;
    int K;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_Complex64_t **B;
    int LDB;
    PLASMA_Complex64_t beta;
    PLASMA_Complex64_t *C;
    int LDC;
    void *fake1;

    quark_unpack_args_14(quark, transA, transB, M, N, K, alpha,
                         A, LDA, B, LDB, beta, C, LDC, fake1);
    cblas_zgemm(
        CblasColMajor,
        (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
        M, N, K,
        CBLAS_SADDR(alpha), A, LDA,
        *B, LDB,
        CBLAS_SADDR(beta), C, LDC);
}
