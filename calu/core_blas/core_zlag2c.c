/**
 *
 * @file core_zlag2c.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 **/
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlag2c = PCORE_zlag2c
#define CORE_zlag2c PCORE_zlag2c
#endif
void CORE_zlag2c(int m, int n,
                 const PLASMA_Complex64_t *A, int lda,
                 PLASMA_Complex32_t *B, int ldb, int *info)
{
    *info = LAPACKE_zlag2c_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}
/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlag2c(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       const PLASMA_Complex64_t *A, int lda,
                       PLASMA_Complex32_t *B, int ldb,
                       PLASMA_sequence *sequence, PLASMA_request *request)
{
    DAG_CORE_LAG2C;
    QUARK_Insert_Task(quark, CORE_zlag2c_quark, task_flags,
        sizeof(int),                        &m,         VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    B,                 OUTPUT,
        sizeof(int),                        &ldb,       VALUE,
        sizeof(PLASMA_sequence*),           &sequence,  VALUE,
        sizeof(PLASMA_request*),            &request,   VALUE,
        0);
}


/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlag2c_quark = PCORE_zlag2c_quark
#define CORE_zlag2c_quark PCORE_zlag2c_quark
#endif
void CORE_zlag2c_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_Complex32_t *B;
    int ldb;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int info;

    quark_unpack_args_8(quark, m, n, A, lda, B, ldb, sequence, request);
    info = LAPACKE_zlag2c_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
    if (sequence->status == PLASMA_SUCCESS && info != 0)
        plasma_sequence_flush(quark, sequence, request, info);
}

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clag2z = PCORE_clag2z
#define CORE_clag2z PCORE_clag2z
#endif
void CORE_clag2z(int m, int n,
                 const PLASMA_Complex32_t *A, int lda,
                 PLASMA_Complex64_t *B, int ldb)
{
    LAPACKE_clag2z_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_clag2z(Quark *quark, Quark_Task_Flags *task_flags,
                       int m, int n, int nb,
                       const PLASMA_Complex32_t *A, int lda,
                       PLASMA_Complex64_t *B, int ldb)
{
    QUARK_Insert_Task(quark, CORE_clag2z_quark, task_flags,
        sizeof(int),                        &m,     VALUE,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex32_t)*nb*nb,    A,             INPUT,
        sizeof(int),                        &lda,   VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    B,             INOUT,
        sizeof(int),                        &ldb,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_clag2z_quark = PCORE_clag2z_quark
#define CORE_clag2z_quark PCORE_clag2z_quark
#endif
void CORE_clag2z_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex32_t *A;
    int lda;
    PLASMA_Complex64_t *B;
    int ldb;

    quark_unpack_args_6(quark, m, n, A, lda, B, ldb);
    LAPACKE_clag2z_work(LAPACK_COL_MAJOR, m, n, A, lda, B, ldb);
}

