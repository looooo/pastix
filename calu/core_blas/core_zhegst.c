/**
 *
 * @file core_zhegst.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"

#undef REAL
#define COMPLEX

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zhegst = PCORE_zhegst
#define CORE_zhegst PCORE_zhegst
#endif
void CORE_zhegst(int itype, PLASMA_enum uplo, int N,
                 PLASMA_Complex64_t *A, int LDA,
                 PLASMA_Complex64_t *B, int LDB, int *INFO)
{
    *INFO = LAPACKE_zhegst_work(
        LAPACK_COL_MAJOR,
        itype,
        lapack_const(uplo),
        N, A, LDA, B, LDB );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zhegst(Quark *quark, Quark_Task_Flags *task_flags,
                       int itype, PLASMA_enum uplo, int n,
                       PLASMA_Complex64_t *A, int lda,
                       PLASMA_Complex64_t *B, int ldb,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
    QUARK_Insert_Task(quark, CORE_zhegst_quark, task_flags,
        sizeof(int),                        &itype,      VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(PLASMA_Complex64_t)*lda*n,    A,             INOUT,
        sizeof(int),                        &lda,       VALUE,
#ifdef COMPLEX
        sizeof(PLASMA_Complex64_t)*ldb*n,    B,             INOUT,
#else
        sizeof(PLASMA_Complex64_t)*ldb*n,    B,             INPUT,
#endif
        sizeof(int),                        &ldb,       VALUE,
        sizeof(PLASMA_sequence*),           &sequence,  VALUE,
        sizeof(PLASMA_request*),            &request,   VALUE,
        sizeof(int),                        &iinfo,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zhegst_quark = PCORE_zhegst_quark
#define CORE_zhegst_quark PCORE_zhegst_quark
#endif
void CORE_zhegst_quark(Quark *quark)
{
    int itype;
    PLASMA_enum uplo;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_Complex64_t *B;
    int ldb;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int iinfo;

    int info;

    quark_unpack_args_10(quark, itype, uplo, n, A, lda, B, ldb, sequence, request, iinfo);
    info = LAPACKE_zhegst_work(
        LAPACK_COL_MAJOR,
        itype,
        lapack_const(uplo),
        n, A, lda, B, ldb);
    if (sequence->status == PLASMA_SUCCESS && info != 0)
      plasma_sequence_flush(quark, sequence, request, iinfo+info);
}
