/**
 *
 * @file core_zpotrf.c
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
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zpotrf = PCORE_zpotrf
#define CORE_zpotrf PCORE_zpotrf
#endif
void CORE_zpotrf(PLASMA_enum uplo, int N, PLASMA_Complex64_t *A, int LDA, int *INFO)
{
    *INFO = LAPACKE_zpotrf_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        N, A, LDA );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zpotrf(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int n, int nb,
                       PLASMA_Complex64_t *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
    DAG_CORE_POTRF;
    QUARK_Insert_Task(quark, CORE_zpotrf_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &n,         VALUE,
        sizeof(PLASMA_Complex64_t)*nb*nb,    A,                 INOUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(PLASMA_sequence*),           &sequence,  VALUE,
        sizeof(PLASMA_request*),            &request,   VALUE,
        sizeof(int),                        &iinfo,     VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zpotrf_quark = PCORE_zpotrf_quark
#define CORE_zpotrf_quark PCORE_zpotrf_quark
#endif
void CORE_zpotrf_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int iinfo;

    int info;

    quark_unpack_args_7(quark, uplo, n, A, lda, sequence, request, iinfo);
    info = LAPACKE_zpotrf_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo),
        n, A, lda);
    if (sequence->status == PLASMA_SUCCESS && info != 0)
      plasma_sequence_flush(quark, sequence, request, iinfo+info);
}
