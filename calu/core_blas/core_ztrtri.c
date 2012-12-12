/**
 *
 * @file core_ztrtri.c
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
#include <lapacke.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_ztrtri = PCORE_ztrtri
#define CORE_ztrtri PCORE_ztrtri
#endif
void CORE_ztrtri(PLASMA_enum uplo, PLASMA_enum diag, int N, PLASMA_Complex64_t *A, int LDA, int *info)
{
    *info = LAPACKE_ztrtri_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo), lapack_const(diag),
        N, A, LDA);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_ztrtri(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, PLASMA_enum diag,
                       int n, int nb,
                       PLASMA_Complex64_t *A, int lda,
                       PLASMA_sequence *sequence, PLASMA_request *request,
                       int iinfo)
{
    QUARK_Insert_Task(
        quark, CORE_ztrtri_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(PLASMA_enum),                &diag,      VALUE,
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
#pragma weak CORE_ztrtri_quark = PCORE_ztrtri_quark
#define CORE_ztrtri_quark PCORE_ztrtri_quark
#endif
void CORE_ztrtri_quark(Quark *quark)
{
    PLASMA_enum uplo;
    PLASMA_enum diag;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    PLASMA_sequence *sequence;
    PLASMA_request *request;
    int iinfo;

    int info;

    quark_unpack_args_8(quark, uplo, diag, N, A, LDA, sequence, request, iinfo);
    info = LAPACKE_ztrtri_work(
        LAPACK_COL_MAJOR,
        lapack_const(uplo), lapack_const(diag),
        N, A, LDA);
    if ((sequence->status == PLASMA_SUCCESS) && (info > 0))
        plasma_sequence_flush(quark, sequence, request, iinfo + info);
}
