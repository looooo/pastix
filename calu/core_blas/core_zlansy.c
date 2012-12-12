/**
 *
 * @file core_zlansy.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
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
#pragma weak CORE_zlansy = PCORE_zlansy
#define CORE_zlansy PCORE_zlansy
#endif
void CORE_zlansy(int norm, PLASMA_enum uplo, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                 double *work, double *normA)
{
    *normA = LAPACKE_zlansy_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm), lapack_const(uplo),
        N, A, LDA, work);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlansy(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, PLASMA_enum uplo, int N,
                       const PLASMA_Complex64_t *A, int LDA, int szeA,
                       int szeW, double *result)
{
    szeW = max(1, szeW);
    DAG_CORE_LANSY;
    QUARK_Insert_Task(quark, CORE_zlansy_quark, task_flags,
        sizeof(PLASMA_enum),                &norm,  VALUE,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(PLASMA_Complex64_t)*szeA,     A,             INPUT,
        sizeof(int),                        &LDA,   VALUE,
        sizeof(double)*szeW,                 NULL,          SCRATCH,
        sizeof(double),                      result,        OUTPUT,
        0);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlansy_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, PLASMA_enum uplo, int N,
                          const PLASMA_Complex64_t *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF)
{
    szeW = max(1, szeW);
    DAG_CORE_LANSY;

    if ( result == fake ) {
        QUARK_Insert_Task(quark, CORE_zlansy_quark, task_flags,
            sizeof(PLASMA_enum),                &norm,  VALUE,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(int),                        &N,     VALUE,
            sizeof(PLASMA_Complex64_t)*szeA,     A,             INPUT,
            sizeof(int),                        &LDA,   VALUE,
            sizeof(double)*szeW,                 NULL,          SCRATCH,
            sizeof(double)*szeF,                 result,        OUTPUT | GATHERV,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_zlansy_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &norm,  VALUE,
            sizeof(PLASMA_enum),                &uplo,  VALUE,
            sizeof(int),                        &N,     VALUE,
            sizeof(PLASMA_Complex64_t)*szeA,     A,             INPUT,
            sizeof(int),                        &LDA,   VALUE,
            sizeof(double)*szeW,                 NULL,          SCRATCH,
            sizeof(double),                      result,        OUTPUT,
            sizeof(double)*szeF,                 fake,          OUTPUT | GATHERV,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlansy_quark = PCORE_zlansy_quark
#define CORE_zlansy_quark PCORE_zlansy_quark
#endif
void CORE_zlansy_quark(Quark *quark)
{
    double *normA;
    int norm;
    PLASMA_enum uplo;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    double *work;

    quark_unpack_args_7(quark, norm, uplo, N, A, LDA, work, normA);
    *normA = LAPACKE_zlansy_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm), lapack_const(uplo),
        N, A, LDA, work);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlansy_f1_quark = PCORE_zlansy_f1_quark
#define CORE_zlansy_f1_quark PCORE_zlansy_f1_quark
#endif
void CORE_zlansy_f1_quark(Quark *quark)
{
    double *normA;
    int norm;
    PLASMA_enum uplo;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    double *work;
    double *fake;

    quark_unpack_args_8(quark, norm, uplo, N, A, LDA, work, normA, fake);
    *normA = LAPACKE_zlansy_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm),
        lapack_const(uplo),
        N, A, LDA, work);
}
