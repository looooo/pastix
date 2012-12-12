/**
 *
 * @file core_zlange.c
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
#pragma weak CORE_zlange = PCORE_zlange
#define CORE_zlange PCORE_zlange
#endif
void CORE_zlange(int norm, int M, int N,
                 const PLASMA_Complex64_t *A, int LDA,
                 double *work, double *normA)
{
    *normA = LAPACKE_zlange_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm),
        M, N, A, LDA, work);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlange(Quark *quark, Quark_Task_Flags *task_flags,
                       int norm, int M, int N,
                       const PLASMA_Complex64_t *A, int LDA, int szeA,
                       int szeW, double *result)
{
    szeW = max(1, szeW);
    DAG_CORE_LANGE;
    QUARK_Insert_Task(quark, CORE_zlange_quark, task_flags,
        sizeof(PLASMA_enum),                &norm,  VALUE,
        sizeof(int),                        &M,     VALUE,
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
void QUARK_CORE_zlange_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum norm, int M, int N,
                          const PLASMA_Complex64_t *A, int LDA, int szeA,
                          int szeW, double *result,
                          double *fake, int szeF)
{
    szeW = max(1, szeW);
    DAG_CORE_LANGE;

    if ( result == fake ) {
        QUARK_Insert_Task(quark, CORE_zlange_quark, task_flags,
            sizeof(PLASMA_enum),                &norm,  VALUE,
            sizeof(int),                        &M,     VALUE,
            sizeof(int),                        &N,     VALUE,
            sizeof(PLASMA_Complex64_t)*szeA,     A,             INPUT,
            sizeof(int),                        &LDA,   VALUE,
            sizeof(double)*szeW,                 NULL,          SCRATCH,
            sizeof(double),                      result,        OUTPUT | GATHERV,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_zlange_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &norm,  VALUE,
            sizeof(int),                        &M,     VALUE,
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
#pragma weak CORE_zlange_quark = PCORE_zlange_quark
#define CORE_zlange_quark PCORE_zlange_quark
#endif
void CORE_zlange_quark(Quark *quark)
{
    double *normA;
    int norm;
    int M;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    double *work;

    quark_unpack_args_7(quark, norm, M, N, A, LDA, work, normA);
    *normA = LAPACKE_zlange_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm),
        M, N, A, LDA, work);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlange_f1_quark = PCORE_zlange_f1_quark
#define CORE_zlange_f1_quark PCORE_zlange_f1_quark
#endif
void CORE_zlange_f1_quark(Quark *quark)
{
    double *normA;
    int norm;
    int M;
    int N;
    PLASMA_Complex64_t *A;
    int LDA;
    double *work;
    double *fake;

    quark_unpack_args_8(quark, norm, M, N, A, LDA, work, normA, fake);
    *normA = LAPACKE_zlange_work(
        LAPACK_COL_MAJOR,
        lapack_const(norm),
        M, N, A, LDA, work);
}

