/**
 *
 * @file core_dzasum.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <cblas.h>
#include <math.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dzasum = PCORE_dzasum
#define CORE_dzasum PCORE_dzasum
#endif
void CORE_dzasum(PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                 const PLASMA_Complex64_t *A, int lda, double *work)
{
    const PLASMA_Complex64_t *tmpA;
    double *tmpW, sum, abs;
    int i,j;

    switch (uplo) {
    case PlasmaUpper:
        for (j = 0; j < N; j++) {
            tmpA = A+(j*lda);
            sum = 0.0;
            for (i = 0; i < j; i++) {
                abs      = cabs(*tmpA);
                sum     += abs;
                work[i] += abs;
                tmpA++;
            }
            work[j] += sum + cabs(*tmpA);
        }
        break;
    case PlasmaLower:
        for (j = 0; j < N; j++) {
            tmpA = A+(j*lda)+j;

            sum = 0.0;
            work[j] += cabs(*tmpA);

            tmpA++;
            for (i = j+1; i < M; i++) {
                abs      = cabs(*tmpA);
                sum     += abs;
                work[i] += abs;
                tmpA++;
            }
            work[j] += sum;
        }
        break;
    case PlasmaUpperLower:
    default:
        if (storev == PlasmaColumnwise) {
            for (j = 0; j < N; j++) {
                /* work[j] += cblas_dzasum(M, &(A[j*lda]), 1); */
                tmpA = A+(j*lda);
                for (i = 0; i < M; i++) {
                    work[j] +=  cabs(*tmpA);
                    tmpA++;
                }
            }
        }
        else {
            for (j = 0; j < N; j++) {
                tmpA = A+(j*lda);
                tmpW = work;
                for (i = 0; i < M; i++) {
                    /* work[i] += cabs( A[j*lda+i] );*/
                    *tmpW += cabs( *tmpA );
                    tmpA++; tmpW++;
                }
            }
        }
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dzasum(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                       const PLASMA_Complex64_t *A, int lda, int szeA,
                       double *work, int szeW)
{
    QUARK_Insert_Task(
        quark, CORE_dzasum_quark, task_flags,
        sizeof(PLASMA_enum),                &storev,    VALUE,
        sizeof(PLASMA_enum),                &uplo,      VALUE,
        sizeof(int),                        &M,         VALUE,
        sizeof(int),                        &N,         VALUE,
        sizeof(PLASMA_Complex64_t)*szeA,     A,                 INPUT,
        sizeof(int),                        &lda,       VALUE,
        sizeof(double)*szeW,                 work,              INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dzasum_quark = PCORE_dzasum_quark
#define CORE_dzasum_quark PCORE_dzasum_quark
#endif
void CORE_dzasum_quark(Quark *quark)
{
    PLASMA_enum storev;
    PLASMA_enum uplo;
    int M;
    int N;
    PLASMA_Complex64_t *A;
    int lda;
    double *work;

    quark_unpack_args_7(quark, storev, uplo, M, N, A, lda, work);
    CORE_dzasum(storev, uplo, M, N, A, lda, work);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_dzasum_f1(Quark *quark, Quark_Task_Flags *task_flags,
                          PLASMA_enum storev, PLASMA_enum uplo, int M, int N,
                          const PLASMA_Complex64_t *A, int lda, int szeA,
                          double *work, int szeW, double *fake, int szeF)
{
    DAG_CORE_ASUM;
    if ( work == fake ) {
        QUARK_Insert_Task(
            quark, CORE_dzasum_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &storev,    VALUE,
            sizeof(PLASMA_enum),                &uplo,      VALUE,
            sizeof(int),                        &M,         VALUE,
            sizeof(int),                        &N,         VALUE,
            sizeof(PLASMA_Complex64_t)*szeA,     A,                 INPUT,
            sizeof(int),                        &lda,       VALUE,
            sizeof(double)*szeW,                 work,              INOUT | GATHERV,
            sizeof(double)*1,                    NULL,              SCRATCH,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_dzasum_f1_quark, task_flags,
            sizeof(PLASMA_enum),                &storev,    VALUE,
            sizeof(PLASMA_enum),                &uplo,      VALUE,
            sizeof(int),                        &M,         VALUE,
            sizeof(int),                        &N,         VALUE,
            sizeof(PLASMA_Complex64_t)*szeA,     A,                 INPUT,
            sizeof(int),                        &lda,       VALUE,
            sizeof(double)*szeW,                 work,              INOUT,
            sizeof(double)*szeF,                 fake,              OUTPUT | GATHERV,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_dzasum_f1_quark = PCORE_dzasum_f1_quark
#define CORE_dzasum_f1_quark PCORE_dzasum_f1_quark
#endif
void CORE_dzasum_f1_quark(Quark *quark)
{
    PLASMA_enum storev;
    PLASMA_enum uplo;
    int M;
    int N;
    PLASMA_Complex64_t *A;
    int lda;
    double *work;
    double *fake;

    quark_unpack_args_8(quark, storev, uplo, M, N, A, lda, work, fake);
    CORE_dzasum(storev, uplo, M, N, A, lda, work);
}
