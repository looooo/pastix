/**
 *
 * @file core_zlaset2.c
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

/***************************************************************************/
/**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zlaset2 - Sets the elements of the matrix A to alpha.
 *  Not LAPACK compliant! Read below.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies which elements of the matrix are to be set
 *          = PlasmaUpper: STRICT Upper part of A is set to alpha;
 *          = PlasmaLower: STRICT Lower part of A is set to alpha;
 *          = PlasmaUpperLower: ALL elements of A are set to alpha.
 *          Not LAPACK Compliant.
 *
 * @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A.  N >= 0.
 *
 * @param[in] alpha
 *         The constant to which the elements are to be set.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, A has been set to alpha accordingly.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaset2 = PCORE_zlaset2
#define CORE_zlaset2 PCORE_zlaset2
#endif
void CORE_zlaset2(PLASMA_enum uplo, int M, int N,
                  PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA)
{
    if (uplo == PlasmaUpper) {
        LAPACKE_zlaset_work(
            LAPACK_COL_MAJOR,
            lapack_const(uplo),
            M, N-1, alpha, alpha, A+LDA, LDA);
    }
    else if (uplo == PlasmaLower) {
        LAPACKE_zlaset_work(
            LAPACK_COL_MAJOR,
            lapack_const(uplo),
            M-1, N, alpha, alpha, A+1, LDA);
    }
    else {
        LAPACKE_zlaset_work(
            LAPACK_COL_MAJOR,
            lapack_const(uplo),
            M, N, alpha, alpha, A, LDA);
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaset2(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_enum uplo, int M, int N,
                       PLASMA_Complex64_t alpha, PLASMA_Complex64_t *A, int LDA)
{
    DAG_CORE_LASET;
    QUARK_Insert_Task(quark, CORE_zlaset2_quark, task_flags,
        sizeof(PLASMA_enum),                &uplo,  VALUE,
        sizeof(int),                        &M,     VALUE,
        sizeof(int),                        &N,     VALUE,
        sizeof(PLASMA_Complex64_t),         &alpha, VALUE,
        sizeof(PLASMA_Complex64_t)*M*N,     A,      OUTPUT,
        sizeof(int),                        &LDA,   VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaset2_quark = PCORE_zlaset2_quark
#define CORE_zlaset2_quark PCORE_zlaset2_quark
#endif
void CORE_zlaset2_quark(Quark *quark)
{
    PLASMA_enum uplo;
    int M;
    int N;
    PLASMA_Complex64_t alpha;
    PLASMA_Complex64_t *A;
    int LDA;

    quark_unpack_args_6(quark, uplo, M, N, alpha, A, LDA);
    CORE_zlaset2(uplo, M, N, alpha, A, LDA);
}
