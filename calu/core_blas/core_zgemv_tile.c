/**
 *
 * @file core_zgemv_tile.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Mark Gates
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <cblas.h>
#include "common.h"

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 * CORE_zgemv_tile performs one of the matrix-vector operations
 *
 *    y := alpha*A*x    + beta*y,   or
 *    y := alpha*A**T*x + beta*y,   or
 *    y := alpha*A**H*x + beta*y,
 *
 * where alpha and beta are scalars, x and y are vectors,
 * and A is an m by n matrix.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *         @arg PlasmaNoTrans:   y := alpha*A*x    + beta*y
 *         @arg PlasmaTrans:     y := alpha*A**T*x + beta*y
 *         @arg PlasmaConjTrans: y := alpha*A**H*x + beta*y
 *
 * @param[in] m
 *         Number of rows of matrix A.
 *
 * @param[in] n
 *         Number of columns of matrix A.
 *
 * @param[in] alpha
 *         Scalar alpha, passed by pointer so it can depend on runtime value.
 *
 * @param[in] A
 *         On entry, m by n matrix A. Dimension (lda,n).
 *
 * @param[in] lda
 *         Leading dimension of array A. lda >= max(1,m).
 *
 * @param[in] x
 *         On entry, vector x.
 *         If trans == PlasmaNoTrans, the n vector x has dimension 1 + (n-1)*abs(incx).
 *         Else, the m vector x has dimension 1 + (m-1)*abs(incx).
 *
 * @param[in] incx
 *         Increment between elements of x. incx must not be zero.
 *
 * @param[in] beta
 *         Scalar beta, passed by pointer so it can depend on runtime value.
 *
 * @param[in,out] y
 *         On entry, vector y. On exit, y  is overwritten by updated vector y.
 *         If trans == PlasmaNoTrans, the m vector y has dimension 1 + (m-1)*abs(incy).
 *         Else, the n vector y has dimension 1 + (n-1)*abs(incy).
 *
 * @param[in] incy
 *         Increment between elements of y. incy must not be zero.
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgemv_tile = PCORE_zgemv_tile
#define CORE_zgemv_tile PCORE_zgemv_tile
#endif
void CORE_zgemv_tile(PLASMA_enum trans,
                     int m, int n,
                     const PLASMA_Complex64_t *alpha, const PLASMA_Complex64_t *A, int lda,
                                                      const PLASMA_Complex64_t *x, int incx,
                     const PLASMA_Complex64_t *beta,        PLASMA_Complex64_t *y, int incy )
{
    cblas_zgemv(
        CblasColMajor,
        (CBLAS_TRANSPOSE)trans,
        m, n,
        CBLAS_SADDR(*alpha), A, lda,
                             x, incx,
        CBLAS_SADDR(*beta),  y, incy );
}

/***************************************************************************//**
 *
 * Version of zgemv for tile storage, to avoid dependency problem when
 * computations are done within the tile. alpha and beta are passed as
 * pointers so they can depend on runtime values.
 *
 * @param[in] Alock
 *          Pointer to tile owning submatrix A.
 *
 * @param[in] xlock
 *          Pointer to tile owning subvector x.
 *
 * @param[in] ylock
 *          Pointer to tile owning subvector y.
 *
 **/
void QUARK_CORE_zgemv_tile(Quark *quark, Quark_Task_Flags *task_flags,
                           PLASMA_enum trans,
                           int m, int n,
                           const PLASMA_Complex64_t *alpha, const PLASMA_Complex64_t *A, int lda,
                                                            const PLASMA_Complex64_t *x, int incx,
                           const PLASMA_Complex64_t *beta,        PLASMA_Complex64_t *y, int incy,
                           const PLASMA_Complex64_t *Alock,
                           const PLASMA_Complex64_t *xlock,
                           const PLASMA_Complex64_t *ylock)
{
    /* Quick return. Bad things happen if sizeof(...)*m*n is zero in QUARK_Insert_Task */
    if ( m == 0 || n == 0 )
        return;

    DAG_SET_PROPERTIES("gemv", "lightslateblue");
    QUARK_Insert_Task(quark, CORE_zgemv_tile_quark, task_flags,
        sizeof(PLASMA_enum),             &trans,  VALUE,
        sizeof(int),                     &m,      VALUE,
        sizeof(int),                     &n,      VALUE,
        sizeof(PLASMA_Complex64_t),      alpha,           INPUT,
        sizeof(PLASMA_Complex64_t)*m*n,  A,               NODEP,          /* input; see Alock */
        sizeof(int),                     &lda,    VALUE,
        sizeof(PLASMA_Complex64_t)*n,    x,               NODEP,          /* input; see xlock */
        sizeof(int),                     &incx,   VALUE,
        sizeof(PLASMA_Complex64_t),      beta,            INPUT,
        sizeof(PLASMA_Complex64_t)*m,    y,                       NODEP,  /* inout; see ylock */
        sizeof(int),                     &incy,   VALUE,
        sizeof(PLASMA_Complex64_t)*m*n,  Alock,           INPUT,
        sizeof(PLASMA_Complex64_t)*n,    xlock,           INPUT,
        sizeof(PLASMA_Complex64_t)*m,    ylock,                   INOUT,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgemv_tile_quark = PCORE_zgemv_tile_quark
#define CORE_zgemv_tile_quark PCORE_zgemv_tile_quark
#endif
void CORE_zgemv_tile_quark(Quark *quark)
{
    PLASMA_enum trans;
    int m, n, lda, incx, incy;
    const PLASMA_Complex64_t *alpha, *beta;
    const PLASMA_Complex64_t *A, *x;
    PLASMA_Complex64_t *y;

    quark_unpack_args_11( quark, trans, m, n, alpha, A, lda, x, incx, beta, y, incy );
    cblas_zgemv(
        CblasColMajor,
        (CBLAS_TRANSPOSE)trans,
        m, n,
        CBLAS_SADDR(*alpha), A, lda,
                             x, incx,
        CBLAS_SADDR(*beta),  y, incy );
}
