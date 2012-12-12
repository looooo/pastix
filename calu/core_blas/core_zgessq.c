/**
 *
 * @file core_zgessq.c
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
#include <math.h>
#include <lapacke.h>
#include "common.h"

#define COMPLEX

#define UPDATE( __nb, __value )                                         \
    __value = fabs(*ptr);                                               \
    if (__value != 0. ){                                                \
        if ( *scale < __value ) {                                       \
            *sumsq = __nb + (*sumsq) * ( *scale / __value ) * ( *scale / __value ); \
            *scale = __value;                                           \
        } else {                                                        \
            *sumsq = *sumsq + __nb * ( __value / *scale ) *  ( __value / *scale ); \
        }                                                               \
    }

/*****************************************************************************
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zgessq returns the values scl and ssq such that
 *
 *    ( scl**2 )*ssq = sum( A( i, j )**2 ) + ( scale**2 )*sumsq,
 *                     i,j
 *
 * with i from 0 to M-1 and j form 0 to N-1. The value of sumsq is
 * assumed to be at least unity and the value of ssq will then satisfy
 *
 *    1.0 .le. ssq .le. ( sumsq + 2*m*n ).
 *
 * scale is assumed to be non-negative and scl returns the value
 *
 *    scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
 *           i
 *
 * scale and sumsq must be supplied in SCALE and SUMSQ respectively.
 * SCALE and SUMSQ are overwritten by scl and ssq respectively.
 *
 * The routine makes only one pass through the tile A.
 * See also LAPACK zlassq.f
 *
 *******************************************************************************
 *
 *  @param[in] M
 *          The number of rows in the tile A.
 *
 *  @param[in] N
 *          The number of columns in the tile A.
 *
 *  @param[in] A
 *          The M-by-N matrix on which to compute the norm.
 *
 *  @param[in] LDA
 *          The leading dimension of the tile A. LDA >= max(1,M).
 *
 *  @param[in,out] scale
 *          On entry, the value  scale  in the equation above.
 *          On exit, scale is overwritten with the value scl.
 *
 *  @param[in,out] sumsq
 *          On entry, the value  sumsq  in the equation above.
 *          On exit, SUMSQ is overwritten with the value ssq.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PLASMA_SUCCESS successful exit
 *          \retval -k, the k-th argument had an illegal value
 *
 */
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgessq = PCORE_zgessq
#define CORE_zgessq PCORE_zgessq
#endif
int CORE_zgessq(int M, int N,
                const PLASMA_Complex64_t *A, int LDA,
                double *scale, double *sumsq)
{
    int i, j;
    double tmp;
    double *ptr;

    for(j=0; j<N; j++) {
        ptr = (double*) ( A + j * LDA );
        for(i=0; i<M; i++, ptr++) {
            tmp = fabs(*ptr);
            UPDATE( 1., tmp );

#ifdef COMPLEX
            ptr++;
            tmp = fabs(*ptr);
            UPDATE( 1., tmp );
#endif
        }
    }
    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgessq_f1( Quark *quark, Quark_Task_Flags *task_flags,
                           int m, int n, const PLASMA_Complex64_t *A, int lda,
                           double *scale, double *sumsq,
                           double *fake, int szeF, int paramF )
{
    DAG_CORE_LASSQ;
    if ( (fake == scale) && (paramF & GATHERV) ) {
        QUARK_Insert_Task(quark, CORE_zgessq_quark, task_flags,
            sizeof(int),                      &m,    VALUE,
            sizeof(int),                      &n,    VALUE,
            sizeof(PLASMA_Complex64_t)*lda*n, A,         INPUT,
            sizeof(int),                      &lda,  VALUE,
            sizeof(double)*1,                 scale,     INOUT | GATHERV,
            sizeof(double)*1,                 sumsq,     INOUT,
            0);
    } else {
        QUARK_Insert_Task(quark, CORE_zgessq_f1_quark, task_flags,
            sizeof(int),                      &m,    VALUE,
            sizeof(int),                      &n,    VALUE,
            sizeof(PLASMA_Complex64_t)*lda*n, A,         INPUT,
            sizeof(int),                      &lda,  VALUE,
            sizeof(double)*1,                 scale,     INOUT,
            sizeof(double)*1,                 sumsq,     INOUT,
            sizeof(double)*szeF,              fake,      paramF,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgessq_quark = PCORE_zgessq_quark
#define CORE_zgessq_quark PCORE_zgessq_quark
#endif
void CORE_zgessq_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    double *scale;
    double *sumsq;

    quark_unpack_args_6( quark, m, n, A, lda, scale, sumsq );
    CORE_zgessq( m, n, A, lda, scale, sumsq );
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgessq_f1_quark = PCORE_zgessq_f1_quark
#define CORE_zgessq_f1_quark PCORE_zgessq_f1_quark
#endif
void CORE_zgessq_f1_quark(Quark *quark)
{
    int m;
    int n;
    PLASMA_Complex64_t *A;
    int lda;
    double *scale;
    double *sumsq;
    double *fake;

    quark_unpack_args_7( quark, m, n, A, lda, scale, sumsq, fake );
    CORE_zgessq( m, n, A, lda, scale, sumsq );
}
