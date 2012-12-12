/**
 *
 * @file core_zgeqp3_norms.c
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
#include <math.h>
#include <cblas.h>
#include "common.h"

#define A(m,n) BLKADDR( A, PLASMA_Complex64_t, m, n )

/* Update scale and sumsq during accumulation of (Aij)^2 */
#define UPDATE( __ptr, __scale, __sumsq, __nb ) \
    {                                                                   \
        double value = fabs(*(__ptr));                                  \
        if ( value != 0. ){                                             \
            if ( *(__scale) < value ) {                                 \
                *(__sumsq) = __nb + (*(__sumsq)) * ( *(__scale) / value ) * ( *(__scale) / value ); \
                *(__scale) = value;                                     \
            } else {                                                    \
                *(__sumsq) = *(__sumsq) + __nb * ( value / *(__scale) ) *  ( value / *(__scale) ); \
            }                                                           \
        }                                                               \
    }

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zgeqp3_norms computes frobenius norm of each column of A. The result is
 *  stored duplicated in norms1 and norms2. The algorithm is based on the
 *  function used in CORE_zgessq to avoid overflow during computation of the
 *  square of each element. See also xlassq functions in LAPACK.
 *
 *******************************************************************************
 *
 *  @param[in] A
 *          PLASMA descriptor of the matrix A.
 *          On entry, the M-by-N matrix described by the descriptor.
 *
 *  @param[out] norms1
 *          Vector of size A.n.
 *          On exit, norms1[j] contains the frobenius norm of the column j from
 *          matrix A.
 *
 *  @param[out] norms2
 *          Vector of size A.n.
 *          On exit, norms2[j] contains the frobenius norm of the column j from
 *          matrix A.
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeqp3_norms = PCORE_zgeqp3_norms
#define CORE_zgeqp3_norms PCORE_zgeqp3_norms
#endif
int CORE_zgeqp3_norms( PLASMA_desc A, double *norms1, double *norms2 )
{
    const PLASMA_Complex64_t *Aij;
    const double *ptr;
    double *sumsq = norms1;
    double *scale = norms2;
    double tmp;
    int m, n, i, j, mb, nb, lda;

    /* Quick return */
    if ( A.n == 0 ) {
        coreblas_error(1, "Illegal value of A.nt");
        return PLASMA_SUCCESS;
    }

    for( j = 0; j < A.n; ++j ) {
        sumsq[j] = 0.;
        scale[j] = 1.;
    }

    for( n = 0; n < A.nt; n++ ) {
        nb  = n == A.nt-1 ? A.n - n * A.nb : A.nb;
        for( m = 0; m < A.mt; m++ ) {

            /* Restore scale and sumsq */
            scale = norms1 + n*A.nb;
            sumsq = norms2 + n*A.nb;

            /* Compute data for tile A(m, n)*/
            mb  = m == A.mt-1 ? A.m - m * A.mb : A.mb;
            Aij = A(m, n);
            lda = BLKLDD( A, m );

            /* Compute frobenius norm per column */
            for(j=0; j<nb; j++) {
                ptr = (double*) ( Aij + j * lda );

                for(i=0; i<mb; i++, ptr++) {
                    UPDATE( ptr, scale, sumsq, 1. );
#ifdef COMPLEX
                    ptr++;
                    UPDATE( ptr, scale, sumsq, 1. );
#endif
                }
                sumsq++;
                scale++;
            }
        }
    }

    /* Set norm1 and norm2 to the computed norm of each column */
    scale = norms1;
    sumsq = norms2;
    for( j=0; j < A.n; j++, scale++, sumsq++ ) {
        tmp = (*scale) * sqrt( *sumsq );
        *sumsq = tmp;
        *scale = tmp;
    }

    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zgeqp3_norms( Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc A,
                              double *norms1, double *norms2 )
{
    Quark_Task *task;
    int m, n;

    DAG_SET_PROPERTIES("norms", "brown");
    task = QUARK_Task_Init( quark, CORE_zgeqp3_norms_quark, task_flags );

    QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_desc),  &A,      VALUE          );
    QUARK_Task_Pack_Arg( quark, task, sizeof(double)*A.nb,  norms1,          OUTPUT );
    QUARK_Task_Pack_Arg( quark, task, sizeof(double)*A.nb,  norms2,          NODEP  );  /* OUTPUT, but implied by norms1 */

    /* depends on block column */
    for( n = 0; n < A.nt; ++n ) {
        for( m = 0; m < A.mt; ++m ) {
            QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_Complex64_t)*A.nb*A.nb, A(m,n), INPUT );
        }
    }

    QUARK_Insert_Task_Packed( quark, task );
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zgeqp3_norms_quark = PCORE_zgeqp3_norms_quark
#define CORE_zgeqp3_norms_quark PCORE_zgeqp3_norms_quark
#endif
void CORE_zgeqp3_norms_quark( Quark *quark )
{
    PLASMA_desc A;
    double *norms1, *norms2;

    quark_unpack_args_3( quark, A, norms1, norms2 );
    CORE_zgeqp3_norms(          A, norms1, norms2 );
}
