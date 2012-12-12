/**
 *
 * @file core_zlacpy_pivot.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.5
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <lapacke.h>
#include "common.h"

#define A(m, n) BLKADDR(descA, PLASMA_Complex64_t, m, n)

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlacpy_pivot = PCORE_zlacpy_pivot
#define CORE_zlacpy_pivot PCORE_zlacpy_pivot
#endif
int CORE_zlacpy_pivot( const PLASMA_desc descA,
                       int k1, int k2, int *ipiv,
                       int *rankin, int *rankout,
                       PLASMA_Complex64_t *A, int lda,
                       int init)
{
    int minMN = min(descA.m, descA.n);
    int i, ip, it, ir, ld, end;
    int *ro = rankout;

    /* Init rankin if first step */
    if ( init ) {
        int val = descA.i;
        for(i=0; i<descA.m; i++, val++) {
            rankin[i] = val;
        }
    }

    /* Generate the rankout */
    ro = rankout;
    ipiv = ipiv;
    for(i=k1-1; i<k2; i++, ro++, ipiv++) {
        *ro = rankin[ *ipiv - 1 ];
        rankin[ *ipiv - 1 ] = rankin[ i ];
    }

    /* Extract the rows */
    ro = rankout;
    for(i=k1-1; i<k2; i++, ro++) {
        ip = (*ro) - descA.i;
        it = ip / descA.mb;
        ir = ip % descA.mb;
        ld = BLKLDD(descA, it);
        cblas_zcopy(descA.n, A(it, 0) + ir, ld,
                             A + i,  lda );
    }

    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlacpy_pivot(Quark *quark, Quark_Task_Flags *task_flags,
                             PLASMA_desc descA,
                             int k1, int k2, int *ipiv,
                             int *rankin, int *rankout,
                             PLASMA_Complex64_t *A, int lda,
                             int pos, int init)
{
    DAG_SET_PROPERTIES( "CPY_PIV"  , "white"   );
    QUARK_Insert_Task(quark, CORE_zlacpy_pivot_quark, task_flags,
        sizeof(PLASMA_desc),                    &descA,         VALUE,
        sizeof(int),                            &k1,            VALUE,
        sizeof(int),                            &k2,            VALUE,
        sizeof(int)*lda,                         ipiv,                INPUT,
        sizeof(int)*lda,                         rankin,              INPUT,
        sizeof(int)*lda,                         rankout,             OUTPUT | GATHERV,
        sizeof(PLASMA_Complex64_t)*lda*descA.nb, A,                   INOUT | GATHERV,
        sizeof(int),                            &lda,           VALUE,
        sizeof(int),                            &pos,           VALUE,
        sizeof(int),                            &init,          VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlacpy_pivot_quark = PCORE_zlacpy_pivot_quark
#define CORE_zlacpy_pivot_quark PCORE_zlacpy_pivot_quark
#endif
void CORE_zlacpy_pivot_quark(Quark *quark)
{
    PLASMA_desc descA;
    PLASMA_Complex64_t *A;
    int lda, pos, k1, k2;
    int *rankin;
    int *rankout;
    int *ipiv;
    int init;

    quark_unpack_args_10(quark, descA, k1, k2, ipiv, rankin, rankout, A, lda, pos, init);
    CORE_zlacpy_pivot(descA, k1, k2, ipiv, rankin, rankout+pos, A+pos, lda, init );
}
