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
    int i, ip, it, ir, ld;
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

#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaepv = PCORE_zlaepv
#define CORE_zlaepv PCORE_zlaepv
#endif
int CORE_zlaepv( const PLASMA_desc descA,
                 int k1, int k2, int *ipiv, int *rankin,
                 PLASMA_Complex64_t *B, int ldb)
{
    int i, ip, it, ir, ld;

    for(i=k1-1; i<k2; i++, piv++) {
        ip = (*piv)-1;
        *piv = rankin[ ip ];
        rankin[ ip ] = rankin[ i ];

        ip = (*piv) - 1 - A.i;
        it = ip / descA.mb;
        ir = ip % descA.mb;
        ld = BLKLDD(descA, it);
        cblas_zcopy(descA.n, A(it, 0) + ir, ld,
                             B + i,         ldb );
    }

    return PLASMA_SUCCESS;
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaepv(Quark *quark, Quark_Task_Flags *task_flags,
                       PLASMA_desc descA,
                       int k1, int k2, int *ipiv, int *rankin,
                       PLASMA_Complex64_t *B, int ldb)
{
    Quark_Task *task;
    int m;

    DAG_SET_PROPERTIES( "CPY_PIV"  , "white"   );
    task = QUARK_Task_Init( quark, CORE_zlaepv_quark, task_flags );
    QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_desc),                     &descA, VALUE);
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),                             &k1,    VALUE);
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),                             &k2,    VALUE);
    QUARK_Task_Pack_Arg( quark, task, sizeof(int)*descA.n,                      ipiv,      INOUT );
    QUARK_Task_Pack_Arg( quark, task, sizeof(int)*descA.m,                      rankin,    INOUT );
    QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_Complex64_t)*ldb*descA.nb,  B,         OUTPUT);
    QUARK_Task_Pack_Arg( quark, task, sizeof(int),                             &ldb,   VALUE);

    for( m = 0; m < descA.mt; m++ ) {
        QUARK_Task_Pack_Arg( quark, task, sizeof(PLASMA_Complex64_t)*descA.mb*descA.nb, A(m,0), INPUT );
    }
    QUARK_Insert_Task_Packed( quark, task );
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaepv_quark = PCORE_zlaepv_quark
#define CORE_zlaepv_quark PCORE_zlaepv_quark
#endif
void CORE_zlaepv_quark(Quark *quark)
{
    PLASMA_desc descA;
    int k1, k2;
    int *ipiv;
    PLASMA_Complex64_t *B;
    int ldb;

    quark_unpack_args_6(quark, descA, k1, k2, ipiv, rankin, B, ldb);
    CORE_zlaepv(descA, k1, k2, ipiv, rankin, B, ldb );
}
