/**
 *
 * @file core_zlaswp.c
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
#include <lapacke.h>
#include "common.h"

#define A(m, n) BLKADDR(descA, PLASMA_Complex64_t, m, n)

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswp = PCORE_zlaswp
#define CORE_zlaswp PCORE_zlaswp
#endif
void CORE_zlaswp(int N, PLASMA_Complex64_t *A, int LDA, int I1, int I2, const int *IPIV, int INC)
{
    LAPACKE_zlaswp_work( LAPACK_COL_MAJOR, N, A, LDA, I1, I2, IPIV, INC );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaswp(Quark *quark, Quark_Task_Flags *task_flags,
                       int n, PLASMA_Complex64_t *A, int lda,
                       int i1,  int i2, const int *ipiv, int inc)
{
    DAG_CORE_LASWP;
    QUARK_Insert_Task(
        quark, CORE_zlaswp_quark, task_flags,
        sizeof(int),                      &n,    VALUE,
        sizeof(PLASMA_Complex64_t)*lda*n,  A,        INOUT | LOCALITY,
        sizeof(int),                      &lda,  VALUE,
        sizeof(int),                      &i1,   VALUE,
        sizeof(int),                      &i2,   VALUE,
        sizeof(int)*n,                     ipiv,     INPUT,
        sizeof(int),                      &inc,  VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswp_quark = PCORE_zlaswp_quark
#define CORE_zlaswp_quark PCORE_zlaswp_quark
#endif
void CORE_zlaswp_quark(Quark *quark)
{
    int n, lda, i1, i2, inc;
    int *ipiv;
    PLASMA_Complex64_t *A;

    quark_unpack_args_7(quark, n, A, lda, i1, i2, ipiv, inc);
    LAPACKE_zlaswp_work(LAPACK_COL_MAJOR, n, A, lda, i1, i2, ipiv, inc );
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaswp_f2(Quark *quark, Quark_Task_Flags *task_flags,
                          int n, PLASMA_Complex64_t *A, int lda,
                          int i1,  int i2, const int *ipiv, int inc,
                          PLASMA_Complex64_t *fake1, int szefake1, int flag1,
                          PLASMA_Complex64_t *fake2, int szefake2, int flag2)
{
    DAG_CORE_LASWP;
    QUARK_Insert_Task(
        quark, CORE_zlaswp_f2_quark, task_flags,
        sizeof(int),                        &n,     VALUE,
        sizeof(PLASMA_Complex64_t)*lda*n,    A,         INOUT | LOCALITY,
        sizeof(int),                        &lda,   VALUE,
        sizeof(int),                        &i1,    VALUE,
        sizeof(int),                        &i2,    VALUE,
        sizeof(int)*n,                       ipiv,      INPUT,
        sizeof(int),                        &inc,   VALUE,
        sizeof(PLASMA_Complex64_t)*szefake1, fake1,     flag1,
        sizeof(PLASMA_Complex64_t)*szefake2, fake2,     flag2,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswp_f2_quark = PCORE_zlaswp_f2_quark
#define CORE_zlaswp_f2_quark PCORE_zlaswp_f2_quark
#endif
void CORE_zlaswp_f2_quark(Quark* quark)
{
    int n, lda, i1, i2, inc;
    int *ipiv;
    PLASMA_Complex64_t *A;
    void *fake1, *fake2;

    quark_unpack_args_9(quark, n, A, lda, i1, i2, ipiv, inc, fake1, fake2);
    LAPACKE_zlaswp_work(LAPACK_COL_MAJOR, n, A, lda, i1, i2, ipiv, inc );
}

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zlaswp_ontile apply the zlaswp function on a matrix stored in
 *  tile layout
 *
 *******************************************************************************
 *
 *  @param[in,out] descA
 *          The descriptor of the matrix A to permute.
 *
 *  @param[in] i1
 *          The first element of IPIV for which a row interchange will
 *          be done.
 *
 *  @param[in] i2
 *          The last element of IPIV for which a row interchange will
 *          be done.
 *
 *  @param[in] ipiv
 *          The pivot indices; Only the element in position i1 to i2
 *          are accessed. The pivot are offset by A.i.
 *
 *  @param[in] inc
 *          The increment between successive values of IPIV.  If IPIV
 *          is negative, the pivots are applied in reverse order.
 *
 *******************************************************************************
 */
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswp_ontile = PCORE_zlaswp_ontile
#define CORE_zlaswp_ontile PCORE_zlaswp_ontile
#endif
int CORE_zlaswp_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc)
{
    int i, j, ip, it;
    PLASMA_Complex64_t *A1;
    int lda1, lda2;

    /* Change i1 to C notation */
    i1--;
    if ( descA.nt > 1 ) {
        coreblas_error(1, "Illegal value of descA.nt");
        return -1;
    }
    if ( i1 < 0 ) {
        coreblas_error(2, "Illegal value of i1");
        return -2;
    }
    if ( (i2 < i1) || (i2 > descA.m) ) {
        coreblas_error(3, "Illegal value of i2");
        return -3;
    }
    if ( ! ( (i2 - i1 - i1%descA.mb -1) < descA.mb ) ) {
        coreblas_error(2, "Illegal value of i1,i2. They have to be part of the same block.");
        return -3;
    }

    it = i1 / descA.mb;
    if (inc > 0) {
        A1 = A(it, 0);
        lda1 = BLKLDD(descA, 0);

        for (j = i1; j < i2; ++j, ipiv+=inc) {
            ip = (*ipiv) - descA.i - 1;
            if ( ip != j )
            {
                it = ip / descA.mb;
                i  = ip % descA.mb;
                lda2 = BLKLDD(descA, it);
                cblas_zswap(descA.n, A1       + j, lda1,
                                     A(it, 0) + i, lda2 );
            }
        }
    }
    else
    {
        A1 = A(it, 0);
        lda1 = BLKLDD(descA, descA.mt-1);

        i1--;
        ipiv = &ipiv[(1-i2)*inc];
        for (j = i2-1; j > i1; --j, ipiv+=inc) {
            ip = (*ipiv) - descA.i - 1;
            if ( ip != j )
            {
                it = ip / descA.mb;
                i  = ip % descA.mb;
                lda2 = BLKLDD(descA, it);
                cblas_zswap(descA.n, A1       + j, lda1,
                                     A(it, 0) + i, lda2 );
            }
        }
    }

    return PLASMA_SUCCESS;
}
/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaswp_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, PLASMA_Complex64_t *Aij,
                              int i1,  int i2, const int *ipiv, int inc, PLASMA_Complex64_t *fakepanel)
{
    DAG_CORE_LASWP;
    if (fakepanel == Aij) {
        QUARK_Insert_Task(
            quark, CORE_zlaswp_ontile_quark, task_flags,
            sizeof(PLASMA_desc),              &descA,     VALUE,
            sizeof(PLASMA_Complex64_t)*1,      Aij,           INOUT | LOCALITY,
            sizeof(int),                      &i1,        VALUE,
            sizeof(int),                      &i2,        VALUE,
            sizeof(int)*(i2-i1+1)*abs(inc),   ipiv,           INPUT,
            sizeof(int),                      &inc,       VALUE,
            sizeof(PLASMA_Complex64_t)*1,      fakepanel,     SCRATCH,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_zlaswp_ontile_quark, task_flags,
            sizeof(PLASMA_desc),              &descA,     VALUE,
            sizeof(PLASMA_Complex64_t)*1,      Aij,           INOUT | LOCALITY,
            sizeof(int),                      &i1,        VALUE,
            sizeof(int),                      &i2,        VALUE,
            sizeof(int)*(i2-i1+1)*abs(inc),   ipiv,           INPUT,
            sizeof(int),                      &inc,       VALUE,
            sizeof(PLASMA_Complex64_t)*1,      fakepanel,     INOUT,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswp_ontile_quark = PCORE_zlaswp_ontile_quark
#define CORE_zlaswp_ontile_quark PCORE_zlaswp_ontile_quark
#endif
void CORE_zlaswp_ontile_quark(Quark *quark)
{
    int i1, i2, inc;
    int *ipiv;
    PLASMA_Complex64_t *A, *fake;
    PLASMA_desc descA;

    quark_unpack_args_7(quark, descA, A, i1, i2, ipiv, inc, fake);
    CORE_zlaswp_ontile(descA, i1, i2, ipiv, inc);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaswp_ontile_f2(Quark *quark, Quark_Task_Flags *task_flags,
                                 PLASMA_desc descA, PLASMA_Complex64_t *Aij,
                                 int i1,  int i2, const int *ipiv, int inc,
                                 PLASMA_Complex64_t *fake1, int szefake1, int flag1,
                                 PLASMA_Complex64_t *fake2, int szefake2, int flag2)
{
    DAG_CORE_LASWP;
    QUARK_Insert_Task(
        quark, CORE_zlaswp_ontile_f2_quark, task_flags,
        sizeof(PLASMA_desc),                &descA, VALUE,
        sizeof(PLASMA_Complex64_t)*1,        Aij,       INOUT | LOCALITY,
        sizeof(int),                        &i1,    VALUE,
        sizeof(int),                        &i2,    VALUE,
        sizeof(int)*(i2-i1+1)*abs(inc),      ipiv,      INPUT,
        sizeof(int),                        &inc,   VALUE,
        sizeof(PLASMA_Complex64_t)*szefake1, fake1, flag1,
        sizeof(PLASMA_Complex64_t)*szefake2, fake2, flag2,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswp_ontile_f2_quark = PCORE_zlaswp_ontile_f2_quark
#define CORE_zlaswp_ontile_f2_quark PCORE_zlaswp_ontile_f2_quark
#endif
void CORE_zlaswp_ontile_f2_quark(Quark *quark)
{
    int i1, i2, inc;
    int *ipiv;
    PLASMA_Complex64_t *A;
    PLASMA_desc descA;
    void *fake1, *fake2;

    quark_unpack_args_8(quark, descA, A, i1, i2, ipiv, inc, fake1, fake2);
    CORE_zlaswp_ontile(descA, i1, i2, ipiv, inc);
}

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zswptr_ontile apply the zlaswp function on a matrix stored in
 *  tile layout, followed by a ztrsm on the first tile of the panel.
 *
 *******************************************************************************
 *
 *  @param[in,out] descA
 *          The descriptor of the matrix A to permute.
 *
 *  @param[in] i1
 *          The first element of IPIV for which a row interchange will
 *          be done.
 *
 *  @param[in] i2
 *          The last element of IPIV for which a row interchange will
 *          be done.
 *
 *  @param[in] ipiv
 *          The pivot indices; Only the element in position i1 to i2
 *          are accessed. The pivot are offset by A.i.
 *
 *  @param[in] inc
 *          The increment between successive values of IPIV.  If IPIV
 *          is negative, the pivots are applied in reverse order.
 *
 *  @param[in] Akk
 *          The triangular matrix Akk. The leading descA.nb-by-descA.nb
 *          lower triangular part of the array Akk contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. The
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in] ldak
 *          The leading dimension of the array Akk. ldak >= max(1,descA.nb).
 *
 *******************************************************************************
 */
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zswptr_ontile = PCORE_zswptr_ontile
#define CORE_zswptr_ontile PCORE_zswptr_ontile
#endif
int CORE_zswptr_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc,
                       const PLASMA_Complex64_t *Akk, int ldak)
{
    PLASMA_Complex64_t zone  = 1.0;
    int lda;
    int m = descA.mt == 1 ? descA.m : descA.mb;

    if ( descA.nt > 1 ) {
        coreblas_error(1, "Illegal value of descA.nt");
        return -1;
    }
    if ( i1 < 1 ) {
        coreblas_error(2, "Illegal value of i1");
        return -2;
    }
    if ( (i2 < i1) || (i2 > m) ) {
        coreblas_error(3, "Illegal value of i2");
        return -3;
    }

    CORE_zlaswp_ontile(descA, i1, i2, ipiv, inc);

    lda = BLKLDD(descA, 0);
    cblas_ztrsm( CblasColMajor, CblasLeft, CblasLower,
                 CblasNoTrans, CblasUnit,
                 m, descA.n, CBLAS_SADDR(zone),
                 Akk,     ldak,
                 A(0, 0), lda );

    return PLASMA_SUCCESS;
}
/***************************************************************************//**
 *
 **/
void QUARK_CORE_zswptr_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, PLASMA_Complex64_t *Aij,
                              int i1,  int i2, const int *ipiv, int inc,
                              const PLASMA_Complex64_t *Akk, int ldak)
{
    DAG_CORE_TRSM;
    QUARK_Insert_Task(
        quark, CORE_zswptr_ontile_quark, task_flags,
        sizeof(PLASMA_desc),              &descA, VALUE,
        sizeof(PLASMA_Complex64_t)*1,      Aij,       INOUT | LOCALITY,
        sizeof(int),                      &i1,    VALUE,
        sizeof(int),                      &i2,    VALUE,
        sizeof(int)*(i2-i1+1)*abs(inc),    ipiv,      INPUT,
        sizeof(int),                      &inc,   VALUE,
        sizeof(PLASMA_Complex64_t)*ldak,   Akk,       INPUT,
        sizeof(int),                      &ldak,  VALUE,
        0);
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zswptr_ontile_quark = PCORE_zswptr_ontile_quark
#define CORE_zswptr_ontile_quark PCORE_zswptr_ontile_quark
#endif
void CORE_zswptr_ontile_quark(Quark *quark)
{
    int i1, i2, inc, ldak;
    int *ipiv;
    PLASMA_Complex64_t *A, *Akk;
    PLASMA_desc descA;

    quark_unpack_args_8(quark, descA, A, i1, i2, ipiv, inc, Akk, ldak);
    CORE_zswptr_ontile(descA, i1, i2, ipiv, inc, Akk, ldak);
}

/***************************************************************************//**
 *
 * @ingroup CORE_PLASMA_Complex64_t
 *
 *  CORE_zlaswpc_ontile apply the zlaswp function on a matrix stored in
 *  tile layout
 *
 *******************************************************************************
 *
 *  @param[in,out] descA
 *          The descriptor of the matrix A to permute.
 *
 *  @param[in] i1
 *          The first element of IPIV for which a column interchange will
 *          be done.
 *
 *  @param[in] i2
 *          The last element of IPIV for which a column interchange will
 *          be done.
 *
 *  @param[in] ipiv
 *          The pivot indices; Only the element in position i1 to i2
 *          are accessed. The pivot are offset by A.i.
 *
 *  @param[in] inc
 *          The increment between successive values of IPIV.  If IPIV
 *          is negative, the pivots are applied in reverse order.
 *
 *******************************************************************************
 */
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswpc_ontile = PCORE_zlaswpc_ontile
#define CORE_zlaswpc_ontile PCORE_zlaswpc_ontile
#endif
int CORE_zlaswpc_ontile(PLASMA_desc descA, int i1, int i2, const int *ipiv, int inc)
{
    int i, j, ip, it;
    PLASMA_Complex64_t *A1;
    int lda;

    if ( descA.mt > 1 ) {
        coreblas_error(1, "Illegal value of descA.mt");
        return -1;
    }
    if ( i1 < 1 ) {
        coreblas_error(2, "Illegal value of i1");
        return -2;
    }
    if ( (i2 < i1) || (i2 > descA.n) ) {
        coreblas_error(3, "Illegal value of i2");
        return -3;
    }
    if ( ! ( (i2 - i1 - i1%descA.nb -1) < descA.nb ) ) {
        coreblas_error(2, "Illegal value of i1,i2. They have to be part of the same block.");
        return -3;
    }

    lda = BLKLDD(descA, 0);

    it = i1 / descA.nb;
    if (inc > 0) {
        A1 = A(0, it);

        for (j = i1-1; j < i2; ++j, ipiv+=inc) {
            ip = (*ipiv) - descA.j - 1;
            if ( ip != j )
            {
                it = ip / descA.nb;
                i  = ip % descA.nb;
                cblas_zswap(descA.m, A1       + j*lda, 1,
                                     A(0, it) + i*lda, 1 );
            }
        }
    }
    else
    {
        A1 = A(0, it);
        i1 -= 2;
        ipiv = &ipiv[(1-i2)*inc];
        for (j = i2-1; j > i1; --j, ipiv+=inc) {
            ip = (*ipiv) - descA.j - 1;
            if ( ip != j )
            {
                it = ip / descA.nb;
                i  = ip % descA.nb;
                cblas_zswap(descA.m, A1       + j*lda, 1,
                                     A(0, it) + i*lda, 1 );
            }
        }
    }

    return PLASMA_SUCCESS;
}
/***************************************************************************//**
 *
 **/
void QUARK_CORE_zlaswpc_ontile(Quark *quark, Quark_Task_Flags *task_flags,
                              PLASMA_desc descA, PLASMA_Complex64_t *Aij,
                              int i1,  int i2, const int *ipiv, int inc, PLASMA_Complex64_t *fakepanel)
{
    DAG_CORE_LASWP;
    if (fakepanel == Aij) {
        QUARK_Insert_Task(
            quark, CORE_zlaswpc_ontile_quark, task_flags,
            sizeof(PLASMA_desc),              &descA,     VALUE,
            sizeof(PLASMA_Complex64_t)*1,      Aij,           INOUT | LOCALITY,
            sizeof(int),                      &i1,        VALUE,
            sizeof(int),                      &i2,        VALUE,
            sizeof(int)*(i2-i1+1)*abs(inc),   ipiv,           INPUT,
            sizeof(int),                      &inc,       VALUE,
            sizeof(PLASMA_Complex64_t)*1,      fakepanel,     SCRATCH,
            0);
    } else {
        QUARK_Insert_Task(
            quark, CORE_zlaswpc_ontile_quark, task_flags,
            sizeof(PLASMA_desc),              &descA,     VALUE,
            sizeof(PLASMA_Complex64_t)*1,      Aij,           INOUT | LOCALITY,
            sizeof(int),                      &i1,        VALUE,
            sizeof(int),                      &i2,        VALUE,
            sizeof(int)*(i2-i1+1)*abs(inc),   ipiv,           INPUT,
            sizeof(int),                      &inc,       VALUE,
            sizeof(PLASMA_Complex64_t)*1,      fakepanel,     INOUT,
            0);
    }
}

/***************************************************************************//**
 *
 **/
#if defined(PLASMA_HAVE_WEAK)
#pragma weak CORE_zlaswpc_ontile_quark = PCORE_zlaswpc_ontile_quark
#define CORE_zlaswpc_ontile_quark PCORE_zlaswpc_ontile_quark
#endif
void CORE_zlaswpc_ontile_quark(Quark *quark)
{
    int i1, i2, inc;
    int *ipiv;
    PLASMA_Complex64_t *A, *fake;
    PLASMA_desc descA;

    quark_unpack_args_7(quark, descA, A, i1, i2, ipiv, inc, fake);
    CORE_zlaswpc_ontile(descA, i1, i2, ipiv, inc);
}

