/**
 *
 * @file pzgeadd.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define B(m,n) BLKADDR(B, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *
 **/
void plasma_pzgeadd(plasma_context_t *plasma)
{
    PLASMA_Complex64_t alpha;
    PLASMA_desc A;
    PLASMA_desc B;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int X, Y;
    int m, n;
    int next_m;
    int next_n;
    int ldam, ldbm;

    plasma_unpack_args_5(alpha, A, B, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    while (m >= A.mt && n < A.nt) {
        n++;
        m = m-A.mt;
    }

    while (n < A.nt) {
        next_m = m;
        next_n = n;

        next_m += PLASMA_SIZE;
        while (next_m >= A.mt && next_n < A.nt) {
            next_n++;
            next_m = next_m-A.mt;
        }

        X = m == A.mt-1 ? A.m-A.mb*m : A.nb;
        Y = n == A.nt-1 ? A.n-A.nb*n : A.nb;
        ldam = BLKLDD(A, m);
        ldbm = BLKLDD(B, m);
        CORE_zgeadd(X, Y, alpha, A(m, n), ldam, B(m, n), ldbm);

        m = next_m;
        n = next_n;
    }
}

/***************************************************************************//**
 *
 **/
void plasma_pzgeadd_quark(PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc B,
                         PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int X, Y;
    int m, n;
    int ldam, ldbm;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for (m = 0; m < A.mt; m++) {
        X = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        ldam = BLKLDD(A, m);
        ldbm = BLKLDD(B, m);

        for (n = 0; n < A.nt; n++) {
            Y = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            QUARK_CORE_zgeadd(
                plasma->quark, &task_flags,
                X, Y, A.mb,
                alpha, A(m, n), ldam,
                       B(m, n), ldbm);
        }
    }
}
