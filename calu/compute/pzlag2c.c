/**
 *
 * @file pzlag2c.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions mixed zc -> ds
 *
 **/
#include "common.h"

#define A(m,n)  BLKADDR(A,  PLASMA_Complex64_t, m, n)
#define B(m,n)  BLKADDR(B,  PLASMA_Complex64_t, m, n)
#define SA(m,n) BLKADDR(SA, PLASMA_Complex32_t, m, n)
#define SB(m,n) BLKADDR(SB, PLASMA_Complex32_t, m, n)
/***************************************************************************//**
 *
 **/
void plasma_pzlag2c(plasma_context_t *plasma)
{
    PLASMA_desc A;
    PLASMA_desc SB;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int X, Y;
    int m, n;
    int next_m;
    int next_n;
    int ldam, ldbm;
    int info = PLASMA_SUCCESS;

    plasma_unpack_args_4(A, SB, sequence, request);
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
        ldbm = BLKLDD(SB, m);
        CORE_zlag2c(X, Y, A(m, n), ldam, SB(m, n), ldbm, &info);

        if (info != 0)
            plasma_request_fail(sequence, request, info);

        m = next_m;
        n = next_n;
    }
}

/***************************************************************************//**
 *
 **/
void plasma_pzlag2c_quark(PLASMA_desc A, PLASMA_desc SB,
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

    for(m = 0; m < A.mt; m++) {
        X = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        ldam = BLKLDD(A, m);
        ldbm = BLKLDD(SB, m);
        for(n = 0; n < A.nt; n++) {
            Y = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            QUARK_CORE_zlag2c(
                plasma->quark, &task_flags,
                X, Y, A.mb,
                A(m, n), ldam,
                SB(m, n), ldbm,
                sequence, request);
        }
    }
}

/***************************************************************************//**
 *
 **/
void plasma_pclag2z(plasma_context_t *plasma)
{
    PLASMA_desc SA;
    PLASMA_desc B;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int X, Y;
    int m, n;
    int ldam, ldbm;
    int next_m;
    int next_n;

    plasma_unpack_args_4(SA, B, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    while (m >= SA.mt && n < SA.nt) {
        n++;
        m = m-SA.mt;
    }

    while (n < SA.nt) {
        next_m = m;
        next_n = n;

        next_m += PLASMA_SIZE;
        while (next_m >= SA.mt && next_n < SA.nt) {
            next_n++;
            next_m = next_m-SA.mt;
        }

        X = m == SA.mt-1 ? SA.m-SA.mb*m : SA.nb;
        Y = n == SA.nt-1 ? SA.n-SA.nb*n : SA.nb;
        ldam = BLKLDD(SA, m);
        ldbm = BLKLDD(B, m);
        CORE_clag2z(X, Y, SA(m, n), ldam, B(m, n), ldbm);

        m = next_m;
        n = next_n;
    }
}

/***************************************************************************//**
 *
 **/
void plasma_pclag2z_quark(PLASMA_desc SA, PLASMA_desc B,
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

    for(m = 0; m < SA.mt; m++) {
        X = m == SA.mt-1 ? SA.m-m*SA.mb : SA.mb;
        ldam = BLKLDD(SA, m);
        ldbm = BLKLDD(B, m);
        for(n = 0; n < SA.nt; n++) {
            Y = n == SA.nt-1 ? SA.n-n*SA.nb : SA.nb;
            QUARK_CORE_clag2z(
                plasma->quark, &task_flags,
                X, Y, SA.mb,
                SA(m, n), ldam,
                B(m, n), ldbm);
        }
    }
}
