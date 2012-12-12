/**
 *
 * @file pztrsmpl.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n)    BLKADDR(A, PLASMA_Complex64_t, m, n)
#define B(m,n)    BLKADDR(B, PLASMA_Complex64_t, m, n)
#define L(m,n)    BLKADDR(L, PLASMA_Complex64_t, m, n)
#define IPIV(m,n) &(IPIV[(int64_t)A.nb*((int64_t)(m)+(int64_t)A.mt*(int64_t)(n))])
/***************************************************************************//**
 *  Parallel forward substitution for tile LU - static scheduling
 **/
void plasma_pztrsmpl(plasma_context_t *plasma)
{
    PLASMA_desc A;
    PLASMA_desc B;
    PLASMA_desc L;
    int *IPIV;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int k, m, n;
    int next_k;
    int next_m;
    int next_n;
    int ldak, ldbk, ldam, ldbm;
    int tempkm, tempnn, tempkmin, tempmm, tempkn;
    int ib;

    plasma_unpack_args_6(A, B, L, IPIV, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;
    ss_init(B.mt, B.nt, -1);

    ib = PLASMA_IB;
    k = 0;
    n = PLASMA_RANK;
    while (n >= B.nt) {
        k++;
        n = n-B.nt;
    }
    m = k;

    while (k < min(A.mt, A.nt) && n < B.nt) {
        next_n = n;
        next_m = m;
        next_k = k;

        next_m++;
        if (next_m == A.mt) {
            next_n += PLASMA_SIZE;
            while (next_n >= B.nt && next_k < min(A.mt, A.nt)) {
                next_k++;
                next_n = next_n-B.nt;
            }
            next_m = next_k;
        }

        tempkm   = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn   = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        tempkmin = k == min(A.mt, A.nt)-1 ? min(A.m, A.n)-k*A.mb : A.mb;
        tempnn   = n == B.nt-1 ? B.n-n*B.nb : B.nb;
        tempmm   = m == A.mt-1 ? A.m-m*A.mb : A.mb;

        ldak = BLKLDD(A, k);
        ldbk = BLKLDD(B, k);
        ldam = BLKLDD(A, m);
        ldbm = BLKLDD(B, m);

        if (m == k) {
            ss_cond_wait(k, n, k-1);
            CORE_zgessm(
                tempkm, tempnn, tempkmin, ib,
                IPIV(k, k),
                A(k, k), ldak,
                B(k, n), ldbk);
            ss_cond_set(k, n, k);
        }
        else {
            ss_cond_wait(m, n, k-1);
            CORE_zssssm(
                A.nb, tempnn, tempmm, tempnn, tempkn, ib,
                B(k, n), ldbk,
                B(m, n), ldbm,
                L(m, k), L.mb,
                A(m, k), ldam,
                IPIV(m, k));
            ss_cond_set(m, n, k);
        }
        n = next_n;
        m = next_m;
        k = next_k;
    }
    ss_finalize();
}

/***************************************************************************//**
 *  Parallel forward substitution for tile LU - dynamic scheduling
 **/
void plasma_pztrsmpl_quark(PLASMA_desc A, PLASMA_desc B, PLASMA_desc L, int *IPIV,
                           PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldak, ldam, ldbk, ldbm;
    int tempkm, tempnn, tempkmin, tempmm, tempkn;
    int ib;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    ib = PLASMA_IB;
    for (k = 0; k < min(A.mt, A.nt); k++) {
        tempkm   = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn   = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        tempkmin = k == min(A.mt, A.nt)-1 ? min(A.m, A.n)-k*A.mb : A.mb;
        ldak = BLKLDD(A, k);
        ldbk = BLKLDD(B, k);
        for (n = 0; n < B.nt; n++) {
            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
            QUARK_CORE_zgessm(
                plasma->quark, &task_flags,
                tempkm, tempnn, tempkmin, ib, L.nb,
                IPIV(k, k),
                A(k, k), ldak,
                B(k, n), ldbk);
        }
        for (m = k+1; m < A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            ldbm = BLKLDD(B, m);
            for (n = 0; n < B.nt; n++) {
                tempnn  = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                QUARK_CORE_zssssm(
                    plasma->quark, &task_flags,
                    A.nb, tempnn, tempmm, tempnn, tempkn, ib, L.nb,
                    B(k, n), ldbk,
                    B(m, n), ldbm,
                    L(m, k), L.mb,
                    A(m, k), ldam,
                    IPIV(m, k));
            }
        }
    }
}
