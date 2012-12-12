/**
 *
 * @file pzunmqr.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define B(m,n) BLKADDR(B, PLASMA_Complex64_t, m, n)
#define T(m,n) BLKADDR(T, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel application of Q using tile V - QR factorization - static scheduling
 **/
void plasma_pzunmqr(plasma_context_t *plasma)
{
    PLASMA_enum side;
    PLASMA_enum trans;
    PLASMA_desc A;
    PLASMA_desc B;
    PLASMA_desc T;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int k, m, n;
    int next_k;
    int next_m;
    int next_n;
    int ldak, ldbk, ldam, ldbm;
    int tempkm, tempnn, tempkmin, tempmm;
    int minMT, minM;
    int ib = PLASMA_IB;
    PLASMA_Complex64_t *work;

    plasma_unpack_args_7(side, trans, A, B, T, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    if (side != PlasmaLeft) {
        plasma_request_fail(sequence, request, PLASMA_ERR_NOT_SUPPORTED);
        return;
    }
    if (trans != PlasmaConjTrans) {
        plasma_request_fail(sequence, request, PLASMA_ERR_NOT_SUPPORTED);
        return;
    }

    work = (PLASMA_Complex64_t*)plasma_private_alloc(plasma, ib*T.nb, T.dtyp);
    ss_init(B.mt, B.nt, -1);

   if (A.m > A.n) {
      minM  = A.n;
      minMT = A.nt;
    } else {
      minM  = A.m;
      minMT = A.mt;
    }

    k = 0;
    n = PLASMA_RANK;
    while (n >= B.nt) {
        k++;
        n = n-B.nt;
    }
    m = k;

    while (k < minMT && n < B.nt) {
        next_n = n;
        next_m = m;
        next_k = k;

        next_m++;
        if (next_m == A.mt) {
            next_n += PLASMA_SIZE;
            while (next_n >= B.nt && next_k < minMT) {
                next_k++;
                next_n = next_n-B.nt;
            }
            next_m = next_k;
        }

        tempkmin = k == minMT-1 ? minM-k*A.nb : A.nb;
        tempkm   = k == B.mt-1 ? B.m-k*B.mb : B.mb;
        tempnn   = n == B.nt-1 ? B.n-n*B.nb : B.nb;
        tempmm   = m == B.mt-1 ? B.m-m*B.mb : B.mb;

        ldak = BLKLDD(A, k);
        ldbk = BLKLDD(B, k);
        ldam = BLKLDD(A, m);
        ldbm = BLKLDD(B, m);

        if (m == k) {
            ss_cond_wait(k, n, k-1);
            CORE_zunmqr(
                side, trans,
                tempkm, tempnn, tempkmin, ib,
                A(k, k), ldak,
                T(k, k), T.mb,
                B(k, n), ldbk,
                work, T.nb);
            ss_cond_set(k, n, k);
        }
        else {
            ss_cond_wait(m, n, k-1);
            CORE_ztsmqr(
                side, trans,
                A.mb, tempnn, tempmm, tempnn, tempkmin, ib,
                B(k, n), ldbk,
                B(m, n), ldbm,
                A(m, k), ldam,
                T(m, k), T.mb,
                work, ib);
            ss_cond_set(m, n, k);
        }
        n = next_n;
        m = next_m;
        k = next_k;
    }
    plasma_private_free(plasma, work);
    ss_finalize();
}

/***************************************************************************//**
 *  Parallel application of Q using tile V - QR factorization - dynamic scheduling
 **/
void plasma_pzunmqr_quark(PLASMA_enum side, PLASMA_enum trans,
                          PLASMA_desc A, PLASMA_desc B, PLASMA_desc T,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldak, ldbk, ldam, ldan, ldbm;
    int tempkm, tempnn, tempkmin, tempmm, tempkn;
    int ib, minMT, minM;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    ib = PLASMA_IB;
    if (A.m > A.n) {
      minM  = A.n;
      minMT = A.nt;
    } else {
      minM  = A.m;
      minMT = A.mt;
    }

    /*
     *  PlasmaLeft / PlasmaConjTrans
     */
    if (side == PlasmaLeft ) {
        if (trans == PlasmaConjTrans) {
            for (k = 0; k < minMT; k++) {
                tempkm   = k == B.mt-1 ? B.m-k*B.mb : B.mb;
                tempkmin = k == minMT-1 ? minM-k*A.nb : A.nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);
                for (n = 0; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    QUARK_CORE_zunmqr(
                        plasma->quark, &task_flags,
                        side, trans,
                        tempkm, tempnn, tempkmin, ib, T.nb,
                        A(k, k), ldak,
                        T(k, k), T.mb,
                        B(k, n), ldbk);
                }
                for (m = k+1; m < B.mt; m++) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_ztsmqr(
                            plasma->quark, &task_flags,
                            side, trans,
                            B.mb, tempnn, tempmm, tempnn, tempkmin, ib, T.nb,
                            B(k, n), ldbk,
                            B(m, n), ldbm,
                            A(m, k), ldam,
                            T(m, k), T.mb);
                    }
                }
            }
        }
        /*
         *  PlasmaLeft / PlasmaNoTrans
         */
        else {
            for (k = minMT-1; k >= 0; k--) {
                tempkm = k == B.mt-1 ? B.m-k*B.mb : B.mb;
                tempkmin = k == minMT-1 ? minM-k*A.nb : A.nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);
                for (m = B.mt-1; m > k; m--) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldam = BLKLDD(A, m);
                    ldbm = BLKLDD(B, m);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_ztsmqr(
                            plasma->quark, &task_flags,
                            side, trans,
                            B.mb, tempnn, tempmm, tempnn, tempkmin, ib, T.nb,
                            B(k, n), ldbk,
                            B(m, n), ldbm,
                            A(m, k), ldam,
                            T(m, k), T.mb);
                    }
                }
                for (n = 0; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    QUARK_CORE_zunmqr(
                        plasma->quark, &task_flags,
                        side, trans,
                        tempkm, tempnn, tempkmin, ib, T.nb,
                        A(k, k), ldak,
                        T(k, k), T.mb,
                        B(k, n), ldbk);
                }
            }
        }
    }
    /*
     *  PlasmaRight / PlasmaConjTrans
     */
    else {
        if (trans == PlasmaConjTrans) {
            for (k = minMT-1; k >= 0; k--) {
                tempkn = k == B.nt-1 ? B.n-k*B.nb : B.nb;
                tempkmin = k == minMT-1 ? minM-k*A.nb : A.nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);
                for (n = B.nt-1; n > k; n--) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    ldan = BLKLDD(A, n);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        QUARK_CORE_ztsmqr(
                            plasma->quark, &task_flags,
                            side, trans,
                            tempmm, B.nb, tempmm, tempnn, tempkmin, ib, T.nb,
                            B(m, k), ldbm,
                            B(m, n), ldbm,
                            A(n, k), ldan,
                            T(n, k), T.mb);
                    }
                }
                for (m = 0; m < B.mt; m++) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldbm = BLKLDD(B, m);
                    QUARK_CORE_zunmqr(
                        plasma->quark, &task_flags,
                        side, trans,
                        tempmm, tempkn, tempkmin, ib, T.nb,
                        A(k, k), ldak,
                        T(k, k), T.mb,
                        B(m, k), ldbm);
                }
            }
        }
        /*
         *  PlasmaRight / PlasmaNoTrans
         */
        else {
            for (k = 0; k < minMT; k++) {
                tempkn   = k == B.nt-1 ? B.n-k*B.nb : B.nb;
                tempkmin = k == minMT-1 ? minM-k*A.nb : A.nb;
                ldak = BLKLDD(A, k); 
                for (m = 0; m < B.mt; m++) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldbm = BLKLDD(B, m);
                    QUARK_CORE_zunmqr(
                        plasma->quark, &task_flags,
                        side, trans,
                        tempmm, tempkn, tempkmin, ib, T.nb,
                        A(k, k), ldak,
                        T(k, k), T.mb,
                        B(m, k), ldbm);
                }
                for (n = k+1; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    ldan = BLKLDD(A, n);
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        QUARK_CORE_ztsmqr(
                            plasma->quark, &task_flags,
                            side, trans,
                            tempmm, B.nb, tempmm, tempnn, tempkmin, ib, T.nb,
                            B(m, k), ldbm,
                            B(m, n), ldbm,
                            A(n, k), ldan,
                            T(n, k), T.mb);
                    }
                }
            }
        }
    }
    
}
