/**
 *
 * @file pzunmlq.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
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
 *  Parallel application of Q using tile V - LQ factorization - static scheduling
 **/
void plasma_pzunmlq(plasma_context_t *plasma)
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
    int ldak, ldbk, ldbm;
    int tempmm, tempnn, tempkm, tempkmin;
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
    ss_init(B.mt, B.nt, min(A.mt, A.nt));

    if (A.m > A.n) {
        minM  = A.n;
        minMT = A.nt;
    } else {
        minM  = A.m;
        minMT = A.mt;
    }

    k = minMT-1;
    n = PLASMA_RANK;
    while (n >= B.nt) {
        k--;
        n = n-B.nt;
    }
    m = B.mt-1;

    while (k >= 0 && n < B.nt) {
        next_n = n;
        next_m = m;
        next_k = k;

        next_m--;
        if (next_m == k-1) {
            next_n += PLASMA_SIZE;
            while (next_n >= B.nt && next_k >= 0) {
                next_k--;
                next_n = next_n-B.nt;
            }
            next_m = B.mt-1;
        }

        tempkmin = k == minMT-1 ? minM-k*A.nb : A.nb;
        tempkm   = k == B.mt-1 ? B.m-k*B.mb : B.mb;
        tempnn   = n == B.nt-1 ? B.n-n*B.nb : B.nb;
        tempmm   = m == B.mt-1 ? B.m-m*B.mb : B.mb;

        ldak = BLKLDD(A, k);
        ldbk = BLKLDD(B, k);
        ldbm = BLKLDD(B, m);

        if (m == k) {
            CORE_zunmlq(
                    side, trans,
                    tempkm, tempnn, tempkmin, ib,
                    A(k, k), ldak,
                    T(k, k), T.mb,
                    B(k, n), ldbk,
                    work, T.nb);
            ss_cond_set(k, n, k);
        }
        else {
            ss_cond_wait(m, n, k+1);
            CORE_ztsmlq(
                    side, trans,
                    A.mb, tempnn, tempmm, tempnn, tempkmin, ib,
                    B(k, n), ldbk,
                    B(m, n), ldbm,
                    A(k, m), ldak,
                    T(k, m), T.mb,
                    work, ib);
            ss_cond_set(m, n, k);
        }
        m = next_m;
        n = next_n;
        k = next_k;
    }
    plasma_private_free(plasma, work);
    ss_finalize();
}

/***************************************************************************//**
 *  Parallel application of Q using tile V - LQ factorization - dynamic scheduling
 **/
void plasma_pzunmlq_quark(PLASMA_enum side, PLASMA_enum trans,
        PLASMA_desc A, PLASMA_desc B, PLASMA_desc T,
        PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldak, ldbk, ldbm;
    int tempmm, tempnn, tempkn, tempkm, tempkmin;
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

    if (side == PlasmaLeft ) {
        if (trans == PlasmaNoTrans) {
            /*
             *  PlasmaLeft / PlasmaNoTrans
             */
            for (k = 0; k < minMT; k++) {
                tempkm   = k == B.mt -1 ? B.m -k*B.mb : B.mb;
                tempkmin = k == minMT-1 ? minM-k*A.nb : A.nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);
                for (n = 0; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    QUARK_CORE_zunmlq(
                            plasma->quark, &task_flags,
                            side, trans,
                            tempkm, tempnn, tempkmin, ib, T.nb,
                            A(k, k), ldak,
                            T(k, k), T.mb,
                            B(k, n), ldbk);
                }
                for (m = k+1; m < B.mt; m++) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldbm = BLKLDD(B, m);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_ztsmlq(
                                plasma->quark, &task_flags,
                                side, trans,
                                B.mb, tempnn, tempmm, tempnn, tempkmin, ib, T.nb,
                                B(k, n), ldbk,
                                B(m, n), ldbm,
                                A(k, m), ldak,
                                T(k, m), T.mb);
                    }
                }
            }
        }
        else {
            /*
             *  PlasmaLeft / PlasmaConjTrans
             */
            for (k = minMT-1; k >= 0; k--) {
                tempkm   = k == B.mt -1 ? B.m -k*B.mb : B.mb;
                tempkmin = k == minMT-1 ? minM-k*A.nb : A.nb;
                ldak = BLKLDD(A, k);
                ldbk = BLKLDD(B, k);
                for (m = B.mt-1; m > k; m--) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldbm = BLKLDD(B, m);
                    for (n = 0; n < B.nt; n++) {
                        tempnn   = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_ztsmlq(
                                plasma->quark, &task_flags,
                                side, trans,
                                B.mb, tempnn, tempmm, tempnn, tempkmin, ib, T.nb,
                                B(k, n), ldbk,
                                B(m, n), ldbm,
                                A(k, m), ldak,
                                T(k, m), T.mb);
                    }
                }
                for (n = 0; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    QUARK_CORE_zunmlq(
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
    else {
        if (trans == PlasmaNoTrans) {
            /*
             *  PlasmaRight / PlasmaNoTrans
             */
            for (k = minMT-1; k >= 0; k--) {
                tempkn   = k == B.nt -1 ? B.n -k*B.nb : B.nb;
                tempkmin = k == minMT-1 ? minM-k*A.nb : A.nb;
                ldak = BLKLDD(A, k);
                for (n = B.nt-1; n > k; n--) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        QUARK_CORE_ztsmlq(
                                plasma->quark, &task_flags,
                                side, trans,
                                tempmm, B.nb, tempmm, tempnn, tempkmin, ib, T.nb,
                                B(m, k), ldbm,
                                B(m, n), ldbm,
                                A(k, n), ldak,
                                T(k, n), T.mb);
                    }
                }
                for (m = 0; m < B.mt; m++) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldbm = BLKLDD(B, m);
                    QUARK_CORE_zunmlq(
                            plasma->quark, &task_flags,
                            side, trans,
                            tempmm, tempkn, tempkmin, ib, T.nb,
                            A(k, k), ldak,
                            T(k, k), T.mb,
                            B(m, k), ldbm);
                }
            }
        }
        else {
            /*
             *  PlasmaRight / PlasmaConjTrans
             */
            for (k = 0; k < minMT; k++) {
                tempkn   = k == B.nt -1 ? B.n -k*B.nb : B.nb;
                tempkmin = k == minMT-1 ? minM-k*A.mb : A.mb;
                ldak = BLKLDD(A, k);
                for (m = 0; m < B.mt; m++) {
                    tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                    ldbm = BLKLDD(B, m);
                    QUARK_CORE_zunmlq(
                            plasma->quark, &task_flags,
                            side, trans,
                            tempmm, tempkn, tempkmin, ib, T.nb,
                            A(k, k), ldak,
                            T(k, k), T.mb,
                            B(m, k), ldbm);
                }
                for (n = k+1; n < B.nt; n++) {
                    tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                    for (m = 0; m < B.mt; m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        QUARK_CORE_ztsmlq(
                                plasma->quark, &task_flags,
                                side, trans,
                                tempmm, B.nb, tempmm, tempnn, tempkmin, ib, T.nb,
                                B(m, k), ldbm,
                                B(m, n), ldbm,
                                A(k, n), ldak,
                                T(k, n), T.mb);
                    }
                }
            }
        }
    }
}
