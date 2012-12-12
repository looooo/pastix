/**
 *
 * @file pzgelqf.c
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

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define T(m,n) BLKADDR(T, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel tile LQ factorization - static scheduling
 **/
void plasma_pzgelqf(plasma_context_t *plasma)
{
    PLASMA_desc A;
    PLASMA_desc T;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int k, m, n;
    int next_k;
    int next_m;
    int next_n;
    int ldak, ldam;
    int tempkm, tempkn, tempmm, tempnn;
    int ib = PLASMA_IB;
    PLASMA_Complex64_t *work, *tau;

    plasma_unpack_args_4(A, T, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;
    work = (PLASMA_Complex64_t*)plasma_private_alloc(plasma, ib*T.nb, T.dtyp);
    tau  = (PLASMA_Complex64_t*)plasma_private_alloc(plasma, A.nb, A.dtyp);
    ss_init(A.mt, A.nt, -1);

    k = 0;
    m = PLASMA_RANK;
    while (m >= A.mt) {
        k++;
        m = m-A.mt+k;
    }
    n = k;

    while (k < min(A.mt, A.nt) && m < A.mt) {
        next_m = m;
        next_n = n;
        next_k = k;

        next_n++;
        if (next_n == A.nt) {
            next_m += PLASMA_SIZE;
            while (next_m >= A.mt && next_k < min(A.nt, A.mt)) {
                next_k++;
                next_m = next_m-A.mt+next_k;
            }
            next_n = next_k;
        }

        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;

        ldak = BLKLDD(A, k);
        ldam = BLKLDD(A, m);

        if (m == k) {
            if (n == k) {
                ss_cond_wait(k, k, k-1);
                CORE_zgelqt(
                    tempkm, tempkn, ib,
                    A(k, k), ldak,
                    T(k, k), T.mb,
                    tau, work);
                ss_cond_set(k, k, k);
            }
            else {
                ss_cond_wait(k, n, k-1);
                CORE_ztslqt(
                    tempkm, tempnn, ib,
                    A(k, k), ldak,
                    A(k, n), ldak,
                    T(k, n), T.mb,
                    tau, work);
                ss_cond_set(k, n, k);
            }
        }
        else {
            if (n == k) {
                ss_cond_wait(k, k, k);
                ss_cond_wait(m, k, k-1);
                CORE_zunmlq(
                    PlasmaRight, PlasmaConjTrans,
                    tempmm, tempkn, tempkn, ib,
                    A(k, k), ldak,
                    T(k, k), T.mb,
                    A(m, k), ldam,
                    work, T.nb);
            }
            else {
                ss_cond_wait(k, n, k);
                ss_cond_wait(m, n, k-1);
                CORE_ztsmlq(
                    PlasmaRight, PlasmaConjTrans,
                    tempmm, A.nb, tempmm, tempnn, A.nb, ib,
                    A(m, k), ldam,
                    A(m, n), ldam,
                    A(k, n), ldak,
                    T(k, n), T.mb,
                    work, T.nb);
                ss_cond_set(m, n, k);
            }
        }
        m = next_m;
        n = next_n;
        k = next_k;
    }
    plasma_private_free(plasma, work);
    plasma_private_free(plasma, tau);
    ss_finalize();
}

/***************************************************************************//**
 *  Parallel tile LQ factorization - dynamic scheduling
 **/
void plasma_pzgelqf_quark(PLASMA_desc A, PLASMA_desc T,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldak, ldam;
    int tempkm, tempkn, tempmm, tempnn;
    int ib;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    ib = PLASMA_IB;
    for (k = 0; k < min(A.mt, A.nt); k++) {
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        ldak = BLKLDD(A, k);
        QUARK_CORE_zgelqt(
            plasma->quark, &task_flags,
            tempkm, tempkn, ib, T.nb,
            A(k, k), ldak,
            T(k, k), T.mb);

        for (m = k+1; m < A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            QUARK_CORE_zunmlq(
                plasma->quark, &task_flags,
                PlasmaRight, PlasmaConjTrans,
                tempmm, tempkn, tempkn, ib, T.nb,
                A(k, k), ldak,
                T(k, k), T.mb,
                A(m, k), ldam);
        }
        for (n = k+1; n < A.nt; n++) {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            QUARK_CORE_ztslqt(
                plasma->quark, &task_flags,
                tempkm, tempnn, ib, T.nb,
                A(k, k), ldak,
                A(k, n), ldak,
                T(k, n), T.mb);

            for (m = k+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                QUARK_CORE_ztsmlq(
                    plasma->quark, &task_flags,
                    PlasmaRight, PlasmaConjTrans,
                    tempmm, A.nb, tempmm, tempnn, A.mb, ib, T.nb,
                    A(m, k), ldam,
                    A(m, n), ldam,
                    A(k, n), ldak,
                    T(k, n), T.mb);
            }
        }
    }
}
