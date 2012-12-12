/**
 *
 * @file pzgetrf_incpiv.c
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
#define L(m,n) BLKADDR(L, PLASMA_Complex64_t, m, n)
#define IPIV(m,n) &(IPIV[(int64_t)A.mb*((int64_t)(m)+(int64_t)A.mt*(int64_t)(n))])
/***************************************************************************//**
 *  Parallel tile LU factorization - static scheduling
 **/
void plasma_pzgetrf_incpiv(plasma_context_t *plasma)
{
    PLASMA_desc A;
    PLASMA_desc L;
    int *IPIV;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int k, m, n;
    int next_k;
    int next_m;
    int next_n;
    int ldak, ldam;
    int info;
    int tempkn, tempkm, tempmm, tempnn;
    int ib = PLASMA_IB;
    PLASMA_Complex64_t *work;

    plasma_unpack_args_5(A, L, IPIV, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;
    work = (PLASMA_Complex64_t*)plasma_private_alloc(plasma, ib*L.nb, L.dtyp);
    ss_init(A.mt, A.nt, -1);

    k = 0;
    n = PLASMA_RANK;
    while (n >= A.nt) {
        k++;
        n = n-A.nt+k;
    }
    m = k;

    while (k < min(A.mt, A.nt) && n < A.nt && !ss_aborted()) {
        next_n = n;
        next_m = m;
        next_k = k;

        next_m++;
        if (next_m == A.mt) {
            next_n += PLASMA_SIZE;
            while (next_n >= A.nt && next_k < min(A.mt, A.nt)) {
                next_k++;
                next_n = next_n-A.nt+next_k;
            }
            next_m = next_k;
        }

        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;

        ldak = BLKLDD(A, k);
        ldam = BLKLDD(A, m);

        if (n == k) {
            if (m == k) {
                ss_cond_wait(k, k, k-1);
                CORE_zgetrf_incpiv(
                    tempkm, tempkn, ib,
                    A(k, k), ldak,
                    IPIV(k, k), &info);
                if (info != 0 && m == A.mt-1) {
                    plasma_request_fail(sequence, request, info + A.nb*k);
                    ss_abort();
                }
                ss_cond_set(k, k, k);
            }
            else {
                ss_cond_wait(m, k, k-1);
                CORE_ztstrf(
                    tempmm, tempkn, ib, A.nb,
                    A(k, k), ldak,
                    A(m, k), ldam,
                    L(m, k), L.mb,
                    IPIV(m, k),
                    work, L.nb, &info);
                if (info != 0 && m == A.mt-1) {
                    plasma_request_fail(sequence, request, info + A.nb*k);
                    ss_abort();
                }
                ss_cond_set(m, k, k);
            }
        }
        else {
            if (m == k) {
                ss_cond_wait(k, k, k);
                ss_cond_wait(k, n, k-1);
                CORE_zgessm(
                    tempkm, tempnn, tempkm, ib,
                    IPIV(k, k),
                    A(k, k), ldak,
                    A(k, n), ldak);
            }
            else {
                ss_cond_wait(m, k, k);
                ss_cond_wait(m, n, k-1);
                CORE_zssssm(
                    A.nb, tempnn, tempmm, tempnn, A.nb, ib,
                    A(k, n), ldak,
                    A(m, n), ldam,
                    L(m, k), L.mb,
                    A(m, k), ldam,
                    IPIV(m, k));
                ss_cond_set(m, n, k);
            }
        }
        n = next_n;
        m = next_m;
        k = next_k;
    }
    plasma_private_free(plasma, work);
    ss_finalize();
}

/***************************************************************************//**
 *  Parallel tile LU factorization - dynamic scheduling
 **/
void plasma_pzgetrf_incpiv_quark(PLASMA_desc A, PLASMA_desc L, int *IPIV,
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
        QUARK_CORE_zgetrf_incpiv(
            plasma->quark, &task_flags,
            tempkm, tempkn, ib, L.nb,
            A(k, k), ldak, IPIV(k, k),
            sequence, request,
            k == A.mt-1, A.nb*k);

        for (n = k+1; n < A.nt; n++) {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            QUARK_CORE_zgessm(
                plasma->quark, &task_flags,
                tempkm, tempnn, tempkm, ib, L.nb,
                IPIV(k, k),
                A(k, k), ldak,
                A(k, n), ldak);
        }
        for (m = k+1; m < A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            QUARK_CORE_ztstrf(
                plasma->quark, &task_flags,
                tempmm, tempkn, ib, L.nb,
                A(k, k), ldak,
                A(m, k), ldam,
                L(m, k), L.mb,
                IPIV(m, k),
                sequence, request,
                m == A.mt-1, A.nb*k);

            for (n = k+1; n < A.nt; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_zssssm(
                    plasma->quark, &task_flags,
                    A.nb, tempnn, tempmm, tempnn, A.nb, ib, L.nb,
                    A(k, n), ldak,
                    A(m, n), ldam,
                    L(m, k), L.mb,
                    A(m, k), ldam,
                    IPIV(m, k));
            }
        }
    }
}
