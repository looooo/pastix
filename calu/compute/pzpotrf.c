/**
 *
 * @file pzpotrf.c
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
/***************************************************************************//**
 *  Parallel tile Cholesky factorization - static scheduling
 **/
void plasma_pzpotrf(plasma_context_t *plasma)
{
    PLASMA_enum uplo;
    PLASMA_desc A;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int k, m, n;
    int next_k;
    int next_m;
    int next_n;
    int ldak, ldam, ldan;
    int info;
    int tempkn, tempmn;

    PLASMA_Complex64_t zone  = (PLASMA_Complex64_t) 1.0;
    PLASMA_Complex64_t mzone = (PLASMA_Complex64_t)-1.0;

    plasma_unpack_args_4(uplo, A, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;
    ss_init(A.nt, A.nt, 0);

    k = 0;
    m = PLASMA_RANK;
    while (m >= A.nt) {
        k++;
        m = m-A.nt+k;
    }
    n = 0;

    while (k < A.nt && m < A.nt && !ss_aborted()) {
        next_n = n;
        next_m = m;
        next_k = k;

        next_n++;
        if (next_n > next_k) {
            next_m += PLASMA_SIZE;
            while (next_m >= A.nt && next_k < A.nt) {
                next_k++;
                next_m = next_m-A.nt+next_k;
            }
            next_n = 0;
        }

        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        tempmn = m == A.nt-1 ? A.n-m*A.nb : A.nb;

        ldak = BLKLDD(A, k);
        ldan = BLKLDD(A, n);
        ldam = BLKLDD(A, m);

        if (m == k) {
            if (n == k) {
                /*
                 *  PlasmaLower
                 */
                if (uplo == PlasmaLower) {
                    CORE_zpotrf(
                        PlasmaLower,
                        tempkn,
                        A(k, k), ldak,
                        &info);
                }
                /*
                 *  PlasmaUpper
                 */
                else {
                    CORE_zpotrf(
                        PlasmaUpper,
                        tempkn,
                        A(k, k), ldak,
                        &info);
                }
                if (info != 0) {
                    plasma_request_fail(sequence, request, info + A.nb*k);
                    ss_abort();
                }
                ss_cond_set(k, k, 1);
            }
            else {
                ss_cond_wait(k, n, 1);
                /*
                 *  PlasmaLower
                 */
                if (uplo == PlasmaLower) {
                    CORE_zherk(
                         PlasmaLower, PlasmaNoTrans,
                         tempkn, A.nb,
                         -1.0, A(k, n), ldak,
                          1.0, A(k, k), ldak);
                }
                /*
                 *  PlasmaUpper
                 */
                else {
                    CORE_zherk(
                         PlasmaUpper, PlasmaConjTrans,
                         tempkn, A.nb,
                         -1.0, A(n, k), ldan,
                          1.0, A(k, k), ldak);
                }
            }
        }
        else {
            if (n == k) {
                ss_cond_wait(k, k, 1);
                /*
                 *  PlasmaLower
                 */
                if (uplo == PlasmaLower) {
                    CORE_ztrsm(
                        PlasmaRight, PlasmaLower, PlasmaConjTrans, PlasmaNonUnit,
                        tempmn, A.nb,
                        zone, A(k, k), ldak,
                              A(m, k), ldam);
                }
                /*
                 *  PlasmaUpper
                 */
                else {
                    CORE_ztrsm(
                        PlasmaLeft, PlasmaUpper, PlasmaConjTrans, PlasmaNonUnit,
                        A.nb, tempmn,
                        zone, A(k, k), ldak,
                              A(k, m), ldak);
                }
                ss_cond_set(m, k, 1);
            }
            else {
                ss_cond_wait(k, n, 1);
                ss_cond_wait(m, n, 1);
                /*
                 *  PlasmaLower
                 */
                if (uplo == PlasmaLower) {
                    CORE_zgemm(
                        PlasmaNoTrans, PlasmaConjTrans,
                        tempmn, A.nb, A.nb,
                        mzone, A(m, n), ldam,
                               A(k, n), ldak,
                         zone, A(m, k), ldam);
                }
                /*
                 *  PlasmaUpper
                 */
                else {
                    CORE_zgemm(
                        PlasmaConjTrans, PlasmaNoTrans,
                        A.nb, tempmn, A.nb,
                        mzone, A(n, k), ldan,
                               A(n, m), ldan,
                         zone, A(k, m), ldak);
                }
            }
        }
        n = next_n;
        m = next_m;
        k = next_k;
    }
    ss_finalize();
}

/***************************************************************************//**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/
void plasma_pzpotrf_quark(PLASMA_enum uplo, PLASMA_desc A,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldak, ldam;
    int tempkm, tempmm;

    PLASMA_Complex64_t zone  = (PLASMA_Complex64_t) 1.0;
    PLASMA_Complex64_t mzone = (PLASMA_Complex64_t)-1.0;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    /*
     *  PlasmaLower
     */
    if (uplo == PlasmaLower) {
        for (k = 0; k < A.mt; k++) {
            tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
            ldak = BLKLDD(A, k);
            QUARK_CORE_zpotrf(
                plasma->quark, &task_flags,
                PlasmaLower, tempkm, A.mb,
                A(k, k), ldak,
                sequence, request, A.nb*k);

            for (m = k+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                QUARK_CORE_ztrsm(
                    plasma->quark, &task_flags,
                    PlasmaRight, PlasmaLower, PlasmaConjTrans, PlasmaNonUnit,
                    tempmm, A.mb, A.mb,
                    zone, A(k, k), ldak,
                          A(m, k), ldam);
            }
            for (m = k+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                QUARK_CORE_zherk(
                    plasma->quark, &task_flags,
                    PlasmaLower, PlasmaNoTrans,
                    tempmm, A.mb, A.mb,
                    -1.0, A(m, k), ldam,
                     1.0, A(m, m), ldam);

                for (n = k+1; n < m; n++) {
                    QUARK_CORE_zgemm(
                        plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaConjTrans,
                        tempmm, A.mb, A.mb, A.mb,
                        mzone, A(m, k), ldam,
                               A(n, k), A.mb,
                        zone,  A(m, n), ldam);
                }
            }
        }
    }
    /*
     *  PlasmaUpper
     */
    else {
        for (k = 0; k < A.nt; k++) {
            tempkm = k == A.nt-1 ? A.n-k*A.nb : A.nb;
            ldak = BLKLDD(A, k);
            QUARK_CORE_zpotrf(
                plasma->quark, &task_flags,
                PlasmaUpper,
                tempkm, A.mb,
                A(k, k), ldak,
                sequence, request, A.nb*k);

            for (m = k+1; m < A.nt; m++) {
                tempmm = m == A.nt-1 ? A.n-m*A.nb : A.nb;
                QUARK_CORE_ztrsm(
                    plasma->quark, &task_flags,
                    PlasmaLeft, PlasmaUpper, PlasmaConjTrans, PlasmaNonUnit,
                    A.nb, tempmm, A.mb,
                    zone, A(k, k), ldak,
                          A(k, m), ldak);
            }
            for (m = k+1; m < A.nt; m++) {
                tempmm = m == A.nt-1 ? A.n-m*A.nb : A.nb;
                ldam = BLKLDD(A, m);
                QUARK_CORE_zherk(
                    plasma->quark, &task_flags,
                    PlasmaUpper, PlasmaConjTrans,
                    tempmm, A.mb, A.mb,
                    -1.0, A(k, m), ldak,
                     1.0, A(m, m), ldam);

                for (n = k+1; n < m; n++) {
                    QUARK_CORE_zgemm(
                        plasma->quark, &task_flags,
                        PlasmaConjTrans, PlasmaNoTrans,
                        A.mb, tempmm, A.mb, A.mb,
                        mzone, A(k, n), ldak,
                               A(k, m), ldak,
                        zone,  A(n, m), A.mb);
                }
            }
        }
    }
}
