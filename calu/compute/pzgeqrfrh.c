/**
 *
 * @file pzgeqrfrh.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n)  BLKADDR(A, PLASMA_Complex64_t, (m), (n))
#define T(m,n)  BLKADDR(T, PLASMA_Complex64_t, (m), (n))
#define T2(m,n) BLKADDR(T, PLASMA_Complex64_t, (m), (n)+A.nt)
/***************************************************************************//**
 *  Parallel tile QR factorization (reduction Householder) - dynamic scheduling
 **/
void plasma_pzgeqrfrh_quark(PLASMA_desc A, PLASMA_desc T, int BS,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int M, RD;
    int ldaM, ldam, ldaMRD;
    int tempkmin, tempkn, tempMm, tempnn, tempmm, tempMRDm;
    int ib;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    ib = PLASMA_IB;

    for (k = 0; k < min(A.mt, A.nt); k++) {
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        for (M = k; M < A.mt; M += BS) {
            tempMm = M == A.mt-1 ? A.m-M*A.mb : A.mb;
            tempkmin = min(tempMm, tempkn);
            ldaM = BLKLDD(A, M);
            QUARK_CORE_zgeqrt(
                plasma->quark, &task_flags,
                tempMm, tempkn, ib, T.nb,
                A(M, k), ldaM,
                T(M, k), T.mb);
            for (n = k+1; n < A.nt; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_zunmqr(
                    plasma->quark, &task_flags,
                    PlasmaLeft, PlasmaConjTrans,
                    tempMm, tempnn, tempkmin, ib, T.nb,
                    A(M, k), ldaM,
                    T(M, k), T.mb,
                    A(M, n), ldaM);
            }
            for (m = M+1; m < min(M+BS, A.mt); m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                QUARK_CORE_ztsqrt(
                    plasma->quark, &task_flags,
                    tempmm, tempkn, ib, T.nb,
                    A(M, k), ldaM,
                    A(m, k), ldam,
                    T(m, k), T.mb);

                for (n = k+1; n < A.nt; n++) {
                    tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                    QUARK_CORE_ztsmqr(
                        plasma->quark, &task_flags,
                        PlasmaLeft, PlasmaConjTrans,
                        A.nb, tempnn, tempmm, tempnn, A.nb, ib, T.nb,
                        A(M, n), ldaM,
                        A(m, n), ldam,
                        A(m, k), ldam,
                        T(m, k), T.mb);
                }
            }
        }
        for (RD = BS; RD < A.mt-k; RD *= 2) {
            for (M = k; M+RD < A.mt; M += 2*RD) {
                tempMRDm = M+RD == A.mt-1 ? A.m-(M+RD)*A.mb : A.mb;
                ldaM   = BLKLDD(A, M   );
                ldaMRD = BLKLDD(A, M+RD);
                QUARK_CORE_zttqrt(
                    plasma->quark, &task_flags,
                    tempMRDm, tempkn, ib, T.nb,
                    A (M   , k), ldaM,
                    A (M+RD, k), ldaMRD,
                    T2(M+RD, k), T.mb);

                for (n = k+1; n < A.nt; n++) {
                    tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                    QUARK_CORE_zttmqr(
                        plasma->quark, &task_flags,
                        PlasmaLeft, PlasmaConjTrans,
                        A.nb, tempnn, tempMRDm, tempnn, A.nb, ib, T.nb,
                        A (M,    n), ldaM,
                        A (M+RD, n), ldaMRD,
                        A (M+RD, k), ldaMRD,
                        T2(M+RD, k), T.mb);
                }
            }
        }
    }
}
