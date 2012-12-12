/**
 *
 * @file pzungqrrh.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Dulceneia Becker
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n)  BLKADDR(A, PLASMA_Complex64_t, (m), (n))
#define Q(m,n)  BLKADDR(Q, PLASMA_Complex64_t, (m), (n))
#define T(m,n)  BLKADDR(T, PLASMA_Complex64_t, (m), (n))
#define T2(m,n) BLKADDR(T, PLASMA_Complex64_t, (m), (n)+(A.nt))
/***************************************************************************//**
 *  Parallel construction of Q using tile V (application to identity;
 * reduction Householder) - dynamic scheduling
 **/
void plasma_pzungqrrh_quark(PLASMA_desc A, PLASMA_desc Q,
                            PLASMA_desc T, int BS,
                            PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int K, M, RD, lastRD;
    int ldaM, ldam, ldaMRD;
    int ldbM, ldbm, ldbMRD;
    int tempkn, tempMm, tempnn, tempmm, tempMRDm, tempkmin;
    int ib;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    ib = PLASMA_IB;
    K = min(A.mt, A.nt);
    for (k = K-1; k >= 0; k--) {
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        lastRD = 0;
        for (RD = BS; RD < A.mt-k; RD *= 2)
            lastRD = RD;
        for (RD = lastRD; RD >= BS; RD /= 2) {
            for (M = k; M+RD < A.mt; M += 2*RD) {
                tempMRDm = M+RD == A.mt-1 ? A.m-(M+RD)*A.mb : A.mb;
                ldbM   = BLKLDD(Q, M   );
                ldbMRD = BLKLDD(Q, M+RD);
                ldaMRD = BLKLDD(A, M+RD);
                for (n = 0; n < Q.nt; n++) {
                    tempnn = n == Q.nt-1 ? Q.n-n*Q.nb : Q.nb;
                    QUARK_CORE_zttmqr(
                        plasma->quark, &task_flags,
                        PlasmaLeft, PlasmaNoTrans,
                        A.nb, tempnn, tempMRDm, tempnn,
                        tempkn, ib, T.nb,
                        Q (M,    n), ldbM,
                        Q (M+RD, n), ldbMRD,
                        A (M+RD, k), ldaMRD,
                        T2(M+RD, k), T.mb); 
                }
            }
        }
        for (M = k; M < A.mt; M += BS) {
            tempMm   = M == A.mt-1 ? A.m-M*A.mb : A.mb;
            tempkmin = min(tempMm, tempkn);
            ldaM = BLKLDD(A, M);
            ldbM = BLKLDD(Q, M);
            for (m = min(M+BS, A.mt)-1; m > M; m--) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldbm = BLKLDD(Q, m);
                ldam = BLKLDD(A, m);

                for (n = 0; n < Q.nt; n++) {
                    tempnn = n == Q.nt-1 ? Q.n-n*Q.nb : Q.nb;
                    QUARK_CORE_ztsmqr(
                        plasma->quark, &task_flags,
                        PlasmaLeft, PlasmaNoTrans,
                        A.nb, tempnn, tempmm, tempnn,
                        tempkn, ib, T.nb,
                        Q(M, n), ldbM,
                        Q(m, n), ldbm,
                        A(m, k), ldam,
                        T(m, k), T.mb);
                }
            }
            for (n = 0; n < Q.nt; n++) {
                tempnn = n == Q.nt-1 ? Q.n-n*Q.nb : Q.nb;
                QUARK_CORE_zunmqr(
                    plasma->quark, &task_flags,
                    PlasmaLeft, PlasmaNoTrans,
                    tempMm, tempnn,
                    tempkmin, ib, T.nb,
                    A(M, k), ldaM,
                    T(M, k), T.mb,
                    Q(M, n), ldbM);
            }
        }
    }
}
