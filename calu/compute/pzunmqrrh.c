/**
 *
 * @file pzunmqrrh.c
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
#define B(m,n)  BLKADDR(B, PLASMA_Complex64_t, (m), (n))
#define T(m,n)  BLKADDR(T, PLASMA_Complex64_t, (m), (n))
#define T2(m,n) BLKADDR(T, PLASMA_Complex64_t, (m), (n)+A.nt)
/***************************************************************************//**
 *  Parallel application of Q using tile V - QR factorization (reduction
 *  Householder) - dynamic scheduling
 **/
void plasma_pzunmqrrh_quark(PLASMA_enum side, PLASMA_enum trans,
        PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS,
        PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int K, M, RD, lastRD;
    int ldaM, ldam, ldan, ldaMRD;
    int ldbM, ldbm, ldbMRD;
    int tempMm, tempkn, tempnn, tempmm, tempMRDm, tempkmin;
    int ib;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    ib = PLASMA_IB;
    K = min(A.mt, A.nt);

    if (side == PlasmaLeft ) {
        if (trans == PlasmaConjTrans) {
            /*
             *  PlasmaLeft / PlasmaConjTrans
             */
            for (k = 0; k < K; k++) {
                tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                for (M = k; M < A.mt; M += BS) {
                    tempMm   = M == A.mt-1 ? A.m-M*A.mb : A.mb;
                    tempkmin = min(tempMm, tempkn);
                    ldaM = BLKLDD(A, M);
                    ldbM = BLKLDD(B, M);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_zunmqr(
                                plasma->quark, &task_flags,
                                side, trans,
                                tempMm, tempnn,
                                tempkmin, ib, T.nb,
                                A(M, k), ldaM,
                                T(M, k), T.mb,
                                B(M, n), ldbM);
                    }
                    for (m = M+1; m < min(M+BS, A.mt); m++) {
                        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                        ldbm = BLKLDD(B, m);
                        ldam = BLKLDD(A, m);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_ztsmqr(
                                    plasma->quark, &task_flags,
                                    side, trans,
                                    A.nb, tempnn, tempmm, tempnn,
                                    tempkn, ib, T.nb,
                                    B(M, n), ldbM,
                                    B(m, n), ldbm,
                                    A(m, k), ldam,
                                    T(m, k), T.mb);
                        }
                    }
                }
                for (RD = BS; RD < A.mt-k; RD *= 2) {
                    for (M = k; M+RD < A.mt; M += 2*RD) {
                        tempMRDm = M+RD == A.mt-1 ? A.m-(M+RD)*A.mb : A.mb;
                        ldbM   = BLKLDD(B, M   );
                        ldbMRD = BLKLDD(B, M+RD);
                        ldaMRD = BLKLDD(A, M+RD);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_zttmqr(
                                    plasma->quark, &task_flags,
                                    side, trans,
                                    A.nb, tempnn, tempMRDm, tempnn,
                                    tempkn, ib, T.nb,
                                    B (M,    n), ldbM,
                                    B (M+RD, n), ldbMRD,
                                    A (M+RD, k), ldaMRD,
                                    T2(M+RD, k), T.mb);
                        }
                    }
                }
            }
        } else {
            /*
             *  PlasmaLeft / PlasmaNoTrans
             */
            for (k = K-1; k >= 0; k--) {
                tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                lastRD = 0;
                for (RD = BS; RD < A.mt-k; RD *= 2)
                    lastRD = RD;
                for (RD = lastRD; RD >= BS; RD /= 2) {
                    for (M = k; M+RD < A.mt; M += 2*RD) {
                        tempMRDm = M+RD == A.mt-1 ? A.m-(M+RD)*A.mb : A.mb;
                        ldbM   = BLKLDD(B, M   );
                        ldbMRD = BLKLDD(B, M+RD);
                        ldaMRD = BLKLDD(A, M+RD);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_zttmqr(
                                    plasma->quark, &task_flags,
                                    side, trans,
                                    A.nb, tempnn, tempMRDm, tempnn,
                                    tempkn, ib, T.nb,
                                    B (M,    n), ldbM,
                                    B (M+RD, n), ldbMRD,
                                    A (M+RD, k), ldaMRD,
                                    T2(M+RD, k), T.mb);
                        }
                    }
                }
                for (M = k; M < A.mt; M += BS) {
                    tempMm   = M == A.mt-1 ? A.m-M*A.mb : A.mb;
                    tempkmin = min(tempMm, tempkn);
                    ldaM = BLKLDD(A, M);
                    ldbM = BLKLDD(B, M);
                    for (m = min(M+BS, A.mt)-1; m > M; m--) {
                        tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                        ldbm = BLKLDD(B, m);
                        ldam = BLKLDD(A, m);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_ztsmqr(
                                    plasma->quark, &task_flags,
                                    side, trans,
                                    A.nb, tempnn, tempmm, tempnn,
                                    tempkn, ib, T.nb,
                                    B(M, n), ldbM,
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
                                tempMm, tempnn,
                                tempkmin, ib, T.nb,
                                A(M, k), ldaM,
                                T(M, k), T.mb,
                                B(M, n), ldbM);
                    }
                }
            }
        }
    } else {
        if (trans == PlasmaConjTrans) {
            /*
             *  PlasmaRight / PlasmaConjTrans
             */
              for (k = K-1; k >= 0; k--) {
                  tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                  lastRD = 0;
                  for (RD = BS; RD < A.mt-k; RD *= 2)
                      lastRD = RD;
                  for (RD = lastRD; RD >= BS; RD /= 2) {
                      for (M = k; M+RD < A.mt; M += 2*RD) {
                          tempMRDm = M+RD == A.mt-1 ? A.m-(M+RD)*A.mb : A.mb;
                          ldaMRD = BLKLDD(A, M+RD);
                          for (m = 0; m < B.mt; m++) {
                              ldbm   = BLKLDD(B, m);
                              tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                              QUARK_CORE_zttmqr(
                                      plasma->quark, &task_flags,
                                      side, trans,
                                      tempmm, B.nb, tempmm, tempMRDm,
                                      tempkn, ib, T.nb,
                                      B (m, M), ldbm,
                                      B (m, M+RD), ldbm,
                                      A (M+RD, k), ldaMRD,
                                      T2(M+RD, k), T.mb);
                          }
                      }
                  }
                  for (M = k; M < A.mt; M += BS) {
                      tempMm   = M == A.mt-1 ? A.m-M*A.mb : A.mb;
                      tempkmin = min(tempMm, tempkn);
                      ldaM = BLKLDD(A, M);
                      ldbM = BLKLDD(B, M);
                      for (n = min(M+BS, A.mt)-1; n > M; n--) {
                          ldan = BLKLDD(A, n);
                          tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                          for (m = 0; m < B.mt; m++) {
                              ldbm = BLKLDD(B, m);
                              tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                              QUARK_CORE_ztsmqr(
                                      plasma->quark, &task_flags,
                                      side, trans,
                                      tempmm, tempMm, tempmm, tempnn,
                                      tempkn, ib, T.nb,
                                      B(m, M), ldbm,
                                      B(m, n), ldbm,
                                      A(n, k), ldan,
                                      T(n, k), T.mb);
                          }
                      }
                      for (m = 0; m < B.mt; m++) {
                          ldbm = BLKLDD(B, m);
                          tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                          QUARK_CORE_zunmqr(
                                  plasma->quark, &task_flags,
                                  side, trans,
                                  tempmm, tempMm,
                                  tempkmin, ib, T.nb,
                                  A(M, k), ldaM,
                                  T(M, k), T.mb,
                                  B(m, M), ldbm);
                      }
                  }
              }
        } else {
            /*
             *  PlasmaRight / PlasmaNoTrans
             */
            for (k = 0; k < K; k++) {
                tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                for (M = k; M < A.mt; M += BS) {
                    tempMm   = M == A.mt-1 ? A.m-M*A.mb : A.mb;
                    tempkmin = min(tempMm, tempkn);
                    ldaM = BLKLDD(A, M);
                    for (m = 0; m < B.mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        QUARK_CORE_zunmqr(
                                plasma->quark, &task_flags,
                                side, trans,
                                tempmm, tempMm,
                                tempkmin, ib, T.nb,
                                A(M, k), ldaM,
                                T(M, k), T.mb,
                                B(m, M), ldbm);
                    }
                    for (n = M+1; n < min(M+BS,  A.mt); n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        ldan = BLKLDD(A, n);
                        for (m = 0; m < B.mt; m++) {
                            tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                            ldbm = BLKLDD(B, m);
                            QUARK_CORE_ztsmqr(
                                    plasma->quark, &task_flags,
                                    side, trans,
                                    tempmm, tempMm, tempmm, tempnn,
                                    tempkn, ib, T.nb,
                                    B(m, M), ldbm,
                                    B(m, n), ldbm,
                                    A(n, k), ldan,
                                    T(n, k), T.mb);
                        }
                    }
                }
                for (RD = BS; RD < A.mt-k; RD *= 2) {
                    for (M = k; M+RD < A.mt; M += 2*RD) {
                        tempMRDm = M+RD == A.mt-1 ? A.m-(M+RD)*A.mb : A.mb;
                        ldaMRD = BLKLDD(A, M+RD);
                        for (m = 0; m < B.mt; m++) {
                            tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                            ldbm   = BLKLDD(B, m);
                            QUARK_CORE_zttmqr(
                                    plasma->quark, &task_flags,
                                    side, trans,
                                    tempmm, B.nb, tempmm, tempMRDm,
                                    tempkn, ib, T.nb,
                                    B (m, M   ), ldbm,
                                    B (m, M+RD), ldbm,
                                    A (M+RD, k), ldaMRD,
                                    T2(M+RD, k), T.mb);
                        }
                    }
                }
            }
        }
    }
}
