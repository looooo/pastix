/**
 *
 * @file pzunmlqrh.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
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
 *  Parallel application of Q using tile V - LQ factorization (reduction
 *  Householder) - dynamic scheduling
 **/
void plasma_pzunmlqrh_quark(PLASMA_enum side, PLASMA_enum trans,
        PLASMA_desc A, PLASMA_desc B, PLASMA_desc T, int BS,
        PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int K, N, RD, lastRD;
    int ldaN, ldak;
    int ldbN, ldbm, ldbNRD;
    int tempNn, tempkm, tempnn, tempmm, tempNRDn, tempkmin;
    int ib;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    ib = PLASMA_IB;
    K = min(A.mt, A.nt);

    if (side == PlasmaLeft ) {
        if (trans == PlasmaNoTrans) {
            /*
             *  PlasmaLeft / PlasmaNoTrans
             */
            for (k = 0; k < K; k++) {
                tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                ldak = BLKLDD(A, k);
                for (N = k; N < A.nt; N += BS) {
                    tempNn   = N == A.nt-1 ? A.n-N*A.nb : A.nb;
                    tempkmin = min(tempkm,tempNn);
                    ldaN = BLKLDD(A, N);
                    ldbN = BLKLDD(B, N);
                    for (n = 0; n < B.nt; n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        QUARK_CORE_zunmlq(
                                plasma->quark, &task_flags,
                                side, trans,
                                tempNn, tempnn,
                                tempkmin, ib, T.nb,
                                A(k, N), ldak,
                                T(k, N), T.mb,
                                B(N, n), ldbN);
                    }
                    for (m = N+1; m < min(N+BS, A.nt); m++) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_ztsmlq(
                                    plasma->quark, &task_flags,
                                    side, trans,
                                    B.nb, tempnn, tempmm, tempnn,
                                    tempkm, ib, T.nb,
                                    B(N, n), ldbN,
                                    B(m, n), ldbm,
                                    A(k, m), ldak,
                                    T(k, m), T.mb);
                        }
                    }
                }
                for (RD = BS; RD < A.nt-k; RD *= 2) {
                    for (N = k; N+RD < A.nt; N += 2*RD) {
                        tempNRDn = N+RD == A.nt-1 ? A.n-(N+RD)*A.nb : A.nb;
                        ldbN   = BLKLDD(B, N   );
                        ldbNRD = BLKLDD(B, N+RD);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_zttmlq(
                                    plasma->quark, &task_flags,
                                    side, trans,
                                    B.mb, tempnn, tempNRDn, tempnn,
                                    tempkm, ib, T.nb,
                                    B (N,    n), ldbN,
                                    B (N+RD, n), ldbNRD,
                                    A (k, N+RD), ldak,
                                    T2(k, N+RD), T.mb);
                        }
                    }
                }
            }
        } else {
            /*
             *  PlasmaLeft / PlasmaConjTrans
             */
            for (k = K-1; k >= 0; k--) {
                tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                ldak = BLKLDD(A, k);
                lastRD = 0;
                for (RD = BS; RD < A.nt-k; RD *= 2)
                    lastRD = RD;
                for (RD = lastRD; RD >= BS; RD /= 2) {
                    for (N = k; N+RD < A.nt; N += 2*RD) {
                        tempNRDn = N+RD == A.nt-1 ? A.n-(N+RD)*A.nb : A.nb;
                        ldbN   = BLKLDD(B, N   );
                        ldbNRD = BLKLDD(B, N+RD);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_zttmlq(
                                    plasma->quark, &task_flags,
                                    side, trans,
                                    B.nb, tempnn, tempNRDn, tempnn,
                                    tempkm, ib, T.nb,
                                    B (N,    n), ldbN,
                                    B (N+RD, n), ldbNRD,
                                    A (k, N+RD), ldak,
                                    T2(k, N+RD), T.mb);
                        }
                    }
                }
                for (N = k; N < A.nt; N += BS) {
                    tempNn   = N == A.nt-1 ? A.n-N*A.nb : A.nb;
                    tempkmin = min(tempkm,tempNn);
                    ldaN = BLKLDD(A, N);
                    ldbN = BLKLDD(B, N);
                    for (m = min(N+BS, A.nt)-1; m > N; m--) {
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        ldbm = BLKLDD(B, m);
                        for (n = 0; n < B.nt; n++) {
                            tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                            QUARK_CORE_ztsmlq(
                                    plasma->quark, &task_flags,
                                    side, trans,
                                    B.mb, tempnn, tempmm, tempnn,
                                    tempkm, ib, T.nb,
                                    B(N, n), ldbN,
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
                                tempNn, tempnn,
                                tempkmin, ib, T.nb,
                                A(k, N), ldak,
                                T(k, N), T.mb,
                                B(N, n), ldbN);
                    }
                }
            }

        }
    } else {
        if (trans == PlasmaNoTrans) {
            /*
             *  PlasmaRight / PlasmaNoTrans
             */
              for (k = K-1; k >= 0; k--) {
                  tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                  ldak = BLKLDD(A, k);
                  lastRD = 0;
                  for (RD = BS; RD < A.nt-k; RD *= 2)
                      lastRD = RD;
                  for (RD = lastRD; RD >= BS; RD /= 2) {
                      for (N = k; N+RD < A.nt; N += 2*RD) {
                          tempNRDn = N+RD == A.nt-1 ? A.n-(N+RD)*A.nb : A.nb;
                          for (m = 0; m < B.mt; m++) {
                              ldbm   = BLKLDD(B, m);
                              tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                              QUARK_CORE_zttmlq(
                                      plasma->quark, &task_flags,
                                      side, trans,
                                      tempmm, B.nb, tempmm, tempNRDn,
                                      tempkm, ib, T.nb,
                                      B (m, N   ), ldbm,
                                      B (m, N+RD), ldbm,
                                      A (k, N+RD), ldak,
                                      T2(k, N+RD), T.mb);
                          }
                      }
                  }
                  for (N = k; N < A.nt; N += BS) {
                      tempNn   = N == A.nt-1 ? A.n-N*A.nb : A.nb;
                      tempkmin = min(tempkm,tempNn);
                      for (n = min(N+BS, A.nt)-1; n > N; n--) {
                          tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                          for (m = 0; m < B.mt; m++) {
                              tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                              ldbm = BLKLDD(B, m);
                              QUARK_CORE_ztsmlq(
                                      plasma->quark, &task_flags,
                                      side, trans,
                                      tempmm, B.nb, tempmm, tempnn,
                                      tempkm, ib, T.nb,
                                      B(m, N), ldbm,
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
                                  tempmm, tempNn,
                                  tempkmin, ib, T.nb,
                                  A(k, N), ldak,
                                  T(k, N), T.mb,
                                  B(m, N), ldbm);
                      }
                  }
              }
        } else {
            /*
             *  PlasmaRight / PlasmaConjTrans
             */
            for (k = 0; k < K; k++) {
                tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
                ldak = BLKLDD(A, k);
                for (N = k; N < A.nt; N += BS) {
                    tempNn = N == A.nt-1 ? A.n-N*A.nb : A.nb;
                    tempkmin = min(tempkm,tempNn);
                    ldaN = BLKLDD(A, N);
                    for (m = 0; m < B.mt; m++) {
                        ldbm = BLKLDD(B, m);
                        tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                        QUARK_CORE_zunmlq(
                                plasma->quark, &task_flags,
                                side, trans,
                                tempmm, tempNn,
                                tempkmin, ib, T.nb,
                                A(k, N), ldaN,
                                T(k, N), T.mb,
                                B(m, N), ldbm);
                    }
                    for (n = N+1; n < min(N+BS, A.nt); n++) {
                        tempnn = n == B.nt-1 ? B.n-n*B.nb : B.nb;
                        for (m = 0; m < B.mt; m++) {
                            tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                            ldbm = BLKLDD(B, m);
                            QUARK_CORE_ztsmlq(
                                    plasma->quark, &task_flags,
                                    side, trans,
                                    tempmm, tempNn, tempmm, tempnn,
                                    tempkm, ib, T.nb,
                                    B(m, N), ldbm,
                                    B(m, n), ldbm,
                                    A(k, n), ldak,
                                    T(k, n), T.mb);
                        }
                    }
                }
                for (RD = BS; RD < A.nt-k; RD *= 2) {
                    for (N = k; N+RD < A.nt; N += 2*RD) {
                        tempNRDn = N+RD == A.nt-1 ? A.n-(N+RD)*A.nb : A.nb;
                        for (m = 0; m < B.mt; m++) {
                            tempmm = m == B.mt-1 ? B.m-m*B.mb : B.mb;
                            ldbm   = BLKLDD(B, m);
                            QUARK_CORE_zttmlq(
                                    plasma->quark, &task_flags,
                                    side, trans,
                                    tempmm, B.nb, tempmm, tempNRDn,
                                    tempkm, ib, T.nb,
                                    B (m, N   ), ldbm,
                                    B (m, N+RD), ldbm,
                                    A (k, N+RD), ldak,
                                    T2(k, N+RD), T.mb);
                        }
                    }
                }
            }
        }
    }
}
