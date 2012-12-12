/**
 *
 * @file pzungqr.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Hatem Ltaief
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define Q(m,n) BLKADDR(Q, PLASMA_Complex64_t, m, n)
#define T(m,n) BLKADDR(T, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel construction of Q using tile V (application to identity) - dynamic scheduling
 **/
void plasma_pzungqr_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldak, ldqk, ldam, ldqm;
    int tempmm, tempnn, tempkmin, tempkm;
    int tempAkm, tempAkn;
    int ib;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    ib = PLASMA_IB;
    for (k = min(A.mt, A.nt)-1; k >= 0; k--) {
        tempAkm  = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempAkn  = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        tempkmin = min( tempAkn, tempAkm );
        tempkm   = k == Q.mt-1 ? Q.m-k*Q.mb : Q.mb;
        ldak = BLKLDD(A, k);
        ldqk = BLKLDD(Q, k);
        for (m = Q.mt - 1; m > k; m--) {
            tempmm = m == Q.mt-1 ? Q.m-m*Q.mb : Q.mb;
            ldam = BLKLDD(A, m);
            ldqm = BLKLDD(Q, m);
            for (n = 0; n < Q.nt; n++) {
                tempnn = n == Q.nt-1 ? Q.n-n*Q.nb : Q.nb;
                QUARK_CORE_ztsmqr(
                    plasma->quark, &task_flags,
                    PlasmaLeft, PlasmaNoTrans,
                    Q.mb, tempnn, tempmm, tempnn, tempAkn, ib, T.nb,
                    Q(k, n), ldqk,
                    Q(m, n), ldqm,
                    A(m, k), ldam,
                    T(m, k), T.mb);
            }
        }
        for (n = 0; n < Q.nt; n++) {
            tempnn = n == Q.nt-1 ? Q.n-n*Q.nb : Q.nb;
            QUARK_CORE_zunmqr(
                plasma->quark, &task_flags,
                PlasmaLeft, PlasmaNoTrans,
                tempkm, tempnn, tempkmin, ib, T.nb,
                A(k, k), ldak,
                T(k, k), T.mb,
                Q(k, n), ldqk);
        }
    }
}
