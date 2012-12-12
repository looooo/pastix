/**
 *
 * @file pzunglq.c
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
void plasma_pzunglq_quark(PLASMA_desc A, PLASMA_desc Q, PLASMA_desc T,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldak, ldqm;
    int tempnn, tempmm, tempkmin, tempkn;
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
        tempkn   = k == Q.nt-1 ? Q.n-k*Q.nb : Q.nb;
        ldak = BLKLDD(A, k);
        for (n = Q.nt-1; n > k; n--) {
            tempnn = n == Q.nt-1 ? Q.n-n*Q.nb : Q.nb;
            for (m = 0; m < Q.mt; m++) {
                tempmm = m == Q.mt-1 ? Q.m-m*Q.mb : Q.mb;
                ldqm = BLKLDD(Q, m);
                QUARK_CORE_ztsmlq(
                    plasma->quark, &task_flags,
                    PlasmaRight, PlasmaNoTrans,
                    tempmm, Q.nb, tempmm, tempnn, tempAkm, ib, T.nb,
                    Q(m, k), ldqm,
                    Q(m, n), ldqm,
                    A(k, n), ldak,
                    T(k, n), T.mb);
            }
        }
        for (m = 0; m < Q.mt; m++) {
            tempmm = m == Q.mt-1 ? Q.m-m*Q.mb : Q.mb;
            ldqm = BLKLDD(Q, m);
            QUARK_CORE_zunmlq(
                plasma->quark, &task_flags,
                PlasmaRight, PlasmaNoTrans,
                tempmm, tempkn, tempkmin, ib, T.nb,
                A(k, k), ldak,
                T(k, k), T.mb,
                Q(m, k), ldqm);
        }
    }
}
