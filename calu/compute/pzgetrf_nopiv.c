/**
 *
 * @file pzgetrf_nopiv.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.5
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
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/
void plasma_pzgetrf_nopiv_quark(PLASMA_desc A,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n, ib;
    int ldak, ldam;
    int tempkm, tempkn, tempmm, tempnn;

    PLASMA_Complex64_t zone  = (PLASMA_Complex64_t) 1.0;
    PLASMA_Complex64_t mzone = (PLASMA_Complex64_t)-1.0;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    ib = PLASMA_IB;

    /*
     *  PlasmaLower
     */
    for (k = 0; k < A.mt; k++) {
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        ldak = BLKLDD(A, k);
        QUARK_CORE_zgetrf_nopiv(
            plasma->quark, &task_flags,
            tempkm, tempkn, ib, A.mb,
            A(k, k), ldak);

        for (m = k+1; m < A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            QUARK_CORE_ztrsm(
                plasma->quark, &task_flags,
                PlasmaRight, PlasmaUpper, PlasmaNoTrans, PlasmaNonUnit,
                tempmm, tempkn, A.mb,
                zone, A(k, k), ldak,
                A(m, k), ldam);
        }
        for (n = k+1; n < A.nt; n++) {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            QUARK_CORE_ztrsm(
                plasma->quark, &task_flags,
                PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaUnit,
                tempkm, tempnn, A.mb,
                zone, A(k, k), ldak,
                A(k, n), ldak);

            for (m = k+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                QUARK_CORE_zgemm(
                    plasma->quark, &task_flags,
                    PlasmaNoTrans, PlasmaNoTrans,
                    tempmm, tempnn, A.mb, A.mb,
                    mzone, A(m, k), ldam,
                    A(k, n), ldak,
                    zone,  A(m, n), ldam);
            }
        }
    }
}
