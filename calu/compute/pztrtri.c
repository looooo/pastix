/**
 *
 * @file pztrtri.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel tile triangular matrix inverse - dynamic scheduling
 **/
void plasma_pztrtri_quark(PLASMA_enum uplo, PLASMA_enum diag, PLASMA_desc A,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n;
    int ldam, ldan;
    int tempkn, tempmm, tempnn;

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
        for (n = 0; n < A.nt; n++) {
            tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            ldan = BLKLDD(A, n);
            for (m = n+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                QUARK_CORE_ztrsm(
                    plasma->quark, &task_flags,
                    PlasmaRight, uplo, PlasmaNoTrans, diag,
                    tempmm, tempnn, A.mb,
                    mzone, A(n, n), ldan,
                           A(m, n), ldam);
            }
            for (m = n+1; m < A.mt; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                ldam = BLKLDD(A, m);
                for (k = 0; k < n; k++) {
                    tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                    QUARK_CORE_zgemm(
                        plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaNoTrans,
                        tempmm, tempkn, tempnn, A.mb,
                        zone, A(m, n), ldam,
                              A(n, k), ldan,
                        zone, A(m, k), ldam);
                }
            }
            for (m = 0; m < n; m++) {
                tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
                QUARK_CORE_ztrsm(
                    plasma->quark, &task_flags,
                    PlasmaLeft, uplo, PlasmaNoTrans, diag,
                    tempnn, tempmm, A.mb,
                    zone, A(n, n), ldan,
                          A(n, m), ldan);
            }
            QUARK_CORE_ztrtri(
                plasma->quark, &task_flags,
                uplo, diag,
                tempnn, A.mb,
                A(n, n), ldan,
                sequence, request, A.nb*n);
        }
    }
    /*
     *  PlasmaUpper
     */
    else {
        for (m = 0; m < A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            for (n = m+1; n < A.nt; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_ztrsm(
                    plasma->quark, &task_flags,
                    PlasmaLeft, uplo, PlasmaNoTrans, diag,
                    tempmm, tempnn, A.mb,
                    mzone, A(m, m), ldam,
                           A(m, n), ldam);
            }
            for (n = 0; n < m; n++) {
                tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                ldan = BLKLDD(A, n);
                for (k = m+1; k < A.nt; k++) {
                    tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
                    QUARK_CORE_zgemm(
                        plasma->quark, &task_flags,
                        PlasmaNoTrans, PlasmaNoTrans,
                        tempnn, tempkn, tempmm, A.mb,
                        zone, A(n, m), ldan,
                              A(m, k), ldam,
                        zone, A(n, k), ldan);
                }
                QUARK_CORE_ztrsm(
                    plasma->quark, &task_flags,
                    PlasmaRight, uplo, PlasmaNoTrans, diag,
                    tempnn, tempmm, A.mb,
                    zone, A(m, m), ldam,
                          A(n, m), ldan);
            }
            QUARK_CORE_ztrtri(
                plasma->quark, &task_flags,
                uplo, diag,
                tempmm, A.mb,
                A(m, m), ldam,
                sequence, request, A.mb*m);
        }
    }
}
