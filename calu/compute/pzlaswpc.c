/**
 *
 * @file pzlaswpc.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define B(m, n) BLKADDR(B, PLASMA_Complex64_t, m, n)
#define IPIV(k) &(IPIV[(int64_t)B.mb*(int64_t)(k)])

/***************************************************************************//**
 *  Parallel tile column interchanges - dynamic scheduling
 **/
void plasma_pzlaswpc_quark(PLASMA_desc B, int *IPIV, int inc,
                           PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int m, n;
    int tempj, tempn, tempmm, tempnn;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    if ( inc > 0 )
    {
        for (n = 0; n < B.nt; n++) {
            tempj = n * B.nb;
            tempn = B.n - tempj;
            tempnn = n == B.nt-1 ? tempn : B.nb;

            for (m = 0; m < B.mt; m++) {
                tempmm = m == B.mt-1 ? B.m - m * B.mb : B.mb;

                QUARK_CORE_zlaswpc_ontile(
                    plasma->quark, &task_flags,
                    plasma_desc_submatrix(B, m*B.mb, tempj, tempmm, tempn),
                    B(m, n), 1, tempnn, IPIV(n), inc, B(m, B.nt-1) );
            }
        }
    }
    else
    {
        for (n = B.nt-1; n > -1; n--) {
            tempj = n * B.nb;
            tempn = B.n - tempj;
            tempnn = n == B.nt-1 ? tempn : B.nb;

            for (m = 0; m < B.mt; m++) {
                tempmm = m == B.mt-1 ? B.m - m * B.mb : B.mb;

                QUARK_CORE_zlaswpc_ontile(
                    plasma->quark, &task_flags,
                    plasma_desc_submatrix(B, m*B.mb, tempj, tempmm, tempn),
                    B(m, n), 1, tempnn, IPIV(n), inc, B(m, 0) );
            }
        }
    }
}
