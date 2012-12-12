/**
 *
 * @file pzlaset2.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel initializztion a 2-D array A to 
 *  ALPHA on the offdiagonals.
 **/
void plasma_pzlaset2_quark(PLASMA_enum uplo, PLASMA_Complex64_t alpha, 
                           PLASMA_desc A,
                           PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int i, j;
    int ldai, ldaj;
    int tempim;
    int tempjm, tempjn;
    int minmn = min(A.mt, A.nt);

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;

    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    if (uplo == PlasmaLower) {
       for (j = 0; j < minmn; j++){
           tempjm = j == A.mt-1 ? A.m-j*A.mb : A.mb;
           tempjn = j == A.nt-1 ? A.n-j*A.nb : A.nb;
           ldaj = BLKLDD(A, j);
           QUARK_CORE_zlaset2(
               plasma->quark, &task_flags,
               PlasmaLower, tempjm, tempjn, alpha,
               A(j, j), ldaj);

           for (i = j+1; i < A.mt; i++){
               tempim = i == A.mt-1 ? A.m-i*A.mb : A.mb;
               ldai = BLKLDD(A, i);
               QUARK_CORE_zlaset2(
                   plasma->quark, &task_flags,
                   PlasmaUpperLower, tempim, tempjn, alpha,
                   A(i, j), ldai);
           }
       }
    }
    else if (uplo == PlasmaUpper) {
       for (j = 1; j < A.nt; j++){
           tempjn = j == A.nt-1 ? A.n-j*A.nb : A.nb;
           for (i = 0; i < min(j, A.mt); i++){
               tempim = i == A.mt-1 ? A.m-i*A.mb : A.mb;
               ldai = BLKLDD(A, i);
               QUARK_CORE_zlaset2(
                   plasma->quark, &task_flags,
                   PlasmaUpperLower, tempim, tempjn, alpha,
                   A(i, j), ldai);
           }
       }
       for (j = 0; j < minmn; j++){
           tempjm = j == A.mt-1 ? A.m-j*A.mb : A.mb;
           tempjn = j == A.nt-1 ? A.n-j*A.nb : A.nb;
           ldaj = BLKLDD(A, j);
           QUARK_CORE_zlaset2(
               plasma->quark, &task_flags,
               PlasmaUpper, tempjm, tempjn, alpha,
               A(j, j), ldaj);
       }
    }
    else {
       for (i = 0; i < A.mt; i++){
           tempim = i == A.mt-1 ? A.m-i*A.mb : A.mb;
           ldai = BLKLDD(A, i);
           for (j = 0; j < A.nt; j++){
               tempjn = j == A.nt-1 ? A.n-j*A.nb : A.nb;
               QUARK_CORE_zlaset2(
                   plasma->quark, &task_flags,
                   PlasmaUpperLower, tempim, tempjn, alpha,
                   A(i, j), ldai);
           }
       }
    } 
}
