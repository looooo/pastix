/**
 *
 * @file pztrsm.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define W(m)   BLKADDR(W, PLASMA_Complex64_t, m, 0)

/***************************************************************************//**
 *  Parallel tile ReVerse triangular solve - dynamic scheduling
 **/
void plasma_pztrsmrv_quark(PLASMA_enum side, PLASMA_enum uplo, PLASMA_enum trans, PLASMA_enum diag,
                           PLASMA_Complex64_t alpha, PLASMA_desc A, PLASMA_desc W, 
                           PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m;
    int ldak, ldam, ldwk, ldwm;
    int tempm, tempkm, tempkn, tempmm;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for (k = A.mt-1; k >= 0; k--) {
        tempm = (k+1)*A.mb;
        tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
        tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;
        ldak = BLKLDD(A, k);
        ldwk = BLKLDD(W, k);
        
        QUARK_CORE_zlacpy(
            plasma->quark, &task_flags,
            PlasmaLower, tempkm, tempkn, A.mb,
            A(k, k), ldak, W(k), ldwk );

        QUARK_CORE_zlaset2(
            plasma->quark, &task_flags,
            PlasmaLower, tempkn, tempkn, 
            0.0, A(k, k), ldak );

        for(m=k+1; m<A.mt; m++) {
            tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            ldwm = BLKLDD(W, m);

            QUARK_CORE_zlacpy(
                plasma->quark, &task_flags,
                PlasmaUpperLower, tempmm, tempkn, A.mb,
                A(m, k), ldam, W(m), ldwm );
            
            QUARK_CORE_zlaset(
                plasma->quark, &task_flags,
                PlasmaUpperLower, tempmm, tempkn,
                0.0, 0.0, A(m, k), ldam );
        }            

        if (k*A.mb+tempkn < A.m) {
           plasma_pzgemm_quark(
               PlasmaNoTrans, PlasmaNoTrans, -alpha,
               plasma_desc_submatrix(A, 0,     tempm,  tempm,     A.n-tempm),
               plasma_desc_submatrix(W, tempm, 0,      W.m-tempm, tempkn   ),
               alpha,
               plasma_desc_submatrix(A, 0,     k*A.mb, tempm,     tempkn   ),
               sequence, request);

           plasma_pzgemm_quark(
               PlasmaNoTrans, PlasmaNoTrans, -alpha,
               plasma_desc_submatrix(A, tempm, tempm,  A.m-tempm, A.n-tempm),
               plasma_desc_submatrix(W, tempm, 0,      W.m-tempm, tempkn   ),
               (PLASMA_Complex64_t)0.,
               plasma_desc_submatrix(A, tempm, k*A.mb, A.m-tempm, tempkn   ),
               sequence, request);
        }

        plasma_pztrsm_quark(
            PlasmaRight, PlasmaLower, 
            PlasmaNoTrans, PlasmaUnit, 
            alpha,
            plasma_desc_submatrix(W, k*W.mb, 0,      tempkm, tempkn),
            plasma_desc_submatrix(A, 0,      k*A.mb, A.m,    tempkn),
            sequence, request);
    }
}
