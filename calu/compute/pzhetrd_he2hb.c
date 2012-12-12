/**
 *
 * @file pzhetrd_he2hb.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define T(m,n) BLKADDR(T, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *  Parallel tile BAND Tridiagonal Reduction - dynamic scheduler
 **/
void plasma_pzhetrd_he2hb_quark(PLASMA_enum uplo, 
                          PLASMA_desc A, PLASMA_desc T,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int k, m, n, i, j;
    int ldak, ldam, ldan, ldaj, ldai;
    int tempkn, tempmm, tempnn, tempjj;
    int ib;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;

    /* CASE NB>N only 1 tile*/
    if(A.mt>A.m)
        return;


    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    ib = PLASMA_IB;
    if (uplo == PlasmaLower) {
       for (k = 0; k < A.nt-1; k++){
           tempkn = k+1 == A.nt-1 ? A.n-(k+1)*A.nb : A.nb;
           ldak = BLKLDD(A, k+1);
           QUARK_CORE_zgeqrt(
               plasma->quark, &task_flags,
               tempkn, A.nb, ib, T.nb,
               A(k+1, k), ldak,
               T(k+1, k), T.mb);

           /* LEFT and RIGHT on the symmetric diagonal block */
           QUARK_CORE_zherfb(
               plasma->quark, &task_flags,
               PlasmaLower,
               tempkn, tempkn, ib, T.nb,
               A(k+1,   k), ldak,
               T(k+1,   k), T.mb,
               A(k+1, k+1), ldak);

           /* RIGHT on the remaining tiles until the bottom */
           for (m = k+2; m < A.mt ; m++) {
               tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
               ldam = BLKLDD(A, m);
               QUARK_CORE_zunmqr(
                   plasma->quark, &task_flags,
                   PlasmaRight, PlasmaNoTrans,
                   tempmm, A.nb, tempkn, ib, T.nb,
                   A(k+1,   k), ldak,
                   T(k+1,   k), T.mb,
                   A(m  , k+1), ldam);
           }

           for (m = k+2; m < A.mt; m++) {
               tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
               ldam = BLKLDD(A, m);
               QUARK_CORE_ztsqrt(
                   plasma->quark, &task_flags,
                   tempmm, A.nb, ib, T.nb,
                   A(k+1, k), ldak,
                   A(m  , k), ldam,
                   T(m  , k), T.mb);

               /* LEFT */
               for (i = k+2; i < m; i++) {
                   ldai = BLKLDD(A, i);
                   QUARK_CORE_ztsmqr_hetra1(
                       plasma->quark, &task_flags,
                       PlasmaLeft, PlasmaConjTrans,
                       A.mb, A.nb, tempmm, A.nb, A.nb, ib, T.nb,
                       A(i, k+1), ldai,
                       A(m,   i), ldam,
                       A(m,   k), ldam,
                       T(m,   k), T.mb);
               }

               /* RIGHT */
               for (j = m+1; j < A.mt ; j++) {
                   tempjj = j == A.mt-1 ? A.m-j*A.mb : A.mb;
                   ldaj = BLKLDD(A, j);
                   QUARK_CORE_ztsmqr(
                       plasma->quark, &task_flags,
                       PlasmaRight, PlasmaNoTrans,
                       tempjj, A.nb, tempjj, tempmm, A.nb, ib, T.nb,
                       A(j, k+1), ldaj,
                       A(j,   m), ldaj,
                       A(m,   k), ldam,
                       T(m,   k), T.mb);
               }
       
               /* LEFT->RIGHT */
               QUARK_CORE_ztsmqr_corner(
                   plasma->quark, &task_flags,
                   A.nb, A.nb, tempmm, A.nb, tempmm, tempmm, A.nb, ib, T.nb,
                   A(k+1, k+1), ldak,
                   A(m  , k+1), ldam,
                   A(m  ,   m), ldam,
                   A(m  ,   k), ldam,
                   T(m  ,   k), T.mb);
           }
       }
    }
    else {
       for (k = 0; k < A.nt-1; k++){
           tempkn = k+1 == A.nt-1 ? A.n-(k+1)*A.nb : A.nb;
           ldak = BLKLDD(A, k+1);
           QUARK_CORE_zgelqt(
               plasma->quark, &task_flags,
               A.nb, tempkn, ib, T.nb,
               A(k, k+1), A.nb,
               T(k, k+1), T.mb);

           /* RIGHT and LEFT on the symmetric diagonal block             */
           QUARK_CORE_zherfb(
               plasma->quark, &task_flags,
               PlasmaUpper,
               tempkn, tempkn, ib, T.nb,
               A(k,   k+1), A.nb,
               T(k,   k+1), T.mb,
               A(k+1, k+1), ldak);

           /* LEFT on the remaining tiles until the left side */
           for (n = k+2; n < A.nt ; n++) {
               tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
               QUARK_CORE_zunmlq(
                   plasma->quark, &task_flags,
                   PlasmaLeft, PlasmaNoTrans,
                   A.nb, tempnn, tempkn, ib, T.nb,
                   A(k,   k+1), A.nb,
                   T(k,   k+1), T.mb,
                   A(k+1,   n), ldak);
           }

           for (n = k+2; n < A.nt; n++) {
               tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
               ldan = BLKLDD(A, n);
               QUARK_CORE_ztslqt(
                   plasma->quark, &task_flags,
                   A.nb, tempnn, ib, T.nb,
                   A(k, k+1), A.nb,
                   A(k,   n), A.nb,
                   T(k,   n), T.mb);

               /* RIGHT */
               for (i = k+2; i < n; i++) {
                   ldai = BLKLDD(A, i);
                   QUARK_CORE_ztsmlq_hetra1(
                       plasma->quark, &task_flags,
                       PlasmaRight, PlasmaConjTrans,
                       A.mb, A.nb, A.nb, tempnn, A.nb, ib, T.nb,
                       A(k+1, i), ldak,
                       A(i,   n), ldai,
                       A(k,   n), A.nb,
                       T(k,   n), T.mb);
               }

               /* LEFT */
               for (j = n+1; j < A.nt ; j++) {
                   tempjj = j == A.nt-1 ? A.n-j*A.nb : A.nb;
                   ldaj = BLKLDD(A, j);
                   QUARK_CORE_ztsmlq(
                       plasma->quark, &task_flags,
                       PlasmaLeft, PlasmaNoTrans,
                       A.nb, tempjj, tempnn, tempjj, A.nb, ib, T.nb,
                       A(k+1, j), ldak,
                       A(n,   j), ldan,
                       A(k,   n), A.nb,
                       T(k,   n), T.mb);
               }
       
               /* RIGHT->LEFT */
               QUARK_CORE_ztsmlq_corner(
                   plasma->quark, &task_flags,
                   A.nb, A.nb, A.nb, tempnn, tempnn, tempnn, A.nb, ib, T.nb,
                   A(k+1, k+1), ldak,
                   A(k+1,   n), ldak,
                   A(n  ,   n), ldan,
                   A(k  ,   n), A.nb,
                   T(k  ,   n), T.mb);
           }
       }
    }
}
