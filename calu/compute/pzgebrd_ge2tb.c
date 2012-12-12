/**
 *
 * @file pzgebrd_ge2tb.c
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

/***************************************************************************//**
 *  Parallel tile BAND Bidiagonal Reduction - static scheduler
 *  Could be optimized by using the algorithms from Trefethen book
 *
 * WARNING: do never call this function because unmqr and unmlq are
 * not implementing all the cases required in static.
 *
 **/
void plasma_pzgebrd_ge2tb(plasma_context_t *plasma)
{
    PLASMA_desc A;
    PLASMA_desc T;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int k;
    int tempkm, tempkn;

    plasma_unpack_args_4(A, T, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    if (A.m >= A.n){
       for (k = 0; k < A.nt; k++) {
           tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
           tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;

           plasma_static_call_4(plasma_pzgeqrf,
               PLASMA_desc, plasma_desc_submatrix(A, k*A.mb, k*A.nb, A.m-k*A.mb, tempkn),
               PLASMA_desc, plasma_desc_submatrix(T, k*T.mb, k*T.nb, T.m-k*T.mb, tempkn),
               PLASMA_sequence*, sequence,
               PLASMA_request*, request);

           plasma_static_call_7(plasma_pzunmqr,
               PLASMA_enum, PlasmaLeft,
               PLASMA_enum, PlasmaConjTrans,
               PLASMA_desc, plasma_desc_submatrix(A, k*A.mb,     k*A.nb, A.m-k*A.mb, tempkn),
               PLASMA_desc, plasma_desc_submatrix(A, k*A.mb, (k+1)*A.nb, A.m-k*A.mb, A.n-(k+1)*A.nb),
               PLASMA_desc, plasma_desc_submatrix(T, k*T.mb,     k*T.nb, T.m-k*T.mb, tempkn),
               PLASMA_sequence*, sequence,
               PLASMA_request*, request);

           if (k+1 < A.nt){
              tempkn = k+1 == A.nt-1 ? A.n-(k+1)*A.nb : A.nb;

              plasma_static_call_4(plasma_pzgelqf,
                  PLASMA_desc, plasma_desc_submatrix(A, k*A.mb, (k+1)*A.nb, tempkm, A.n-(k+1)*A.nb),
                  PLASMA_desc, plasma_desc_submatrix(T, k*T.mb, (k+1)*T.nb, tempkm, T.n-(k+1)*T.nb),
                  PLASMA_sequence*, sequence,
                  PLASMA_request*, request);

              plasma_static_call_7(plasma_pzunmlq,
                  PLASMA_enum, PlasmaRight,
                  PLASMA_enum, PlasmaConjTrans,
                  PLASMA_desc, plasma_desc_submatrix(A,     k*A.mb, (k+1)*A.nb, tempkm,         A.n-(k+1)*A.nb),
                  PLASMA_desc, plasma_desc_submatrix(A, (k+1)*A.mb, (k+1)*A.nb, A.m-(k+1)*A.mb, A.n-(k+1)*A.nb),
                  PLASMA_desc, plasma_desc_submatrix(T,     k*T.mb, (k+1)*T.nb, tempkm,         T.n-(k+1)*T.nb),
                  PLASMA_sequence*, sequence,
                  PLASMA_request*, request);
           }
       }
    }
    else{
       for (k = 0; k < A.mt; k++) {
           tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
           tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;

           plasma_static_call_4(plasma_pzgelqf,
               PLASMA_desc, plasma_desc_submatrix(A, k*A.mb, k*A.nb, tempkm, A.n-k*A.nb),
               PLASMA_desc, plasma_desc_submatrix(T, k*T.mb, k*T.nb, tempkm, T.n-k*T.nb),
               PLASMA_sequence*, sequence,
               PLASMA_request*, request);

           plasma_static_call_7(plasma_pzunmlq,
               PLASMA_enum, PlasmaRight,
               PLASMA_enum, PlasmaConjTrans,
               PLASMA_desc, plasma_desc_submatrix(A,     k*A.mb, k*A.nb, tempkm,         A.n-k*A.nb),
               PLASMA_desc, plasma_desc_submatrix(A, (k+1)*A.mb, k*A.nb, A.m-(k+1)*A.mb, A.n-k*A.nb),
               PLASMA_desc, plasma_desc_submatrix(T,     k*T.mb, k*T.nb, tempkm,         T.n-k*T.nb),
               PLASMA_sequence*, sequence,
               PLASMA_request*, request);

           if (k+1 < A.mt){
              tempkm = k+1 == A.mt-1 ? A.m-(k+1)*A.mb : A.mb;
              tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;

              plasma_static_call_4(plasma_pzgeqrf,
                  PLASMA_desc, plasma_desc_submatrix(A, (k+1)*A.mb, k*A.nb, A.m-(k+1)*A.mb, tempkn),
                  PLASMA_desc, plasma_desc_submatrix(T, (k+1)*T.mb, k*T.nb, T.m-(k+1)*T.mb, tempkn),
                  PLASMA_sequence*, sequence,
                  PLASMA_request*, request);

              plasma_static_call_7(plasma_pzunmqr,
                  PLASMA_enum, PlasmaLeft,
                  PLASMA_enum, PlasmaConjTrans,
                  PLASMA_desc, plasma_desc_submatrix(A, (k+1)*A.mb,     k*A.nb, A.m-(k+1)*A.mb, tempkn),
                  PLASMA_desc, plasma_desc_submatrix(A, (k+1)*A.mb, (k+1)*A.nb, A.m-(k+1)*A.mb, A.n-(k+1)*A.nb),
                  PLASMA_desc, plasma_desc_submatrix(T, (k+1)*T.mb,     k*T.nb, T.m-(k+1)*T.mb, tempkn),
                  PLASMA_sequence*, sequence,
                  PLASMA_request*, request);
           }
       }
    }
}

/***************************************************************************//**
 *  Parallel tile BAND Bidiagonal Reduction - dynamic scheduler
 *  Could be optimized by using the algorithms from Trefethen book
 **/
void plasma_pzgebrd_ge2tb_quark(PLASMA_desc A, PLASMA_desc T,
                                PLASMA_sequence *sequence, PLASMA_request *request)
{
    int k;
    int tempkm, tempkn;

    if (A.m >= A.n){
       for (k = 0; k < A.nt; k++) {
           tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
           tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;

           plasma_pzgeqrf_quark(
               plasma_desc_submatrix(A, k*A.mb, k*A.nb, A.m-k*A.mb, tempkn),
               plasma_desc_submatrix(T, k*T.mb, k*T.nb, T.m-k*T.mb, tempkn),
               sequence, request);

           plasma_pzunmqr_quark(
               PlasmaLeft,
               PlasmaConjTrans,
               plasma_desc_submatrix(A, k*A.mb,     k*A.nb, A.m-k*A.mb, tempkn),
               plasma_desc_submatrix(A, k*A.mb, (k+1)*A.nb, A.m-k*A.mb, A.n-(k+1)*A.nb),
               plasma_desc_submatrix(T, k*T.mb,     k*T.nb, T.m-k*T.mb, tempkn),
               sequence, request);

           if (k+1 < A.nt){
              tempkn = k+1 == A.nt-1 ? A.n-(k+1)*A.nb : A.nb;

              plasma_pzgelqf_quark(
                  plasma_desc_submatrix(A, k*A.mb, (k+1)*A.nb, tempkm, A.n-(k+1)*A.nb),
                  plasma_desc_submatrix(T, k*T.mb, (k+1)*T.nb, T.mb,   T.n-(k+1)*T.nb),
                  sequence, request);

              plasma_pzunmlq_quark(
                  PlasmaRight, PlasmaConjTrans,
                  plasma_desc_submatrix(A,     k*A.mb, (k+1)*A.nb, tempkm,         A.n-(k+1)*A.nb),
                  plasma_desc_submatrix(A, (k+1)*A.mb, (k+1)*A.nb, A.m-(k+1)*A.mb, A.n-(k+1)*A.nb),
                  plasma_desc_submatrix(T,     k*T.mb, (k+1)*T.nb, T.mb,           T.n-(k+1)*T.nb),
                  sequence, request);
           }
       }
    }
    else{
       for (k = 0; k < A.mt; k++) {
           tempkm = k == A.mt-1 ? A.m-k*A.mb : A.mb;
           tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;

           plasma_pzgelqf_quark(
               plasma_desc_submatrix(A, k*A.mb, k*A.nb, tempkm, A.n-k*A.nb),
               plasma_desc_submatrix(T, k*T.mb, k*T.nb, T.mb,   T.n-k*T.nb),
               sequence, request);

           plasma_pzunmlq_quark(
               PlasmaRight, PlasmaConjTrans,
               plasma_desc_submatrix(A,     k*A.mb, k*A.nb, tempkm,         A.n-k*A.nb),
               plasma_desc_submatrix(A, (k+1)*A.mb, k*A.nb, A.m-(k+1)*A.mb, A.n-k*A.nb),
               plasma_desc_submatrix(T,     k*T.mb, k*T.nb, T.mb,           T.n-k*T.nb),
               sequence, request);

           if (k+1 < A.mt){
              tempkm = k+1 == A.mt-1 ? A.m-(k+1)*A.mb : A.mb;
              tempkn = k == A.nt-1 ? A.n-k*A.nb : A.nb;

              plasma_pzgeqrf_quark(
                   plasma_desc_submatrix(A, (k+1)*A.mb, k*A.nb, A.m-(k+1)*A.mb, tempkn),
                   plasma_desc_submatrix(T, (k+1)*T.mb, k*T.nb, T.m-(k+1)*T.mb, tempkn),
                   sequence, request);

              plasma_pzunmqr_quark(
                  PlasmaLeft, PlasmaConjTrans,
                  plasma_desc_submatrix(A, (k+1)*A.mb,     k*A.nb, A.m-(k+1)*A.mb, tempkn),
                  plasma_desc_submatrix(A, (k+1)*A.mb, (k+1)*A.nb, A.m-(k+1)*A.mb, A.n-(k+1)*A.nb),
                  plasma_desc_submatrix(T, (k+1)*T.mb,     k*T.nb, T.m-(k+1)*T.mb, tempkn),
                  sequence, request);
           }
       }
    }
}
