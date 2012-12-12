/**
 *
 * @file pzlacpy.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Emmanuel Agullo
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"

#define A(m,n) BLKADDR(A, PLASMA_Complex64_t, m, n)
#define B(m,n) BLKADDR(B, PLASMA_Complex64_t, m, n)
/***************************************************************************//**
 *
 **/
void plasma_pzlacpy(plasma_context_t *plasma)
{
    PLASMA_enum uplo;
    PLASMA_desc A;
    PLASMA_desc B;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    int X, Y;
    int m, n;
    int next_m;
    int next_n;
    int ldam, ldbm;

    plasma_unpack_args_5(uplo, A, B, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    switch (uplo) {
    /*
     *  PlasmaUpper
     */
    case PlasmaUpper:
        m = 0;
        n = PLASMA_RANK;
        while (n >= A.nt) {
            m++;
            n = n - A.nt + m;
        }

        while (m < A.mt) {
            next_m = m;
            next_n = n;

            next_n += PLASMA_SIZE;
            while (next_n >= A.nt && next_m < A.mt) {
                next_m++;
                next_n = next_n - A.nt + next_m;
            }

            X = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            Y = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            ldam = BLKLDD(A, m);
            ldbm = BLKLDD(B, m);
            CORE_zlacpy(
                m == n ? uplo : PlasmaUpperLower,
                X, Y,
                A(m, n), ldam,
                B(m, n), ldbm);

            n = next_n;
            m = next_m;
        }
        break;
    /*
     *  PlasmaLower
     */
    case PlasmaLower:
        n = 0;
        m = PLASMA_RANK;
        while (m >= A.mt) {
            n++;
            m = m - A.mt + n;
        }

        while (n < A.nt) {
            next_m = m;
            next_n = n;

            next_m += PLASMA_SIZE;
            while (next_m >= A.mt && next_n < A.nt) {
                next_n++;
                next_m = next_m - A.mt + next_n;
            }

            X = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            Y = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            ldam = BLKLDD(A, m);
            ldbm = BLKLDD(B, m);
            CORE_zlacpy(
                m == n ? uplo : PlasmaUpperLower,
                X, Y,
                A(m, n), ldam,
                B(m, n), ldbm);

            n = next_n;
            m = next_m;
        }
        break;
    /*
     *  PlasmaUpperLower
     */
    case PlasmaUpperLower:
    default:
        n = 0;
        m = PLASMA_RANK;
        while (m >= A.mt) {
            n++;
            m = m - A.mt;
        }

        while (n < A.nt) {
            next_m = m;
            next_n = n;

            next_m += PLASMA_SIZE;
            while (next_m >= A.mt && next_n < A.nt) {
                next_n++;
                next_m = next_m - A.mt;
            }

            X = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            Y = n == A.nt-1 ? A.n-n*A.nb : A.nb;
            ldam = BLKLDD(A, m);
            ldbm = BLKLDD(B, m);
            CORE_zlacpy(
                PlasmaUpperLower,
                X, Y,
                A(m, n), ldam,
                B(m, n), ldbm);

            n = next_n;
            m = next_m;
        }
        break;
    }
}
/***************************************************************************//**
 *
 **/
void plasma_pzlacpy_quark(PLASMA_enum uplo, PLASMA_desc A, PLASMA_desc B,
                          PLASMA_sequence *sequence, PLASMA_request *request)
{
    plasma_context_t *plasma;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    int X, Y;
    int m, n;
    int ldam, ldbm;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    switch (uplo) {
    /*
     *  PlasmaUpper
     */
    case PlasmaUpper:
        for (m = 0; m < A.mt; m++) {
            X = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            ldbm = BLKLDD(B, m);
            if (m < A.nt) {
                Y = m == A.nt-1 ? A.n-m*A.nb : A.nb;
                QUARK_CORE_zlacpy(
                    plasma->quark, &task_flags,
                    PlasmaUpper,
                    X, Y, A.mb,
                    A(m, m), ldam,
                    B(m, m), ldbm);
            }
            for (n = m+1; n < A.nt; n++) {
                Y = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_zlacpy(
                    plasma->quark, &task_flags,
                    PlasmaUpperLower,
                    X, Y, A.mb,
                    A(m, n), ldam,
                    B(m, n), ldbm);
            }
        }
        break;
    /*
     *  PlasmaLower
     */
    case PlasmaLower:
        for (m = 0; m < A.mt; m++) {
            X = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            ldbm = BLKLDD(B, m);
            if (m < A.nt) {
                Y = m == A.nt-1 ? A.n-m*A.nb : A.nb;
                QUARK_CORE_zlacpy(
                    plasma->quark, &task_flags,
                    PlasmaLower,
                    X, Y, A.mb,
                    A(m, m), ldam,
                    B(m, m), ldbm);
            }
            for (n = 0; n < min(m, A.nt); n++) {
                Y = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_zlacpy(
                    plasma->quark, &task_flags,
                    PlasmaUpperLower,
                    X, Y, A.mb,
                    A(m, n), ldam,
                    B(m, n), ldbm);
            }
        }
        break;
    /*
     *  PlasmaUpperLower
     */
    case PlasmaUpperLower:
    default:
        for (m = 0; m < A.mt; m++) {
            X = m == A.mt-1 ? A.m-m*A.mb : A.mb;
            ldam = BLKLDD(A, m);
            ldbm = BLKLDD(B, m);
            for (n = 0; n < A.nt; n++) {
                Y = n == A.nt-1 ? A.n-n*A.nb : A.nb;
                QUARK_CORE_zlacpy(
                    plasma->quark, &task_flags,
                    PlasmaUpperLower,
                    X, Y, A.mb,
                    A(m, n), ldam,
                    B(m, n), ldbm);
            }
        }
    }
}
