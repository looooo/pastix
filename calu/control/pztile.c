/**
 *
 * @file ztile.c
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "auxiliary.h"
#include "tile.h"
#include "quark.h"

#define AF77(m, n) &(Af77[ ((int64_t)A.nb*(int64_t)lda*(int64_t)(n)) + (int64_t)(A.mb*(m)) ])
#define ABDL(m, n) BLKADDR(A, PLASMA_Complex64_t, m, n)

void CORE_ztile_zero_quark(Quark* quark);

/***************************************************************************//**
 *  Conversion from LAPACK F77 matrix layout to tile layout - static scheduling
 **/
void plasma_pzlapack_to_tile(plasma_context_t *plasma)
{
    PLASMA_Complex64_t *Af77;
    int lda;
    PLASMA_desc A;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    PLASMA_Complex64_t *f77;
    PLASMA_Complex64_t *bdl;

    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;
    int next_m;
    int next_n;

    plasma_unpack_args_5(Af77, lda, A, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    while (m >= A.mt && n < A.nt) {
        n++;
        m = m-A.mt;
    }

    while (n < A.nt) {
        next_m = m;
        next_n = n;

        next_m += PLASMA_SIZE;
        while (next_m >= A.mt && next_n < A.nt) {
            next_n++;
            next_m = next_m-A.mt;
        }

        X1 = n == 0 ? A.j%A.nb : 0;
        X2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
        Y1 = m == 0 ? A.i%A.mb : 0;
        Y2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;

        f77 = AF77(m, n);
        bdl = ABDL(m, n);
        ldt = BLKLDD(A, m);
        CORE_zlacpy(
            PlasmaUpperLower, (Y2-Y1), (X2-X1),
            &(f77[X1*lda+Y1]), lda, 
            &(bdl[X1*lda+Y1]), ldt);
        
        m = next_m;
        n = next_n;
    }
}

/***************************************************************************//**
 *  Conversion from LAPACK F77 matrix layout to tile layout - dynamic scheduling
 **/
void plasma_pzlapack_to_tile_quark(PLASMA_Complex64_t *Af77, int lda, PLASMA_desc A,
                                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_Complex64_t *f77;
    PLASMA_Complex64_t *bdl;
    plasma_context_t *plasma;
    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for (m = 0; m < A.mt; m++)
    {
        ldt = BLKLDD(A, m);
        for (n = 0; n < A.nt; n++)
        {
            X1 = n == 0 ? A.j%A.nb : 0;
            Y1 = m == 0 ? A.i%A.mb : 0;
            X2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
            Y2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;

            f77 = AF77(m, n);
            bdl = ABDL(m, n);
            QUARK_CORE_zlacpy(
                plasma->quark, &task_flags,
                PlasmaUpperLower, (Y2-Y1), (X2-X1), A.mb,
                &(f77[X1*lda+Y1]), lda, 
                &(bdl[X1*lda+Y1]), ldt);
        }
    }
}

/***************************************************************************//**
 *  Conversion from LAPACK F77 matrix layout to tile layout - static scheduling
 **/
void plasma_pztile_to_lapack(plasma_context_t *plasma)
{
    PLASMA_desc A;
    PLASMA_Complex64_t *Af77;
    int lda;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    PLASMA_Complex64_t *f77;
    PLASMA_Complex64_t *bdl;

    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;
    int next_m;
    int next_n;

    plasma_unpack_args_5(A, Af77, lda, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    while (m >= A.mt && n < A.nt) {
        n++;
        m = m-A.mt;
    }

    while (n < A.nt) {
        next_m = m;
        next_n = n;

        next_m += PLASMA_SIZE;
        while (next_m >= A.mt && next_n < A.nt) {
            next_n++;
            next_m = next_m-A.mt;
        }

        X1 = n == 0 ? A.j%A.nb : 0;
        Y1 = m == 0 ? A.i%A.mb : 0;
        X2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
        Y2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;

        f77 = AF77(m, n);
        bdl = ABDL(m, n);
        ldt = BLKLDD(A, m);
        CORE_zlacpy(
            PlasmaUpperLower, (Y2-Y1), (X2-X1),
            &(bdl[X1*lda+Y1]), ldt,
            &(f77[X1*lda+Y1]), lda);

        m = next_m;
        n = next_n;
    }
}

/***************************************************************************//**
 *  Conversion from LAPACK F77 matrix layout to tile layout - dynamic scheduling
 **/
void plasma_pztile_to_lapack_quark(PLASMA_desc A, PLASMA_Complex64_t *Af77, int lda,
                                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_Complex64_t *f77;
    PLASMA_Complex64_t *bdl;
    plasma_context_t *plasma;
    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for (m = 0; m < A.mt; m++)
    {
        ldt = BLKLDD(A, m);
        for (n = 0; n < A.nt; n++)
        {
            X1 = n == 0 ? A.j%A.nb : 0;
            Y1 = m == 0 ? A.i%A.mb : 0;
            X2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
            Y2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;

            f77 = AF77(m, n);
            bdl = ABDL(m, n);
            QUARK_CORE_zlacpy(
                plasma->quark, &task_flags,
                PlasmaUpperLower, (Y2-Y1), (X2-X1), A.mb,
                &(bdl[X1*lda+Y1]), ldt,
                &(f77[X1*lda+Y1]), lda);
        }
    }
}

/***************************************************************************//**
 *  Zeroes a submatrix in tile layout - static scheduling
 **/
void plasma_pztile_zero(plasma_context_t *plasma)
{
    PLASMA_desc A;
    PLASMA_sequence *sequence;
    PLASMA_request *request;

    PLASMA_Complex64_t *bdl;
    int x, y;
    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;
    int next_m;
    int next_n;

    plasma_unpack_args_3(A, sequence, request);
    if (sequence->status != PLASMA_SUCCESS)
        return;

    n = 0;
    m = PLASMA_RANK;
    while (m >= A.mt && n < A.nt) {
        n++;
        m = m-A.mt;
    }

    while (n < A.nt) {
        next_m = m;
        next_n = n;

        next_m += PLASMA_SIZE;
        while (next_m >= A.mt && next_n < A.nt) {
            next_n++;
            next_m = next_m-A.mt;
        }

        X1 = n == 0 ? A.j%A.nb : 0;
        Y1 = m == 0 ? A.i%A.mb : 0;
        X2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
        Y2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;

        bdl = ABDL(m, n);
        ldt = BLKLDD(A, m);
        for (x = X1; x < X2; x++)
            for (y = Y1; y < Y2; y++)
                bdl[ldt*x+y] = 0.0;

        m = next_m;
        n = next_n;
    }
}

/***************************************************************************//**
 *  Zeroes a submatrix in tile layout - dynamic scheduling
 **/
void plasma_pztile_zero_quark(PLASMA_desc A, PLASMA_sequence *sequence, PLASMA_request *request)
{
    PLASMA_Complex64_t *bdl;
    plasma_context_t *plasma;
    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;

    plasma = plasma_context_self();
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for (m = 0; m < A.mt; m++)
    {
        ldt = BLKLDD(A, m);
        for (n = 0; n < A.nt; n++)
        {
            X1 = n == 0 ? A.j%A.nb : 0;
            Y1 = m == 0 ? A.i%A.mb : 0;
            X2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
            Y2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;

            bdl = ABDL(m, n);
            QUARK_Insert_Task(plasma->quark, CORE_ztile_zero_quark, &task_flags,
                sizeof(int),                       &X1,  VALUE,
                sizeof(int),                       &X2,  VALUE,
                sizeof(int),                       &Y1,  VALUE,
                sizeof(int),                       &Y2,  VALUE,
                sizeof(PLASMA_Complex64_t)*A.bsiz, bdl,      OUTPUT | LOCALITY,
                sizeof(int),                       &ldt, VALUE,
                0);
        }
    }
}

/***************************************************************************//**
 *
 **/
void CORE_ztile_zero_quark(Quark* quark)
{
    int X1;
    int X2;
    int Y1;
    int Y2;
    PLASMA_Complex64_t *A;
    int lda;

    int x, y;

    quark_unpack_args_6(quark, X1, X2, Y1, Y2, A, lda);

    for (x = X1; x < X2; x++)
        for (y = Y1; y < Y2; y++)
            A[lda*x+y] = 0.0;
}
