/**
 * 
 * @file core_alloc.c
 *
 *  PLASMA core_blas kernel
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/
#include <stdlib.h>
#include "common.h"

/***************************************************************************//**
 *
 **/
void CORE_free_quark(Quark *quark)
{
    void *A;

    quark_unpack_args_1(quark, A);
    if (A != NULL)
        free(A);
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_free(Quark *quark, Quark_Task_Flags *task_flags, void *A, int szeA)
{
    QUARK_Insert_Task(
        quark, CORE_free_quark, task_flags,
        szeA, A, INOUT,
        0);
}

void CORE_foo_quark(Quark *quark) {
    void *A;
    quark_unpack_args_1(quark, A);
}

void CORE_foo2_quark(Quark *quark) {
    void *A, *B;
    quark_unpack_args_2(quark, A, B);
}

/***************************************************************************//**
 *
 **/
void CORE_pivot_update_quark(Quark *quark)
{
    int i, m, n, offset, init;
    int *indices;
    int *ipiv;

    quark_unpack_args_6(quark, m, n, indices, ipiv, offset, init);

    if ( init ) {
        for(i=0; i<m; i++) {
            indices[i] = offset+i;
        }
    }
    for(i=0; i<n; i++) {
        ipiv[i] = indices[ ipiv[i]-1 ]+1;
    }
}

/***************************************************************************//**
 *
 **/
void QUARK_CORE_pivot_update(Quark *quark, Quark_Task_Flags *task_flags, 
                             int m, int n, int *indices, int *ipiv, 
                             int offset, int init)
{
    DAG_SET_PROPERTIES( "PIV_UP"  , "white"   );
    QUARK_Insert_Task(quark, CORE_pivot_update_quark, task_flags,
        sizeof(int),  &m,       VALUE,
        sizeof(int),  &n,       VALUE,
        sizeof(int)*m, indices,     INPUT,
        sizeof(int)*n, ipiv,        INOUT,
        sizeof(int),  &offset,  VALUE,
        sizeof(int),  &init,    VALUE,
        0);
}
