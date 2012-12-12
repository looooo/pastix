/**
 *
 * @file descriptor.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_DESCRIPTOR_H_
#define _PLASMA_DESCRIPTOR_H_

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Internal routines
 **/
inline static void *plasma_geteltaddr( const PLASMA_desc *A, int m, int n, int eltsize);
inline static void *plasma_getaddr(PLASMA_desc A, int m, int n);
PLASMA_desc plasma_desc_init(PLASMA_enum dtyp, int mb, int nb, int bsiz, int lm, int ln, int i, int j, int m, int n);
PLASMA_desc plasma_desc_submatrix(PLASMA_desc descA, int i, int j, int m, int n);
int plasma_desc_check(PLASMA_desc *desc);
int plasma_desc_mat_alloc(PLASMA_desc *desc);
int plasma_desc_mat_free(PLASMA_desc *desc);

/***************************************************************************//**
 *  Internal function to return adress of block (m,n)
 **/
inline static void *plasma_getaddr(PLASMA_desc A, int m, int n)
{
    size_t mm = m+A.i/A.mb;
    size_t nn = n+A.j/A.nb;
    size_t eltsize = plasma_element_size(A.dtyp);
    size_t offset = 0;

    if (mm < A.lm1) {
        if (nn < A.ln1)
            offset = A.bsiz*(mm+A.lm1*nn);
        else
            offset = A.A12 + (A.mb*(A.ln%A.nb)*mm);
    }
    else {
        if (nn < A.ln1)
            offset = A.A21 + ((A.lm%A.mb)*A.nb*nn);
        else
            offset = A.A22;
    }

    return (void*)((intptr_t)A.mat + (offset*eltsize) );
}

/***************************************************************************//**
 *  Internal function to return adress of element A(m,n)
 **/
inline static void *plasma_geteltaddr( const PLASMA_desc *A, int m, int n, int eltsize)
{
    size_t mm = m/A->mb;
    size_t nn = n/A->nb;
    size_t offset = 0;

    if (mm < A->lm1) {
        if (nn < A->ln1)
            offset = A->bsiz*(mm+A->lm1*nn) + m%A->mb + A->mb*(n%A->nb);
        else
            offset = A->A12 + (A->mb*(A->ln%A->nb)*mm) + m%A->mb + A->mb*(n%A->nb);
    }
    else {
        if (nn < A->ln1)
            offset = A->A21 + ((A->lm%A->mb)*A->nb*nn) + m%A->mb + (A->lm%A->mb)*(n%A->nb);
        else
            offset = A->A22 + m%A->mb  + (A->lm%A->mb)*(n%A->nb);
    }
    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

/***************************************************************************//**
 *  User routines
 **/
int PLASMA_Desc_Create(PLASMA_desc **desc, void *mat, PLASMA_enum dtyp, int mb, int nb, int bsiz, int lm, int ln, int i, int j, int m, int n);
int PLASMA_Desc_Destroy(PLASMA_desc **desc);

#ifdef __cplusplus
}
#endif

#endif
