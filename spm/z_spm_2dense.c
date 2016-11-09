/**
 *
 * @file z_spm_2dense.c
 *
 * Convert a sparse matrix into a dense matrix.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Alban Bellot
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 *
 **/
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "pastix.h"
#include "common.h"
#include "spm.h"
#include "z_spm.h"

pastix_complex64_t *
z_spmCSC2dense( const pastix_spm_t *spm )
{
    pastix_int_t i, j, lda, baseval;
    pastix_complex64_t *A, *valptr;
    pastix_int_t *colptr, *rowptr;

    assert( spm->fmttype == PastixCSC );
    assert( spm->flttype == PastixComplex64 );

    lda = spm->gNexp;
    A = (pastix_complex64_t*)malloc(lda * lda * sizeof(pastix_complex64_t));
    memset( A, 0, lda * lda * sizeof(pastix_complex64_t));

    baseval = spmFindBase( spm );
    i = 0; j = 0;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (pastix_complex64_t*)(spm->values);

    /**
     * Constant degree of fredom of 1
     */
    if ( spm->dof == 1 ) {
        switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
        case PastixHermitian:
            for(i=0; i<spm->n; i++, colptr++)
            {
                for(j=colptr[0]; j<colptr[1]; j++, rowptr++, valptr++)
                {
                    if( i == (*rowptr-baseval) ) {
                        /* Make sure the matrix is hermitian */
                        A[ i * lda + (*rowptr - baseval) ] = creal(*valptr) + I * 0.;
                    }
                    else {
                        A[ i * lda + (*rowptr - baseval) ] = *valptr;
                        A[ (*rowptr - baseval) * lda + i ] = conj(*valptr);
                    }
                }
            }
            break;
#endif
        case PastixSymmetric:
            for(i=0; i<spm->n; i++, colptr++)
            {
                for(j=colptr[0]; j<colptr[1]; j++, rowptr++, valptr++)
                {
                    A[ i * lda + (*rowptr - baseval) ] = *valptr;
                    A[ (*rowptr - baseval) * lda + i ] = *valptr;
                }
            }
            break;
        case PastixGeneral:
        default:
            for(i=0; i<spm->n; i++, colptr++)
            {
                for(j=colptr[0]; j<colptr[1]; j++, rowptr++, valptr++)
                {
                    A[ i * lda + (*rowptr - baseval) ] = *valptr;
                }
            }
        }
    }
    /**
     * General degree of freedom (constant or variable)
     */
    else {
        pastix_int_t  k, ii, jj, dofi, dofj, col, row;
        pastix_int_t *dofs = spm->dofs;

        switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
        case PastixHermitian:
            for(i=0; i<spm->n; i++, colptr++)
            {
                dofi = ( spm->dof > 1 ) ? spm->dof : dofs[i+1] - dofs[i];
                col = dofs[i];

                for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
                {
                    j = (*rowptr - baseval);
                    dofj = ( spm->dof > 1 ) ? spm->dof : dofs[j+1] - dofs[j];
                    row = dofs[j];

                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            A[ (col + ii) * lda + (row + jj) ] = *valptr;
                            A[ (row + jj) * lda + (col + ii) ] = conj(*valptr);
                        }
                    }
                }
            }
            break;
#endif
        case PastixSymmetric:
            for(i=0; i<spm->n; i++, colptr++)
            {
                dofi = ( spm->dof > 1 ) ? spm->dof : dofs[i+1] - dofs[i];
                col = dofs[i];

                for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
                {
                    j = (*rowptr - baseval);
                    dofj = ( spm->dof > 1 ) ? spm->dof : dofs[j+1] - dofs[j];
                    row = dofs[j];

                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            A[ (col + ii) * lda + (row + jj) ] = *valptr;
                            A[ (row + jj) * lda + (col + ii) ] = *valptr;
                        }
                    }
                }
            }
            break;
        case PastixGeneral:
        default:
            for(i=0; i<spm->n; i++, colptr++)
            {
                dofi = ( spm->dof > 1 ) ? spm->dof : dofs[i+1] - dofs[i];
                col = dofs[i];

                for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
                {
                    j = (*rowptr - baseval);
                    dofj = ( spm->dof > 1 ) ? spm->dof : dofs[j+1] - dofs[j];
                    row = dofs[j];

                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            A[ (col + ii) * lda + (row + jj) ] = *valptr;
                        }
                    }
                }
            }
        }
    }
    return A;
}

pastix_complex64_t *
z_spmCSR2dense( const pastix_spm_t *spm )
{
    pastix_int_t i, j, lda, baseval;
    pastix_complex64_t *A, *valptr;
    pastix_int_t *colptr, *rowptr;

    assert( spm->fmttype == PastixCSR );
    assert( spm->flttype == PastixComplex64 );

    lda = spm->gNexp;
    A = (pastix_complex64_t*)malloc(lda * lda * sizeof(pastix_complex64_t));
    memset( A, 0, lda * lda * sizeof(pastix_complex64_t));

    baseval = spmFindBase( spm );
    i = 0; j = 0;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (pastix_complex64_t*)(spm->values);

    /**
     * Constant degree of fredom of 1
     */
    if ( spm->dof == 1 ) {
        switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
        case PastixHermitian:
            for(i=0; i<spm->n; i++, rowptr++)
            {
                for(j=rowptr[0]; j<rowptr[1]; j++, colptr++, valptr++)
                {
                    if( i == (*colptr-baseval) ) {
                        /* Make sure the matrix is hermitian */
                        A[ i * lda + i ] = creal(*valptr) + I * 0.;
                    }
                    else {
                        A[ i * lda + (*colptr - baseval) ] = conj(*valptr);
                        A[ (*colptr - baseval) * lda + i ] = *valptr;
                    }
                }
            }
            break;
#endif
        case PastixSymmetric:
            for(i=0; i<spm->n; i++, rowptr++)
            {
                for(j=rowptr[0]; j<rowptr[1]; j++, colptr++, valptr++)
                {
                    A[ i * lda + (*colptr - baseval) ] = *valptr;
                    A[ (*colptr - baseval) * lda + i ] = *valptr;
                }
            }
            break;
        case PastixGeneral:
        default:
            for(i=0; i<spm->n; i++, rowptr++)
            {
                for(j=rowptr[0]; j<rowptr[1]; j++, colptr++, valptr++)
                {
                    A[ i * lda + (*colptr - baseval) ] = *valptr;
                }
            }
        }
    }
    /**
     * General degree of freedom (constant or variable)
     */
    else {
        pastix_int_t  k, ii, jj, dofi, dofj, col, row;
        pastix_int_t *dofs = spm->dofs;

        switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
        case PastixHermitian:
            for(i=0; i<spm->n; i++, rowptr++)
            {
                dofi = ( spm->dof > 1 ) ? spm->dof : dofs[i+1] - dofs[i];
                row = dofs[i];

                for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
                {
                    j = (*colptr - baseval);
                    dofj = ( spm->dof > 1 ) ? spm->dof : dofs[j+1] - dofs[j];
                    col = dofs[j];

                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            A[ (col + ii) * lda + (row + jj) ] = *valptr;
                            A[ (row + jj) * lda + (col + ii) ] = conj(*valptr);
                        }
                    }
                }
            }
            break;
#endif
        case PastixSymmetric:
            for(i=0; i<spm->n; i++, rowptr++)
            {
                dofi = ( spm->dof > 1 ) ? spm->dof : dofs[i+1] - dofs[i];
                row = dofs[i];

                for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
                {
                    j = (*colptr - baseval);
                    dofj = ( spm->dof > 1 ) ? spm->dof : dofs[j+1] - dofs[j];
                    col = dofs[j];

                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            A[ (col + ii) * lda + (row + jj) ] = *valptr;
                            A[ (row + jj) * lda + (col + ii) ] = *valptr;
                        }
                    }
                }
            }
            break;
        case PastixGeneral:
        default:
            for(i=0; i<spm->n; i++, rowptr++)
            {
                dofi = ( spm->dof > 1 ) ? spm->dof : dofs[i+1] - dofs[i];
                row = dofs[i];

                for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
                {
                    j = (*colptr - baseval);
                    dofj = ( spm->dof > 1 ) ? spm->dof : dofs[j+1] - dofs[j];
                    col = dofs[j];

                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            A[ (col + ii) * lda + (row + jj) ] = *valptr;
                        }
                    }
                }
            }
        }
    }
    return A;
}

pastix_complex64_t *
z_spmIJV2dense( const pastix_spm_t *spm )
{
    pastix_int_t i, j, k, lda, baseval;
    pastix_complex64_t *A, *valptr;
    pastix_int_t *colptr, *rowptr;

    assert( spm->fmttype == PastixIJV );
    assert( spm->flttype == PastixComplex64 );

    lda = spm->gNexp;
    A = (pastix_complex64_t*)malloc(lda * lda * sizeof(pastix_complex64_t));
    memset( A, 0, lda * lda * sizeof(pastix_complex64_t));

    baseval = spmFindBase( spm );
    i = 0; j = 0;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (pastix_complex64_t*)(spm->values);

    /**
     * Constant degree of fredom of 1
     */
    if ( spm->dof == 1 ) {
        switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
        case PastixHermitian:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++, valptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;

                if( i == j ) {
                    /* Make sure the matrix is hermitian */
                    A[ i * lda + i ] = creal(*valptr) + I * 0.;
                }
                else {
                    A[ j * lda + i ] = *valptr;
                    A[ i * lda + j ] = conj(*valptr);
                }
            }
            break;
#endif
        case PastixSymmetric:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++, valptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;

                A[ j * lda + i ] = *valptr;
                A[ i * lda + j ] = *valptr;
            }
            break;
        case PastixGeneral:
        default:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++, valptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;

                A[ j * lda + i ] = *valptr;
            }
        }
    }
    /**
     * General degree of freedom (constant or variable)
     */
    else {
        pastix_int_t  k, ii, jj, dofi, dofj, col, row;
        pastix_int_t *dofs = spm->dofs;

        switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
        case PastixHermitian:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;

                dofi = ( spm->dof > 1 ) ? spm->dof : dofs[i+1] - dofs[i];
                row = dofs[i];

                dofj = ( spm->dof > 1 ) ? spm->dof : dofs[j+1] - dofs[j];
                col = dofs[j];

                for(jj=0; jj<dofj; jj++)
                {
                    for(ii=0; ii<dofi; ii++, valptr++)
                    {
                        if( col+jj == row+ii ) {
                            /* Make sure the matrix is hermitian */
                            A[ (col+jj) * lda + (row+ii) ] = creal(*valptr) + I * 0.;
                        }
                        else {
                            A[ (col+jj) * lda + (row+ii) ] = *valptr;
                            A[ (row+ii) * lda + (col+jj) ] = conj(*valptr);
                        }
                    }
                }
            }
            break;
#endif
        case PastixSymmetric:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;

                dofi = ( spm->dof > 1 ) ? spm->dof : dofs[i+1] - dofs[i];
                row = dofs[i];

                dofj = ( spm->dof > 1 ) ? spm->dof : dofs[j+1] - dofs[j];
                col = dofs[j];

                for(jj=0; jj<dofj; jj++)
                {
                    for(ii=0; ii<dofi; ii++, valptr++)
                    {
                        A[ (col+jj) * lda + (row+ii) ] = *valptr;
                        A[ (row+ii) * lda + (col+jj) ] = *valptr;
                    }
                }
            }
            break;

        case PastixGeneral:
        default:
            for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;

                dofi = ( spm->dof > 1 ) ? spm->dof : dofs[i+1] - dofs[i];
                row = dofs[i];

                dofj = ( spm->dof > 1 ) ? spm->dof : dofs[j+1] - dofs[j];
                col = dofs[j];

                for(jj=0; jj<dofj; jj++)
                {
                    for(ii=0; ii<dofi; ii++, valptr++)
                    {
                        A[ (col+jj) * lda + (row+ii) ] = *valptr;
                    }
                }
            }
        }
    }
    return A;
}

pastix_complex64_t *
z_spm2dense( const pastix_spm_t *spm )
{
    switch (spm->fmttype) {
    case PastixCSC:
        return z_spmCSC2dense( spm );
    case PastixCSR:
        return z_spmCSR2dense( spm );
    case PastixIJV:
        return z_spmIJV2dense( spm );
    }
    return NULL;
}

void
z_spmDensePrint( pastix_int_t m, pastix_int_t n, pastix_complex64_t *A, pastix_int_t lda )
{
    pastix_int_t i, j;

    for(i=0; i<m; i++)
    {
        for(j=0; j<n; j++)
        {
            if ( cabs( A[ lda *j + i] ) != 0. ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                fprintf( stderr, "%ld %ld (%e, %e)\n",
                         i, j, creal(A[lda * j + i]), cimag(A[lda * j + i]) );
#else
                fprintf( stderr, "%ld %ld %e\n",
                         i, j, A[lda * j + i] );
#endif
            }
        }
    }
    return;
}

