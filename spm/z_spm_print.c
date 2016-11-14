/**
 *
 * @file z_spm_print.c
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

void
z_spmCSCPrint( FILE *f, const pastix_spm_t *spm )
{
    pastix_int_t i, j, baseval;
    pastix_int_t k, ii, jj, dofi, dofj, col, row;
    pastix_complex64_t *valptr;
    pastix_int_t *colptr, *rowptr, *dofs;

    assert( (spm->fmttype == PastixCSC) || (spm->fmttype == PastixCSR) );
    assert( spm->flttype == PastixComplex64 );

    baseval = spmFindBase( spm );
    i = 0; j = 0;

    if ( spm->fmttype == PastixCSC ) {
        colptr = spm->colptr;
        rowptr = spm->rowptr;
    }
    else {
        colptr = spm->rowptr;
        rowptr = spm->colptr;
    }

    valptr = (pastix_complex64_t*)(spm->values);
    dofs   = spm->dofs;

    switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PastixHermitian:
        for(j=0; j<spm->n; j++, colptr++)
        {
            dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
            col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j] - baseval;

            for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
            {
                i = (*rowptr - baseval);
                dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
                row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

                if ( spm->layout == PastixColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                            }
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                            }
                        }
                    }
                }
            }
        }
        break;
#endif
    case PastixSymmetric:
        for(j=0; j<spm->n; j++, colptr++)
        {
            dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
            col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j] - baseval;

            for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
            {
                i = (*rowptr - baseval);
                dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
                row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

                if ( spm->layout == PastixColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, *valptr );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, *valptr );
                            }
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, *valptr );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, *valptr );
                            }
                        }
                    }
                }
            }
        }
        break;
    case PastixGeneral:
    default:
        for(j=0; j<spm->n; j++, colptr++)
        {
            dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
            col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j] - baseval;

            for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
            {
                i = (*rowptr - baseval);
                dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
                row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

                if ( spm->layout == PastixColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                        }
                    }
                }
            }
        }
    }
    return;
}

void
z_spmCSRPrint( FILE *f, const pastix_spm_t *spm )
{
    pastix_int_t i, j, baseval;
    pastix_int_t k, ii, jj, dofi, dofj, col, row;
    pastix_complex64_t *valptr;
    pastix_int_t *colptr, *rowptr, *dofs;

    assert( spm->fmttype == PastixCSR );
    assert( spm->flttype == PastixComplex64 );

    baseval = spmFindBase( spm );
    i = 0; j = 0;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (pastix_complex64_t*)(spm->values);
    dofs   = spm->dofs;

    switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PastixHermitian:
        for(i=0; i<spm->n; i++, rowptr++)
        {
            dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
            row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

            for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
            {
                j = (*colptr - baseval);
                dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
                col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j] - baseval;

                if ( spm->layout == PastixColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                            }
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                            }
                        }
                    }
                }
            }
        }
        break;
#endif
    case PastixSymmetric:
        for(i=0; i<spm->n; i++, rowptr++)
        {
            dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
            row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

            for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
            {
                j = (*colptr - baseval);
                dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
                col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j] - baseval;

                if ( spm->layout == PastixColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, *valptr );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, *valptr );
                            }
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            if ( row == col ) {
                                if (row+ii >= col+jj) {
                                    z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                    if (row+ii > col+jj) {
                                        z_spmPrintElt( f, col + jj, row + ii, *valptr );
                                    }
                                }
                            }
                            else {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                z_spmPrintElt( f, col + jj, row + ii, *valptr );
                            }
                        }
                    }
                }
            }
        }
        break;
    case PastixGeneral:
    default:
        for(i=0; i<spm->n; i++, rowptr++)
        {
            dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
            row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i] - baseval;

            for(k=rowptr[0]; k<rowptr[1]; k++, colptr++)
            {
                j = (*colptr - baseval);
                dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
                col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j] - baseval;

                if ( spm->layout == PastixColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                        }
                    }
                }
            }
        }
    }
    return;
}

void
z_spmIJVPrint( FILE *f, const pastix_spm_t *spm )
{
    pastix_int_t i, j, baseval;
    pastix_int_t k, ii, jj, dofi, dofj, col, row;
    pastix_complex64_t *valptr;
    pastix_int_t *colptr, *rowptr, *dofs;

    assert( spm->fmttype == PastixIJV );
    assert( spm->flttype == PastixComplex64 );

    baseval = spmFindBase( spm );
    i = 0; j = 0;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (pastix_complex64_t*)(spm->values);
    dofs   = spm->dofs;

    switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PastixHermitian:
        for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
        {
            i = *rowptr - baseval;
            j = *colptr - baseval;

            if ( spm->dof > 0 ) {
                dofi = spm->dof;
                row  = spm->dof * i;
                dofj = spm->dof;
                col  = spm->dof * j;
            }
            else {
                dofi = dofs[i+1] - dofs[i];
                row  = dofs[i] - baseval;
                dofj = dofs[j+1] - dofs[j];
                col  = dofs[j] - baseval;
            }

            if ( spm->layout == PastixColMajor ) {
                for(jj=0; jj<dofj; jj++)
                {
                    for(ii=0; ii<dofi; ii++, valptr++)
                    {
                        if ( row == col ) {
                            if (row+ii >= col+jj) {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                if (row+ii > col+jj) {
                                    z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                                }
                            }
                        }
                        else {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                            z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                        }
                    }
                }
            }
            else {
                for(ii=0; ii<dofi; ii++)
                {
                    for(jj=0; jj<dofj; jj++, valptr++)
                    {
                        if ( row == col ) {
                            if (row+ii >= col+jj) {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                if (row+ii > col+jj) {
                                    z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                                }
                            }
                        }
                        else {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                            z_spmPrintElt( f, col + jj, row + ii, conj(*valptr) );
                        }
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

            if ( spm->dof > 0 ) {
                dofi = spm->dof;
                row  = spm->dof * i;
                dofj = spm->dof;
                col  = spm->dof * j;
            }
            else {
                dofi = dofs[i+1] - dofs[i];
                row  = dofs[i] - baseval;
                dofj = dofs[j+1] - dofs[j];
                col  = dofs[j] - baseval;
            }

            if ( spm->layout == PastixColMajor ) {
                for(jj=0; jj<dofj; jj++)
                {
                    for(ii=0; ii<dofi; ii++, valptr++)
                    {
                        if ( row == col ) {
                            if (row+ii >= col+jj) {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                if (row+ii > col+jj) {
                                    z_spmPrintElt( f, col + jj, row + ii, *valptr );
                                }
                            }
                        }
                        else {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                            z_spmPrintElt( f, col + jj, row + ii, *valptr );
                        }
                    }
                }
            }
            else {
                for(ii=0; ii<dofi; ii++)
                {
                    for(jj=0; jj<dofj; jj++, valptr++)
                    {
                        if ( row == col ) {
                            if (row+ii >= col+jj) {
                                z_spmPrintElt( f, row + ii, col + jj, *valptr );
                                if (row+ii > col+jj) {
                                    z_spmPrintElt( f, col + jj, row + ii, *valptr );
                                }
                            }
                        }
                        else {
                            z_spmPrintElt( f, row + ii, col + jj, *valptr );
                            z_spmPrintElt( f, col + jj, row + ii, *valptr );
                        }
                    }
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

            if ( spm->dof > 0 ) {
                dofi = spm->dof;
                row  = spm->dof * i;
                dofj = spm->dof;
                col  = spm->dof * j;
            }
            else {
                dofi = dofs[i+1] - dofs[i];
                row  = dofs[i] - baseval;
                dofj = dofs[j+1] - dofs[j];
                col  = dofs[j] - baseval;
            }

            if ( spm->layout == PastixColMajor ) {
                for(jj=0; jj<dofj; jj++)
                {
                    for(ii=0; ii<dofi; ii++, valptr++)
                    {
                        z_spmPrintElt( f, row + ii, col + jj, *valptr );
                    }
                }
            }
            else {
                for(ii=0; ii<dofi; ii++)
                {
                    for(jj=0; jj<dofj; jj++, valptr++)
                    {
                        z_spmPrintElt( f, row + ii, col + jj, *valptr );
                    }
                }
            }
        }
    }
    return;
}

void
z_spmPrint( FILE *f, const pastix_spm_t *spm )
{
    switch (spm->fmttype) {
    case PastixCSC:
        z_spmCSCPrint( f, spm );
        break;
    case PastixCSR:
        z_spmCSRPrint( f, spm );
        break;
    case PastixIJV:
        z_spmIJVPrint( f, spm );
    }
    return;
}
