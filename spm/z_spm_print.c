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
z_spmPrint( FILE *f, const pastix_spm_t *spm )
{
    pastix_int_t i, j, baseval;
    pastix_int_t k, ii, jj, dofi, dofj, col, row;
    pastix_complex64_t *valptr;
    pastix_int_t *colptr, *rowptr, *dofs;

    assert( spm->fmttype == PastixCSC );
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
        for(j=0; j<spm->n; j++, colptr++)
        {
            dofj = ( spm->dof > 0 ) ? spm->dof     : dofs[j+1] - dofs[j];
            col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j];

            for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
            {
                i = (*rowptr - baseval);
                dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
                row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i];

                if ( spm->layout == PastixColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            if (row+ii >= col+jj) {
                                fprintf( f, "%ld %ld (%e, %e)\n",
                                         row + ii, col + jj, creal(*valptr), cimag(*valptr) );
                                if (row+ii > col+jj) {
                                    fprintf( f, "%ld %ld (%e, %e)\n",
                                             col + jj, row + ii, creal(conj(*valptr)), cimag(conj(*valptr)) );
                                }
                            }
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            if (row+ii >= col+jj) {
                                fprintf( f, "%ld %ld (%e, %e)\n",
                                         row + ii, col + jj, creal(*valptr), cimag(*valptr) );
                                if (row+ii > col+jj) {
                                    fprintf( f, "%ld %ld (%e, %e)\n",
                                             col + jj, row + ii, creal(conj(*valptr)), cimag(conj(*valptr)) );
                                }
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
            col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j];

            for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
            {
                i = (*rowptr - baseval);
                dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
                row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i];

                if ( spm->layout == PastixColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
                            if (row+ii >= col+jj) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                                fprintf( f, "%ld %ld (%e, %e)\n",
                                         row + ii, col + jj, creal(*valptr), cimag(*valptr) );
                                if (row+ii > col+jj) {
                                    fprintf( f, "%ld %ld (%e, %e)\n",
                                             col + jj, row + ii, creal(*valptr), cimag(*valptr) );
                                }
#else
                                fprintf( f, "%ld %ld %e\n",
                                         row + ii, col + jj, *valptr );
                                if (row+ii > col+jj) {
                                    fprintf( f, "%ld %ld %e\n",
                                             col + jj, row + ii, *valptr );
                                }
#endif
                            }
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
                            if (row+ii >= col+jj) {
#if defined(PRECISION_z) || defined(PRECISION_c)
                                fprintf( f, "%ld %ld (%e, %e)\n",
                                         row + ii, col + jj, creal(*valptr), cimag(*valptr) );
                                if (row+ii > col+jj) {
                                    fprintf( f, "%ld %ld (%e, %e)\n",
                                             col + jj, row + ii, creal(*valptr), cimag(*valptr) );
                                }
#else
                                fprintf( f, "%ld %ld %e\n",
                                         row + ii, col + jj, *valptr );
                                if (row+ii > col+jj) {
                                    fprintf( f, "%ld %ld %e\n",
                                             col + jj, row + ii, *valptr );
                                }
#endif
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
            col  = ( spm->dof > 0 ) ? spm->dof * j : dofs[j];

            for(k=colptr[0]; k<colptr[1]; k++, rowptr++)
            {
                i = (*rowptr - baseval);
                dofi = ( spm->dof > 0 ) ? spm->dof     : dofs[i+1] - dofs[i];
                row  = ( spm->dof > 0 ) ? spm->dof * i : dofs[i];

                if ( spm->layout == PastixColMajor ) {
                    for(jj=0; jj<dofj; jj++)
                    {
                        for(ii=0; ii<dofi; ii++, valptr++)
                        {
#if defined(PRECISION_z) || defined(PRECISION_c)
                            fprintf( f, "%ld %ld (%e, %e)\n",
                                     row + ii, col + jj, creal(*valptr), cimag(*valptr) );
#else
                            fprintf( f, "%ld %ld %e\n",
                                     row + ii, col + jj, *valptr );
#endif
                        }
                    }
                }
                else {
                    for(ii=0; ii<dofi; ii++)
                    {
                        for(jj=0; jj<dofj; jj++, valptr++)
                        {
#if defined(PRECISION_z) || defined(PRECISION_c)
                            fprintf( f, "%ld %ld (%e, %e)\n",
                                     row + ii, col + jj, creal(*valptr), cimag(*valptr) );
#else
                            fprintf( f, "%ld %ld %e\n",
                                     row + ii, col + jj, *valptr );
#endif
                        }
                    }
                }
            }
        }
    }
    return;
}

