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
z_spmPrint( const pastix_spm_t *spm )
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
                        fprintf( stderr, "%ld %ld (%e, %e)\n",
                                 row + jj, col + ii, creal(*valptr), cimag(*valptr) );
                        if (i != j) {
                            fprintf( stderr, "%ld %ld (%e, %e)\n",
                                     col + ii, row + jj, creal(conj(*valptr)), cimag(conj(*valptr)) );
                        }
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
#if defined(PRECISION_z) || defined(PRECISION_c)
                        fprintf( stderr, "%ld %ld (%e, %e)\n",
                                 row + jj, col + ii, creal(*valptr), cimag(*valptr) );
                        if (i != j) {
                            fprintf( stderr, "%ld %ld (%e, %e)\n",
                                     col + ii, row + jj, creal(*valptr), cimag(*valptr) );
                        }
#else
                        fprintf( stderr, "%ld %ld %e\n",
                                 row + jj, col + ii, *valptr );
                        if (i != j) {
                            fprintf( stderr, "%ld %ld %e\n",
                                     col + ii, row + jj, *valptr );
                        }
#endif
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
#if defined(PRECISION_z) || defined(PRECISION_c)
                        fprintf( stderr, "%ld %ld (%e, %e)\n",
                                 row + jj, col + ii, creal(*valptr), cimag(*valptr) );
#else
                        fprintf( stderr, "%ld %ld %e\n",
                                 row + jj, col + ii, *valptr );
#endif
                    }
                }
            }
        }
    }
    return;
}

