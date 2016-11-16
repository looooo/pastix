/**
 *
 * @file z_spm_dofs2flat.c
 *
 *  PaStiX spm routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
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

/**
 *******************************************************************************
 *
 * @ingroup spm_internal
 *
 * z_spmDofs2Flat - Convert a sparse matrix with dofs into a sparse matrix without dofs.
 *
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The sparse matrix which needs to be converted
 *
 *******************************************************************************
 **/
void
z_spmDofs2Flat(pastix_spm_t *spm)
{
    pastix_int_t i, j, k, ii, jj, dofi, dofj, col, row, baseval, cpt;
    baseval = spmFindBase( spm );

    pastix_int_t       *new_col  = calloc(spm->nexp+1,sizeof(pastix_int_t));
    pastix_int_t       *new_row;
    pastix_int_t       *dofs     = spm->dofs;

    pastix_complex64_t *new_vals;
    pastix_complex64_t *vals     = (pastix_complex64_t*)spm->values;

    switch(spm->mtxtype)
    {
    case PastixGeneral:
        new_row  = malloc(sizeof(pastix_int_t)*spm->nnzexp);
        new_vals = malloc(sizeof(pastix_complex64_t)*spm->nnzexp);
        for(i=0;i<spm->n ; i++)
        {
            col  = ( spm->dof > 0 ) ? i        : dofs[i];
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            for(k=spm->colptr[i]-baseval; k<spm->colptr[i+1]-baseval; k++)
            {
                j = spm->rowptr[k]-baseval;
                row  = ( spm->dof > 0 ) ? j        : dofs[j];
                dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
                for(ii=0; ii<dofi; ii++)
                {
                    new_col[col+ii+1] +=  dofj;
                }
            }
        }

        for(i=0; i<spm->nexp; i++)
        {
            new_col[i+1]+=new_col[i];
        }

        cpt = 0;
        for(i=0; i < spm->n;i++)
        {
            col  = ( spm->dof > 0 ) ? i        : dofs[i];
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            for(k=spm->colptr[i]-baseval ; k<spm->colptr[i+1]-baseval ;k++)
            {
                j = spm->rowptr[k] - baseval;
                row  = ( spm->dof > 0 ) ? j        : dofs[j];
                dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
                for(ii=0;ii < dofi; ii++)
                {
                    for(jj=0;jj < dofj ; jj++)
                    {
                        new_vals[new_col[col+ii]] = vals[cpt];
                        new_row[new_col[col+ii]]  = row + jj + baseval;
                        new_col[col+ii]++;
                        cpt++;
                    }
                }
            }
        }

        {
            int tmp;
            int tmp1 = 0;
            for(i=0; i<spm->nexp; i++)
            {
                tmp = new_col[i];
                new_col[i] = tmp1+baseval;
                tmp1 = tmp;
            }
            new_col[i] += baseval;
        }
        spm->gN   = spm->gNexp;
        spm->n    = spm->nexp;
        spm->gnnz = spm->gnnzexp;
        spm->nnz  = spm->nnzexp;

        spm->dof      = 1;
        spm->dofs     = NULL;
        spm->layout   = PastixColMajor;

        spm->colptr   = new_col;
        spm->rowptr   = new_row;
        //spm->loc2glob = NULL; // ?
        spm->values   = new_vals;
        break;

    case PastixSymmetric:
        for(i=0;i<spm->n ; i++)
        {
            col  = ( spm->dof > 0 ) ? i        : dofs[i];
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            for(k=spm->colptr[i]-baseval; k<spm->colptr[i+1]-baseval; k++)
            {
                j = spm->rowptr[k]-baseval;
                row  = ( spm->dof > 0 ) ? j        : dofs[j];
                dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
                for(ii=0; ii<dofi; ii++)
                {
                    for(jj=0; jj<dofj; jj++)
                    {
                        if( i != j )
                            new_col[col+ii+1] +=  1;
                        else
                            if(ii <= jj )
                                new_col[col+ii+1] += 1;
                    }
                }
            }
        }
        for(i=0; i<spm->nexp; i++)
        {
            new_col[i+1] += new_col[i];
        }
        pastix_int_t nnz = new_col[spm->nexp];
        new_row  = malloc(sizeof(pastix_int_t)*nnz);
        new_vals = malloc(sizeof(pastix_complex64_t)*nnz);

        cpt = 0;
        for(i=0; i < spm->n;i++)
        {
            col  = ( spm->dof > 0 ) ? i        : dofs[i];
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            for(k=spm->colptr[i]-baseval ; k<spm->colptr[i+1]-baseval ;k++)
            {
                j = spm->rowptr[k] - baseval;
                row  = ( spm->dof > 0 ) ? j        : dofs[j];
                dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
                for(ii=0;ii < dofi; ii++)
                {
                    for(jj=0;jj < dofj ; jj++)
                    {
                        if( i == j )
                        {
                            if ( ii <= jj )
                            {
                                /* diagonal dominant for spd matrix
                                if( ii == jj)
                                    new_vals[new_col[col+ii]] = 2*vals[cpt];
                                 else
                                */
                                new_vals[new_col[col+ii]] = vals[cpt];
                                new_row[new_col[col+ii]]  = row + jj + baseval;
                                new_col[col+ii]++;
                            }
                        }
                        else
                        {
                            new_vals[new_col[col+ii]] = vals[cpt];
                            new_row[new_col[col+ii]]  = row + jj + baseval;
                            new_col[col+ii]++;

                        }
                        cpt++;
                    }
                }
            }
        }
        {
            int tmp;
            int tmp1 = 0;
            for(i=0; i<spm->nexp; i++)
            {
                tmp = new_col[i];
                new_col[i] = tmp1+baseval;
                tmp1 = tmp;
            }
            new_col[i] += baseval;
        }
        spm->gN   = spm->gNexp;
        spm->n    = spm->nexp;
        spm->gnnz = nnz;
        spm->nnz  = nnz;
        spm->gnnzexp = nnz;
        spm->nnzexp  = nnz;

        spm->dof      = 1;
        spm->dofs     = NULL;
        spm->layout   = PastixColMajor;

        spm->colptr   = new_col;
        spm->rowptr   = new_row;
        //spm->loc2glob = NULL; //
        spm->values   = new_vals;

        break;
    }
    assert(spm->loc2glob == NULL);//to do
}
