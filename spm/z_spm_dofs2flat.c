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


pastix_spm_t *
dofs2flat(pastix_spm_t *spm)
{
    pastix_spm_t* new_spm=malloc(sizeof(pastix_spm_t));
    spmInit(new_spm);

    new_spm->mtxtype = spm->mtxtype;
    new_spm->flttype = spm->flttype;
    new_spm->fmttype = spm->fmttype;

    new_spm->gN   = spm->gNexp;
    new_spm->n    = spm->nexp;
    new_spm->gnnz = spm->gnnzexp;
    new_spm->nnz  = spm->nnzexp;

    new_spm->gNexp   = spm->gNexp;
    new_spm->nexp    = spm->nexp;
    new_spm->gnnzexp = spm->gnnzexp;
    new_spm->nnzexp  = spm->nnzexp;

    new_spm->dof       = 1;
    new_spm->dofs      = NULL;
    new_spm->colmajor  = 1;

    pastix_int_t i, j, k, ii, jj, dofi, dofj, col, row, baseval;
    baseval = spmFindBase( spm );

    pastix_int_t       *new_col  = calloc(spm->nexp+1,sizeof(pastix_int_t));
    pastix_int_t       *new_row  = malloc(sizeof(pastix_int_t)*spm->nnzexp);
    pastix_int_t       *dofs     = spm->dofs;

    pastix_complex64_t *new_vals = malloc(sizeof(pastix_int_t)*spm->nnzexp);
    pastix_complex64_t *vals     = (pastix_complex64_t*)spm->values;

    for(i=0;i<spm->n ; i++)
    {
        col = ( spm->dof > 0 ) ? i : dofs[i];
        for(k=spm->colptr[i]-baseval; k<spm->colptr[i+1]-baseval; k++)
	{
            j = spm->rowptr[k]-baseval;
            row  = ( spm->dof > 0 ) ? j        : dofs[j];
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
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

    int cpt = 0;
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
                    cpt++;
                    new_col[col+ii]++;
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

    new_spm->colptr   = new_col;
    new_spm->rowptr   = new_row;
    new_spm->loc2glob = NULL; // ?
    new_spm->values   = new_vals;

    return new_spm;
}
