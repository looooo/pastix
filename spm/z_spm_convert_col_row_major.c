/**
 *
 * @file z_spm_convert_col_row_major.c
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
 * @precisions normal z -> c d s p
 **/
#include "common.h"
#include "spm.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * z_spmConvertColMaj2RowMaj - convert a matrix in Column Major format to a
 * matrix in Row Major format.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The colmaj matrix at enter,
 *          the rowmaj matrix at exit.
 *
 *******************************************************************************/
void
z_spmConvertColMaj2RowMaj(pastix_spm_t *spm)
{
    assert(spm->colmajor > 0);
    assert(spm->dof != 1);
    assert(spm->dofs);

    spm->colmajor = -1;

    pastix_int_t dofi, dofj, ii, jj, i, j, k, baseval;
    pastix_int_t cpt=0;
    pastix_int_t *dofs = spm->dofs;
    pastix_int_t *tmp;
    pastix_complex64_t *oavals, *navals;

    baseval = spmFindBase( spm );

    oavals = (pastix_complex64_t*)spm->values;
    navals = malloc(sizeof(pastix_complex64_t)*spm->nnzexp);

    switch(spm->fmttype)
    {
    case PastixCSC:
        for(i=0; i<spm->n; i++)
        {
            //dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            dofi = dofs[i+1] - dofs[i];
            for(k=spm->colptr[i]; k<spm->colptr[i+1]; k++)
            {
                j = spm->rowptr[k-baseval] - baseval;
                //dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
                dofj = dofs[j+1] - dofs[j];
                for(ii=0; ii<dofi; ii++)
                {
                    for(jj=0; jj<dofj; jj++)
                    {
                        navals[cpt + jj * dofi + ii] = *oavals;
                        oavals++;
                    }
                }
                cpt += dofi * dofj;
            }
        }
        spm->values = navals;
        break;

    case PastixCSR:
        tmp           = spm->rowptr;
        spm->rowptr   = spm->colptr;
        spm->colptr   = tmp;
        spm->fmttype  = PastixCSC;

        z_spmConvertRowMaj2ColMaj(spm);

        spm->colmajor = -1;
        tmp          = spm->rowptr;
        spm->rowptr  = spm->colptr;
        spm->colptr  = tmp;
        spm->fmttype = PastixCSR;
        break;

    case PastixIJV:
        for(k=0; k<spm->nnz; k++)
        {
            j = spm->rowptr[k]-baseval;
            i = spm->colptr[k]-baseval;
            //dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            //dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
            dofi = dofs[i+1] - dofs[i];
            dofj = dofs[j+1] - dofs[j];
            for(ii=0; ii<dofi; ii++)
            {
                for(jj=0; jj<dofj; jj++)
                {
                    navals[cpt + jj * dofi + ii] = *oavals;
                    oavals++;
                }
            }
            cpt += dofi * dofj;
        }
        spm->values = navals;
    default:
        break;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * z_spmConvertColMaj2RowMaj - convert a matrix in Row Major format to a matrix
 * in Column Major format.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The rowmaj matrix at enter,
 *          the colmaj matrix at exit.
 *
 *******************************************************************************/
void
z_spmConvertRowMaj2ColMaj(pastix_spm_t *spm)
{
    assert(spm->colmajor < 0);
    assert(spm->dof != 1);
    assert(spm->dofs);

    spm->colmajor = 1;
    pastix_int_t dofi, dofj, ii, jj, i, j, k, baseval;
    pastix_int_t cpt=0;
    pastix_int_t *dofs = spm->dofs;
    pastix_int_t *tmp;
    pastix_complex64_t *oavals, *navals;

    baseval = spmFindBase( spm );

    oavals = (pastix_complex64_t*)spm->values;
    navals = malloc(sizeof(pastix_complex64_t)*spm->nnzexp);

    switch(spm->fmttype)
    {
    case PastixCSC :
        for(i=0; i<spm->n; i++)
        {
            //dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            dofi = dofs[i+1] - dofs[i];
            for(k=spm->colptr[i]; k<spm->colptr[i+1]; k++)
            {
                j = spm->rowptr[k-baseval]-baseval;
                //dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
                dofj = dofs[j+1] - dofs[j];
                for(jj=0; jj<dofj; jj++)
                {
                    for(ii=0; ii<dofi; ii++)
                    {
                        navals[cpt + ii * dofj + jj] = *oavals;
                        oavals++;
                    }
                }
                cpt += dofi * dofj;
            }
        }
        spm->values = navals;
        break;

    case PastixCSR :
        tmp          = spm->rowptr;
        spm->rowptr  = spm->colptr;
        spm->colptr  = tmp;
        spm->fmttype = PastixCSC;

        z_spmConvertColMaj2RowMaj(spm);

        spm->colmajor = 1;
        tmp          = spm->rowptr;
        spm->rowptr  = spm->colptr;
        spm->colptr  = tmp;
        spm->fmttype = PastixCSR;
        break;

    case PastixIJV:
        for(k=0; k<spm->nnz; k++)
        {
            j = spm->rowptr[k]-baseval;
            i = spm->colptr[k]-baseval;
            //dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            //dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
            dofi = dofs[i+1] - dofs[i];
            dofj = dofs[j+1] - dofs[j];
                for(jj=0; jj<dofj; jj++)
            {
                for(ii=0; ii<dofi; ii++)
                {
                    navals[cpt + ii * dofj + jj] = *oavals;
                    oavals++;
                }
            }
            cpt += dofi * dofj;
        }
        spm->values = navals;
    default:
        break;
    }
}
