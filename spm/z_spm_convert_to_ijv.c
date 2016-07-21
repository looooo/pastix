/**
 *
 * @file z_spm_convert_to_ijv.c
 *
 *  PaStiX spm routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
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
 * z_spmConvertCSC2IJV - convert a matrix in CSC format to a matrix in IJV
 * format.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The csc matrix at enter,
 *          the ijv matrix at exit.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertCSC2IJV( pastix_spm_t *spm )
{
    pastix_int_t *col_ijv, *colptr;
    pastix_int_t i, j, baseval, nnz;

    /*
     * Check the baseval
     */
    baseval = spmFindBase( spm );
    nnz = spm->nnz;
    spm->fmttype = PastixIJV;

    col_ijv = malloc(nnz*sizeof(pastix_int_t));
    assert( col_ijv );

    colptr = col_ijv;
    for(i=0; i<spm->n; i++)
    {
        for(j=spm->colptr[i]; j<spm->colptr[i+1]; j++)
        {
            *colptr = i+baseval; colptr++;
        }
    }

    /* Transpose values in row major format */
    if( !spm->colmajor ) //A test
    {
        int k,ii,jj,dofi,dofj;
        int cpt=0;
        pastix_complex64_t* oavals = (pastix_complex64_t*)spm->values;
        pastix_complex64_t* navals = malloc(sizeof(pastix_complex64_t)*spm->nnzexp);
        pastix_int_t* dofs=spm->dofs;
        for(k=0; k<spm->nnz; k++)
        {
            j = spm->rowptr[k]-baseval;
            i = col_ijv[k]-baseval;
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
            for(ii=0; ii<dofi; ii++)
            {
                for(jj=0; jj<dofj; jj++)
                {
                    navals[cpt+jj*dofi+ii]=*oavals;
                    oavals++;
                }
            }
            cpt+=dofi*dofj;
        }
        spm->values=navals;
    }

    memFree_null(spm->colptr);
    spm->colptr = col_ijv;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * z_spmConvertCSR2IJV - convert a matrix in CSR format to a matrix in IJV
 * format.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The csr matrix at enter,
 *          the ijv matrix at exit.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertCSR2IJV( pastix_spm_t *spm )
{
    pastix_int_t *row_ijv, *rowptr;
    pastix_int_t i, j, baseval, nnz;

    /*
     * Check the baseval
     */
    baseval = spmFindBase( spm );
    nnz = spm->nnz;
    spm->fmttype = PastixIJV;

    row_ijv = malloc(nnz*sizeof(pastix_int_t));
    assert( row_ijv );

    rowptr = row_ijv;
    for(i=0; i<spm->n; i++)
    {
        for(j=spm->rowptr[i]; j<spm->rowptr[i+1]; j++)
        {
            *rowptr = i+baseval; rowptr++;
        }
    }

    /* Transpose values in column major format */
    if( spm->colmajor )
    {
        int k,ii,jj,dofi,dofj;
        int cpt=0;
        pastix_complex64_t* oavals = (pastix_complex64_t*)spm->values;
        pastix_complex64_t* navals = malloc(sizeof(pastix_complex64_t)*spm->nnzexp);
        pastix_int_t* dofs=spm->dofs;
        for(k=0; k<spm->nnz; k++)
        {
            i = spm->colptr[k]-baseval;
            j = row_ijv[k]-baseval;
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
            for(jj=0; jj<dofj; jj++)
            {
                for(ii=0; ii<dofi; ii++)
                {
                    navals[cpt+ii*dofj+jj]=*oavals;
                    oavals++;
                }
            }
            cpt+=dofi*dofj;
        }
        spm->values=navals;
    }

    memFree_null(spm->rowptr);
    spm->rowptr = row_ijv;

    return PASTIX_SUCCESS;
}
