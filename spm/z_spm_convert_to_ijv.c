/**
 *
 * @file z_spm_convert_to_ijv.c
 *
 * SParse Matrix package conversion routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
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
 * @ingroup spm_dev_convert
 *
 * @brief convert a matrix in CSC format to a matrix in IJV format.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The csc matrix at enter,
 *          the ijv matrix at exit.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
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

    memFree_null(spm->colptr);
    spm->colptr = col_ijv;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_convert
 *
 * @brief convert a matrix in CSR format to a matrix in IJV format.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The csr matrix at enter,
 *          the ijv matrix at exit.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
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

    memFree_null(spm->rowptr);
    spm->rowptr = row_ijv;

    return PASTIX_SUCCESS;
}
