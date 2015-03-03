/**
 *
 * @file z_spm_convert_to_csr.c
 *
 *  PaStiX csc routines
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
#include "csc.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * z_spmConvertCSC2CSR - convert a matrix in CSC format to a matrix in CSR
 * format.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The csc matrix at enter,
 *          the csr matrix at exit.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertCSC2CSR( pastix_csc_t *spm )
{
    pastix_int_t *tmp;
    pastix_int_t  result;

    switch( spm->mtxtype ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PastixHermitian:
    {
        /* Similar to PastixSymmetric case with conjugate of the values */
        pastix_complex64_t *valptr = spm->avals;
        pastix_int_t i;

        for(i=0; i<spm->nnz; i++, valptr++){
            *valptr = conj( *valptr );
        }
    }
#endif
    case PastixSymmetric:
    {
        pastix_int_t *tmp;

        /* Just need to swap the pointers */
        tmp         = spm->rows;
        spm->rows   = spm->colptr;
        spm->colptr = tmp;

        return PASTIX_SUCCESS;
    }
    break;

    case PastixGeneral:
    default:
    {

        /* transpose spm in CSC to trans(spm) in CSR */
        tmp         = spm->rows;
        spm->rows   = spm->colptr;
        spm->colptr = tmp;
        spm->fmttype = PastixCSR;
        
        /* convert trans(spm) in CSR to trans(spm) in CSC */
        result = z_spmConvertCSR2CSC( spm );

        /* transpose trans(spm) in CSC to obtain spm in CSR */
        tmp         = spm->rows;
        spm->rows   = spm->colptr;
        spm->colptr = tmp;
        spm->fmttype = PastixCSR;
    }
    }

    return result;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * z_spmConvertIJV2CSR - convert a matrix in IJV format to a matrix in CSR
 * format.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The ijv matrix at enter,
 *          the csr matrix at exit.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertIJV2CSR( pastix_csc_t *spm )
{
#if !defined(PRECISION_p)
    pastix_complex64_t *navals = NULL;
    pastix_complex64_t *oavals = NULL;
#endif
    pastix_int_t       *spmptx, *otmp;
    pastix_int_t i, j, tmp, baseval, total;
    pastix_csc_t oldspm;

    /* Backup the input */
    memcpy( &oldspm, spm, sizeof(pastix_csc_t) );

    /*
     * Check the baseval, we consider that arrays are sorted by columns or rows
     */
    baseval = pastix_imin( *(oldspm.colptr), *(oldspm.rows) );

    /* Compute the new rows (should be called rowptr) */
    spm->rows = (pastix_int_t *) calloc(spm->n+1,sizeof(pastix_int_t));

    /* Compute the number of edges per row */
    spmptx = spm->rows - baseval;
    otmp   = oldspm.rows;
    for (i=0; i<spm->nnz; i++, otmp++)
    {
        spmptx[ *otmp ] ++;
    }

    /* Compute the indexes in C numbering for the following sort */
    total = 0;
    spmptx = spm->rows;
    for (i=0; i<(spm->n+1); i++, spmptx++)
    {
        tmp = *spmptx;
        *spmptx = total;
        total += tmp;
    }
    assert( total == spm->nnz );

    /* Sort the colptr and avals arrays by rows */
    spm->colptr  = malloc(spm->nnz * sizeof(pastix_int_t));

#if defined(PRECISION_p)
    spm->avals = NULL;
#else
    spm->avals = malloc(spm->nnz * sizeof(pastix_complex64_t));
    navals = (pastix_complex64_t*)(spm->avals);
    oavals = (pastix_complex64_t*)(oldspm.avals);
#endif

    for (j=0; j<spm->nnz; j++)
    {
        i = oldspm.rows[j] - baseval;

        spm->colptr[ spm->rows[i] ] = oldspm.colptr[j];

#if !defined(PRECISION_p)
        navals[ spm->rows[i] ] = oavals[j];
#endif
        (spm->rows[i])++;

        assert( spm->rows[i] <= spm->rows[i+1] );
    }

    /* Rebuild the rows (rowptr) with the correct baseval */
    tmp = spm->rows[0];
    spm->rows[0] = baseval;

    spmptx = spm->rows + 1;
    for (i=1; i<(spm->n+1); i++, spmptx++)
    {
        total = *spmptx;
        *spmptx = tmp + baseval;
        tmp = total;
    }
    assert( spm->rows[ spm->n ] == (spm->nnz+baseval) );

    free( oldspm.rows );
    free( oldspm.colptr );

    if (oldspm.avals != NULL)
        free( oldspm.avals );

    spm->fmttype = PastixCSR;

    return PASTIX_SUCCESS;
}