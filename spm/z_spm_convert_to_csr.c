/**
 *
 * @file z_spm_convert_to_csr.c
 *
 *  PaStiX spm routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @authot Alban Bellot
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
 * @ingroup pastix_spm_internal
 *
 * z_spmConvertCSC2CSR - convert a matrix in CSC format to a matrix in CSR
 * format. If the matrix is PastixSymmetric or PastixHermitian, then the
 * transpose or respectively the conjugate is returned.
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
z_spmConvertCSC2CSR( pastix_spm_t *spm )
{
    pastix_int_t *tmp;
    pastix_int_t  result;
    pastix_int_t  type = spm->mtxtype;

    if(spm->dof != 1)
        spm->mtxtype = PastixGeneral;

    switch( spm->mtxtype ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PastixHermitian:
    {
        /* Similar to PastixSymmetric case with conjugate of the values */
        pastix_complex64_t *valptr = spm->values;
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
        tmp          = spm->rowptr;
        spm->rowptr  = spm->colptr;
        spm->colptr  = tmp;
        spm->fmttype = PastixCSR;
        return PASTIX_SUCCESS;
    }
    break;

    case PastixGeneral:
    default:
    {
        /* Transpose the spm in CSC to trans(spm) in CSR */
        tmp          = spm->rowptr;
        spm->rowptr  = spm->colptr;
        spm->colptr  = tmp;
        spm->fmttype = PastixCSR;

        /* Convert trans(spm) in CSR to trans(spm) in CSC */
        result = z_spmConvertCSR2CSC( spm );

        /* Transpose trans(spm) in CSC to obtain the spm in CSR */
        tmp          = spm->rowptr;
        spm->rowptr  = spm->colptr;
        spm->colptr  = tmp;
        spm->fmttype = PastixCSR;
    }
    }

    if(spm-> dof != 1)
        spm->mtxtype = type;
    spm->colmajor = -1;

    return result;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_internal
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
z_spmConvertIJV2CSR( pastix_spm_t *spm )
{
#if !defined(PRECISION_p)
    pastix_complex64_t *navals = NULL;
    pastix_complex64_t *oavals = NULL;
#endif
    pastix_int_t       *spmptx, *otmp;
    pastix_int_t i, j, tmp, baseval, total;
    pastix_int_t row, col, dofi, dofj;

    pastix_int_t *node     = calloc(spm->nnz+1,sizeof(pastix_int_t));
    pastix_int_t *old_node = calloc(spm->nnz+1,sizeof(pastix_int_t));
    pastix_int_t *dofs     = spm->dofs;

    pastix_spm_t oldspm;

    /*
     * Check the baseval, we consider that arrays are sorted by columns or rows
     */
    baseval = spmFindBase( spm );

#if !defined(PRECISION_p)
    pastix_int_t ii, jj, k;
    /* Transpose values in row major format */
    if( spm->colmajor )
    {
        pastix_int_t cpt=0;
        pastix_int_t* dofs = spm->dofs;

        oavals = (pastix_complex64_t*)spm->values;
        navals = malloc(sizeof(pastix_complex64_t)*spm->nnzexp);

        for(k=0; k<spm->nnz; k++)
        {
            j = spm->rowptr[k]-baseval;
            i = spm->colptr[k]-baseval;
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
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
    }
#endif

    /* Backup the input */
    memcpy( &oldspm, spm, sizeof(pastix_spm_t) );

    /* Compute the new rowptr */
    spm->rowptr = (pastix_int_t *) calloc(spm->n+1,sizeof(pastix_int_t));

    /* Compute the number of edges per row */
    spmptx = spm->rowptr - baseval;
    otmp   = oldspm.rowptr;
    for (i=0; i<spm->nnz; i++, otmp++)
    {
        spmptx[ *otmp ] ++;
    }

    /* Compute the indexes in C numbering for the following sort */
    total = 0;
    spmptx = spm->rowptr;
    for (i=0; i<(spm->n+1); i++, spmptx++)
    {
        tmp = *spmptx;
        *spmptx = total;
        total += tmp;
    }
    assert( total == spm->nnz );

    /* Sort the colptr and avals arrays by rows */
    spm->colptr  = malloc(spm->nnz * sizeof(pastix_int_t));

    for(i=0; i<spm->nnz; i++)
    {
        row = oldspm.rowptr[i]-baseval;
        col = oldspm.colptr[i]-baseval;
        dofi = ( spm->dof > 0 ) ? spm->dof : dofs[col+1] - dofs[col];
        dofj = ( spm->dof > 0 ) ? spm->dof : dofs[row+1] - dofs[row];
        old_node[i+1] = dofi * dofj;
        node[spm->rowptr[row]+1] += dofi * dofj;
        spm->rowptr[row]++;
    }

    for(i=0;i<spm->nnz;i++)
    {
        node[i+1]+=node[i];
        old_node[i+1]+=old_node[i];
    }

    /* Restore the rowptr indexes */
    {
        pastix_int_t tmp, tmp2;
        tmp = spm->rowptr[0];
        spm->rowptr[0] = 0;
        for (j=0; j<spm->n; j++) {
            tmp2 = spm->rowptr[j+1];
            spm->rowptr[j+1] = tmp;
            tmp = tmp2;
        }
    }

#if defined(PRECISION_p)
    spm->values = NULL;
#else
    spm->values = malloc(spm->nnzexp * sizeof(pastix_complex64_t));
    navals = (pastix_complex64_t*)(spm->values);
    oavals = (pastix_complex64_t*)(oldspm.values);
#endif

    for (j=0; j<spm->nnz; j++)
    {
        i = oldspm.rowptr[j] - baseval;

        spm->colptr[ spm->rowptr[i] ] = oldspm.colptr[j];

#if !defined(PRECISION_p)
        dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
        dofj = ( spm->dof > 0 ) ? spm->dof : dofs[oldspm.colptr[j]-baseval+1] - dofs[oldspm.colptr[j]-baseval];
        for(ii=0; ii<dofi; ii++)
        {
            for(jj=0; jj<dofj; jj++)
            {
                navals[node[spm->rowptr[i]]] = oavals[ old_node[j]];
                old_node[j]++;
                node[spm->rowptr[i]]++;
            }
        }
#endif
        (spm->rowptr[i])++;

        assert( spm->rowptr[i] <= spm->rowptr[i+1] );
    }

    /* Rebuild the rows (rowptr) with the correct baseval */
    tmp = spm->rowptr[0];
    spm->rowptr[0] = baseval;

    spmptx = spm->rowptr + 1;
    for (i=1; i<(spm->n+1); i++, spmptx++)
    {
        total = *spmptx;
        *spmptx = tmp + baseval;
        tmp = total;
    }
    assert( spm->rowptr[ spm->n ] == (spm->nnz+baseval) );

    spmExit( &oldspm );

    spm->colmajor = -1;
    spm->fmttype = PastixCSR;

    return PASTIX_SUCCESS;
}

