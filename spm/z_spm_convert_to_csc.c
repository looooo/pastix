/**
 *
 * @file z_spm_convert_to_csc.c
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
 * z_spmConvertIJV2CSC - convert a matrix in IJV format to a matrix in CSC
 * format.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The ijv matrix at enter,
 *          the csc matrix at exit.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertIJV2CSC( pastix_spm_t *spm )
{
#if !defined(PRECISION_p)
    pastix_complex64_t *navals = NULL;
    pastix_complex64_t *oavals = NULL;
#endif
    pastix_int_t       *spmptx, *otmp;
    pastix_int_t i, j, tmp, baseval, total;
    pastix_spm_t oldspm;

    /* Backup the input */
    memcpy( &oldspm, spm, sizeof(pastix_spm_t) );

    /*
     * Check the baseval, we consider that arrays are sorted by columns or rows
     */
    baseval = spmFindBase( spm );

    /* Compute the new colptr */
    spm->colptr = (pastix_int_t *) calloc(spm->n+1,sizeof(pastix_int_t));

    /* Compute the number of edges per row */
    spmptx = spm->colptr - baseval;
    otmp   = oldspm.colptr;
    for (i=0; i<spm->nnz; i++, otmp++)
    {
        spmptx[ *otmp ] ++;
    }

    /* Compute the indexes in C numbering for the following sort */
    total = 0;
    spmptx = spm->colptr;
    for (i=0; i<(spm->n+1); i++, spmptx++)
    {
        tmp = *spmptx;
        *spmptx = total;
        total += tmp;
    }
    assert( total == spm->nnz );

    /* Sort the rows and avals arrays by column */
    spm->rowptr = malloc(spm->nnz * sizeof(pastix_int_t));
    pastix_int_t *node     = calloc(spm->nnz+1,sizeof(pastix_int_t));
    pastix_int_t *old_node = calloc(spm->nnz+1,sizeof(pastix_int_t));
    pastix_int_t *dofs     = spm->dofs;
    pastix_int_t row, col, dofi, dofj;

    for(i=0; i<spm->nnz; i++)
    {
        row = oldspm.rowptr[i]-baseval;
        col = oldspm.colptr[i]-baseval;
        dofi = ( spm->dof > 0 ) ? spm->dof : dofs[col+1] - dofs[col];
        dofj = ( spm->dof > 0 ) ? spm->dof : dofs[row+1] - dofs[row];
        old_node[i+1] = dofi * dofj;
        node[spm->colptr[col]+1] += dofi * dofj;
        spm->colptr[col]++;
    }

    for(i=0;i<spm->nnz;i++)
    {
        node[i+1]+=node[i];
        old_node[i+1]+=old_node[i];
    }

    /* Restore the colptr indexes */
    {
        pastix_int_t tmp, tmp2;
        tmp = spm->colptr[0];
        spm->colptr[0] = 0;
        for (j=0; j<spm->n; j++) {
            tmp2 = spm->colptr[j+1];
            spm->colptr[j+1] = tmp;
            tmp = tmp2;
        }
    }

#if defined(PRECISION_p)
    spm->values = NULL;
#else
    pastix_int_t ii, jj;
    spm->values = malloc(spm->nnzexp * sizeof(pastix_complex64_t));
    navals = (pastix_complex64_t*)(spm->values);
    oavals = (pastix_complex64_t*)(oldspm.values);
#endif

    for (j=0; j<spm->nnz; j++)
    {
        i = oldspm.colptr[j] - baseval;
        spm->rowptr[ spm->colptr[i] ] = oldspm.rowptr[j];

#if !defined(PRECISION_p)
        dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
        dofj = ( spm->dof > 0 ) ? spm->dof : dofs[oldspm.rowptr[j]-baseval+1] - dofs[oldspm.rowptr[j]-baseval];
        for(ii=0; ii<dofi; ii++)
        {
            for(jj=0; jj<dofj; jj++)
            {
                navals[node[spm->colptr[i]]] = oavals[ old_node[j]];
                old_node[j]++;
                node[spm->colptr[i]]++;
            }
        }
#endif
        (spm->colptr[i])++;

        assert( spm->colptr[i] <= spm->colptr[i+1] );
    }

    /* Rebuild the colptr with the correct baseval */
    tmp = spm->colptr[0];
    spm->colptr[0] = baseval;
    spmptx = spm->colptr + 1;
    for (i=1; i<(spm->n+1); i++, spmptx++)
    {
        total = *spmptx;
        *spmptx = tmp + baseval;
        tmp = total;
    }
    assert( spm->colptr[ spm->n ] == (spm->nnz+baseval) );

    spmExit( &oldspm );

    spm->fmttype = PastixCSC;
    spm->colmajor = 1;
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * z_spmConvertCSR2CSC - convert a matrix in CSR format to a matrix in CSC
 * format. If the matrix is PastixSymmetric or PastixHermitian, then the
 * transpose or respectively the conjugate is returned.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The csr matrix at enter,
 *          the csc matrix at exit.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertCSR2CSC( pastix_spm_t *spm )
{
    assert( spm->loc2glob == NULL );
    assert( spm->fmttype == PastixCSR );

    spm->fmttype = PastixCSC;

    pastix_int_t type = spm->mtxtype;
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
        spm->fmttype = PastixCSC;

        return PASTIX_SUCCESS;
    }
    break;

    case PastixGeneral:
    default:
    {
        pastix_int_t       *row_csc;
        pastix_int_t       *col_csc;
        pastix_int_t       *node;
        pastix_int_t       *dofs;

#if !defined(PRECISION_p)
        pastix_complex64_t *val_csc;
        pastix_complex64_t *valptr = (pastix_complex64_t*)(spm->values);
        pastix_int_t ii, jj;
        int    cpt = 0;
#endif
        pastix_int_t  i, j, k, col, row, nnz, baseval;
        pastix_int_t dofi, dofj;

        baseval = spmFindBase( spm );
        nnz = spm->nnz;

        row_csc = malloc(nnz * sizeof(pastix_int_t));
        col_csc = calloc(spm->n+1,sizeof(pastix_int_t));
        node    = calloc(spm->nnz+1,sizeof(pastix_int_t));
        dofs    = spm->dofs;

        assert( row_csc );
        assert( col_csc );
        assert( node );
        assert( ( (dofs) && !(spm->dof > 0) ) ||
                ( !(dofs) && (spm->dof > 0) ) ); // (dofs) xor (spm->dof > 0)

#if !defined(PRECISION_p)
        val_csc = malloc(spm->nnzexp*sizeof(pastix_complex64_t));
        assert( val_csc );
#endif
        /* Count the number of elements per column */
        for (j=0; j<nnz; j++) {
            col = spm->colptr[j] - baseval;
            assert(col < spm->n );
            col_csc[ col+1 ] ++;
        }

        /* Compute the index of each column */
        col_csc[0] = 0;
        for (j=0; j<spm->n; j++){
            col_csc[j+1] += col_csc[j];
        }

        for(i=0; i<spm->n; i++)
        {
            for(k=spm->rowptr[i]; k<spm->rowptr[i+1]; k++)
            {
                j = spm->colptr[k-baseval] - baseval;
                row =  ( spm->dof > 0 ) ? j        : dofs[j];
                col =  ( spm->dof > 0)  ? i        : dofs[i];
                dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
                dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];
                node[col_csc[j]+1] += dofi * dofj;
                col_csc[j]++;
            }
        }

        for(i=0;i <spm->nnz; i++)
            node[i+1] += node[i];

        /* Restore the colptr indexes */
        {
            pastix_int_t tmp, tmp2;
            tmp = col_csc[0];
            col_csc[0] = 0;
            for (j=0; j<spm->n; j++) {
                tmp2 = col_csc[j+1];
                col_csc[j+1] = tmp;
                tmp = tmp2;
            }
        }

        assert( (col_csc[spm->gN]) == nnz );

        for (row=0; row<spm->n; row++) {
            pastix_int_t fcol = spm->rowptr[row  ] - baseval;
            pastix_int_t lcol = spm->rowptr[row+1] - baseval;

            for (k=fcol; k<lcol; k++) {
                col = spm->colptr[k] - baseval;
                j = col_csc[col];
                row_csc[j] = row + baseval;
#if !defined(PRECISION_p)
                dofi = ( spm->dof > 0 ) ? spm->dof : dofs[col+1] - dofs[col];
                dofj = ( spm->dof > 0 ) ? spm->dof : dofs[row+1] - dofs[row];
                for(jj=0; jj < dofj ; jj++)
                {
                    for(ii=0; ii < dofi ; ii++)
                    {
                        val_csc[node[j] + ii * dofj + jj] = valptr[cpt];
                        cpt++;
                    }
                }
#endif
                col_csc[col]++;
            }
        }

        /* Restore the colptr indexes */
        {
            pastix_int_t tmp, tmp2;
            tmp = col_csc[0];
            col_csc[0] = baseval;
            for (j=0; j<spm->n; j++)
            {
                tmp2 = col_csc[j+1];
                col_csc[j+1] = tmp + baseval;
                tmp = tmp2;
            }
        }
        spmExit( spm );
        spm->colptr = col_csc;
        spm->rowptr = row_csc;
#if !defined(PRECISION_p)
        spm->values = val_csc;
#else
        spm->values = NULL;
#endif
    }
    }

    if(spm-> dof != 1)
        spm->mtxtype = type;
    spm->colmajor = 1;
    return PASTIX_SUCCESS;
}
