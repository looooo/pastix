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


int
z_spmConvertCSC2CSR( int ofmttype, pastix_csc_t *spm )
{
    pastix_int_t       *row_csr;
    pastix_int_t       *col_csr;
#if !defined(PRECISION_p)
    pastix_complex64_t *val_csr;
    pastix_complex64_t  val;
    pastix_complex64_t *valptr;
#endif
    pastix_int_t       *count;
    pastix_int_t j, k, col, row, nnz, baseval;

    baseval = pastix_imin( *(spm->colptr), *(spm->rows) );
    nnz=spm->colptr[spm->gN]-baseval;
    spm->fmttype=PastixCSR;

    row_csr = calloc(spm->gN+1,sizeof(pastix_int_t));
    col_csr = malloc(nnz*sizeof(pastix_int_t));
    count = calloc(spm->gN,sizeof(pastix_int_t));

    assert( row_csr );
    assert( col_csr );
    assert( count_csr );

#if !defined(PRECISION_p)
    val_csr = malloc(nnz*sizeof(pastix_complex64_t));
    assert( val_csr );
#endif

    for (j=0;j<=nnz;j++){
        row_csr[spm->rows[j]]+=1;
    }
    row_csr[0]=0;
    for (j=1;j<=spm->gN;j++){
        row_csr[j]+=row_csr[j-1];
    }
    row_csr[0]=baseval;

    assert( row_csr[spm->gN] == nnz+1 );

    for (col=1;col<=spm->gN;col++){
        for (k=0;k<(spm->colptr[col-baseval+1]-spm->colptr[col-baseval]);k++){
            row=spm->rows[spm->colptr[col-baseval]-baseval+k];
            col_csr[row_csr[row-baseval]-baseval+count[row-baseval]]=col;
#if !defined(PRECISION_p)
            valptr=spm->avals+spm->colptr[col-baseval]-baseval+k;
            val=*valptr;
            val_csr[row_csr[row-baseval]-baseval+count[row-baseval]]=val;
#endif
            count[row-baseval]+=1;
        }
    }

    memFree_null(count);
    memFree_null(spm->colptr);
    memFree_null(spm->rows);
    spm->colptr=col_csr;
    spm->rows  =row_csr;
#if !defined(PRECISION_p)
    memFree_null(spm->avals);
    spm->avals =val_csr;
#endif

    return PASTIX_SUCCESS;
}

int
z_spmConvertIJV2CSR( int ofmttype, pastix_csc_t *spm )
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