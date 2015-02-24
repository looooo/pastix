/**
 *
 * @file z_spm_convert_to_csc.c
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
z_spmConvertIJV2CSC( int ofmttype, pastix_csc_t *spm )
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
    spm->rows  = malloc(spm->nnz * sizeof(pastix_int_t));

#if defined(PRECISION_p)
    spm->avals = NULL;
#else
    spm->avals = malloc(spm->nnz * sizeof(pastix_complex64_t));
    navals = (pastix_complex64_t*)(spm->avals);
    oavals = (pastix_complex64_t*)(oldspm.avals);
#endif

    for (j=0; j<spm->nnz; j++)
    {
        i = oldspm.colptr[j] - baseval;

        spm->rows[ spm->colptr[i] ] = oldspm.rows[j];

#if !defined(PRECISION_p)
        navals[ spm->colptr[i] ] = oavals[j];
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

    free( oldspm.colptr );
    free( oldspm.rows );

    if (oldspm.avals != NULL)
        free( oldspm.avals );

    spm->fmttype = PastixCSC;

    return PASTIX_SUCCESS;
}

int
z_spmConvertCSR2CSC( int ofmttype, pastix_csc_t *spm )
{
    pastix_int_t       *row_csc;
    pastix_int_t       *col_csc;
#if !defined(PRECISION_p)
    pastix_complex64_t *val_csc;
    pastix_complex64_t  val;
    pastix_complex64_t *valptr;
#endif
    pastix_int_t       *count;
    pastix_int_t j, k, col, row, nnz, baseval;

    baseval = pastix_imin( *(spm->colptr), *(spm->rows) );
    nnz=spm->rows[spm->gN]-baseval;
    spm->fmttype=PastixCSC;

    row_csc = malloc(nnz*sizeof(pastix_int_t));
    col_csc = calloc(spm->gN+1,sizeof(pastix_int_t));
    count = calloc(spm->gN,sizeof(pastix_int_t));

    assert( row_csc );
    assert( col_csc );
    assert( count_csc );

#if !defined(PRECISION_p)
    val_csc = malloc(nnz*sizeof(pastix_complex64_t));
    assert( val_csc );
#endif

    for (j=0;j<=nnz;j++){
        col_csc[spm->colptr[j]]+=1;
    }
    col_csc[0]=0;
    for (j=1;j<=spm->gN;j++){
        col_csc[j]+=col_csc[j-1];
    }
    col_csc[0]=baseval;

    assert( col_csc[spm->gN] == nnz+1 );

    for (row=1;row<=spm->gN;row++){
        for (k=0;k<spm->rows[row-baseval+1]-spm->rows[row-baseval];k++){
            col=spm->colptr[spm->rows[row-baseval]-baseval+k];
            row_csc[col_csc[col-baseval]-baseval+count[col-baseval]]=row;
#if !defined(PRECISION_p)
            valptr=spm->avals+spm->rows[row-baseval]-baseval+k;
            val=*valptr;
            val_csc[col_csc[col-baseval]-baseval+count[col-baseval]]=val;
#endif
            count[col-baseval]+=1;
        }
    }

    memFree_null(count);
    memFree_null(spm->colptr);
    memFree_null(spm->rows);
    spm->colptr=col_csc;
    spm->rows  =row_csc;
#if !defined(PRECISION_p)
    memFree_null(spm->avals);
    spm->avals =val_csc;
#endif

    return PASTIX_SUCCESS;
}