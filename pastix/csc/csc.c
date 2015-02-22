/**
 *
 * @file csc.c
 *
 *  PaStiX csc routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "csc.h"

int
spmConvertIJV( int ofmttype, pastix_csc_t *spm )
{
    pastix_csc_t oldspm;

    /* Backup the input */
    memcpy( oldspm, spm, sizeof(pastix_csc_t));

    switch( ofmttype ) {
    case PastixCSC:
        return PASTIX_SUCCESS;

    case PastixCSR:
        //return spmConvertCSR( ofmttype, ospm );
        break;

    case PastixIJV:
    {
        pastix_int_t *spmptx, *otmp;
        pastix_int_t i, tmp, baseval, total;

        /* Check the baseval */
        baseval = 999999999999;
        otmp = oldspm->colptr;
        for (i=0; i<spm->nnz; i++, otmp++)
        {
            baseval = pastix_imin( baseval, *otmp );
        }

        /* Compute the new colptr */
        spm->colptr = (pastix_int_t *) calloc(spm->n+1,sizeof(pastix_int_t));

        /* Compute the number of edges per row */
        spmptx = spm->colptr - baseval;
        otmp   = oldspm->colptr;
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
        spm->rows  = (pastix_int_t *)malloc(spm->nnz * sizeof(pastix_int_t));
        spm->avals = (double *)      malloc(spm->nnz * sizeof(double));

        for (j=0; j<spm->nnz; j++)
        {
            i = oldspm->colptr[j] - baseval;

            spm->rows[  spm->colptr[i] ] = oldspm->rows[j];
            spm->avals[ spm->colptr[i] ] = oldspm->avals[j];

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

        free( oldspm->colptr );
        free( oldspm->rows );
        free( oldspm->avals );
        spm->fmttype = PastixIJV;
    }
    break;
    default:
        fprintf("spmConvertCSC: output format unknow\n");
    }

    return PASTIX_SUCCESS;
}

int
spmConvert( int ofmttype, pastix_csc_t *ospm )
{
    if (ospm->fmttype == fmttype) {
        return PASTIX_SUCCESS;
    }
    else {
        switch( ospm->fmttype ) {
        case PastixCSC:
            //return spmConvertCSC( ofmttype, ospm );
            break;

        case PastixCSR:
            //return spmConvertCSR( ofmttype, ospm );
            break;

        case PastixIJV:
            return spmConvertIJV( ofmttype, ospm );

        default:
            return PASTIX_ERR_BADPARAMETER;
        }
    }
    return PASTIX_SUCCESS;
}
