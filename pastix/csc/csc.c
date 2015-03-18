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

static int (*conversionTable[3][3][6])(pastix_csc_t*) = {
    /* From CSC */
    {{ NULL, NULL, NULL, NULL, NULL, NULL },
     { p_spmConvertCSC2CSR,
       NULL,
       s_spmConvertCSC2CSR,
       d_spmConvertCSC2CSR,
       c_spmConvertCSC2CSR,
       z_spmConvertCSC2CSR },
     { p_spmConvertCSC2IJV,
       NULL,
       s_spmConvertCSC2IJV,
       d_spmConvertCSC2IJV,
       c_spmConvertCSC2IJV,
       z_spmConvertCSC2IJV }},
    /* From CSR */
    {{ p_spmConvertCSR2CSC,
       NULL,
       s_spmConvertCSR2CSC,
       d_spmConvertCSR2CSC,
       c_spmConvertCSR2CSC,
       z_spmConvertCSR2CSC },
     { NULL, NULL, NULL, NULL, NULL, NULL },
     { p_spmConvertCSR2IJV,
       NULL,
       s_spmConvertCSR2IJV,
       d_spmConvertCSR2IJV,
       c_spmConvertCSR2IJV,
       z_spmConvertCSR2IJV }},
    /* From IJV */
    {{ p_spmConvertIJV2CSC,
       NULL,
       s_spmConvertIJV2CSC,
       d_spmConvertIJV2CSC,
       c_spmConvertIJV2CSC,
       z_spmConvertIJV2CSC },
     { p_spmConvertIJV2CSR,
       NULL,
       s_spmConvertIJV2CSR,
       d_spmConvertIJV2CSR,
       c_spmConvertIJV2CSR,
       z_spmConvertIJV2CSR },
     { NULL, NULL, NULL, NULL, NULL, NULL }}
};



int
spmConvert( int ofmttype, pastix_csc_t *ospm )
{
    if ( conversionTable[ospm->fmttype][ofmttype][ospm->flttype] ) {
        return conversionTable[ospm->fmttype][ofmttype][ospm->flttype]( ospm );
    }
    else {
        return PASTIX_SUCCESS;
    }
}

pastix_int_t
spmFindBase( pastix_csc_t *spm )
{

    pastix_int_t i, *tmp, baseval;

    /*
     * Check the baseval, we consider that arrays are sorted by columns or rows
     */
    baseval = pastix_imin( *(spm->colptr), *(spm->rowptr) );
    /*
     * if not:
     */
    if ( ( baseval != 0 ) &&
         ( baseval != 1 ) )
    {
        baseval = spm->n;
        tmp = spm->colptr;
        for(i=0; i<spm->nnz; i++, tmp++){
            baseval = pastix_imin( *tmp, baseval );
        }
    }

    return baseval;
}

void
spmExit( pastix_csc_t *spm )
{
    if(spm->colptr != NULL)
        memFree_null(spm->colptr);
    if(spm->rowptr != NULL)
        memFree_null(spm->rowptr);
    if(spm->loc2glob != NULL)
        memFree_null(spm->loc2glob);
    if(spm->values != NULL)
        memFree_null(spm->values);
}
