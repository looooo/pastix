/**
 *
 * @file z_spm_dof_extend.c
 *
 * SParse Matrix package multi-dof matrix expanser.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Alban Bellot
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 **/
#include "common.h"
#include "spm.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_dev
 *
 * @brief Extend a multi-dof sparse matrix to a single dof sparse matrix.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The sparse matrix to extend.
 *
 *******************************************************************************/
void
z_spmDofExtend(pastix_spm_t *spm)
{
    pastix_int_t        i, j, k, ii, jj, dofi, dofj, baseval;
    pastix_int_t       *colptr, *rowptr, *dofs;
    pastix_complex64_t *newval, *oldval, *oldvalptr;

    oldval = oldvalptr = (pastix_complex64_t*)(spm->values);
    newval = spm->values = malloc( spm->nnzexp * sizeof(pastix_complex64_t) );

    baseval = spmFindBase( spm );
    colptr = spm->colptr;
    rowptr = spm->rowptr;
    dofs   = spm->dofs;

    switch(spm->fmttype)
    {
    case PastixCSR:
        /* Swap pointers to call CSC */
        colptr = spm->rowptr;
        rowptr = spm->colptr;

    case PastixCSC:
        /**
         * Loop on col
         */
        for(j=0; j<spm->n; j++, colptr++)
        {
            dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];

            /**
             * Loop on rows
             */
            for(k=colptr[0]; k<colptr[1]; k++, rowptr++, oldval++)
            {
                i = *rowptr - baseval;
                dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];

                for(jj=0; jj<dofj; jj++)
                {
                    for(ii=0; ii<dofi; ii++, newval++)
                    {
                        *newval = *oldval;
                    }
                }
            }
        }
        break;
    /* case PastixCSR: */
    /*     /\** */
    /*      * Loop on row */
    /*      *\/ */
    /*     for(i=0; i<spm->n; i++, rowptr++) */
    /*     { */
    /*         dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i]; */

    /*         /\** */
    /*          * Loop on cols */
    /*          *\/ */
    /*         for(k=rowptr[0]; k<rowptr[1]; k++, colptr++, oldval++) */
    /*         { */
    /*             j = *colptr - baseval; */
    /*             dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j]; */

    /*             for(jj=0; jj<dofj; jj++) */
    /*             { */
    /*                 for(ii=0; ii<dofi; ii++, newval++) */
    /*                 { */
    /*                     *newval = *oldval; */
    /*                 } */
    /*             } */
    /*         } */
    /*     } */
    /*     break; */
    case PastixIJV:
        /**
         * Loop on coordinates
         */
        for(k=0; k<spm->nnz; k++, rowptr++, colptr++, oldval++)
        {
            i = *rowptr - baseval;
            j = *colptr - baseval;
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];
            dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];

            for(jj=0; jj<dofj; jj++)
            {
                for(ii=0; ii<dofi; ii++, newval++)
                {
                    *newval = *oldval;
                }
            }
        }
        break;
    }

    free(oldvalptr);
    return;
}
