/**
 *
 * @file z_spm_dof_extend.c
 *
 *  PaStiX spm routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Alban Bellot
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 **/
#include "common.h"
#include "spm.h"
#include "z_spm.h"

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
    case PastixCSR:
        /**
         * Loop on row
         */
        for(i=0; i<spm->n; i++, rowptr++)
        {
            dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i];

            /**
             * Loop on cols
             */
            for(k=rowptr[0]; k<rowptr[1]; k++, colptr++, oldval++)
            {
                j = *colptr - baseval;
                dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j];

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
