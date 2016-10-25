/**
 *
 * @file z_spm_expand.c
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
 * @precisions normal z -> c d s p
 **/
#include "common.h"
#include "spm.h"
#include "z_spm.h"

/**
 * TODO: This function is incorrect
 */
int
z_spmExpand(pastix_spm_t *spm)
{
    pastix_int_t i, col, row, cpt, dofj, dofi, baseval;
    pastix_complex64_t *oldvalptr;
    pastix_complex64_t *newvalptr;

    if (1) {
        return PASTIX_ERR_NOTIMPLEMENTED;
    }

    baseval = spmFindBase( spm );

    oldvalptr = (pastix_complex64_t*)spm->values;
    spm->values = malloc( spm->nnzexp * sizeof(pastix_complex64_t) );
    newvalptr = (pastix_complex64_t*)spm->values;

    cpt = 0;
    dofi = spm->dof;
    dofj = spm->dof;

    for( col=0; col<spm->n; col++)
    {
        if ( spm->dof <= 0 ) {
            dofi = spm->dofs[col+1] - spm->dofs[col];
        }

        for( row=spm->colptr[col]-baseval; row<spm->colptr[col+1]-baseval; row++)
        {
            if ( spm->dof <= 0 ) {
                dofj = spm->dofs[spm->rowptr[row]-baseval+1] - spm->dofs[spm->rowptr[row]-baseval];
            }

            for( i=0; i<dofi*dofj; i++)
            {
                newvalptr[cpt] = oldvalptr[row] / ((i/dofj) + (i%dofj) + 2); // Col major
                cpt++;
            }
        }
    }

    free( oldvalptr );
}
