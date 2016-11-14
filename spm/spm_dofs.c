/**
 *
 * @file spm_dofs.c
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
 **/
#include "common.h"
#include "spm.h"

#include "z_spm.h"
#include "c_spm.h"
#include "d_spm.h"
#include "s_spm.h"
#include "p_spm.h"

pastix_spm_t *
spmDofExtend( const int type,
              const int dof,
              const pastix_spm_t *spm )
{
    pastix_spm_t *newspm;

    /* Quick return */
    if ( dof == 1 )
        return (pastix_spm_t *)spm;

    if ( spm->dof != 1 ) {
        pastix_error_print( "Cannot extend spm including dofs already\n" );
        return (pastix_spm_t *)spm;
    }

    newspm = spmCopy( spm );

    /**
     * Generate constant dof
     */
    if (type == 0) {
        newspm->dof     = dof;
        newspm->nexp    = spm->n  * dof;
        newspm->gNexp   = spm->gN * dof;
        newspm->nnzexp  = spm->nnz  * dof * dof;
        newspm->gnnzexp = spm->gnnz * dof * dof;
    }
    else {
        pastix_int_t i, j, k, dofi, dofj, baseval;
        pastix_int_t *dofptr, *colptr, *rowptr;

        baseval = spmFindBase( spm );

        newspm->dof  = -1;
        newspm->dofs = malloc( (spm->n+1) * sizeof(pastix_int_t) );
        dofptr = newspm->dofs;

        /**
         * Initialize the dofs array where the degree of freedom of vertex i is
         * dof[i+1] - dof[i]
         */
        *dofptr = baseval;
        for(i=0; i<spm->n; i++, dofptr++) {
            dofi = 1 + ( rand() % dof );
            dofptr[1] = dofptr[0] + dofi;
        }

        newspm->nexp  = *dofptr - baseval;
        newspm->gNexp = newspm->nexp;

        /**
         * Count the number of non zeroes
         */
        newspm->nnzexp = 0;
        colptr = newspm->colptr;
        rowptr = newspm->rowptr;
        dofptr = newspm->dofs;

        switch(spm->fmttype)
        {
        case PastixCSR:
            /* Swap pointers to call CSC */
            colptr = newspm->rowptr;
            rowptr = newspm->colptr;

        case PastixCSC:
            for(j=0; j<newspm->n; j++, colptr++) {
                dofj = dofptr[j+1] - dofptr[j];

                for(k=colptr[0]; k<colptr[1]; k++, rowptr++) {
                    i = *rowptr - baseval;
                    dofi = dofptr[i+1] - dofptr[i];

                    newspm->nnzexp += dofi * dofj;
                }
            }
            break;
        case PastixIJV:
            for(k=0; k<newspm->nnz; k++, rowptr++, colptr++)
            {
                i = *rowptr - baseval;
                j = *colptr - baseval;
                dofi = dofptr[i+1] - dofptr[i];
                dofj = dofptr[j+1] - dofptr[j];

                newspm->nnzexp += dofi * dofj;
            }
        }
        newspm->gnnzexp = newspm->nnzexp;
    }

    switch (spm->flttype) {
    case PastixFloat:
        s_spmDofExtend( newspm );
        break;

    case PastixDouble:
        d_spmDofExtend( newspm );
        break;

    case PastixComplex32:
        c_spmDofExtend( newspm );
        break;

    case PastixComplex64:
        z_spmDofExtend( newspm );
        break;

    case PastixPattern:
        ;
    }

    return newspm;
}
