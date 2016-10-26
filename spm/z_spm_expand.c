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
pastix_spm_t *
z_spmExpand(const pastix_spm_t *spm)
{
    pastix_spm_t *newspm;
    pastix_int_t i, col, row, cpt, dofj, dofi, baseval;
#if !defined(PRECISION_p)
    pastix_complex64_t *oldvalptr;
    pastix_complex64_t *newvalptr;
#endif

    if (spm->dof == 1) {
        return spm;
    }

    if (1) {
        return NULL;
    }

    
/*     baseval = spmFindBase( spm ); */

/* #if !defined(PRECISION_p) */
/*     oldvalptr = (pastix_complex64_t*)spm->values; */
/*     spm->values = malloc( spm->nnzexp * sizeof(pastix_complex64_t) ); */
/*     newvalptr = (pastix_complex64_t*)spm->values; */
/* #endif */

/*     cpt = 0; */
/*     dofi = spm->dof; */
/*     dofj = spm->dof; */

/*     for( col=0; col<spm->n; col++) */
/*     { */
/*         if ( spm->dof <= 0 ) { */
/*             dofi = spm->dofs[col+1] - spm->dofs[col]; */
/*         } */

/*         for( row=spm->colptr[col]-baseval; row<spm->colptr[col+1]-baseval; row++) */
/*         { */
/*             if ( spm->dof <= 0 ) { */
/*                 dofj = spm->dofs[spm->rowptr[row]-baseval+1] - spm->dofs[spm->rowptr[row]-baseval]; */
/*             } */

/*             for( i=0; i<dofi*dofj; i++) */
/*             { */
/* #if !defined(PRECISION_p) */
/*                 newvalptr[cpt] = oldvalptr[row] / ((i/dofj) + (i%dofj) + 2); // Col major */
/*                 cpt++; */
/* #endif */
/*             } */
/*         } */
/*     } */

/* #if !defined(PRECISION_p) */
/*     free( oldvalptr ); */
/* #endif */
}
