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
    pastix_spm_t       *newspm;
    pastix_int_t        i, j, k, ii, jj, dofi, dofj, col, row, baseval, cpt;
    pastix_int_t       *newcol, *newrow, *oldcol, *oldrow, *dofs;
#if !defined(PRECISION_p)
    pastix_complex64_t *newval, *oldval;
#endif

    if ( spm->dof == 1 ) {
        return (pastix_spm_t*)spm;
    }

    if ( spm->layout != PastixColMajor ) {
        pastix_error_print( "Unsupported layout\n" );
        return NULL;
    }

    newspm = malloc( sizeof(pastix_spm_t) );
    spmInit( newspm );

    baseval = spmFindBase( spm );
    oldcol = spm->colptr;
    oldrow = spm->rowptr;
    dofs   = spm->dofs;
#if !defined(PRECISION_p)
    oldval = (pastix_complex64_t*)(spm->values);
#endif

    switch(spm->fmttype)
    {
    case PastixCSC:
        newspm->colptr = newcol = malloc(sizeof(pastix_int_t)*spm->nexp);
        newspm->rowptr = newrow = malloc(sizeof(pastix_int_t)*spm->nnzexp);
#if !defined(PRECISION_p)
        newspm->values  = newval = malloc(sizeof(pastix_complex64_t)*spm->nnzexp);
#endif

        *newcol = 0; newcol++;
        *newrow = 0;
        /**
         * Loop on col
         */
        for(i=0; i<spm->n; i++)
        {
            if ( spm->dof > 0 ) {
                col  = spm->dof * i;
                dofi = spm->dof;
            }
            else {
                col  = dofs[i];
                dofi = dofs[i+1] - dofs[i];
            }

            for(ii=0; ii<dofi; ii++, newcol++)
            {
                /**
                 * Loop on rows
                 */
                for(k=oldcol[i]; k<oldcol[i+1]; k++)
                {
                    j = oldrow[k-baseval]-baseval;
                    if ( spm->dof > 0 ) {
                        row  = spm->dof * j;
                        dofj = spm->dof;
                    }
                    else {
                        row  = dofs[j];
                        dofj = dofs[j+1] - dofs[j];
                    }

                    for(jj=0; jj<dofj; jj++, newrow++)
                    {
                        (*newcol)++;
                        (*newrow) = row + jj + baseval;

#if !defined(PRECISION_p)
                        if ( (spm->mtxtype != PastixGeneral) &&
                             (row + jj < col + ii) )
                        {
                            (*newval) = oldval[ cpt ];
                            newval++;
                        }
                        cpt++;
#endif
                    }
                }
                (*newcol) += baseval;
            }
        }
        break;
    case PastixCSR:
    case PastixIJV:
        free( newspm );
        return NULL;
    }
        /*     for(i=0; i<spm->nexp; i++) */
    /*     { */
    /*         new_col[i+1]+=new_col[i]; */
    /*     } */

    /*     cpt = 0; */
    /*     for(i=0; i < spm->n;i++) */
    /*     { */
    /*         col  = ( spm->dof > 0 ) ? i        : dofs[i]; */
    /*         dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i]; */
    /*         for(k=spm->colptr[i]-baseval ; k<spm->colptr[i+1]-baseval ;k++) */
    /*         { */
    /*             j = spm->rowptr[k] - baseval; */
    /*             row  = ( spm->dof > 0 ) ? j        : dofs[j]; */
    /*             dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j]; */
    /*             for(ii=0;ii < dofi; ii++) */
    /*             { */
    /*                 for(jj=0;jj < dofj ; jj++) */
    /*                 { */
    /*                     new_vals[new_col[col+ii]] = vals[cpt]; */
    /*                     new_row[new_col[col+ii]]  = row + jj + baseval; */
    /*                     new_col[col+ii]++; */
    /*                     cpt++; */
    /*                 } */
    /*             } */
    /*         } */
    /*     } */

    /*     { */
    /*         int tmp; */
    /*         int tmp1 = 0; */
    /*         for(i=0; i<spm->nexp; i++) */
    /*         { */
    /*             tmp = new_col[i]; */
    /*             new_col[i] = tmp1+baseval; */
    /*             tmp1 = tmp; */
    /*         } */
    /*         new_col[i] += baseval; */
    /*     } */
    /*     spm->gN   = spm->gNexp; */
    /*     spm->n    = spm->nexp; */
    /*     spm->gnnz = spm->gnnzexp; */
    /*     spm->nnz  = spm->nnzexp; */

    /*     spm->dof      = 1; */
    /*     spm->dofs     = NULL; */
    /*     spm->layout   = PastixColMajor; */

    /*     spm->colptr   = new_col; */
    /*     spm->rowptr   = new_row; */
    /*     //spm->loc2glob = NULL; // ? */
    /*     spm->values   = new_vals; */
    /*     break; */

    /* case PastixSymmetric: */
    /*     for(i=0;i<spm->n ; i++) */
    /*     { */
    /*         col  = ( spm->dof > 0 ) ? i        : dofs[i]; */
    /*         dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i]; */
    /*         for(k=spm->colptr[i]-baseval; k<spm->colptr[i+1]-baseval; k++) */
    /*         { */
    /*             j = spm->rowptr[k]-baseval; */
    /*             row  = ( spm->dof > 0 ) ? j        : dofs[j]; */
    /*             dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j]; */
    /*             for(ii=0; ii<dofi; ii++) */
    /*             { */
    /*                 for(jj=0; jj<dofj; jj++) */
    /*                 { */
    /*                     if( i != j ) */
    /*                         new_col[col+ii+1] +=  1; */
    /*                     else */
    /*                         if(ii <= jj ) */
    /*                             new_col[col+ii+1] += 1; */
    /*                 } */
    /*             } */
    /*         } */
    /*     } */
    /*     for(i=0; i<spm->nexp; i++) */
    /*     { */
    /*         new_col[i+1] += new_col[i]; */
    /*     } */
    /*     pastix_int_t nnz = new_col[spm->nexp]; */
    /*     new_row  = malloc(sizeof(pastix_int_t)*nnz); */
    /*     new_vals = malloc(sizeof(pastix_complex64_t)*nnz); */

    /*     cpt = 0; */
    /*     for(i=0; i < spm->n;i++) */
    /*     { */
    /*         col  = ( spm->dof > 0 ) ? i        : dofs[i]; */
    /*         dofi = ( spm->dof > 0 ) ? spm->dof : dofs[i+1] - dofs[i]; */
    /*         for(k=spm->colptr[i]-baseval ; k<spm->colptr[i+1]-baseval ;k++) */
    /*         { */
    /*             j = spm->rowptr[k] - baseval; */
    /*             row  = ( spm->dof > 0 ) ? j        : dofs[j]; */
    /*             dofj = ( spm->dof > 0 ) ? spm->dof : dofs[j+1] - dofs[j]; */
    /*             for(ii=0;ii < dofi; ii++) */
    /*             { */
    /*                 for(jj=0;jj < dofj ; jj++) */
    /*                 { */
    /*                     if( i == j ) */
    /*                     { */
    /*                         if ( ii <= jj ) */
    /*                         { */
    /*                             /\* diagonal dominant for spd matrix */
    /*                             if( ii == jj) */
    /*                                 new_vals[new_col[col+ii]] = 2*vals[cpt]; */
    /*                              else */
    /*                             *\/ */
    /*                             new_vals[new_col[col+ii]] = vals[cpt]; */
    /*                             new_row[new_col[col+ii]]  = row + jj + baseval; */
    /*                             new_col[col+ii]++; */
    /*                         } */
    /*                     } */
    /*                     else */
    /*                     { */
    /*                         new_vals[new_col[col+ii]] = vals[cpt]; */
    /*                         new_row[new_col[col+ii]]  = row + jj + baseval; */
    /*                         new_col[col+ii]++; */

    /*                     } */
    /*                     cpt++; */
    /*                 } */
    /*             } */
            /* } */
    /* } */

    newspm->gN      = spm->gNexp;
    newspm->n       = spm->nexp;
    newspm->gnnz    = spm->gnnzexp;
    newspm->nnz     = spm->nnzexp;
    newspm->gnnzexp = spm->gnnzexp;
    newspm->nnzexp  = spm->nnzexp;

    newspm->dof      = 1;
    newspm->dofs     = NULL;
    newspm->layout   = PastixColMajor;

    assert(spm->loc2glob == NULL);//to do

    (void)col; (void)cpt;
    return newspm;
}
