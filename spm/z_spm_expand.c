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
 * @author Theophile Terraz
 * @author Alban Bellot
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s p
 **/
#include "common.h"
#include "spm.h"
#include "z_spm.h"


void z_extandCSC(pastix_spm_t* spm)
{
  pastix_int_t col, row, i, cpt, dofj, dofi, baseval;
  cpt=0;

  pastix_complex64_t *valptr = (pastix_complex64_t*)spm->values;
  pastix_complex64_t *newval = calloc(spm->nnzexp,sizeof(pastix_complex64_t));

  baseval=spmFindBase( spm );

  for( col=0; col<spm->n; col++)
  {
      for( row=spm->colptr[col]-baseval; row<spm->colptr[col+1]-baseval; row++)
      {
          dofi = ( spm->dof > 0 ) ? spm->dof : spm->dofs[col+1] - spm->dofs[col];
          dofj = ( spm->dof > 0 ) ? spm->dof : spm->dofs[spm->rowptr[row]-baseval+1] - spm->dofs[spm->rowptr[row]-baseval];
	  for( i=0; i<dofi*dofj; i++)
          {
              newval[cpt] = valptr[row] / ((i/dofj) + (i%dofj) + 2); // Col major
	      cpt++;
          }
      }
  }
  spm->values=newval;
}
