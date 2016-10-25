/**
 *
 * @file spm_expand.c
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

void print_tab_int(pastix_int_t* tab, pastix_int_t size)
{
    int i;
    for(i=0;i<size;i++)
        printf("%d ",tab[i]);
    printf("\n");
}

void print_tab_cpx(double* tab, pastix_int_t size)
{
    int i;
    for(i=0;i<size;i++)
        printf("%f ",tab[i]);
    printf("\n");
}


void init_rand()
{
  srand(time(NULL));
}

int randminmax(pastix_int_t min, pastix_int_t max)
{
  return min + rand()%(max-min+1);
}

pastix_int_t* computeCorrespondIndex(const pastix_int_t *dofs, pastix_int_t n)
{
    pastix_int_t* dofs_coef=malloc(sizeof(pastix_int_t)*(n+1));
    dofs_coef[0]=0;
    int i;
    for(i=1;i<n+1;i++)
    {
        dofs_coef[i]=dofs_coef[i-1]+dofs[i-1];
    }
    return dofs_coef;
}

void genDof(pastix_int_t step,pastix_int_t* dofs,pastix_int_t size)
{
  int i;
  for(i=0;i<size;i++)
      dofs[i]=i%step+1;
}

void genRandDof(pastix_int_t min,pastix_int_t max,pastix_int_t* dofs,pastix_int_t size)
{
  int i;
  for(i=0;i<size;i++)
    dofs[i]=randminmax(min,max);
}

pastix_int_t* computeEltPerCol(const pastix_spm_t *spm)
{
  pastix_int_t col,row,baseval;

  pastix_int_t* coef=malloc(sizeof(pastix_int_t)*spm->n);
  int i;
  for(i=0;i<spm->n;i++)
      coef[i]=0;
  baseval=spmFindBase( spm );
  row=0;
  for( col=0; col < spm->n; col++ )
    {
      for( row=spm->colptr[col]-baseval; row<spm->colptr[col+1]-baseval; row++)
	{
            coef[col] += spm->dofs[col] * spm->dofs[spm->rowptr[row]-baseval];
	}
    }
  return coef;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_internal
 *
 * dofCst - Expand a CSC matrix without dofs into a CSC matrix with constant dofs
 *
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm to expand
 *
 *******************************************************************************/
void
dofCst(pastix_spm_t* spm,pastix_int_t dof)
{
    int i;
    pastix_int_t* dofs = malloc((spm->n+1)*sizeof(pastix_int_t));

    spm->dof      = dof;
    spm->gNexp    = spm->n*dof;
    spm->nexp     = spm->n*dof;
    spm->nnzexp   = spm->n*dof*dof;
    spm->gnnzexp  = spm->n*dof*dof;
    spm->dofs     = dofs;

    for(i=0;i<spm->n+1;i++)
        dofs[i]=dof*i;

    spmExpand(spm);
}


/**
 *******************************************************************************
 *
 * @ingroup spm_internal
 *
 * dofVar - Expand a CSC matrix without dofs into a CSC matrix with variable dofs
 *
 *
 *******************************************************************************
 *
 * @param[in] spm
 *           The spm to expand
 *
 *******************************************************************************/
void
dofVar(pastix_spm_t* spm)
{
    pastix_int_t *dofs=malloc(sizeof(pastix_int_t)*spm->n);
    pastix_int_t i,nnzexp,nexp;

    spm->gNexp    = spm->gN;
    spm->nexp     = spm->n;
    spm->gnnzexp  = spm->gnnz;
    spm->nnzexp   = spm->nnz;
    spm->dof      = -1;

    genDof(3,dofs,spm->n); //cyclique, pas de 3
    //genRandDof(/*min*/1,/*max*/5 ,spm->dofs,spm->n); //random dofs, entre 1 et 5

    spm->dofs = dofs;
    pastix_int_t *coef = computeEltPerCol(spm); // free coef

    nnzexp = 0;
    nexp   = 0;

    for(i=0;i<spm->n;i++)
    {
        nnzexp += coef[i];
        nexp   += spm->dofs[i];
    }

    spm->nnzexp  = nnzexp;
    spm->gnnzexp = nnzexp;
    spm->gNexp   = nexp;
    spm->nexp    = nexp;

    pastix_int_t* tab = computeCorrespondIndex(dofs,spm->n);

    spm->dofs   = tab;
    spmExpand(spm);

    /*
    print_tab_int(spm->colptr,spm->n);
    print_tab_int(spm->rowptr,spm->nnz);
    print_tab_cpx(spm->values,spm->nnzexp);
    print_tab_int(spm->dofs,spm->n+1);
     */
}


