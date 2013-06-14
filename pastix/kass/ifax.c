#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "common.h"
#include "symbol.h"
#include "sparRow.h"
#include "sort_row.h"
#include "SF_level.h"
#include "ifax.h"


void ifax(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_int_t levelk, pastix_int_t  cblknbr, pastix_int_t *rangtab, pastix_int_t *perm, pastix_int_t *iperm, SymbolMatrix *symbmtx)
{
  pastix_int_t i, j, k;
  pastix_int_t ii, ind, cblk;
  pastix_int_t nnzL;
  pastix_int_t *cblkflag  = NULL;
  pastix_int_t *tmpj      = NULL;
  pastix_int_t *node2cblk = NULL;
  csptr qmat;
  csptr P;

  ASSERT(n == rangtab[cblknbr],MOD_BLEND);
  fprintf(stderr, "\n NUMBER OF NNZ in A %ld \n", (long)(ia[n]));
  fprintf(stderr, "\n NUMBER OF CBLK  in rangtab %ld \n", (long)cblknbr);
  

  MALLOC_INTERN(qmat, 1, struct SparRow);
  initCS(qmat, cblknbr);
  

  MALLOC_INTERN(cblkflag,  cblknbr, pastix_int_t);
  MALLOC_INTERN(node2cblk, n,       pastix_int_t);
  MALLOC_INTERN(tmpj,      cblknbr, pastix_int_t);
  
  /** Fill the node2cblk vector **/
  for(k=0;k<cblknbr;k++)
    for(i=rangtab[k];i<rangtab[k+1];i++)
      node2cblk[i] = k;



  /**** Build the quotient graph ****/
  /*** NOTE : We build the symmetrized matrix !!***/
  bzero(cblkflag, sizeof(pastix_int_t)*cblknbr);
  bzero(tmpj, sizeof(pastix_int_t)*cblknbr);
  for(k=0;k<cblknbr;k++)
    {
      cblkflag[k] = 1;
      tmpj[0] = k; /** The diagonal block is always present **/
                   /** With this the function still works with null diagonal term **/
      ind = 1;

      for(ii=rangtab[k];ii<rangtab[k+1];ii++)
	{
	  i = iperm[ii];
	  for(j=ia[i];j<ia[i+1];j++)
	    {
	      cblk = node2cblk[ perm[ja[j]] ];

	      ASSERT(cblk >= 0 && cblk < cblknbr,MOD_BLEND);

	      if(cblkflag[cblk] == 0)
		{
		  cblkflag[cblk] = 1;
		  tmpj[ind++] = cblk;
		}
	    }
	}

      for(ii=0;ii<ind;ii++)
	cblkflag[ tmpj[ii] ] = 0;


      qmat->nnzrow[k] = ind;
      MALLOC_INTERN(qmat->ja[k], ind, pastix_int_t);
      memcpy(qmat->ja[k], tmpj, sizeof(pastix_int_t)*ind);
      qmat->ma[k] = NULL;
    }


  /*** Reorder the matrix ***/
  sort_row(qmat);

  fprintf(stdout, "NUMBER of CBLK in A %ld \n", (long)cblknbr);
  fprintf(stdout, " NUMBER of BLOCK in A %ld \n", (long)CSnnz(qmat));

  /*** Compute the ILU(k) pattern of the quotient matrix ***/
  MALLOC_INTERN(P, 1, SparMat);
  initCS(P, cblknbr);
  
  fprintf(stderr, "Process the Symbolic block Factorization for a fill level of %ld ...\n", (long)levelk);
  /*nnzL = SF_level(0, qmat, levelk, P);
    nnzL = SF_level(1, qmat, levelk, P);*/

  nnzL = SF_level(2, qmat, levelk, P);

  /*nnzL = SF_GSurow(mat, levelk, P);*/
  /*fprintf(stdout, "nnzA %ld nnzP %ld ; RATIO = %g \n", (long)( CSnnz(mat) + n)/2, (long)CSnnz(P), (double)CSnnz(P)/(double)( (CSnnz(mat)+n)/2.0 ));*/
  fprintf(stderr, "Number of Block in L (not yet recompacted) is %ld (==%ld) \n", (long)nnzL, (long)(CSnnz(P)*2 - P->n));
  
  /** Sort the rows of the symbolic matrix */
  sort_row(P);


  /**********************************/
  /*** Compute the symbol matrix ****/
  /**********************************/
  symbmtx->baseval = 0;
  symbmtx->cblknbr = cblknbr;
  symbmtx->bloknbr = CSnnz(P);
  symbmtx->nodenbr = n;

  MALLOC_INTERN(symbmtx->cblktab, cblknbr+1,        SymbolCblk);
  MALLOC_INTERN(symbmtx->bloktab, symbmtx->bloknbr, SymbolBlok);
  
  ind = 0;
  for(k=0;k<cblknbr;k++)
    {
      symbmtx->cblktab[k].fcolnum = rangtab[k];
      symbmtx->cblktab[k].lcolnum = rangtab[k+1]-1;
      symbmtx->cblktab[k].bloknum = ind;
      for(i=0;i<P->nnzrow[k];i++)
	{
	  j = P->ja[k][i];

	  ASSERT(j >= k,MOD_BLEND);
	  ASSERT(j < cblknbr,MOD_BLEND);

	  symbmtx->bloktab[ind].frownum = rangtab[j];
	  symbmtx->bloktab[ind].lrownum = rangtab[j+1]-1;
	  symbmtx->bloktab[ind].cblknum = j;
	  symbmtx->bloktab[ind].levfval = 0;
	  ind++;
	}
    }
  /*  virtual cblk to avoid side effect in the loops on cblk bloks */
  symbmtx->cblktab[cblknbr].fcolnum = symbmtx->cblktab[cblknbr-1].lcolnum+1;
  symbmtx->cblktab[cblknbr].lcolnum = symbmtx->cblktab[cblknbr-1].lcolnum+1;
  symbmtx->cblktab[cblknbr].bloknum = ind;

  
  
  cleanCS(P);
  cleanCS(qmat);
  memFree_null(qmat);
  memFree_null(node2cblk);
  memFree_null(tmpj);
  memFree_null(cblkflag);

}
