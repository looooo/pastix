#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "common.h"
/* #include "symbol.h" */

#include "sparRow.h"


pastix_int_t initCS(csptr amat, pastix_int_t len)
{
/*---------------------------------------------------------------------- 
| Initialize SparRow structs.
|----------------------------------------------------------------------
| on entry: 
|========== 
| ( amat )  =  Pointer to a SparRow struct.
|     len   =  size of matrix
|
| On return:
|===========
|
|  amat->n
|      ->*nnzrow
|      ->**ja
|      ->**ma
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   amat->n = len;
   amat->inarow = 0;
   if(len > 0)
     {
       MALLOC_INTERN(amat->nnzrow, len, pastix_int_t);
       MALLOC_INTERN(amat->ja,     len, pastix_int_t *);
       MALLOC_INTERN(amat->ma,     len, double *);

       bzero( amat->nnzrow, sizeof(pastix_int_t)*len);
       bzero( amat->ja, sizeof(pastix_int_t *)*len);
       bzero( amat->ma, sizeof(double *)*len);
       amat->jatab = NULL;
       amat->matab = NULL;
     }
   else
     {
       amat->nnzrow = NULL;
       amat->ja = NULL;
       amat->ma = NULL;
       amat->jatab = NULL;
       amat->matab = NULL;
     }

   return 0;
}


pastix_int_t cleanCS(csptr amat)
{
/*----------------------------------------------------------------------
| Free up memory allocated for SparRow structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a SparRow struct.
|     len   =  size of matrix
|--------------------------------------------------------------------*/
   /*	*/
  pastix_int_t i; 
  if (amat == NULL) return 0;
  /*if (amat->n < 1) return 0;*/

  if(amat->inarow == 0)
    {
      if(amat->n > 0) /** The sparse part has not been deallocated **/
	{
	  for (i=0; i<amat->n; i++) {
              //	    if (amat->nnzrow[i] > 0) {
	      if (amat->ma[i]) memFree_null(amat->ma[i]);
	      if (amat->ja[i]) memFree_null(amat->ja[i]);
              //	    }
	  }	
	}
    }
  else
    {
      if(amat->inarow != 2 && amat->jatab != NULL)
	{
	  memFree_null(amat->matab);
	  memFree_null(amat->jatab);
	}
      amat->inarow = 0;

    }

  if (amat->ma) memFree_null(amat->ma);
  if (amat->ja) memFree_null(amat->ja);
  if (amat->nnzrow) memFree_null(amat->nnzrow); 
     

  return 0;
}


pastix_int_t CSnnz(csptr mat)
{
  /*-----------------------------------/
  / Return number of non zero entries  /
  / in a matrix in SparRow format      /
  /-----------------------------------*/
  pastix_int_t nnz;
  pastix_int_t i;

  nnz = 0;
  for(i=0;i<mat->n;i++)
    nnz += mat->nnzrow[i];
      
  return nnz;
}

pastix_int_t CS_RowPerm(csptr mat, pastix_int_t *perm)
{
/*----------------------------------------------------------------------
|
| This subroutine permutes the rows of a matrix in SparRow format. 
| rperm  computes B = P A  where P is a permutation matrix.  
| The permutation P is defined through the array perm: for each j, 
| perm[j] represents the destination row number of row number j. 
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (amat) = a matrix stored in SparRow format.
|
|
| on return:
| ----------
| (amat) = P A stored in SparRow format.
|
| INTeger value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/
  pastix_int_t **addj = NULL;
  pastix_int_t * nnz  = NULL;
  pastix_int_t i, size;
  double **addm = NULL;

   size = mat->n;
   if(size == 0)
     return 0;

   MALLOC_INTERN(addj, size, pastix_int_t *);
   MALLOC_INTERN(addm, size, double *);
   MALLOC_INTERN(nnz,  size, pastix_int_t);

   for (i=0; i<size; i++) {
      addj[perm[i]] = mat->ja[i];
      addm[perm[i]] = mat->ma[i];
      nnz[perm[i]] = mat->nnzrow[i];
   }
   for (i=0; i<size; i++) {
      mat->ja[i] = addj[i];
      mat->ma[i] = addm[i];
      mat->nnzrow[i] = nnz[i];
   }
   memFree(addj);
   memFree(addm);
   memFree(nnz);
   return 0;
}
/*------- end of rperm ------------------------------------------------- 
|---------------------------------------------------------------------*/

pastix_int_t CS_ColPerm(csptr mat, pastix_int_t *perm) 
{
/*----------------------------------------------------------------------
|
| This subroutine permutes the columns of a matrix in SparRow format.
| cperm computes B = A P, where P is a permutation matrix.
| that maps column j INTo column perm(j), i.e., on return 
| The permutation P is defined through the array perm: for each j, 
| perm[j] represents the destination column number of column number j. 
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (mat) = a matrix stored in SparRow format.
|
|
| on return:
| ----------
| (mat) = A P stored in SparRow format.
|
|---------------------------------------------------------------------*/
   pastix_int_t i, j,  size, *aja;
   pastix_int_t *newj = NULL;
   pastix_int_t maxrow;
   size = mat->n;
   if(size == 0)
     return 0;
   
   maxrow = 0;
   for(i=0;i<mat->n;i++)
     if(mat->nnzrow[i]>maxrow)
       maxrow = mat->nnzrow[i];
   
   if(maxrow == 0)
     return 0;

   MALLOC_INTERN(newj, maxrow, pastix_int_t);

   for (i=0; i<size; i++) {
      aja = mat->ja[i];
      for (j=0; j<mat->nnzrow[i]; j++)
	newj[j] = perm[aja[j]];
  
      for (j=0; j<mat->nnzrow[i]; j++)
	 aja[j] = newj[j];
      /*mat->ja[i] = aja;*/
   }
   memFree(newj);
   return 0;
}
/*------- end of cperm ------------------------------------------------- 
|---------------------------------------------------------------------*/
pastix_int_t CS_Perm(csptr mat, pastix_int_t *perm) 
{
/*----------------------------------------------------------------------
|
| This subroutine permutes the rows and columns of a matrix in 
| SparRow format.  dperm computes B = P^T A P, where P is a permutation 
| matrix.
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (amat) = a matrix stored in SparRow format.
|
|
| on return:
| ----------
| (amat) = P^T A P stored in SparRow format.
|
| INTeger value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/
   if (CS_RowPerm(mat, perm)) return 1;
   if (CS_ColPerm(mat, perm)) return 1;
   return 0;
}

