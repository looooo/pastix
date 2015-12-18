/**
 *  PaStiX CSC management routines.
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "z_csc_utils.h"

/* /\* */
/*    Function: cmp_colrow */

/*    Used for qsort to sort arrays of pastix_int_t following their first element. */

/*    Returns the difference between the first element of *p1* and */
/*    the first element of *p2* */

/*    Parameters: */
/*      p1 - the first array to compare */
/*      p2 - the second array to compare */
/* *\/ */
/* int */
/* cmp_colrow(const void *p1, const void *p2) */
/* { */
/*   return ((* (pastix_int_t * const *) p1) - (*(pastix_int_t * const *) p2)); */
/* } */

/*
  Function: z_csc_symgraph


  Modify the CSC to a symetric graph one.
  Don't use it on a lower symetric CSC
  it would give you all the CSC upper + lower.

  External function.

  Parameters:
    n     - Number of columns/vertices
    ia    - Starting index of each column in *ja* and *a*
    ja    - Row index of each element
    a     - Value of each element,can be NULL
    newn  - New number of column
    newia - Starting index of each column in *ja* and *a*
    newja - Row index of each element
    newa  - Value of each element,can be NULL

 */
int z_csc_symgraph( pastix_int_t               n,
                    const pastix_int_t        *ia,
                    const pastix_int_t        *ja,
                    const pastix_complex64_t  *a,
                    pastix_int_t              *newn,
                    pastix_int_t             **newia,
                    pastix_int_t             **newja,
                    pastix_complex64_t       **newa)
{
  return z_csc_symgraph_int(n, ia, ja, a, newn, newia, newja, newa, API_NO);
}

/*
  Function: csc_symgraph_int


  Modify the CSC to a symetric graph one.
  Don't use it on a lower symetric CSC
  it would give you all the CSC upper + lower.

  Parameters:
    n           - Number of columns/vertices
    ia          - Starting index of each column in *ja* and *a*
    ja          - Row index of each element
    a           - Value of each element,can be NULL
    newn        - New number of column
    newia       - Starting index of each column in *ja* and *a*
    newja       - Row index of each element
    newa        - Value of each element,can be NULL
    malloc_flag - flag to indicate if function call is intern to pastix or extern.
 */
int z_csc_symgraph_int (pastix_int_t               n,
                        const pastix_int_t        *ia,
                        const pastix_int_t        *ja,
                        const pastix_complex64_t  *a,
                        pastix_int_t              *newn,
                        pastix_int_t             **newia,
                        pastix_int_t             **newja,
                        pastix_complex64_t       **newa,
                        int                        malloc_flag)
{
  pastix_int_t * nbrEltCol = NULL; /* nbrEltCol[i] = Number of elt to add in column i */
  pastix_int_t * cia       = NULL; /* ia of diff between good CSC and bad CSC */
  pastix_int_t * cja       = NULL; /* ja of diff between good CSC and bad CSC */
  pastix_int_t   nbr2add;          /* Number of elt to add */
  pastix_int_t   itercol, iterrow, iterrow2; /* iterators */
  pastix_int_t   l = ia[n] -1;
  pastix_int_t   newl;

  /* Ncol=Nrow don't need change */
  *newn = n;

  MALLOC_INTERN(nbrEltCol, n, pastix_int_t);
  /* !! Need check for malloc */

  /* Init nbrEltCol */
  for (itercol=0; itercol<n; itercol++)
    {
      nbrEltCol[itercol]=0;
    }

  /* Compute number of element by col to add for correcting the CSC */
  for (itercol=0; itercol<n; itercol++)
    {
      for (iterrow=ia[itercol]-1; iterrow<ia[itercol+1]-1; iterrow++)
        {
          if (ja[iterrow] != (itercol+1))
            {
              /* Not diagonal elt */
              /* So we have a (i,j) and we are looking for a (j,i) elt */
              /* i = itercol+1, j=ja[iterrow] */
              int rowidx=ja[iterrow]-1;
              int flag=0;

              for (iterrow2=ia[rowidx]-1; iterrow2<ia[rowidx+1]-1; iterrow2++)
                {
                  if (ja[iterrow2] == itercol+1)
                    {
                      /* Ok we found (j,i) so stop this madness */
                      flag = 1;
                      break;
                    }
                }

              if (flag==0)
                {
                  /* We never find (j,i) so increase nbrEltCol[j] */
                  (nbrEltCol[ja[iterrow]-1])++;
                }
            }
        }
    }

  /* Compute number of element to add */
  /* And cia=ia part of csc of element to add */
  /* kind of a diff between the corrected one and the original CSC */
  MALLOC_INTERN(cia, n+1, pastix_int_t);
  /* !! Need checking good alloc) */
  nbr2add=0;
  for (itercol=0;itercol<n;itercol++)
    {
      cia[itercol]=nbr2add;
      nbr2add += nbrEltCol[itercol];
    }
  cia[n]=nbr2add;
  /*fprintf(stderr, "nbr of elt to add %ld\n", nbr2add);*/

  if (nbr2add != 0)
    {
      /* Build cja */
      /* like cia, cja is ja part of diff CSC */
      MALLOC_INTERN(cja, nbr2add, pastix_int_t);
      /* !! again we need check of memAlloc */

      /* We walkthrough again the csc */
      for (itercol=0;itercol<n;itercol++)
        {
          for (iterrow=ia[itercol]-1;iterrow<ia[itercol+1]-1;iterrow++)
            {
              if (ja[iterrow] != itercol+1)
                {
                  /* we find (i,j) need to find (j,i) */
                  int rowidx=ja[iterrow]-1;
                  int flag=0;

                  for (iterrow2=ia[rowidx]-1;iterrow2<ia[rowidx+1]-1;iterrow2++)
                    {
                      if (ja[iterrow2] == itercol+1)
                        {
                          /* find (j,i) */
                          flag=1;
                          break;
                        }
                    }

                  if (flag==0)
                    {
                      /* We don't find (j,i) so put in diff CSC (cia,cja,0) */
                      pastix_int_t index=ja[iterrow]-1;
                      /* cia[index] = index to put in cja the elt */
                      cja[cia[index]] = itercol+1;
                      (cia[index])++;
                    }
                }
            }
        }

      /* Restore cia */
      cia[0]=0;
      for (itercol=0;itercol<n;itercol++)
        {
          cia[itercol+1]=cia[itercol]+nbrEltCol[itercol];
        }

      memFree_null(nbrEltCol);

      /* Build corrected csc */
      newl = l+nbr2add;

      if (malloc_flag == API_NO)
        {
          /* Ici on a des malloc car le free est externe */
          MALLOC_EXTERN(*newia, n+1, pastix_int_t);
          MALLOC_EXTERN(*newja, newl, pastix_int_t);
          if (a != NULL)
            MALLOC_EXTERN(*newa, newl, pastix_complex64_t);
        }
      else
        {
          /* Ici on a des memAlloc car le free est interne */
          MALLOC_INTERN(*newia, n+1, pastix_int_t);
          MALLOC_INTERN(*newja, newl, pastix_int_t);
          if (a != NULL)
            MALLOC_INTERN(*newa, newl, pastix_complex64_t);
        }
      iterrow2 = 0; /* iterator of the CSC diff */
      for (itercol=0; itercol<n; itercol++)
        {
          (*newia)[itercol] = ia[itercol]+iterrow2;
          for (iterrow=ia[itercol]-1;iterrow<ia[itercol+1]-1;iterrow++)
            {
              /* we add the new elt with respect of order in row */
              while ((iterrow2<cia[itercol+1]) &&
                     (ja[iterrow] > cja[iterrow2]))
                {
                  /* we have elt(s) to add with a row lower than ja[iterrow] */
                  (*newja)[iterrow+iterrow2]=cja[iterrow2];
                  if (a != NULL)
                    (*newa)[iterrow+iterrow2]=0.;
                  iterrow2++;
                }

              /* Put the elt from the origin CSC */
              (*newja)[iterrow+iterrow2] = ja[iterrow];
              if (a != NULL)
                (*newa)[iterrow+iterrow2] = a[iterrow];
            }

          /* Since we put elt with a row lower than elt in the origin CSC */
          /* We could have some elt to add after the last elt in the column */
          while(iterrow2<cia[itercol+1])
            {
              (*newja)[iterrow+iterrow2]=cja[iterrow2];
              if (a != NULL)
                (*newa)[iterrow+iterrow2]=0.;
              iterrow2++;
            }
        }

      (*newia)[n]=ia[n]+iterrow2;
      memFree_null(cja);

    }
  else
    {
      /* No correction to do */
      memFree_null(nbrEltCol);
      newl = l;
      if (malloc_flag == API_NO)
        {
          /* ici on a des mallocs car le free est externe */
          MALLOC_EXTERN(*newia, n+1, pastix_int_t);
          MALLOC_EXTERN(*newja, l, pastix_int_t);
          if (a != NULL)
            MALLOC_EXTERN(*newa, l, pastix_complex64_t);
        }
      else
        {
          /* ici on a des memAllocs car le free est interne */
          MALLOC_INTERN(*newia, n+1, pastix_int_t);
          MALLOC_INTERN(*newja, l,   pastix_int_t);
          if (a != NULL)
            MALLOC_INTERN(*newa, l, pastix_complex64_t);
        }
      memcpy((*newia), ia, (n+1)*sizeof(pastix_int_t));
      memcpy((*newja), ja, l * sizeof(pastix_int_t));
      if (a != NULL)
        memcpy((*newa) , a , l * sizeof(pastix_complex64_t));
    }
  memFree_null(cia);

  return EXIT_SUCCESS;
}


/**
    Function: z_csc_noDiag

    Supress diagonal term.
    After this call, *ja* can be reallocated to *ia[n] -1*.

    Parameters:
      n  - size of the matrix.
      ia - Index in *ja* and *a* of the first element of each column
      ja - row of each element
      a  - value of each element, can be set to NULL

    Returns:
      ia and ja tabulars modified.
*/
void z_csc_noDiag(pastix_int_t        baseval,
                  pastix_int_t        n,
                  pastix_int_t       *ia,
                  pastix_int_t       *ja,
                  pastix_complex64_t *a)
{
  pastix_int_t i, j;
  pastix_int_t indj;
  pastix_int_t *old_ia = NULL;

  MALLOC_INTERN(old_ia, n+1, pastix_int_t);
  memcpy(old_ia, ia, sizeof(pastix_int_t)*(n+1));

  ASSERT(ia[0]==baseval,MOD_SOPALIN);

  indj = 0;
  /*fprintf(stdout, "NNZ with diag = %ld \n", ia[n]);*/

  for(i=0;i<n;i++)
    {
      /* ia[i] = number of column already counted */
      ia[i] = indj+baseval;
      /* for each row number in each column i */
      for(j=old_ia[i];j<old_ia[i+1];j++)
        /* if element is not diagonal
           we add it in ja and we count it */
        if(ja[j-baseval] != i+baseval)
          {
            ja[indj] = ja[j-baseval];
            if (a != NULL)
              a[indj] = a [j -baseval];
            indj++;
          }
    }
  ia[n] = indj+baseval;

  assert( ia[n] <= old_ia[n] );

  /*fprintf(stdout, "NNZ without diag = %ld \n", ia[n]);*/
  memFree_null(old_ia);
}

/*
  Function: z_csc_check_doubles

  Check if the csc contains doubles and if correct if asked

  Assumes that the CSC is sorted.

  Assumes that the CSC is Fortran numeroted (base 1)

  Parameters:
    n         - Size of the matrix.
    colptr    - Index in *rows* and *values* of the first element of each column
    rows      - row of each element
    values    - value of each element (Can be NULL)
    dof       - Number of degrees of freedom
    flag      - Indicate if user wants correction (<API_BOOLEAN>)
    flagalloc - indicate if allocation on CSC uses internal malloc.

  Returns:
    API_YES - If the matrix contained no double or was successfully corrected.
    API_NO  - Otherwise.
*/
int z_csc_check_doubles(pastix_int_t         n,
                        pastix_int_t        *colptr,
                        pastix_int_t       **rowptr,
                        pastix_complex64_t **values,
                        int                  dof,
                        int                  flag,
                        int                  flagalloc)
{
  pastix_int_t     i,j,k,d;
  int     doubles = 0;
  pastix_int_t   * tmprows = NULL;
  pastix_complex64_t * tmpvals = NULL;
  pastix_int_t     index = 0;
  pastix_int_t     lastindex = 0;

  ASSERT(values == NULL || dof > 0, MOD_SOPALIN);
  ASSERT(flag == API_NO || flag == API_YES, MOD_SOPALIN);
  ASSERT(colptr[0] == 1, MOD_SOPALIN);
  ASSERT(n >= 0, MOD_SOPALIN);

  for (i = 0; i < n; i++)
    {
      for (j = colptr[i]-1; j < colptr[i+1]-1; j = k)
        {
          (*rowptr)[index]   = (*rowptr)[j];
          if (values != NULL)
            for (d = 0; d < dof*dof; d++)
              (*values)[index*dof*dof+d] = (*values)[j*dof*dof+d];

          k = j+1;
          while (k < colptr[i+1]-1 && (*rowptr)[j] == (*rowptr)[k])
            {
              if (flag == API_NO)
                return API_NO;
              if (values != NULL)
                for (d = 0; d < dof*dof; d++)
                  (*values)[index*dof*dof+d] += (*values)[k*dof*dof+d];
              doubles++;
              k++;
            }
          index++;
        }

      colptr[i] = lastindex+1;
      lastindex = index;
    }
  if (flag == API_NO)
    return API_YES;
  ASSERT(index == colptr[n]-1-doubles, MOD_SOPALIN);
  if (doubles > 0)
    {
      colptr[n] = lastindex+1;
      if (flagalloc == API_NO)
        {
          MALLOC_EXTERN(tmprows, lastindex, pastix_int_t);
          if (values != NULL)
            MALLOC_EXTERN(tmpvals, lastindex*dof*dof, pastix_complex64_t);
        }
      else
        {
          MALLOC_INTERN(tmprows, lastindex, pastix_int_t);
          if (values != NULL)
            MALLOC_INTERN(tmpvals, lastindex*dof*dof, pastix_complex64_t);
        }

      memcpy(tmprows, *rowptr,   lastindex*sizeof(pastix_int_t));
      if (values != NULL)
        memcpy(tmpvals, *values, lastindex*dof*dof*sizeof(pastix_complex64_t));
      if (flagalloc == API_NO)
        {
          free(*rowptr);
          if (values != NULL)
            free(*values);
        }
      else
        {
          memFree_null(*rowptr);
          if (values != NULL)
            memFree_null(*values);
        }
      *rowptr   = tmprows;
      if (values != NULL)
        *values = tmpvals;
    }
  return API_YES;

}

/*
  Function: z_csc_checksym

    Check if the CSC graph is symetric.

    For all local column C,

    For all row R in the column C,

    We look in column R if we have the row number C.

    If we can correct we had missing non zeros.

    Assumes that the CSC is Fortran numbered (1 based).

    Assumes that the matrix is sorted.

  Parameters:
    n        - Number of local columns
    colptr   - Starting index of each columns in *ja*
    rows     - Row of each element.
    values   - Value of each element.
    correct  - Flag indicating if we can correct the symmetry.
    alloc    - indicate if allocation on CSC uses internal malloc.
    dof      - Number of degrees of freedom.
*/
int z_csc_checksym(pastix_int_t         n,
                   pastix_int_t        *colptr,
                   pastix_int_t       **rows,
                   pastix_complex64_t **values,
                   int                  correct,
                   int                  alloc,
                   int                  dof)
{
  pastix_int_t            i,j,k,l,d;
  pastix_int_t            index1;
  pastix_int_t            index2;
  int            found;
  pastix_int_t            toaddsize;
  pastix_int_t            toaddcnt;
  pastix_int_t         *  toadd      = NULL;
  pastix_int_t         *  tmpcolptr  = NULL;
  pastix_int_t         *  tmprows    = NULL;
  pastix_complex64_t       *  tmpvals    = NULL;

  /* For all local column C,
     For all row R in the column C,

     If the row number R correspond to a local column,
     We look in column R if we have the row number C.

     Else,
   */
  toaddcnt  = 0;
  toaddsize = 0;
  for (i = 0; i < n; i++)
    {
      for (j = (colptr)[i]-1; j < (colptr)[i+1]-1; j++)
        {
          if ((*rows)[j] != i+1)
            {
              /* not in diagonal */
              k = (*rows)[j];
              found = 0;
              for (l = (colptr)[k-1]-1; l < (colptr)[k-1+1]-1; l++)
                {
                  if (i+1 == (*rows)[l])
                    {
                      found = 1;
                      break;
                    }
                  if (i+1 < (*rows)[l])
                    {
                      /* The CSC is sorted */
                      found = 0;
                      break;
                    }
                }
              if (found == 0)
                {
                  if (correct == API_NO)
                    return EXIT_FAILURE;
                  else
                    {
                      if (toaddsize == 0)
                        {
                          toaddsize = n/2;
                          MALLOC_INTERN(toadd, 2*toaddsize, pastix_int_t);
                        }
                      if (toaddcnt >= toaddsize)
                        {
                          toaddsize += toaddsize/2 + 1;
                          if (NULL ==
                              (toadd =
                               (pastix_int_t*)memRealloc(toadd,
                                                2*toaddsize*sizeof(pastix_int_t))))
                              MALLOC_ERROR("toadd");
                        }
                      toadd[2*toaddcnt]     = (*rows)[j];
                      toadd[2*toaddcnt + 1] = i+1;
                      /* fprintf(stdout, "Adding %ld, %ld\n", (long)(i+1), (long)(*rows)[j]); */
                      toaddcnt++;
                    }
                }
            }
        }
    }

  if (toaddcnt > 0)
    {

      intSort2asc1(toadd, toaddcnt);
      /* Correct here is API_YES, otherwise we would have return EXIT_FAILURE
         Or toaddcnt == 0*/
      MALLOC_INTERN(tmpcolptr, n + 1, pastix_int_t);
      if (alloc == API_NO)
        {
          MALLOC_EXTERN(tmprows, colptr[n]-1 + toaddcnt, pastix_int_t);
          if (values != NULL)
            {
              MALLOC_EXTERN(tmpvals, colptr[n]-1 + toaddcnt, pastix_complex64_t);
            }
        }
      else
        {
          MALLOC_INTERN(tmprows, colptr[n]-1 + toaddcnt, pastix_int_t);
          if (values != NULL)
            {
              MALLOC_INTERN(tmpvals, colptr[n]-1 + toaddcnt, pastix_complex64_t);
            }
        }
      /* Build tmpcolptr

         tmpcolptr[i+1] will contain the number of element of
         the column i
       */
      index1 = 0;
      index2 = 0;
      for (i = 0; i <  n; i++)
        {
          tmpcolptr[i] = index2+1;
          for (j = colptr[i]-1; j < colptr[i+1]-1; j++)
            {
              if (index1 < toaddcnt &&
                  (toadd[2*index1] == i+1) &&
                  (toadd[2*index1+1] < (*rows)[j]))
                {
                  tmprows[index2] = toadd[2*index1+1];
                  if (values != NULL)
                    {
                      for (d = 0; d < dof*dof ; d++)
                        tmpvals[index2*dof*dof+d] = 0.0;
                    }
                  index1++;
                  j--; /* hack do not increment j this step of the loop */
                }
              else
                {
                  tmprows[index2] = (*rows)[j];
                  if (values != NULL)
                    {
                      for (d = 0; d < dof*dof ; d++)
                        tmpvals[index2*dof*dof+d] = (*values)[j*dof*dof+d];
                    }
                }
              index2++;
            }

          while(index1 < toaddcnt && toadd[2*index1] == i+1)
            {
              tmprows[index2] = toadd[2*index1+1];
              if (values != NULL)
                {
                  for (d = 0; d < dof*dof ; d++)
                    tmpvals[index2*dof*dof+d] = 0.0;
                }
              index1++;
              index2++;
            }
        }
      tmpcolptr[n] = index2+1;
      ASSERT((tmpcolptr[n] - 1) == (colptr[n] - 1 + toaddcnt), MOD_SOPALIN);

      memcpy(colptr, tmpcolptr, (n+1)*sizeof(pastix_int_t));
      memFree_null(tmpcolptr);
      memFree_null(toadd);
      if (alloc == API_NO)
        {
          free(*rows);
          if (values != NULL)
            {
              free(*values);
            }
        }
      else
        {
          memFree_null(*rows);
          if (values != NULL)
            {
              memFree_null(*values);
            }
        }
      *rows   = tmprows;
      if (values != NULL)
        {
          *values = tmpvals;
        }
    }
  return EXIT_SUCCESS;
}


/*
  Function: z_csc_colPerm

  Performs column permutation on a CSC

  Parameters:
    n     - Size of the matrix.
    ia    - Index of first element of each column in *ia* and *a*
    ja    - Rows of non zeros of the matrix.
    a     - Values of non zeros of the matrix.
    cperm - Permutation to perform
*/
void z_csc_colPerm(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_complex64_t *a, pastix_int_t *cperm)
{
  pastix_int_t i, k;
  pastix_int_t   *newja = NULL;
  pastix_int_t   *newia = NULL;
  pastix_complex64_t *newa  = NULL;
  int numflag, numflag2;

  numflag = ia[0];
  numflag2 = 1;
  for(i=0;i<n;i++)
    if(cperm[i] == 0)
      {
        numflag2 = 0;
        break;
      }

  if(numflag2 != numflag)
    {
      errorPrint("CSC_colPerm: rperm not in same numbering than the CSC.");
      exit(-1);
    }


  if(numflag == 1)
    {
      z_csc_Fnum2Cnum(ja, ia, n);
      for(i=0;i<n;i++)
        cperm[i]--;
    }

  MALLOC_INTERN(newia, n+1,   pastix_int_t);
  MALLOC_INTERN(newja, ia[n], pastix_int_t);
  MALLOC_INTERN(newa,  ia[n], pastix_complex64_t);


  newia[0] = 0;
  for(i=0;i<n;i++)
    {
#ifdef DEBUG_KASS
      ASSERT(cperm[i]>=0 && cperm[i] < n, MOD_KASS);
#endif
      newia[cperm[i]+1] = ia[i+1]-ia[i];
    }

#ifdef DEBUG_KASS
  for(i=1;i<=n;i++)
    ASSERT(newia[i] >0, MOD_KASS);
#endif

  for(i=1;i<=n;i++)
    newia[i] += newia[i-1];

#ifdef DEBUG_KASS
  ASSERT(newia[n] == ia[n], MOD_KASS);
#endif


  for(i=0;i<n;i++)
    {
      k = cperm[i];
#ifdef DEBUG_KASS
      ASSERT(newia[k+1]-newia[k] == ia[i+1]-ia[i], MOD_KASS);
#endif
      memcpy(newja + newia[k], ja + ia[i], sizeof(pastix_int_t)*(ia[i+1]-ia[i]));
      memcpy(newa + newia[k], a + ia[i], sizeof(pastix_complex64_t)*(ia[i+1]-ia[i]));

    }

  memcpy(ia, newia, sizeof(pastix_int_t)*(n+1));
  memcpy(ja, newja, sizeof(pastix_int_t)*ia[n]);
  memcpy(a, newa, sizeof(pastix_complex64_t)*ia[n]);

  memFree(newia);
  memFree(newja);
  memFree(newa);
  if(numflag == 1)
    {
      z_csc_Cnum2Fnum(ja, ia, n);
      for(i=0;i<n;i++)
        cperm[i]++;
    }
}


/*
  Function: z_csc_colScale

  Moved from kass, only used in MC64
*/
void z_csc_colScale(pastix_int_t        n,
                    pastix_int_t       *ia,
                    pastix_int_t       *ja,
                    pastix_complex64_t *a,
                    pastix_complex64_t *dcol)
{
  pastix_int_t i, j;
  int numflag;
  pastix_complex64_t d;
  numflag = ia[0];

  if(numflag == 1)
    z_csc_Fnum2Cnum(ja, ia, n);

  for(i=0;i<n;i++)
    {
      d = dcol[i];
      for(j=ia[i];j<ia[i+1];j++)
        {
          /***@@@ OIMBE DSCAL **/
          a[j] *= d;
        }
    }

  if(numflag == 1)
    z_csc_Cnum2Fnum(ja, ia, n);
}

/*
  Function: z_csc_rowScale

  Moved from kass, only used in MC64
*/
void z_csc_rowScale(pastix_int_t        n,
                    pastix_int_t       *ia,
                    pastix_int_t       *ja,
                    pastix_complex64_t *a,
                    pastix_complex64_t *drow)
{
  pastix_int_t i, j;
  int numflag;
  numflag = ia[0];

  if(numflag == 1)
    z_csc_Fnum2Cnum(ja, ia, n);

  for(i=0;i<n;i++)
    {
      for(j=ia[i];j<ia[i+1];j++)
        {
#ifdef DEBUG_KASS
          ASSERT(ja[j]>0 && ja[j] <n, MOD_KASS);
#endif
          a[j] *= drow[ja[j]];
        }
    }

  if(numflag == 1)
    z_csc_Cnum2Fnum(ja, ia, n);
}

/*
 * z_csc_sort:
 *
 * Sort CSC columns
 *
 * Parameters:
 *   n  - Number of columns
 *   ia - Index of first element of each column in *ia*.
 *   ja - Rows of each non zeros.
 *   a  - Values of each non zeros.
*/
void z_csc_sort(pastix_int_t        n,
                pastix_int_t       *ia,
                pastix_int_t       *ja,
                pastix_complex64_t *a,
                pastix_int_t        ndof)
{
  pastix_int_t i;
  int numflag;
  pastix_int_t ndof2 = ndof * ndof;
  void * sortptr[3];
  numflag = ia[0];
  if(numflag == 1)
    z_csc_Fnum2Cnum(ja, ia, n);
  if (a != NULL)
    {

      for(i=0;i<n;i++)
        {
          sortptr[0] = &ja[ia[i]];
          sortptr[1] = &a[ia[i]*ndof2];
          sortptr[2] = &ndof2;
          z_qsortIntFloatAsc(sortptr, ia[i+1] - ia[i]);
        }

    }
  else
    {
      for(i=0;i<n;i++)
        intSort1asc1(&ja[ia[i]], ia[i+1] - ia[i]);

    }
  if(numflag == 1)
    z_csc_Cnum2Fnum(ja, ia, n);
}

/*
  Function: z_csc_Fnum2Cnum

  Convert CSC numbering from fortran numbering to C numbering.

  Parameters:
    ja - Rows of each element.
    ia - First index of each column in *ja*
    n  - Number of columns
*/
void z_csc_Fnum2Cnum(pastix_int_t *ja,
                     pastix_int_t *ia,
                     pastix_int_t  n)
{
  pastix_int_t i, j;
  for(i=0;i<=n;i++)
    ia[i]--;

  for(i=0;i<n;i++)
    for(j=ia[i];j<ia[i+1];j++)
      ja[j]--;

}

/*
  Function: z_csc_Cnum2Fnum

  Convert CSC numbering from C numbering to Fortran numbering.

  Parameters:
    ja - Rows of each element.
    ia - First index of each column in *ja*
    n  - Number of columns
*/
void z_csc_Cnum2Fnum(pastix_int_t *ja,
                     pastix_int_t *ia,
                     pastix_int_t  n)
{
  pastix_int_t i, j;

  for(i=0;i<n;i++)
    for(j=ia[i];j<ia[i+1];j++)
      ja[j]++;

  for(i=0;i<=n;i++)
    ia[i]++;
}
/*
  Function: z_csc_buildZerosAndNonZerosGraphs

  Separate a graph in two graphs, following
  wether the diagonal term of a column is null or not.

  Parameters:
    n, colptr, rows, values  - The initial CSC
    n_nz, colptr_nz, rows_nz - The graph of the non-null diagonal part.
    n_z, colptr_z, rows_z    - The graph of the null diagonal part.
    perm                     - Permutation to go from the first graph to
                               the one composed of the two graph concatenated.
    revperm                  - Reverse permutation tabular.
    criteria                 - Value beside which a number is said null.
*/
int z_csc_buildZerosAndNonZerosGraphs(pastix_int_t         n,
                                      pastix_int_t        *colptr,
                                      pastix_int_t        *rows,
                                      pastix_complex64_t  *values,
                                      pastix_int_t        *n_nz,
                                      pastix_int_t       **colptr_nz,
                                      pastix_int_t       **rows_nz,
                                      pastix_int_t        *n_z,
                                      pastix_int_t       **colptr_z,
                                      pastix_int_t       **rows_z,
                                      pastix_int_t        *perm,
                                      pastix_int_t        *revperm,
                                      double               criteria)
{
  pastix_int_t  itercol;
  pastix_int_t  iterrow;

  pastix_int_t  ncoefszeros  = 0;
  pastix_int_t  ncoefsnzeros = 0;
  pastix_int_t  itercol_nz   = 0;
  pastix_int_t  itercol_z    = 0;
  int  seen;
  pastix_int_t  cntrows;

  for (itercol = 0; itercol <n; itercol++)
    {
      seen = 0;
      for (iterrow = colptr[itercol]-1; iterrow < colptr[itercol+1]-1; iterrow++)
        {
          if (itercol == rows[iterrow] -1 )
            {
              if (ABS_FLOAT(values[iterrow]) < criteria)
                {
                  (*n_z) ++;
                  ncoefszeros += colptr[itercol+1] - colptr[itercol];
                  seen = 1;
                }
              else
                {
                  (*n_nz) ++;
                  ncoefsnzeros += colptr[itercol+1] - colptr[itercol];
                  seen = 1;
                }
              break;
            }
        }
      if (colptr[itercol+1] == colptr[itercol]) /* empty column */
        {
          (*n_z)++;
          seen = 1;
        }
      if (seen == 0)/*  column without diag*/
        {
          (*n_z)++;
          ncoefszeros += colptr[itercol+1] - colptr[itercol];
        }
    }
  fprintf(stdout, "n_z %ld\n",  (long)*n_z);
  fprintf(stdout, "n_nz %ld\n", (long)*n_nz);
  ASSERT(*n_z+*n_nz == n, MOD_SOPALIN);
  if (*n_z == 0 || *n_nz == 0)
    return PASTIX_SUCCESS;

  for (itercol = 0; itercol <n; itercol++)
    {
      seen = 0;
      for (iterrow = colptr[itercol]-1; iterrow < colptr[itercol+1]-1; iterrow++)
        {
          if (itercol == rows[iterrow] -1 )
            {
              if (ABS_FLOAT(values[iterrow]) < criteria)
                {
                  perm[itercol] = (*n_nz) + itercol_z + 1;
                  itercol_z++;
                  seen = 1;
                }
              else
                {
                  perm[itercol] = itercol_nz + 1;
                  itercol_nz++;
                  seen = 1;
                }
            }
        }
      if (colptr[itercol] == colptr[itercol+1])
        { /* empty column */
          perm[itercol] = (*n_nz) + itercol_z + 1;
          itercol_z++;
          seen =1;
        }
      if (seen == 0)/*  column without diag*/
        {
          perm[itercol] = (*n_nz) + itercol_z + 1;
          itercol_z++;
        }
    }

  ASSERT(itercol_nz == *n_nz, MOD_SOPALIN);
  ASSERT(itercol_z  == *n_z, MOD_SOPALIN);
  for(itercol = 0; itercol < n; itercol++)
    revperm[perm[itercol]-1] = itercol + 1;

  MALLOC_INTERN(*colptr_nz, *n_nz + 1, pastix_int_t);
  MALLOC_INTERN(*rows_nz, ncoefsnzeros, pastix_int_t);
  cntrows = 0;
  for (itercol = 0; itercol <*n_nz; itercol++)
    {
      (*colptr_nz)[itercol] = cntrows + 1;
      for (iterrow = colptr[revperm[itercol]-1]-1; iterrow < colptr[revperm[itercol]-1+1]-1; iterrow++)
        {
          if (perm[rows[iterrow]-1] - 1 < *n_nz )
            {
              (*rows_nz)[cntrows] = perm[rows[iterrow]-1];
              cntrows++;
            }

        }

    }
  (*colptr_nz)[*n_nz] = cntrows+1;

  MALLOC_INTERN(*colptr_z, *n_z + 1, pastix_int_t);
  MALLOC_INTERN(*rows_z, ncoefszeros, pastix_int_t);
  cntrows = 0;
  for (itercol = 0; itercol <*n_z; itercol++)
    {
      (*colptr_z)[itercol] = cntrows + 1;
      for (iterrow = colptr[revperm[itercol+*n_nz]-1]-1; iterrow < colptr[revperm[itercol+*n_nz]-1+1]-1; iterrow++)
        {
          if (perm[rows[iterrow]-1]  > *n_nz )
            {
              (*rows_z)[cntrows] = perm[rows[iterrow]-1] - *n_nz;
              cntrows++;
            }

        }

    }
  (*colptr_z)[*n_z] = cntrows+1;

  return PASTIX_SUCCESS;
}

/*
  Function: z_csc_isolate

  Isolate a list of unknowns at the end of the CSC.

  Parameters:
    n            - Number of columns.
    colptr       - Index of first element of each column in *ia*.
    rows         - Rows of each non zeros.
    n_isolate    - Number of unknow to isolate.
    isolate_list - List of unknown to isolate.
    perm         - permutation tabular.
    revperm      - reverse permutation tabular.
*/
int z_csc_isolate(pastix_int_t  n,
                  pastix_int_t *colptr,
                  pastix_int_t *rows,
                  pastix_int_t  n_isolate,
                  pastix_int_t *isolate_list,
                  pastix_int_t *perm,
                  pastix_int_t *revperm)
{
  pastix_int_t  itercol;
  pastix_int_t  iterrow;

  pastix_int_t  iter_isolate  = 0;
  pastix_int_t  iter_non_isolate  = 0;
  pastix_int_t *tmpcolptr = NULL;
  pastix_int_t *tmprows   = NULL;

  if (n_isolate == 0)
    {
      errorPrintW("No schur complement\n");
      return PASTIX_SUCCESS;
    }
  intSort1asc1(isolate_list, n_isolate);

  for (itercol = 0; itercol <n; itercol++)
    {
      if (iter_isolate < n_isolate &&
          itercol == isolate_list[iter_isolate]-1)
        {
          revperm[n-n_isolate+iter_isolate] = itercol;
          iter_isolate++;
        }
      else
        {
          revperm[iter_non_isolate] = itercol;
          iter_non_isolate++;
        }
    }
  ASSERT(iter_non_isolate == n - n_isolate, MOD_SOPALIN);
  ASSERT(iter_isolate == n_isolate, MOD_SOPALIN);

  MALLOC_INTERN(tmpcolptr, n - n_isolate + 1,       pastix_int_t);
  memset(tmpcolptr, 0, (n - n_isolate + 1)*sizeof(pastix_int_t));

  for(itercol = 0; itercol < n; itercol++)
    perm[revperm[itercol]] = itercol;

  for(itercol = 0; itercol < n; itercol++)
    {
      ASSERT(perm[itercol] < n, MOD_SOPALIN);
      ASSERT(perm[itercol] > -1, MOD_SOPALIN);
    }

  tmpcolptr[0] = 1;
  for (itercol = 0; itercol <n; itercol++)
    {
      if (perm[itercol] < n - n_isolate)
        {
          for (iterrow = colptr[itercol]-1; iterrow < colptr[itercol+1]-1; iterrow ++)
            {
              /* Count edges in each column of the new graph */
              if (perm[rows[iterrow]-1] < n-n_isolate)
                {
                  tmpcolptr[perm[itercol]+1]++;
                }
            }
        }
    }

  for (itercol = 0; itercol <n - n_isolate; itercol++)
    tmpcolptr[itercol+1] += tmpcolptr[itercol];

  MALLOC_INTERN(tmprows,   tmpcolptr[n- n_isolate]-1, pastix_int_t);
  for (itercol = 0; itercol <n; itercol++)
    {
      if (perm[itercol] < n - n_isolate)
        {
          for (iterrow = colptr[itercol]-1; iterrow < colptr[itercol+1]-1; iterrow ++)
            {
              /* Count edges in each column of the new graph */
              if (perm[rows[iterrow]-1] < n-n_isolate)
                {
                  tmprows[tmpcolptr[perm[itercol]]-1] = perm[rows[iterrow]-1]+1;
                  tmpcolptr[perm[itercol]]++;
                }
            }
        }
    }


  /* restore tmpcolptr */

  for (itercol = 1; itercol <n - n_isolate ; itercol++)
    {
      tmpcolptr[n - n_isolate - itercol] = tmpcolptr[n - n_isolate - itercol - 1];
    }
  tmpcolptr[0] = 1;



  ASSERT(colptr[n] >= tmpcolptr[n-n_isolate], MOD_SOPALIN);
  memcpy(colptr, tmpcolptr, (n-n_isolate + 1)*sizeof(pastix_int_t));
  memcpy(rows,   tmprows,   (colptr[n-n_isolate]-1)*sizeof(pastix_int_t));

  memFree_null(tmpcolptr);
  memFree_null(tmprows);

  return PASTIX_SUCCESS;
}

/*
  Function: csc_save

  Save a csc on disk.

  Parameters:
    n       - number of columns
    colptr  - First cscd starting index of each column in *ja* and *a*
    rows    - Row of each element in first CSCD
    values  - value of each cscd in first CSCD (can be NULL)
    dof     - Number of degrees of freedom
    outfile - Output stream.

  Return:
    PASTIX_SUCCESS

*/
int z_csc_save(pastix_int_t        n,
               pastix_int_t       *colptr,
               pastix_int_t       *rows,
               pastix_complex64_t *values,
               int                 dof,
               FILE               *outfile)
{

  pastix_int_t i;
  fprintf(outfile, "%ld %ld %d\n", (long)n, (long)dof, ((values == NULL)?0:1));
  /* Copie du colptr */
  for (i=0; i<n+1; i++)
    {
      fprintf(outfile, "%ld ", (long)colptr[i]);
      if (i%4 == 3) fprintf(outfile, "\n");
    }
  if ((i-1)%4 !=3) fprintf(outfile, "\n");

  /* Copie de JA */
  for (i=0; i<colptr[n]-1; i++)
    {
      fprintf(outfile, "%ld ", (long)rows[i]);
      if (i%4 == 3) fprintf(outfile, "\n");
    }
  if ((i-1)%4 !=3) fprintf(outfile, "\n");

  /* Copie de Avals */

  if (values != NULL)
    {
      for (i=0; i<(colptr[n]-1)*dof*dof; i++)
        {
#ifdef TYPE_COMPLEX
          fprintf(outfile, "%lg %lg ", (double)(creal(values[i])), (double)(cimag(values[i])));
#else
          fprintf(outfile, "%lg ", (double)(values[i]));
#endif
          if (i%4 == 3) fprintf(outfile, "\n");
        }
      if ((i-1)%4 !=3) fprintf(outfile, "\n");
    }
  return PASTIX_SUCCESS;

}

/*
  Function: csc_load

  Load a csc from disk.

  Fill *n*, *colptr*, *rows*, *values* and *dof* from *infile*.

  Parameters:
    n       - number of columns
    colptr  - First cscd starting index of each column in *ja* and *a*
    rows    - Row of each element in first CSCD
    values  - value of each cscd in first CSCD (can be NULL)
    dof     - Number of degrees of freedom
    outfile - Output stream.

  Return:
    PASTIX_SUCCESS

*/
int z_csc_load(pastix_int_t        *n,
               pastix_int_t       **colptr,
               pastix_int_t       **rows,
               pastix_complex64_t **values,
               int                 *dof,
               FILE                *infile)
{
  int  hasval;
  long tmp1, tmp2, tmp3, tmp4;
  double tmpflt1, tmpflt2, tmpflt3, tmpflt4;
#ifdef TYPE_COMPLEX
  double tmpflt5, tmpflt6, tmpflt7, tmpflt8;
#endif
  pastix_int_t i;
  if (3 != fscanf(infile, "%ld %ld %d\n", &tmp1, &tmp2, &hasval)){
    errorPrint("CSCD badly formated");
    return EXIT_FAILURE;
  }
  *n   = (pastix_int_t)tmp1;
  *dof = (int)tmp2;
  /* Copie de IA */
  *colptr = NULL;
  MALLOC_INTERN(*colptr, *n+1, pastix_int_t);
  for (i=0; i<*n+1+1-4; i+=4)
    {
      if (4 != fscanf(infile, "%ld %ld %ld %ld", &tmp1, &tmp2, &tmp3, &tmp4)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*colptr)[i  ] = tmp1;
      (*colptr)[i+1] = tmp2;
      (*colptr)[i+2] = tmp3;
      (*colptr)[i+3] = tmp4;
    }
  switch (*n +1 - i)
    {
    case 3:
      if (3 != fscanf(infile, "%ld %ld %ld", &tmp1, &tmp2, &tmp3)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*colptr)[i  ] = tmp1;
      (*colptr)[i+1] = tmp2;
      (*colptr)[i+2] = tmp3;
      break;
    case 2:
      if (2 != fscanf(infile, "%ld %ld", &tmp1, &tmp2)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*colptr)[i  ] = tmp1;
      (*colptr)[i+1] = tmp2;
      break;
    case 1:
      if (1 != fscanf(infile, "%ld", &tmp1)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*colptr)[i  ] = tmp1;
      break;
    }


  /* Copie de JA */
  (*rows) = NULL;
  MALLOC_INTERN(*rows, (*colptr)[*n]-(*colptr)[0], pastix_int_t);
  for (i=0; i< (*colptr)[*n]-(*colptr)[0]+1-4; i+=4)
    {
      if (4 != fscanf(infile, "%ld %ld %ld %ld", &tmp1, &tmp2, &tmp3, &tmp4)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*rows)[i  ] = tmp1;
      (*rows)[i+1] = tmp2;
      (*rows)[i+2] = tmp3;
      (*rows)[i+3] = tmp4;
    }

  switch ( (*colptr)[*n]-(*colptr)[0] - i)
    {
    case 3:
      if (3 != fscanf(infile, "%ld %ld %ld", &tmp1, &tmp2, &tmp3)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*rows)[i  ] = tmp1;
      (*rows)[i+1] = tmp2;
      (*rows)[i+2] = tmp3;
      break;
    case 2:
      if (2 != fscanf(infile, "%ld %ld", &tmp1, &tmp2)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*rows)[i  ] = tmp1;
      (*rows)[i+1] = tmp2;
      break;
    case 1:
      if (1 != fscanf(infile, "%ld", &tmp1)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*rows)[i  ] = tmp1;
      break;
    }

  /* Copie de Avals */
  if (hasval)
    {
      (*values) = NULL;

      MALLOC_INTERN(*values,  (*colptr)[*n]-(*colptr)[0], pastix_complex64_t);

      for (i=0; i< (*colptr)[*n]-(*colptr)[0]+1-4; i+=4)
        {
#ifdef TYPE_COMPLEX
          if (8 != fscanf(infile, "%lg %lg %lg %lg %lg %lg %lg %lg",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4,
                          &tmpflt5, &tmpflt6, &tmpflt7, &tmpflt8)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (pastix_complex64_t)(tmpflt1 + I * tmpflt2);
          (*values)[i+1] = (pastix_complex64_t)(tmpflt3 + I * tmpflt4);
          (*values)[i+2] = (pastix_complex64_t)(tmpflt5 + I * tmpflt6);
          (*values)[i+3] = (pastix_complex64_t)(tmpflt7 + I * tmpflt8);
#else
          if (4 != fscanf(infile, "%lg %lg %lg %lg",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (pastix_complex64_t)tmpflt1;
          (*values)[i+1] = (pastix_complex64_t)tmpflt2;
          (*values)[i+2] = (pastix_complex64_t)tmpflt3;
          (*values)[i+3] = (pastix_complex64_t)tmpflt4;
#endif
        }
      switch ( (*colptr)[*n]-(*colptr)[0] - i )
        {
        case 3:
#ifdef TYPE_COMPLEX
          if (6 != fscanf(infile, "%lg %lg %lg %lg %lg %lg",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4,
                          &tmpflt5, &tmpflt6)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (pastix_complex64_t)(tmpflt1 + I * tmpflt2);
          (*values)[i+1] = (pastix_complex64_t)(tmpflt3 + I * tmpflt4);
          (*values)[i+2] = (pastix_complex64_t)(tmpflt5 + I * tmpflt6);
#else
          if (3 != fscanf(infile, "%lg %lg %lg",
                          &tmpflt1, &tmpflt2, &tmpflt3)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (pastix_complex64_t)tmpflt1;
          (*values)[i+1] = (pastix_complex64_t)tmpflt2;
          (*values)[i+2] = (pastix_complex64_t)tmpflt3;
#endif
          break;
        case 2:
#ifdef TYPE_COMPLEX
          if (4 != fscanf(infile, "%lg %lg %lg %lg",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (pastix_complex64_t)(tmpflt1 + I * tmpflt2);
          (*values)[i+1] = (pastix_complex64_t)(tmpflt3 + I * tmpflt4);
#else
          if (2 != fscanf(infile, "%lg %lg",
                          &tmpflt1, &tmpflt2)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (pastix_complex64_t)tmpflt1;
          (*values)[i+1] = (pastix_complex64_t)tmpflt2;
#endif
          break;
        case 1:
#ifdef TYPE_COMPLEX
          if (2 != fscanf(infile, "%lg %lg",
                          &tmpflt1, &tmpflt2)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (pastix_complex64_t)(tmpflt1 + I * tmpflt2);
#else
          if (1 != fscanf(infile, "%lg",
                          &tmpflt1)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (pastix_complex64_t)tmpflt1;
#endif
          break;
        }
    }
  return PASTIX_SUCCESS;
}
