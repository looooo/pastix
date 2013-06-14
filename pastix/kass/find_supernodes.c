/************************************************************/
/**                                                        **/
/**   NAME       : find_supernodes.c                       **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 10/02/2006      **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "common.h"
#include "find_supernodes.h"

/** local function **/
void  post_order(pastix_int_t n, pastix_int_t *father, pastix_int_t *T,  pastix_int_t *perm, pastix_int_t *iperm);
void compute_subtree_size(pastix_int_t n, pastix_int_t *father, pastix_int_t *perm, pastix_int_t *iperm, pastix_int_t *T);
void compute_elimination_tree(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_int_t *perm, pastix_int_t *iperm, pastix_int_t *father);

void  find_supernodes(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_int_t *perm, pastix_int_t *iperm, pastix_int_t *snodenbr, pastix_int_t *snodetab, pastix_int_t *treetab)
{
  /********************************************************************************************/
  /* This function computes the supernodal partition of a reordered matrix.                   */
  /* The permutation of the matrix given on entry is modified to obtain a                     */
  /* postorder of the elimination tree: this does not affect the fill-in                      */
  /* properties of the initial ordering.                                                      */
  /* The matrix pattern is assumed to be symmetric                                            */
  /* NOTE :                                                                                   */
  /* This function can take on entry the lower triangular part or the whole matrix A:         */
  /* this does not change the results                                                         */
  /*                                                                                          */
  /* On entry:                                                                                */
  /* n, ia, ja : the initial matrix A in CSC format                                           */
  /*             in C numbering                                                               */
  /*  perm : the permutation vector                                                           */
  /* iperm : the inverse permutation vector                                                   */
  /* snodetab : allocated to contain at most n integers                                       */
  /* treetab  : allocated to contain at most n integers                                       */
  /*            it can be set to NULL in which case the return value is treetab = NULL        */
  /* On return :                                                                              */
  /* perm, iperm : postorder of the node in the elimination tree deduced from the initial     */
  /*               ordering                                                                   */
  /* snodenbr  : number of supernodes found                                                   */
  /* snodetab  : snodetab[i] is the beginning in the new ordering of the ith supernodes       */
  /* treetab   : treetab[s] is the number of the father of supernode i on the supernodal      */
  /*             elimination tree                                                             */
  /********************************************************************************************/

  pastix_int_t *father     = NULL; /** father[i] is the father of node i in he elimination tree of A **/
  pastix_int_t *T          = NULL; /** T[j] is the number of node in the subtree rooted in node j in
                              the elimination tree of A **/
  pastix_int_t *S          = NULL; /** S[i] is the number of sons for node i in the elimination tree **/
  pastix_int_t *isleaf     = NULL;
  pastix_int_t *prev_rownz = NULL;
  pastix_int_t i, j, k;
  pastix_int_t pi, pj;
  pastix_int_t dad;


  MALLOC_INTERN(T,          n, pastix_int_t);
  MALLOC_INTERN(S,          n, pastix_int_t);
  MALLOC_INTERN(father,     n, pastix_int_t);
  MALLOC_INTERN(isleaf,     n, pastix_int_t);
  MALLOC_INTERN(prev_rownz, n, pastix_int_t);
#ifdef DEBUG_BLEND
  assert(ia[0] == 0);
#endif


#ifdef DEBUG_BLEND
  /** Check the permutation vector **/
  for(i=0;i<n;i++)
    {
      assert(perm[i] >= 0);
      assert(perm[i] < n);
    }

  bzero(S, sizeof(pastix_int_t)*n);
  for(i=0;i<n;i++)
    S[perm[i]]++;

  k = 0;
  for(i=0;i<n;i++)
    if(S[i] != 1)
      k++;
  if(k>0)
    errorPrint("perm array is not valid, number of error =  %ld", (long)k);
  assert(k==0);
#endif


  /*** Compute the elimination tree of A ***/
  compute_elimination_tree(n, ia, ja, perm, iperm, father);


  /*** Compute the postorder of the elimination tree ***/
  /*** This operation modifies perm and iperm ***/
  post_order(n, father, T, perm, iperm);

  /*** Compute the number of descendant of each node i in the elimination tree ***/
  compute_subtree_size(n, father, perm, iperm, T);


  bzero(isleaf, sizeof(pastix_int_t)*n);
  bzero(prev_rownz, sizeof(pastix_int_t)*n);


  for(j=0;j<n;j++)
    {
      pj = iperm[j];
      for(i=ia[pj];i<ia[pj+1];i++)
        {
          pi = perm[ja[i]];
          if(pi > j)
            {
              k = prev_rownz[pi];
              if(k < j - T[pj]+1 )
                isleaf[j] = 1;

              prev_rownz[pi] = j;
            }
        }
    }


  /*** Compute the number of sons of each node in the elimination tree ***/
  bzero(S, sizeof(pastix_int_t)*n);
  for(i=0;i<n;i++)
    if(father[i] != i)
      S[father[i]]++;

  for(i=0;i<n;i++)
    if(S[i] != 1)
      isleaf[perm[i]] = 1;

  (*snodenbr) = 0;
  for(i=0;i<n;i++)
    if(isleaf[i] == 1)
      {
        snodetab[(*snodenbr)] = i;
        (*snodenbr)++;
      }
  snodetab[(*snodenbr)] = n;


  if(treetab != NULL)
    {
      /*** Node to supernodet conversion vector ***/
      for(i=0;i<(*snodenbr);i++)
        for(j=snodetab[i];j<snodetab[i+1];j++)
          S[j] = i;

      /*** Fill the treetab info ***/
      for(i=0;i<(*snodenbr);i++)
        {
          k=(*snodenbr);
          for(j=snodetab[i];j<snodetab[i+1];j++)
            {
              dad = S[perm[father[iperm[j]]]];
              if( dad < k && dad > i)
                k = dad;
            }
          treetab[i] = k;
          if(k==(*snodenbr))
            {
              /*fprintf(stdout, "THIS IS A ROOT %d \n", i);*/
              treetab[i] = i; /** This is a root **/
            }
#ifdef DEBUG_BLEND
          assert(treetab[i] >= i);
#endif
        }


    }

  memFree(prev_rownz);
  memFree(isleaf);
  memFree(father);
  memFree(S);
  memFree(T);






}


void  post_order(pastix_int_t n, pastix_int_t *father, pastix_int_t *T,  pastix_int_t *perm, pastix_int_t *iperm)
{
  /********************************************************************************/
  /* This function compute the post order of the elimination tree given on entry  */
  /* On entry:                                                                    */
  /* n : number of nodes                                                          */
  /* father : father[i] is the node number of the father of node i                */
  /*          if node i is a root then father[i] = i                              */
  /* T : a temporary vector of size n                                             */
  /* perm, iperm : ordering of the matrix (value optional: if set then the post   */
  /*                            ordering try to keep the initial ordering as much */
  /*                            as possible)                                      */
  /* On return :                                                                  */
  /* perm, iperm : permutation and inverse permutation vector of the postorder    */
  /********************************************************************************/
  pastix_int_t i;
  pastix_int_t j, k, t;


  /*** First compute the number of node in the subtree rooted in node i ***/
  compute_subtree_size(n, father, perm, iperm, T);

  /** When multiple roots are present we have to compute the start index of each root **/
  t=0;
  for(k=0;k<n;k++)
    {
      i = iperm[k];
      if(father[i] == i)
        {
          /** This is a root **/
          j = T[i];
          T[i] += t;
          t += j;
        }
    }

#ifdef DEBUG_BLEND
  for(i=0;i<n;i++)
    assert(T[i] <= T[father[i]]);
#endif



 for(k=n-1;k>=0;k--)
   {
     i = iperm[k];
     perm[i] = T[father[i]]; /** We MUST HAVE father[i] == i for a root ! **/
     T[father[i]]-= T[i];
     T[i] = perm[i]-1;
#ifdef DEBUG_BLEND
     assert(perm[father[i]] >= perm[i]);
#endif
   }



  /** We need to retrieve 1 for the C numbering compliance **/
  for(i=0;i<n;i++)
    perm[i]--;

#ifdef DEBUG_BLEND
  /** Check the permutation vector **/
  for(i=0;i<n;i++)
    {
      /*fprintf(stderr, "(%d = %d) ", i, perm[i]);*/
      assert(perm[i] >= 0);
      assert(perm[i] < n);
    }

  bzero(iperm, sizeof(pastix_int_t)*n);
  for(i=0;i<n;i++)
    iperm[perm[i]]++;

  k = 0;
  for(i=0;i<n;i++)
    if(iperm[i] != 1)
      k++;
  if(k>0)
    errorPrint("Number of errors in perm vector in postorder %ld", (long)k);
  assert(k==0);
#endif


  /** Compute the iperm vector **/
  for(i=0;i<n;i++)
    iperm[perm[i]] = i;



}



void compute_subtree_size(pastix_int_t n, pastix_int_t *father, pastix_int_t *perm, pastix_int_t *iperm, pastix_int_t *T)
{
  /********************************************/
  /*  Compute the size of each subtree given  */
  /*  the number of the father of each node   */
  /********************************************/
  pastix_int_t k, i;
  (void)perm;

  /*** OIMBE pas la peine d'utiliser un tas; il suffit de parcourir iperm pour assurer
       de toujours traiter un fils avant son pere ***/

  bzero(T, sizeof(pastix_int_t)*n);

  for(k=0;k<n;k++)
    {
      i = iperm[k];
      T[i]++;
      if(i!=father[i])
        T[father[i]] += T[i];
    }



#ifdef DEBUG_BLEND
 {
   pastix_int_t sum;
   sum = 0;
   for(i=0;i<n;i++)
     if(father[i] == i)
       sum += T[i];

   if(sum != n)
     errorPrint("compute_subtree_size: sum of the subtree = %ld n = %ld", (long)sum, (long)n);
   assert(n == sum);
 }
#endif

}




void compute_elimination_tree(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_int_t *perm, pastix_int_t *iperm, pastix_int_t *father)
{
  /******************************************************************************/
  /* Compute the elimination tree of a matrix A (without computing the symbolic */
  /* factorization) associated with a reordering of the matrix                  */
  /* On entry:                                                                  */
  /* n, ia, ja : the adjacency graph of the matrix (symmetrized)                */
  /* perm : a permutation vector of the matrix                                  */
  /* iperm : the inverse permutation vector                                     */
  /* On return:                                                                 */
  /* father : father[i] = father of node i on the eliminination tree            */
  /* If node i is a root then father[i] = i                                     */
  /* NOTE : father is allocated at a size of n interger on input                */
  /******************************************************************************/
  pastix_int_t i, j, k;
  pastix_int_t node;
  pastix_int_t vroot;

  /** Optim **/
  pastix_int_t flag, ind;
  pastix_int_t *jrev = NULL;
  pastix_int_t *jf   = NULL;


  MALLOC_INTERN(jrev, n, pastix_int_t);
  for(i=0;i<n;i++)
    jrev[i] = -1;
  MALLOC_INTERN(jf, n, pastix_int_t);
  bzero(jf, sizeof(pastix_int_t)*n);

  for(i=0;i<n;i++)
    father[i] = -1;

  for(i=0;i<n;i++)
    {
      ind = 0;
      node = iperm[i];
      for(j=ia[node];j<ia[node+1];j++)
        {

          k = ja[j];
          if(perm[k] < perm[node])
            {
              flag = 1;
              vroot = k;
              while(father[vroot] != -1 && father[vroot] != node)
                {
                  if(jrev[vroot] >= 0)
                    {
                      flag = 0;
                      break;
                    }
                  jrev[vroot] = ind;
                  jf[ind] = vroot;
                  ind++;

                  vroot = father[vroot];
                }
              if(flag == 1)
                father[vroot] = node;
            }
        }
      /** reinit jrev **/
      for(j=0;j<ind;j++)
        jrev[jf[j]]=-1;
    }

  memFree_null(jrev);
  memFree_null(jf);

  for(i=0;i<n;i++)
    if(father[i] == -1)
      father[i]=i;

#ifdef DEBUG_KASS
  /*** Check to see if a father has a lower rank in the permutation array than one of its sons ***/
  for(i=0;i<n;i++)
    {
      if(perm[i] > perm[father[i]])
        {
          fprintf(stderr, "Node %ld perm=%ld Father %ld perm=%ld \n", (long)i, (long)perm[i], (long)father[i], (long)perm[father[i]]);
          assert(perm[i] <= perm[father[i]]);
        }
    }
#endif
}
