/************************************************************/
/**                                                        **/
/**   NAME       : amalgamate.c                            **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 15/08/2006      **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "common.h"
#include "sparRow.h"
#include "SF_Direct.h"


void UnionSet(pastix_int_t *set1,
              pastix_int_t  n1,
              pastix_int_t *set2,
              pastix_int_t  n2,
              pastix_int_t *set,
              pastix_int_t *n);

pastix_int_t SF_Direct(csptr A,
                             pastix_int_t  cblknbr,
                       const pastix_int_t *rangtab,
                             pastix_int_t *treetab,
                       csptr P)
{
    /********************************************************/
    /* This function computes the direct factor nnz pattern */
    /* of a matrix A given the supernode partition          */
    /********************************************************/
    pastix_int_t i,j,k, nnz;
    pastix_int_t ind, nnznbr, father;
    pastix_int_t *tmpj      = NULL;
    pastix_int_t *tmp       = NULL;
    pastix_int_t *tmp2      = NULL;
    pastix_int_t *ja        = NULL;
    pastix_int_t *node2cblk = NULL;

    MALLOC_INTERN(tmpj,      A->n, pastix_int_t);
    MALLOC_INTERN(tmp,       A->n, pastix_int_t);
    MALLOC_INTERN(node2cblk, A->n, pastix_int_t);

    for(k=0;k<cblknbr;k++)
        for(i=rangtab[k];i<rangtab[k+1];i++)
            node2cblk[i] = k;

    /** Compute the nnz structure of each supernode in A **/
    for(k=0;k<cblknbr;k++)
    {
        ind = rangtab[k];

        /** Put the diagonal elements (A does not contains them) **/
        j = 0;
        for(i=rangtab[k];i<rangtab[k+1];i++)
            tmpj[j++] = i;
        nnznbr = j;

        for(i=rangtab[k];i<rangtab[k+1];i++)
        {
            j = 0;
            while(j<A->nnzrow[i] && A->ja[i][j] <= i)
                j++;

            /** Realise la fusion de 2 listes triees croissantes **/
            UnionSet(tmpj, nnznbr, A->ja[i]+j, A->nnzrow[i]-j, tmp, &ind);

            /** echange tmpj et le resultat de la fusion **/
            nnznbr = ind;
            tmp2 = tmpj;
            tmpj = tmp;
            tmp  = tmp2;
        }
#ifdef DEBUG_KASS
        ind = 0;
        for(j=rangtab[k]; j < rangtab[k+1];j++)
            ASSERT(tmpj[ind++] == j, MOD_KASS);
#endif

#ifdef DEBUG_KASS
        ASSERT(nnznbr > 0, MOD_KASS);
#endif
        P->nnzrow[k] = nnznbr;
        MALLOC_INTERN(P->ja[k], nnznbr, pastix_int_t);
        memcpy(P->ja[k], tmpj, sizeof(pastix_int_t)*nnznbr);
        P->ma[k]=NULL;
    }

    P->n = cblknbr;

    /** Compute the symbolic factorization **/
    for(k=0;k<cblknbr;k++)
    {
        /*father = treetab[k];*/
        i = 0;
        ja = P->ja[k];
        while(i < P->nnzrow[k] && node2cblk[ja[i]] <= k)
            i++;
        if(i<P->nnzrow[k])
            father = node2cblk[ja[i]];
        else
            father = -1;
        treetab[k] = father;

        if(father != k && father > 0)
        {
            /*i = 0;
             while(i < P->nnzrow[k] && P->ja[k][i] <= rangtab[father])
             i++;
             if(i == P->nnzrow[k])
             continue;*/

            UnionSet(P->ja[k]+i, P->nnzrow[k]-i, P->ja[father], P->nnzrow[father], tmpj, &nnznbr);

            memFree(P->ja[father]);
            MALLOC_INTERN(P->ja[father], nnznbr, pastix_int_t);
            memcpy(P->ja[father], tmpj, sizeof(pastix_int_t)*nnznbr);
            P->nnzrow[father] = nnznbr;
        }
    }


#ifdef DEBUG_KASS
    for(i=0;i<cblknbr;i++)
    {
        pastix_int_t j;
        ASSERT(P->nnzrow[i] >= rangtab[i+1]-rangtab[i], MOD_KASS);
        k = 0;
        for(j=rangtab[i]; j < rangtab[i+1];j++)
            ASSERT(P->ja[i][k++] == j, MOD_KASS);
    }

    /** Check that all terms of A are in the pattern **/
    for(k=0;k<cblknbr;k++)
    {
        /** Put the diagonal elements (A does not contains them) **/
        for(i=rangtab[k];i<rangtab[k+1];i++)
        {
            j = 0;
            while(j<A->nnzrow[i] && A->ja[i][j] < i)
                j++;

            for(ind = j;ind < A->nnzrow[i];ind++)
                assert(A->ja[i][ind] >= i);
            for(ind = j+1;ind < A->nnzrow[i];ind++)
                assert(A->ja[i][ind] > A->ja[i][ind-1]);


            UnionSet(P->ja[k], P->nnzrow[k], A->ja[i]+j, A->nnzrow[i]-j, tmp, &ind);
            if(ind > P->nnzrow[k])
                fprintf(stderr, "k=%ld [%ld %ld]  i=%ld ind %ld nnz %ld \n",
                        (long)k, (long)rangtab[k], (long)rangtab[k+1], (long)i, (long)ind, (long)P->nnzrow[k]);

            ASSERT(ind <= P->nnzrow[k], MOD_KASS);
        }
    }
#endif

    memFree(node2cblk);
    memFree(tmpj);
    memFree(tmp);

    nnz = 0;
    for(i=0; i<P->n; i++)
    {
        pastix_int_t ncol;
        ncol = rangtab[i+1]-rangtab[i];
        nnz += (ncol*(ncol+1))/2;
#ifdef DEBUG_KASS
        ASSERT(P->nnzrow[i] >= ncol, MOD_KASS);
        if(P->nnzrow[i] >= n)
            fprintf(stderr,"P->nnzrow[%ld] = %ld \n", (long)i, (long)P->nnzrow[i]);
        ASSERT(P->nnzrow[i] < n, MOD_KASS);
#endif
        nnz += (P->nnzrow[i]-ncol)*ncol;
    }
#ifdef DEBUG_KASS
    print_one("NNZL = %ld \n", (long)nnz);
#endif

    return nnz;
}


void UnionSet(pastix_int_t *set1, pastix_int_t n1, pastix_int_t *set2, pastix_int_t n2, pastix_int_t *set, pastix_int_t *n)
{
    /********************************************************/
    /* Compute the union of two sorted set                  */
    /* set must have a big enough size to contain the union */
    /* i.e n1+n2                                            */
    /********************************************************/
    pastix_int_t ind, ind1, ind2;

    ind = 0;
    ind1 = 0;
    ind2 = 0;

#ifdef DEBUG_KASS
    {
        pastix_int_t i;
        for(i=0;i<n1-1;i++)
            ASSERT(set1[i] < set1[i+1], MOD_KASS);
        for(i=0;i<n2-1;i++)
            ASSERT(set2[i] < set2[i+1], MOD_KASS);
    }
#endif

    while(ind1 < n1 && ind2 < n2)
    {
        if(set1[ind1] == set2[ind2])
        {
            set[ind] = set1[ind1];
            ind++;
            ind1++;
            ind2++;
            continue;
        }

        if(set1[ind1] < set2[ind2])
        {
            set[ind] = set1[ind1];
            ind++;
            ind1++;
            continue;
        }

        if(set1[ind1] > set2[ind2])
        {
            set[ind] = set2[ind2];
            ind++;
            ind2++;
            continue;
        }
    }

    while(ind1 < n1)
    {
        set[ind] = set1[ind1];
        ind++;
        ind1++;
    }
    while(ind2 < n2)
    {
        set[ind] = set2[ind2];
        ind++;
        ind2++;
    }
#ifdef DEBUG_KASS
    ASSERT(ind <= ind1 +ind2, MOD_KASS);
    ASSERT(ind >= MAX(ind1, ind2), MOD_KASS);
#endif


    *n = ind;
}
