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
#include "graph.h"
#include "kass.h"

pastix_int_t
SF_Direct(const kass_csr_t   *graphA,
                pastix_int_t  cblknbr,
          const pastix_int_t *rangtab,
                pastix_int_t *treetab,
                kass_csr_t   *graphL)
{
    /********************************************************/
    /* This function computes the direct factor nnz pattern */
    /* of a matrix A given the supernode partition          */
    /********************************************************/
    pastix_int_t i,j,k, nnz;
    pastix_int_t nnznbr, father;
    pastix_int_t *tmpj      = NULL;
    pastix_int_t *tmp       = NULL;
    pastix_int_t *tmp2      = NULL;
    pastix_int_t *ja        = NULL;
    pastix_int_t *node2cblk = NULL;

    MALLOC_INTERN(tmpj,      graphA->n, pastix_int_t);
    MALLOC_INTERN(tmp,       graphA->n, pastix_int_t);
    MALLOC_INTERN(node2cblk, graphA->n, pastix_int_t);

    for(k=0;k<cblknbr;k++)
        for(i=rangtab[k];i<rangtab[k+1];i++)
            node2cblk[i] = k;

    /* Compute the nnz structure of each supernode in A */
    kass_csrInit( cblknbr, graphL );
    for(k=0; k<cblknbr; k++)
    {
        /* Put the diagonal elements (In case A does not contains them) */
        j = 0;
        for(i=rangtab[k]; i<rangtab[k+1]; i++)
            tmpj[j++] = i;
        nnznbr = j;

        for(i=rangtab[k]; i<rangtab[k+1]; i++)
        {
            j = 0;

            /* We take only the elements greater than i */
            while( (j < graphA->nnz[i]) && (graphA->rows[i][j] <= i) )
                j++;

            /* Merge the actual list with the edges of the ith vertex */
            nnznbr = pastix_intset_union( nnznbr,           tmpj,
                                          graphA->nnz[i]-j, graphA->rows[i] + j,
                                          tmp);

            /* Swap tmpj and the merged set tmp */
            tmp2 = tmpj;
            tmpj = tmp;
            tmp  = tmp2;
        }

#ifdef DEBUG_KASS
        ind = 0;
        for(j=rangtab[k]; j < rangtab[k+1];j++)
            ASSERT(tmpj[ind++] == j, MOD_KASS);
        ASSERT(nnznbr > 0, MOD_KASS);
#endif

        /* Update graphL */
        graphL->nnz[k] = nnznbr;

        MALLOC_INTERN(graphL->rows[k], nnznbr, pastix_int_t);
        memcpy(graphL->rows[k], tmpj, sizeof(pastix_int_t)*nnznbr);
    }

    /** Compute the symbolic factorization **/
    for(k=0;k<cblknbr;k++)
    {
        /*father = treetab[k];*/
        i = 0;
        ja = graphL->rows[k];
        while( (i < graphL->nnz[k]) && (node2cblk[ja[i]] <= k) )
            i++;

        if( i < graphL->nnz[k] )
            father = node2cblk[ja[i]];
        else
            father = -1;
        treetab[k] = father;

        /* Merge son's nodes into father's list */
        if( (father != k) && (father > 0) )
        {
            nnznbr = pastix_intset_union( graphL->nnz[k] - i,  graphL->rows[k] + i,
                                          graphL->nnz[father], graphL->rows[father],
                                          tmpj );

            memFree( graphL->rows[father] );
            MALLOC_INTERN( graphL->rows[father], nnznbr, pastix_int_t);
            memcpy( graphL->rows[father], tmpj, sizeof(pastix_int_t)*nnznbr );
            graphL->nnz[father] = nnznbr;
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

            ind = pastix_intset_union( P->nnzrow[k],   P->ja[k],
                                       A->nnzrow[i]-j, A->ja[i]+j,
                                       tmp );

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
    for(i=0; i<cblknbr; i++)
    {
        pastix_int_t ncol, nrow;
        ncol = rangtab[i+1]-rangtab[i];
        nrow = graphL->nnz[i];

        nnz += (ncol*(ncol+1))/2;
        nnz += (ncol*(nrow-ncol));

#ifdef DEBUG_KASS
        ASSERT(P->nnzrow[i] >= ncol, MOD_KASS);
        if(P->nnzrow[i] >= n)
            fprintf(stderr,"P->nnzrow[%ld] = %ld \n", (long)i, (long)P->nnzrow[i]);
        ASSERT(P->nnzrow[i] < n, MOD_KASS);
#endif
    }
    return nnz;
}
