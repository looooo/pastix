/**
 *
 * @file amalgamate.c
 *
 *  PaStiX symbolic factorization routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * This file contains routines to amalgamate a given graph with a fillin ratio.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "kass.h"
#include "queue.h"

/** 2 percents **/
#define RAT_CBLK 0.02
#define INFINI 1e9

double cblk_time_fact (pastix_int_t n, pastix_int_t *ja, pastix_int_t colnbr);
double cblk_time_solve(pastix_int_t n, pastix_int_t *ja, pastix_int_t colnbr);

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbfact_internal
 *
 * amalgamate_get_sonslist - This function returns the list of sons of a given
 * node.
 *
 *******************************************************************************
 *
 * @param[in] node
 *          Node from which the list is wanted.
 *
 * @param[in] sonindex
 *          Integer array of the indexes of the first son of each node in sontab array.
 *
 * @param[in] sontab
 *          Integer array which contains the list of sons for each node.
 *
 * @param[in] colweight
 *          Integer array of the weight of each node.
 *
 * @param[out] list
 *          Integer array that can contains the list of sons, should be of size
 *          the number of unknowns to prevent overflow.
 *          On exit, contains the list of sons of node.
 *
 *******************************************************************************
 *
 * @return
 *          \retval The number of sons of node.
 *
 *******************************************************************************/
static inline pastix_int_t
amalgamate_get_sonslist(       pastix_int_t  node,
                         const pastix_int_t *sonindex,
                         const pastix_int_t *sontab,
                         const pastix_int_t *colweight,
                               pastix_int_t *list )
{
    pastix_int_t i, s;
    pastix_int_t ind;
    ind = 0;
    for(i=sonindex[node];i<sonindex[node+1];i++)
    {
        s = sontab[i];
        if(colweight[s] <= 0)
        {
            ind += amalgamate_get_sonslist(s, sonindex, sontab, colweight, list+ind);
        }
        else
            list[ind++] = s;
    }
    return ind;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbfact_internal
 *
 * amalgamate_merge_cost - This function computes the cost of merging two nodes
 * a and b together. The additional cost is returned.
 *
 *******************************************************************************
 *
 * @param[in] a
 *          First node to merge. a >= 0.
 *
 * @param[in] b
 *          Second node to merge. b >= 0.
 *
 * @param[in] P
 *          The graph that contains a and b nodes with their associated rows.
 *
 * @param[in] colweight
 *          Integer array of size P->n with the weight of each node.
 *
 *******************************************************************************
 *
 * @return
 *          \retval The additional non zero entries required if both nodes are
 *          merged together.
 *
 *******************************************************************************/
static inline pastix_int_t
amalgamate_merge_cost(       pastix_int_t  a,
                             pastix_int_t  b,
                       const kass_csr_t   *P,
                       const pastix_int_t *colweight )
{
    pastix_int_t ia, na, ib, nb;
    pastix_int_t *ja, *jb;
    pastix_int_t cost, costa, costb;

    ja    = P->rows[a];
    jb    = P->rows[b];
    na    = P->nnz[a];
    nb    = P->nnz[b];
    costa = colweight[a];
    costb = colweight[b];

    ia = ib = 0;

    /* The diagonal elements of row a does not create fill-in */
    while( (ia < na) && (ja[ia] < jb[0]) )
        ia++;

    cost = 0;
    while((ia < na) && (ib < nb))
    {
        if(ja[ia] < jb[ib])
        {
            cost += costb;
            ia++;
            continue;
        }
        if(ja[ia] > jb[ib])
        {
            cost += costa;
            ib++;
            continue;
        }

        assert( ja[ia] == jb[ib] );
        ia++;
        ib++;
    }

    cost += costb * ( na - ia );
    cost += costa * ( nb - ib );

    assert(cost >= 0);
    return cost;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbfact_internal
 *
 * amalgamate_merge_gain - This function computes the computational gain of
 * merging two nodes a and b together according to the cost function
 * cblktime. The cost variation is returned.
 *
 *******************************************************************************
 *
 * @param[in] a
 *          First node to merge. a >= 0.
 *
 * @param[in] b
 *          Second node to merge. b >= 0.
 *
 * @param[in] P
 *          The graph that contains a and b nodes with their associated rows.
 *
 * @param[in] colweight
 *          Integer array of size P->n with the weight of each node.
 *
 * @param[in,out] tmp
 *          Workspace array to compute the union of the rows associated to both nodes.
 *
 * @param[in] cblktime
 *          Function that returns the computational cost of a cblk.
 *
 *******************************************************************************
 *
 * @return
 *          \retval The cost variation between merging the two nodes together vs
 *          keeping them separated according to the cblktime function.
 *
 *******************************************************************************/
static inline double
amalgamate_merge_gain(       pastix_int_t  a,
                             pastix_int_t  b,
                       const kass_csr_t   *P,
                       const pastix_int_t *colweight,
                             pastix_int_t *tmp,
                      double (*cblktime)(pastix_int_t, pastix_int_t *, pastix_int_t))
{
    double costa, costb, costm;
    pastix_int_t nm;

    costa = cblktime(P->nnz[a], P->rows[a], colweight[a]);
    costb = cblktime(P->nnz[b], P->rows[b], colweight[b]);

    nm = pastix_intset_union( P->nnz[a], P->rows[a],
                              P->nnz[b], P->rows[b],
                              tmp);

    costm = cblktime(nm, tmp, colweight[a] + colweight[b]);

    return costm - costa - costb;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbfact_internal
 *
 * amalgamate_merge_col - This function merge two nodes together withing the
 * given graph.
 *
 *******************************************************************************
 *
 * @param[in] a
 *          First node to merge. a >= 0.
 *
 * @param[in] b
 *          Second node to merge. b >= 0.
 *
 * @param[in,out] P
 *          The graph that contains a and b nodes with their associated rows.
 *          On exit, a is merged into b node.
 *
 * @param[in,out] tmp
 *          Workspace array to compute the union of the rows associated to both
 *          nodes.
 *
 *******************************************************************************/
static inline void
amalgamate_merge_col(pastix_int_t  a,
                     pastix_int_t  b,
                     kass_csr_t   *P,
                     pastix_int_t *tmp )
{
    pastix_int_t n;

    n = pastix_intset_union( P->nnz[a], P->rows[a],
                             P->nnz[b], P->rows[b],
                             tmp);

    if(P->rows[a] != NULL)
    {
        P->nnz[a] = 0;
        memFree(P->rows[a]);
    }

    if(P->rows[b] != NULL) {
        memFree(P->rows[b]);
    }

    P->nnz[b] = n;

    MALLOC_INTERN(P->rows[b], n, pastix_int_t);
    memcpy(P->rows[b], tmp, n*sizeof(pastix_int_t));
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbfact
 *
 * amalgamate - This function takes a supernode graph and performs amalgamation
 * of the columns until a fill tolerance is reached. Two level of fill tolerance
 * are available, the first one, rat_cblk, defines a fill tolerance based on the
 * graph structure only, and the second one, rat_blas, allows to keep amalgmate
 * the structure until the higher ratio while it improves the BLAS performance
 * based on the internal cost model.
 *
 *******************************************************************************
 *
 * @param[in] rat_cblk
 *          Percentage of fill-in allowed using structure only amalgamation. Must be >= 0.
 *
 * @param[in] rat_blas
 *          Percentage of fill-in allowed using structure BLAS performance
 *          model. rat_blas >= rat_cblk.
 *          The ratio is always on the number of non zero entries of the
 *          pattern, but the criteria to merge columns is based only on reducing
 *          the cost of computations. The cost function used is a model for
 *          Cholesky factorization or solve in double precision, other
 *          factorizations and/or precisions are not implemented. They are
 *          considered to be similar up to a factor, and so doesn't affect the
 *          amalgamation.
 *          The solve cost is used for ILU(k) factorization and the
 *          factorization cost for direct factorization.
 *
 * @param[in] nnzL
 *          The number of non zero entries in the matrix described by graphL,
 *          usually computed through kassFactLevel() or kassFactDirect().
 *
 * @param[in,out] graphL
 *          The graph of the matrix L to amalgamate. On exit, the compressed
 *          graph is returned (not the quotient graph !!)
 *
 * @param[in] snodenbr
 *          The number of supernodes, if an initial supernode partition is
 *          given, -1 otherwise.
 *
 * @param[in] snodetab
 *          Integer array of size snodenbr+1, that contains the initial
 *          supernode partition, NULL otherwise.
 *
 * @param[in,out] treetab
 *          Integer array of size: snodenbr if snodetab is provided, graphL->n otherwise.
 *          This array contains the father of each supernode/node in the given graph.
 *          On exit, the array is updated to reflect the new amalgamated graph.
 *
 * @param[out] cblknbr
 *          The number of supernodes found after amalgamation process.
 *
 * @param[out] rangtab
 *          Integer array of size cblknbr+1.
 *          On exit, contains the new supernode partition found with
 *          amalgamation process.
 *
 * @param[out] nodetab
 *          Integer array of size cblknbr.
 *          On exit, contains the inverse permutation associated with the new
 *          amalgamated graph.
 *
 *******************************************************************************/
void
amalgamate(double rat_cblk, double rat_blas,
                 kass_csr_t    *graphL,
                 pastix_int_t   nnzL,
                 pastix_int_t   snodenbr,
           const pastix_int_t  *snodetab,
                 pastix_int_t  *treetab,
                 pastix_int_t  *cblknbr,
                 pastix_int_t **rangtab,
                 pastix_int_t  *nodetab,
                 MPI_Comm       pastix_comm )
{
    double (*cblktime)(pastix_int_t, pastix_int_t *, pastix_int_t);

    pastix_int_t *nnzadd    = NULL;
    pastix_int_t *sonindex  = NULL;
    pastix_int_t *sontab    = NULL;
    pastix_int_t *tmp       = NULL;
    pastix_int_t *tmp2      = NULL;
    pastix_int_t *newnum    = NULL;
    pastix_int_t *colweight = NULL;

    pastix_int_t  i, j, k, ind;
    pastix_int_t  father;
    pastix_int_t  n, nn, nbcblk_merged;
    pastix_int_t  fillwhile, fillcblk, fillblas, fill;

    double *gain   = NULL;
    double  key;
    Queue heap;
    int blas_gain = 0;
    int procnum;

    (void)pastix_comm;

    MPI_Comm_rank(pastix_comm, &procnum);

    /* If snodetab is NULL, we are doing incomplete factorization */
    if ( snodetab == NULL )
        cblktime = cblk_time_solve;
    else
        cblktime = cblk_time_fact;

    /*
     * If rat is negative, then the amalgamation is computed first with nnz
     * structure up to rat_cblk ratio and then an extra amalgamation is performed
     * as long as it improves the BLAS performance, and stay under the (-rat)
     * ratio.
     */
    if(rat_blas < 0)
    {
        rat_blas = -rat_blas;
        rat_cblk = MIN(rat_blas, RAT_CBLK);
    }
    else
    {
        rat_blas = rat_blas;
        rat_cblk = rat_cblk;
    }

    /*
     * We always start by the structural amalgamation and then we start the
     * almalgamation leaded by the reduction of the blas cost
     */
    blas_gain = 0;

    /*
     * Compute the initial weight of each cblk
     */
    n = graphL->n;
    MALLOC_INTERN(colweight, n, pastix_int_t);

    if(snodetab != NULL)
    {
        assert( graphL->n == snodenbr );
        for(i=0; i<n; i++) {
            colweight[i] = snodetab[i+1] - snodetab[i];
            assert( (graphL->nnz[i] >= colweight[i]) &&
                    (colweight[i] > 0) );

#if defined(PASTIX_DEBUG_SYMBOL)
            {
                pastix_int_t j;
                /* Check that the first elements of graphL are those from the supernode (the diagonal elements) */
                for(j=0;j<n;j++)
                {
                    pastix_int_t l, k = 0;
                    for(l=snodetab[j]; l<snodetab[j+1]; l++) {
                        assert(graphL->rows[j][k++] == l);
                    }
                }
            }
#endif
        }

        /** Number of unknowns **/
        nn = snodetab[snodenbr];
    }
    else
    {
        nn = n;
        for(i=0; i<n; i++)
            colweight[i] = 1;
    }
    assert( nn >= n );

    /**********************************************/
    /*** Compute the maximum extra fill allowed ***/
    /**********************************************/
    /*
     * - fillblas is the limit of nnz added during the whole amalgamation
     *   process
     * - fillcblk is the limit of nnz added under the amalgamation process that
     *   tries to reduce the number of supernode
     * - Between fillcblk and fillblas the amalgamation tries to reduce the time
     *   of the solve or factorization time
     */
    fillcblk  = (pastix_int_t)lround((double)nnzL * rat_cblk);
    fillblas  = (pastix_int_t)lround((double)nnzL * rat_blas);
    fillwhile = fillcblk;

#if defined(PASTIX_DEBUG_SYMBOL)
    {
        double key = 0.0;

        pastix_print(procnum, 0, "fillcblk %ld fillmax %ld \n", (long)fillcblk, (long)fillblas);

        for(i=0;i<graphL->n;i++)
            key += cblktime(graphL->nnz[i],
                            graphL->rows[i],
                            colweight[i]);
        pastix_print(procnum, 0, "COST of the NON AMALGAMATED MATRIX = %e\n", key);
    }
#endif

    MALLOC_INTERN(tmp,   nn, pastix_int_t);

    /**********************************************/
    /* Compute the list of sons of each supernode */
    /**********************************************/
    {
        pastix_int_t nbsons = 0;

        /* Compute the number of sons of each cblk */
        memset( tmp, 0, n * sizeof(pastix_int_t));
        for(i=0;i<n-1;i++) {
            if(treetab[i] >= 0) { /** IF THIS SNODE IS NOT A ROOT **/
                tmp[treetab[i]]++;
                nbsons++;
            }
        }
        assert( nbsons < n );

        MALLOC_INTERN(sonindex, n+1,    pastix_int_t);
        MALLOC_INTERN(sontab,   nbsons, pastix_int_t);

        /* Create the sonindex array */
        sonindex[0] = 0;
        for(i=0;i<n;i++) {
            sonindex[i+1] = sonindex[i] + tmp[i];
        }

        /* Fill the sontab array */
        memcpy( tmp, sonindex, n*sizeof(pastix_int_t) );
        for(i=0;i<n-1;i++)
        {
            pastix_int_t father = treetab[i];
            assert( (father != i) && (father < n) );
            if ( father >= 0 ) /** IF THIS SNODE IS NOT A ROOT **/
            {
                sontab[ tmp[father] ] = i;
                tmp[father]++;
                assert( tmp[father] <= sonindex[father+1] );
            }
        }
    }

    /****************************************************************/
    /* Compute the fill to merge a column of graphL with its father */
    /****************************************************************/
    MALLOC_INTERN(nnzadd, n, pastix_int_t);
    MALLOC_INTERN(gain,   n, double);
    for(i=0;i<n;i++)
    {
        father = treetab[i];
        if(father == -1 || father == i)
        {
            nnzadd[i] = INFINI;
            gain[i] = INFINI;
            continue;
        }

        nnzadd[i] = amalgamate_merge_cost(i, father, graphL, colweight);
        gain[i]   = (double)(nnzadd[i]);
    }

    /*****************************************************/
    /** Merge all the columns so it doesn't add fill-in **/
    /*****************************************************/
    nbcblk_merged = 0;
    for(i=0; i<n; i++)
    {
        assert(colweight[i] > 0);

        if( (colweight[i] != 0 ) &&
            (nnzadd[i]    == 0 ) )
        {
            father = treetab[i];
            assert((father > 0) && (father != i));

            nbcblk_merged++;

            /* Merge the node i and its father **/
            amalgamate_merge_col( i, father, graphL, tmp );

            /* Update cost of each node */
            colweight[father] += colweight[i];
            colweight[i]       = 0;

            /* Update the nnzadd and gain array for the father */
            k = treetab[father];
            if((k != -1) && (k != father))
            {
                nnzadd[father] = amalgamate_merge_cost(father, k, graphL, colweight);
                gain[father]   = (double)(nnzadd[father]);
            }

            /* Update the list of sons of i now */
            ind = amalgamate_get_sonslist(i, sonindex, sontab, colweight, tmp);
            for(j=0;j<ind;j++)
            {
                k = tmp[j];
                assert( k != father );
                treetab[k] = father;
                nnzadd[k]  = amalgamate_merge_cost(k, father, graphL, colweight);
                gain[k]    = (double)(nnzadd[k]);
            }
        }
    }

    pastix_print(procnum, 0, "Number of cblk after amalgamation initialization = %ld (reduced by %ld)\n",
                 (long)(n-nbcblk_merged), (long)nbcblk_merged);

    /* Initialize the first round with column sorted by nnzadd */
    queueInit(&heap, (n-nbcblk_merged));
    for(i=0;i<n;i++) {
        if( (colweight[i] >  0) &&
            (treetab[i]   >  0) &&
            (treetab[i]   != i) )
        {
            queueAdd(&heap, i, gain[i]);
        }
    }

    /*********************************************/
    /* Merge supernodes until we reach fillwhile */
    /*********************************************/
  restart:

    /* Merge supernodes untill we reach the fillwhile limit */
    fill = 0.0;
    while( (queueSize(&heap) > 0) &&
           (fill < fillwhile) )
    {
        i = queueGet2(&heap, &key, NULL);

        /* If we pop an old version of i, let's skip it */
        if( (gain[i] != key) || (colweight[i] <= 0) )
            continue;

        /* check if we reached the fillwhile */
        if((fill + nnzadd[i]) > fillwhile)
            continue;
        else
            fill += nnzadd[i];

        nbcblk_merged++;
        father = treetab[i];
        assert((father > 0) && (father != i));
        assert( colweight[father] > 0 );

        /* Merge the node i and its father **/
        amalgamate_merge_col( i, father, graphL, tmp );

        /* Update cost of each node */
        colweight[father] += colweight[i];
        colweight[i]       = 0;

        /* Update the nnzadd and gain array for the father */
        k = treetab[father];
        if((k != -1) && (k != father))
        {
            nnzadd[father] = amalgamate_merge_cost(father, k, graphL, colweight);
            if(blas_gain == 1) {
                gain[father] = amalgamate_merge_gain(father, k, graphL, colweight, tmp2, cblktime) / nnzadd[father];
                if(gain[father] <= 0)
                    queueAdd(&heap, father, gain[father]);
            }
            else {
                gain[father] = (double)(nnzadd[father]);
                queueAdd(&heap, father, gain[father]);
            }
        }

        /* Update the list of sons of i now */
        ind = amalgamate_get_sonslist(i, sonindex, sontab, colweight, tmp);
        for(j=0;j<ind;j++)
        {
            k = tmp[j];
            assert( k != father );
            treetab[k] = father;
            nnzadd[k] = amalgamate_merge_cost(k, father, graphL, colweight);

            if(blas_gain == 1) {
                gain[k] = amalgamate_merge_gain(k, father, graphL, colweight, tmp2, cblktime) / nnzadd[k];
                if(gain[k] <= 0)
                    queueAdd(&heap, k, gain[k]);
            }
            else {
                gain[k] = (double)(nnzadd[k]);
                queueAdd(&heap, k, gain[k]);
            }
        }
    }

    pastix_print(procnum, 0, "Number of cblk after amalgamation phase = %ld (reduced by %ld)\n",
                 (long)(n-nbcblk_merged), (long)nbcblk_merged);

    queueExit(&heap);
    if(fillwhile < fillblas)
    {
        /* Now the gain of amalgamation is based on the BLAS model */
        assert(blas_gain == 0);

        blas_gain = 1;
        fillwhile = fillblas;

        /* Recompute the gain using BLAS model to restart the second round */
        queueInit(&heap, (n-nbcblk_merged));
        for(i=0;i<n;i++)
        {
            father = treetab[i];
            if(father == -1 || father == i)
            {
                gain[i] = INFINI;
                continue;
            }
            gain[i] = amalgamate_merge_gain(i, father, graphL, colweight, tmp, cblktime)/nnzadd[i];

            if( (colweight[i] >  0) &&
                (treetab[i]   >  0) &&
                (treetab[i]   != i) &&
                (gain[i] <= 0.) )
            {
                queueAdd(&heap, i, gain[i]);
            }
        }

        /* Allocate second temporary workspace for amalgamate_merge_gain */
        MALLOC_INTERN(tmp2, nn, pastix_int_t);
        goto restart;
    }

    if (tmp2 != NULL) memFree(tmp2);
    memFree(nnzadd);
    memFree(gain);

    /********************************/
    /* Compute the new partition    */
    /********************************/

    /* Create the new numering array of the comptacted supernodes (stored in
     * tmp) and compute the size of each supernode */
    newnum = tmp;
    MALLOC_INTERN(*rangtab, n-nbcblk_merged+1, pastix_int_t);
    memset(newnum,   0,                   n * sizeof(pastix_int_t));
    memset(*rangtab, 0, (n-nbcblk_merged+1) * sizeof(pastix_int_t));

    k = 0;
    for(i=0;i<n;i++) {
        if(colweight[i] > 0) {
            newnum[i] = k++;
            (*rangtab)[ k ] = colweight[i];
        }
    }
    assert( k == (n-nbcblk_merged) );
    *cblknbr = k;

    /* Create the rangtab */
    for(i=0; i<k; i++)
        (*rangtab)[i+1] += (*rangtab)[i];

    /* Create the invp vector to adapt to the new partition */
    for(i=0;i<n;i++)
    {
        father = i;
        if(colweight[i] <= 0)
        {
            /* find the cblk this node is in */
            while(colweight[father] <= 0)
            {
                father = treetab[father];
                assert(father > 0);
            }
        }

        if(snodetab != NULL)
            for(j=snodetab[i];j<snodetab[i+1];j++)
                nodetab[(*rangtab)[newnum[father]]++] = j;
        else
            nodetab[(*rangtab)[newnum[father]]++] = i;
    }

    /* Reset rangtab to its real value */
    for(i=*cblknbr; i>0; i--)
        (*rangtab)[i] = (*rangtab)[i-1];
    (*rangtab)[0] = 0;

    /* Compact the graph to remove merged nodes */
    kass_csrCompact( graphL );
    assert( graphL->n == *cblknbr );

    /* Apply the new permutation to graphL */
    {
        /* Create the perm vector */
        for(i=0;i<nn;i++)
            tmp[nodetab[i]] = i;

        for(i=0;i<graphL->n;i++)
        {
            pastix_int_t *ja;
            ja = graphL->rows[i];
            for(j=0; j<graphL->nnz[i]; j++)
                ja[j] = tmp[ja[j]];
        }
    }

#if defined(PASTIX_DEBUG_SYMBOL)
    {
        double key = 0.0;

        for(i=0;i<graphL->n;i++)
            key += cblktime(graphL->nnz[i],
                            graphL->rows[i],
                            colweight[i]);
        pastix_print(procnum, 0, "COST of the AMALGAMATED MATRIX = %e\n", key);
    }
#endif

    memFree(tmp);
    memFree(sonindex);
    memFree(sontab);
    memFree(colweight);
}

