
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
#include "kass.h"
#include "queue.h"
#include "perf.h"
#include "sparRow.h"
#include "amalgamate.h"

/*#define BLAS_GAIN*/    /** Amalgamation use the best ratio time/nnzadd
                         to merge cblk **/
/*#define BLAS_SOLVE*/   /** Amalgamation seeks to optimize the triangular
                         solve **/

#ifdef BLAS_SOLVE
#define CBLKTIME  cblk_time_fact
#else
#define CBLKTIME  cblk_time_solve
#endif

/** 2 percents **/
#define RAT_CBLK 0.02

#define INFINI 10e6

#define print_one(fmt, ...)    if( procnum == 0) fprintf(stdout, fmt, ##__VA_ARGS__)

double cblk_time_fact (pastix_int_t n, pastix_int_t *ja, pastix_int_t colnbr);
double cblk_time_solve(pastix_int_t n, pastix_int_t *ja, pastix_int_t colnbr);
double merge_gain(pastix_int_t a, pastix_int_t b, csptr P, pastix_int_t *colweight, pastix_int_t *tmp);
void   merge_col (pastix_int_t a, pastix_int_t b, csptr P);
void   get_son(pastix_int_t node, pastix_int_t *sonindex, pastix_int_t *sontab, pastix_int_t *colweight, pastix_int_t *ns, pastix_int_t *list);

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

static inline double
amalgamate_merge_gain(       pastix_int_t  a,
                             pastix_int_t  b,
                       const kass_csr_t   *P,
                       const pastix_int_t *colweight,
                             pastix_int_t *tmp)
{
    double costa, costb, costm;
    pastix_int_t nm;

    costa = CBLKTIME(P->nnz[a], P->rows[a], colweight[a]);
    costb = CBLKTIME(P->nnz[b], P->rows[b], colweight[b]);

    nm = pastix_intset_union( P->nnz[a], P->rows[a],
                              P->nnz[b], P->rows[b],
                              tmp);

    costm = CBLKTIME(nm, tmp, colweight[a] + colweight[b]);

    assert( costm <= (costa + costb) );
    return costm - costa - costb;
}

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

void amalgamate2(double rat,
                 kass_csr_t    *graphL,
                 pastix_int_t   nnzL,
                 pastix_int_t   snodenbr,
                 pastix_int_t  *snodetab,
                 pastix_int_t  *treetab,
                 pastix_int_t  *cblknbr,
                 pastix_int_t **rangtab,
                 pastix_int_t  *nodetab,
                 MPI_Comm pastix_comm )
{
    /********************************************************************/
    /* amagamate takes a supernode graph (P,colwgt) and performs some   */
    /* amalgation of the column until a fill tolerance (rat) is reached */
    /* It returns the partition of columns obtained                     */
    /* On return                                                        */
    /* P  is updated to form the column compressed graph (not the
     quotient graph !!) */
    /* treetab is also updated to fit with the new column compressed    */
    /*   graph
     */
    /* cblknbr, rangtab is the supernode partition of the initial
     matrix P */
    /********************************************************************/

    pastix_int_t *nnzadd    = NULL;
    pastix_int_t *sonindex  = NULL;

    pastix_int_t *sontab    = NULL;
    pastix_int_t *tmp       = NULL;
    pastix_int_t *tmp2      = NULL;
    pastix_int_t *newnum    = NULL;
    pastix_int_t i,j,k,ind;
    pastix_int_t father;
    pastix_int_t cblknum;
    pastix_int_t *colweight = NULL;
    double *gain   = NULL;
    pastix_int_t toto;
    pastix_int_t n, nn;
    double rat_cblk;
    int blas_gain = 0;

    double key;
    Queue heap;
    long fillmax, fill;
    long fillcblk, fillwhile;

    int procnum;
    (void)pastix_comm;

    MPI_Comm_rank(pastix_comm, &procnum);

    /* TODO: need some explanation */
    if(rat < 0)
    {
        rat = -rat;
        rat_cblk = MIN(rat, RAT_CBLK);
    }
    else
        rat_cblk = rat;

    /*
     * We always begin by amalgation leaded by the
     * reduction of the number of cblk (until rat_cblk)
     */
    blas_gain = 0;

    /*** Set the weight of each column ***/
    n = graphL->n;
    MALLOC_INTERN(colweight, n, pastix_int_t);

    if(snodetab != NULL)
    {
        assert( graphL->n == snodenbr );
        for(i=0; i<n; i++)
            colweight[i] = snodetab[i+1] - snodetab[i];
        nn = snodetab[snodenbr]; /** Number of unknowns **/
    }
    else
    {
        nn = n;
        for(i=0; i<n; i++)
            colweight[i] = 1;
    }

#ifdef DEBUG_KASS
    if(snodetab != NULL)
    {
        for(i=0;i<n;i++)
        {
            assert( graphL->nnz[i] >= (snodetab[i+1]-snodetab[i]) )
            k = 0;
            for(j=snodetab[i]; j<snodetab[i+1]; j++)
                assert(graphL->rows[i][k++] == j, MOD_KASS);
        }
    }
    for(i=0;i<n;i++)
        ASSERT(colweight[i] >= 0, MOD_KASS);
#endif


    /**********************************************/
    /*** Compute the maximum extra fill allowed ***/
    /**********************************************/
    /** fillmax is the limit of nnz added du to the whole amalgamation
     process **/
    /** fillcblk (<fillmax) is the limit under which the amalgamation
     tries to reduce the number of supernode **/
    /** between fillcblk and fillmax the amalgamation tries to reduce
     the time of the solve or factorization time **/

    fillmax   = (long)nnzL * (long)rat;
    fillcblk  = (long)nnzL * (long)rat_cblk;
    fillwhile = fillcblk;

#ifdef DEBUG_KASS
    {
        double key = 0.0;

        print_one("fillcblk %ld fillmax %ld \n", (long)fillcblk, (long)fillmax);

        for(i=0;i<graphL->n;i++)
            key += CBLKTIME(graphL->nnz[i],
                            graphL->row[i],
                            colweight[i]);
        fprintf(stderr, "COST of the NON AMALGAMTED MATRIX = %g \n", key);
    }
#endif

    MALLOC_INTERN(nnzadd,   n, pastix_int_t);
    MALLOC_INTERN(gain,     n, double);
    MALLOC_INTERN(tmp,  nn, pastix_int_t);
    MALLOC_INTERN(tmp2, nn, pastix_int_t);

    /**********************************************/
    /*** Compute the son list of each supernode ***/
    /**********************************************/
    /*** Compute the number of sons of each cblk ***/
    bzero(tmp, sizeof(pastix_int_t)*n);
    for(i=0;i<n-1;i++)
        if(treetab[i] >= 0) /** IF THIS SNODE IS NOT A ROOT **/
            tmp[treetab[i]]++;

    ind = 0;
    for(i=0;i<n;i++)
        ind += tmp[i];

    MALLOC_INTERN(sonindex, n+1, pastix_int_t);
    MALLOC_INTERN(sontab,   ind, pastix_int_t);
    ind = 0;
    for(i=0;i<n;i++)
    {
        sonindex[i] = ind;
        ind += tmp[i];
    }
    sonindex[n] = ind;

    bzero(tmp, sizeof(pastix_int_t)*n);

    for(i=0;i<n-1;i++)
    {
        cblknum = treetab[i];
        if(cblknum >= 0) /** IF THIS SNODE IS NOT A ROOT **/
        {
            sontab[sonindex[cblknum]+tmp[cblknum]] = i;
            tmp[cblknum]++;
        }
    }


    /***********************************************************/
    /* Compute the fill to merge a column of P with its father */
    /***********************************************************/
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

        if(blas_gain == 1)
            gain[i] = amalgamate_merge_gain(i, father, graphL, colweight, tmp2)/nnzadd[i];
        else
            gain[i] = (double)(nnzadd[i]);
    }

    /*****************************************************/
    /** Merge all the columns so it doesn't add fill-in **/
    /*****************************************************/
    toto = 0;
    for(i=0;i<n;i++)
    {

#ifdef DEBUG_KASS
        ASSERT(colweight[i] > 0, MOD_KASS);
#endif
        if(colweight[i] != 0 && nnzadd[i] == 0 && gain[i] <= 0)
        {
            father = treetab[i];

            toto++;
#ifdef DEBUG_KASS
            ASSERT(father > 0 && father != i, MOD_KASS);
#endif
            /** We merge the snode i and its father **/
            amalgamate_merge_col( i, father, graphL, tmp );

            /*fprintf(stderr, "MERGE %ld (%ld) and %ld (%ld) \n", (long)i, (long)colweight[i],
             (long)father, (long)colweight[father]); */

            colweight[father] += colweight[i];
            colweight[i] = 0; /** mark this node as does not exist **/

            /**  we update nnzadd for the father **/
            k = treetab[father];
            if(k != -1 && k != father)
            {
                nnzadd[father] = amalgamate_merge_cost(father, k, graphL, colweight);

                if(blas_gain == 1)
                    gain[father] = amalgamate_merge_gain(father, k, graphL, colweight, tmp2)/nnzadd[father];
                else
                    gain[father] = (double)(nnzadd[father]);
            }

            /** We update the sons of i now **/
            get_son(i, sonindex, sontab, colweight, &ind, tmp);
            for(j=0;j<ind;j++)
            {
                k = tmp[j];
                treetab[k] = father;
                nnzadd[k] = amalgamate_merge_cost(k, father, graphL, colweight);


                /** @@@ OIMBE inutile : blas_gain vaut 0 au depart maintenant **/
                if(blas_gain == 1)
                    gain[k] = amalgamate_merge_gain(k, father, graphL, colweight, tmp2)/nnzadd[k];
                else
                    gain[k] = (double)(nnzadd[k]);

#ifdef DEBUG_KASS
                ASSERT(nnzadd[k] > 0, MOD_KASS);
#endif
            }

        }
    }
#ifdef DEBUG_KASS
    print_one("Apres amalgamation init cblk = %ld \n", (long)(n-toto));
#endif
    /*** Put in a sort heap the column sorted by their nnzadd ***/
  debut:
    queueInit(&heap, n);
    for(i=0;i<n;i++)
        if(colweight[i] > 0 && (treetab[i] > 0 && treetab[i] != i))
        {
#ifdef DEBUG_KASS
            if(blas_gain != 1)
            {
                if(nnzadd[i] <= 0)
                    fprintf(stderr, "nnzadd[%ld] = %ld \n", (long)i, (long)nnzadd[i]);
                ASSERT(nnzadd[i]>0, MOD_KASS);
            }
#endif

            /*#ifdef BLAS_GAIN*/
            if(blas_gain == 1)
            {
                if(gain[i] <= 0)
                    queueAdd(&heap, i, gain[i]);
            }
            /*#else*/
            else
                queueAdd(&heap, i, gain[i]);
            /*#endif*/
        }

    /*******************************************/
    /* Merge supernodes until we reach fillmax */
    /*******************************************/
    /*** Merge supernodes untill we reach the fillmax limit ****/
    fill = 0.0;


    while(queueSize(&heap)>0 && fill < fillwhile)
    {
        i = queueGet2(&heap, &key, NULL);


        /*if(nnzadd[i] != (pastix_int_t)key || colweight[i] <= 0)*/
        if(gain[i] != key || colweight[i] <= 0)
            continue;

        if(fill + nnzadd[i] > fillmax)
            break;
        else
            fill += nnzadd[i];




        toto++;
        father = treetab[i];
#ifdef DEBUG_KASS
        ASSERT(father > 0 && father != i, MOD_KASS);
        ASSERT(colweight[father]>0, MOD_KASS);
#endif

        /*fprintf(stderr, "Merge col %ld and %ld gain = %g \n", (long)i, (long)father, gain[i]);*/

        /*
         * We merge the snode i and its father and we update treetab, nnzadd,
         * colweight
         */
        amalgamate_merge_col(i, father, graphL, tmp);
        colweight[father] += colweight[i];
        colweight[i]      = 0;

        /**  we update nnzadd for the father **/
        k = treetab[father];
        if(k != -1 && k != father)
        {
            nnzadd[father] = amalgamate_merge_cost(father, k, graphL, colweight);

            if(blas_gain == 1)
                gain[father] = amalgamate_merge_gain(father, k, graphL, colweight, tmp2)/nnzadd[father];
            else
                gain[father] = (double)(nnzadd[father]);

            /*queueAdd(&heap, father, (double) nnzadd[father]);*/
            /*#ifdef BLAS_GAIN*/
            if(blas_gain == 1)
            {
                if(gain[father] <= 0)
                    queueAdd(&heap, father, gain[father]);
            }
            else
                /*#else*/
                queueAdd(&heap, father, gain[father]);
            /*#endif*/
        }

        /** We update the sons of i now **/
        get_son(i, sonindex, sontab, colweight, &ind, tmp);
        for(j=0;j<ind;j++)
        {
            k = tmp[j];
            treetab[k] = father;
            nnzadd[k] = amalgamate_merge_cost(k, father, graphL, colweight);
            /*#ifdef BLAS_GAIN*/
            if(blas_gain == 1)
            {
                gain[k] = amalgamate_merge_gain(k, father, graphL, colweight, tmp2)/nnzadd[k];
                if(gain[k] <= 0)
                    queueAdd(&heap, k, gain[k]);
            }
            /*#else*/
            else
            {
                gain[k] = (double)(nnzadd[k]);
                queueAdd(&heap, k, gain[k]);
            }
            /*#endif*/
        }
    }
#ifdef DEBUG_KASS
    print_one("After amalg phase cblk = %ld fillwhile %ld fillmax %ld \n", (long)toto, (long)fillwhile, (long)fillmax);
#endif
    if(fillwhile < fillmax)
    {

        fillwhile = fillmax;

        /** Now the gain of amalgamation is based on the BLAS model **/
        queueExit(&heap);
#ifdef DEBUG_KASS
        ASSERT(blas_gain == 0, MOD_KASS);
#endif
        blas_gain = 1;

        /** Recompute the gain using BLAS model **/
        for(i=0;i<n;i++)
        {
            father = treetab[i];
            if(father == -1 || father == i)
            {
                gain[i] = INFINI;
                continue;
            }
            gain[i] = amalgamate_merge_gain(i, father, graphL, colweight, tmp2)/nnzadd[i];
        }


        goto debut;
    }


    memFree(nnzadd);
    memFree(gain);
    queueExit(&heap);

    /*fprintf(stderr, "FINAL cblk = %ld \n", (long)(n-toto));*/



    /********************************/
    /* Compute the new partition    */
    /********************************/

    /** Count the number of supernodes **/

    /** tmp will be the newnum of node i in the rangtab **/
    newnum = tmp;
    bzero(newnum, sizeof(pastix_int_t)*n);
    k = 0;
    for(i=0;i<n;i++)
        if(colweight[i] > 0)
            newnum[i] = k++;
    *cblknbr = k;
#ifdef DEBUG_KASS
    print_one("Number of cblk after amal = %ld \n", (long)*cblknbr);
#endif

    MALLOC_INTERN(*rangtab, k+1, pastix_int_t);
    bzero(*rangtab, sizeof(pastix_int_t)*(k+1));

    for(i=0;i<n;i++)
        if(colweight[i] > 0)
            (*rangtab)[newnum[i] + 1]+= colweight[i];

    for(i=1;i<= (*cblknbr);i++)
        (*rangtab)[i] += (*rangtab)[i-1];


    for(i=0;i<n;i++)
    {
        if(colweight[i] > 0)
        {
            if(snodetab != NULL)
                for(j=snodetab[i];j<snodetab[i+1];j++)
                    nodetab[(*rangtab)[newnum[i]]++] = j;
            else
                nodetab[(*rangtab)[newnum[i]]++] = i;
        }
        else
        {
            /** find the cblk this node is in **/
            father = i;
            while(colweight[father] <= 0)
            {
                father = treetab[father];
#ifdef DEBUG_KASS
                ASSERT(father > 0, MOD_KASS);
#endif
            }
            if(snodetab != NULL)
                for(j=snodetab[i];j<snodetab[i+1];j++)
                    nodetab[(*rangtab)[newnum[father]]++] = j;
            else
                nodetab[(*rangtab)[newnum[father]]++] = i;
        }
    }

    /** reset rangtab to its real value **/
    for(i=*cblknbr; i>0; i--)
        (*rangtab)[i] = (*rangtab)[i-1];
    (*rangtab)[0] = 0;





#ifdef DEBUG_KASS
    /*for(i=0;i<*cblknbr+1;i++)
     fprintf(stderr, "rangtab[%ld] = %ld \n", (long)i, (long)(*rangtab)[i]);
     exit(0);*/
    for(i=0;i<n;i++)
        if(colweight[i] > 0)
            ASSERT(colweight[i] == (*rangtab)[newnum[i]+1]-(*rangtab)[newnum[i]], MOD_KASS);


    /** check the iperm vector (nodetab) **/
    {
        pastix_int_t *flag;
        fprintf(stderr, "Cblknbr = %ld NN = %ld \n", *cblknbr, (long)nn);
        ASSERT( (*rangtab)[*cblknbr] == nn, MOD_KASS);

        MALLOC_INTERN(flag, nn, pastix_int_t);


        bzero(flag, sizeof(pastix_int_t)*nn);
        for(i=0;i<nn;i++)
        {
            ASSERT(nodetab[i] >= 0 && nodetab[i] < nn, MOD_KASS);
            flag[nodetab[i]]++;
        }
        for(i=0;i<nn;i++)
        {
            if(flag[nodetab[i]] != 1)
                fprintf(stderr, "(Nodetab[%ld] = %ld falg = %ld ) ", (long)i, (long)nodetab[i], (long)flag[nodetab[i]]);
            ASSERT(flag[nodetab[i]] == 1, MOD_KASS);
        }

        memFree(flag);
    }
#endif

    /* Compact the graph to remove merged nodes */
    kass_csrCompact( graphL );
    assert( graphL->n == *cblknbr );

    /*** Apply the new permutation to P ****/
    /** tmp is the perm vector **/
    {
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

#ifdef DEBUG_KASS
    /** Check some things about P **/
    for(i=0;i<P->n;i++)
    {
        if(P->nnzrow[i] < (*rangtab)[i+1]-(*rangtab)[i])
            for(j=0;j< P->nnzrow[i];j++)
                fprintf(stderr, "ja = %ld, rang = %ld \n", (long)P->ja[i][j],(long)(*rangtab)[i]+j);

        ASSERT(P->nnzrow[i] >= (*rangtab)[i+1]-(*rangtab)[i], MOD_KASS);
        for(j=0;j< (*rangtab)[i+1]-(*rangtab)[i];j++)
        {
            if(P->ja[i][j] != (*rangtab)[i]+j)
                fprintf(stderr, "Cblk %ld j %ld ja %ld rangtab[%ld]=%ld rangtab[%ld] = %ld \n",
                        i, j, P->ja[i][j], i, (*rangtab)[i], i+1, (*rangtab)[i+1]);
            ASSERT(P->ja[i][j] == (*rangtab)[i]+j, MOD_KASS);
        }

        /** The matrix should be still sorted **/
        for(j=1;j<P->nnzrow[i];j++)
            ASSERT(P->ja[i][j] > P->ja[i][j-1], MOD_KASS);

    }

    /*for(i=0;i<nn;i++)
     fprintf(stderr, "%ld ", nodetab[i]);
     fprintf(stderr, "\n");*/

    key = 0.0;
    for(i=0;i<P->n;i++)
        key += CBLKTIME(P->nnzrow[i], P->ja[i], (*rangtab)[i+1]-(*rangtab)[i]);
    fprintf(stderr, "COST of the AMALGAMATED MATRIX = %g \n", key);
#endif

    memFree(tmp);
    memFree(tmp2);
    memFree(sonindex);
    memFree(sontab);
    memFree(colweight);
}

