#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "common.h"
#include "queue.h"
#include "bulles.h"

#define NBLEAVES btree->leavesnbr
#define NBBUBMAX btree->nodemax

/*
 Function: bubbleInitTree

 Initialize a tree of bubbles

 Parameters:
 btree        - Pointer to bubble tree to allocate and initialize.
 nbleaves     - Number of leaves in the tree.

 Returns:

 */
void bubbleInitTree(BubbleTree * btree, int nbleaves)
{

    int i,n;
    int nblevel=2;
    int nbbubblesmax;

    n=nbleaves;
    while(n >>= 1)
        nblevel++;

    nbbubblesmax = (1 << nblevel)-1;
    nbbubblesmax *= 2;
    print_debug(DBG_BUBBLES, "Bulles : nbleaves = %d, nblevel = %d et nbbubblesmax = %d\n",
                nbleaves, nblevel, nbbubblesmax);

    btree->leavesnbr = nbleaves;
    btree->nodenbr   = nbleaves;
    btree->nodemax   = nbbubblesmax;
    MALLOC_INTERN(btree->sonstab, nbbubblesmax, int);
    MALLOC_INTERN(btree->nodetab, nbbubblesmax, BubbleTreeNode);

    for(i=0; i<NBLEAVES; i++)
    {
        btree->sonstab[i] = -1;
        btree->nodetab[i].fathnum   = -1;
        btree->nodetab[i].sonsnbr   = 0;
        btree->nodetab[i].fsonnum   = 0;
        btree->nodetab[i].fcandnum  = i;
        btree->nodetab[i].lcandnum  = i;
        btree->nodetab[i].localcost = 0.0;
        btree->nodetab[i].costlevel = 0.0;
        btree->nodetab[i].treelevel = 0;
        btree->nodetab[i].nbtasks   = 0;
        btree->nodetab[i].priomin   = INTVALMAX;
        btree->nodetab[i].priomax   = 0;
        btree->nodetab[i].treelevel = 0;
        MALLOC_INTERN(btree->nodetab[i].taskheap, 1, Queue);
        queueInit(btree->nodetab[i].taskheap, 1000);
    }

    for(i=NBLEAVES; i<NBBUBMAX; i++)
    {
        btree->sonstab[i] = -1;
        btree->nodetab[i].fathnum   = -1;
        btree->nodetab[i].sonsnbr   = 0;
        btree->nodetab[i].fsonnum   = 0;
        btree->nodetab[i].fcandnum  = -1;
        btree->nodetab[i].lcandnum  = -1;
        btree->nodetab[i].localcost = 0.0;
        btree->nodetab[i].costlevel = 0.0;
        btree->nodetab[i].treelevel = 0;
        btree->nodetab[i].nbtasks   = 0;
        btree->nodetab[i].priomin   = INTVALMAX;
        btree->nodetab[i].priomax   = 0;
        MALLOC_INTERN(btree->nodetab[i].taskheap, 1, Queue);
        queueInit(btree->nodetab[i].taskheap, 100);
    }
}

/*
 Function: bubbleFree

 Free the BubbleTree structure.

 Parameters:
 btree        - Pointer to the bubble tree to be free, NULL is returned.

 Returns:

 */
void bubbleFree(BubbleTree *btree)
{
    int i;

    for (i=0; i<NBBUBMAX; i++)
    {
        queueExit(btree->nodetab[i].taskheap);
        memFree_null(btree->nodetab[i].taskheap);
    }

    memFree_null(btree->nodetab);
    memFree_null(btree->sonstab);
}

/*
 Function: bubbleAdd

 Add a task in the BubbleTree structure defined by the set of candidates
 processors (first and last) and by the cost and the deep of this task in the
 dependence's tree.

 Parameters:
 btree        - Pointer to Bubble tree in which we need to add some task.
 fcandnum     - Indice of the first candidate for the task to add.
 lcandnum     - Indice of the last candidate for the task to add.
 costlevel    - Cost of the path from the root to this task.
 treelevel    - Depth of the task in the task's tree.

 Returns:
 The indice of the bubble corresponding to the set of candidates.
 Exit if the search or the add of the bubble fails.

 */
int
bubbleAdd(BubbleTree * btree,
          pastix_int_t fathnum,
          pastix_int_t fcandnum,
          pastix_int_t lcandnum,
          double       localcost,
          double       costlevel,
          pastix_int_t treelevel)
{
    int i;

    /* It's a leaf */
    if( fcandnum == lcandnum ) {

#ifdef BUBBLE_DEBUG
        errorPrintW("bubbleAdd : Ajout bulle processeur %d, costlevel %f, treelevel %d",
                    fcandnum, costlevel, treelevel);
#endif
        if ((btree->nodetab[fcandnum].fathnum == -1) ||
            (btree->nodetab[fcandnum].fathnum == fathnum))
        {
            btree->nodetab[fcandnum].fathnum    = fathnum;
            btree->nodetab[fcandnum].localcost += localcost;
            btree->nodetab[fcandnum].costlevel  = pastix_imin(btree->nodetab[fcandnum].costlevel, -costlevel);
            btree->nodetab[fcandnum].treelevel  = pastix_imin(btree->nodetab[fcandnum].treelevel, -treelevel);
            btree->nodetab[fcandnum].nbtasks ++;
            return fcandnum;
        }

        /* If we are in crossproc mode, a leaf can be duplicated in multiple branches,
         * so we handle other occurences as new nodes */
    }

    for (i=NBLEAVES; i<NBBUBMAX; i++)
    {
        /* The node already exists */
        if ((btree->nodetab[i].fcandnum == fcandnum) &&
            (btree->nodetab[i].lcandnum == lcandnum) &&
            (btree->nodetab[i].fathnum == fathnum || fathnum == -1))
        {
            btree->nodetab[i].localcost += localcost;
            btree->nodetab[i].costlevel  = pastix_imin(btree->nodetab[i].costlevel, -costlevel);
            btree->nodetab[i].treelevel  = pastix_imin(btree->nodetab[i].treelevel, -treelevel);
            btree->nodetab[i].nbtasks ++;
            return i;
        }
        /* The node doesn't exist, it is created */
        else if (btree->nodetab[i].fcandnum == -1)
        {
            btree->nodetab[i].fathnum   =  fathnum;
            btree->nodetab[i].fcandnum  =  fcandnum;
            btree->nodetab[i].lcandnum  =  lcandnum;
            btree->nodetab[i].localcost =  localcost;
            btree->nodetab[i].costlevel = -costlevel;
            btree->nodetab[i].treelevel = -treelevel;
            btree->nodetab[i].nbtasks ++;
            btree->nodenbr++;
            return i;
        }
    }

    /* Don't never be here */
    errorPrint("La bulle [%ld %ld] n'a pas pu etre ajoutee\n",
               (long)fcandnum, (long)lcandnum);
    assert(0);
    return -1;
}


/*
 Function: bubblebuildtree

 Construct the tree by filling nodetab and sonstab tables thanks to
 the bubbles added whith bubbleAdd. This function must be call only
 one time for one bubble tree. And only after the add of all bubbles.

 Parameters:
 btree        - Pointer to Bubble tree
 verbose      - If verbose, the tree is printed

 Returns:
 void
 */
void bubbleBuildTree(const BubbleTree * btree){

    int i, j;
    int bubsize;
    int fprocnum, lprocnum, fathnum;
    int maxlevel    = 0;
    int *sonsnbrtmp = NULL;

    for (i=0; i<btree->nodenbr; i++)
    {
        bubsize  = NBLEAVES;
        fprocnum = btree->nodetab[i].fcandnum;
        lprocnum = btree->nodetab[i].lcandnum;

        if (btree->nodetab[i].fathnum == -1) {
            for(j=NBLEAVES; j<btree->nodenbr; j++)
            {
                /* The father node is the smallest bubble which contains it */
                if ((i != j) && (fprocnum >= btree->nodetab[j].fcandnum)
                    && (lprocnum <= btree->nodetab[j].lcandnum)
                    && (bubsize > (btree->nodetab[j].lcandnum - btree->nodetab[j].fcandnum)))
                {
                    bubsize  = btree->nodetab[j].lcandnum - btree->nodetab[j].fcandnum;
                    btree->nodetab[i].fathnum = j;
                }
            }
        }

        if (btree->nodetab[i].fathnum != -1)
            btree->nodetab[btree->nodetab[i].fathnum].sonsnbr++;
    }

    for(i=NBLEAVES; i<btree->nodenbr; i++)
        maxlevel = MAX(maxlevel, btree->nodetab[i].treelevel);
    for(i=NBLEAVES; i<btree->nodenbr; i++)
        btree->nodetab[i].treelevel = maxlevel + 1 - btree->nodetab[i].treelevel;

    /* Compute the number of son and indice of first son */
    MALLOC_INTERN(sonsnbrtmp, btree->nodenbr, int);
    sonsnbrtmp[0] = 0;
    for (i=1; i<btree->nodenbr; i++)
    {
        btree->nodetab[i].fsonnum = btree->nodetab[i-1].fsonnum + btree->nodetab[i-1].sonsnbr;
        sonsnbrtmp[i] = 0;
    }

    /* Fill in sonstab */
    for (i=0; i<btree->nodenbr; i++)
    {
        fathnum = btree->nodetab[i].fathnum;
        if (fathnum != -1)
        {
            btree->sonstab[btree->nodetab[fathnum].fsonnum + sonsnbrtmp[fathnum]] = i;
            sonsnbrtmp[fathnum]++;
            ASSERTDBG(sonsnbrtmp[fathnum] <= btree->nodetab[fathnum].sonsnbr, MOD_BLEND);
        }
    }
    memFree_null(sonsnbrtmp);

    return;
}


/*
 Function: bubbleprint

 Generate a dot file for graphviz with the bubble tree.

 Parameters:
 btree        - Pointer to Bubble tree
 bcost        - Cost in time of all tasks included in each bubble

 Returns:
 void
 */
void bubblePrint(const BubbleTree * btree,
                 FILE *out)
{
    pastix_int_t i;
    pastix_int_t fprocnum, lprocnum;
    double percent, totalcost;
    double *subtrees_costs   = (double*)malloc( btree->nodenbr * sizeof(double) );
    int    *subtrees_nbtasks = (int*)   malloc( btree->nodenbr * sizeof(int) );
    memset( subtrees_costs,   0, btree->nodenbr * sizeof(double) );
    memset( subtrees_nbtasks, 0, btree->nodenbr * sizeof(int) );

    /* const double     * bcost, */
    /* double totalcost, */

    fprintf(out,"digraph G {\n"
            "\tcolor=white\n"
            "\trankdir=BT;\n");

    totalcost = 0.;
    for (i=0; i<btree->nodenbr; i++) {
        double localcost = btree->nodetab[i].localcost;
        int    localtask = btree->nodetab[i].nbtasks;
        int    father    = btree->nodetab[i].fathnum;

        subtrees_costs[i] += localcost;
        totalcost         += localcost;
        subtrees_nbtasks[i] += localtask;

        while ( father != -1 ) {
            subtrees_costs[father]   += localcost;
            subtrees_nbtasks[father] += localtask;
            father = btree->nodetab[father].fathnum;
        }
    }

    for (i=0; i<btree->nodenbr; i++)
    {
        int    father = btree->nodetab[i].fathnum;
        double subptotal  = (subtrees_costs[i] / totalcost) * 100.;
        double subpparent = (subtrees_costs[i] / (subtrees_costs[father] - btree->nodetab[father].localcost)) * 100.;

        fprocnum = btree->nodetab[i].fcandnum;
        lprocnum = btree->nodetab[i].lcandnum;

        percent = subtrees_costs[i] / totalcost;

        if (btree->nodetab[i].fathnum == -1)
            fprintf(out, "\t%ld [style=filled, label=\"%ld-%ld\\n#tsk: %d (%d)\\n"
                    "Local Time: %lf, %3.0lf%%\\nSubtree Time: %lf, %3.0lf%% (%3.0lf%%)\",color=\"/set39/%d\"]\n",
                    (long)i, (long)fprocnum, (long)lprocnum,
                    btree->nodetab[i].nbtasks,
                    subtrees_nbtasks[i],
                    btree->nodetab[i].localcost,
                    btree->nodetab[i].localcost / totalcost,
                    subtrees_costs[i], subptotal, subpparent,
                    ( percent > ( (double)( lprocnum-fprocnum+1 ) / (double)(btree->leavesnbr) ) ) ? 4 : 7 );
        else
            fprintf(out, "\t%ld -> %ld\n\t%ld [style=filled, label=\"%ld-%ld\\n#tsk : %d (%d)\\n"
                    "Local Time: %lf, %3.0lf%%\\nSubtree Time: %lf, %3.0lf%% (%3.0lf%%)\",color=\"/set39/%d\"]\n",
                    (long)i, (long)btree->nodetab[i].fathnum,
                    (long)i, (long)fprocnum, (long)lprocnum,
                    btree->nodetab[i].nbtasks,
                    subtrees_nbtasks[i],
                    btree->nodetab[i].localcost,
                    (btree->nodetab[i].localcost / totalcost) * 100.,
                    subtrees_costs[i], subptotal, subpparent,
                    ( percent > ( (double)( lprocnum-fprocnum+1 ) / (double)(btree->leavesnbr) ) ) ? 4 : 7 );
    }

    fprintf(out, "}\n");

    free(subtrees_costs);
}
