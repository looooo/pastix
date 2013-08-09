/**
 *
 * @file elimin_tree.c
 *
 *  PaStiX analyse routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains basic functions to manipulate elimination tree structure.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "symbol.h"
#include "elimin.h"
#include "cost.h"
#include "cand.h"
#include "dof.h"
#include "extrastruct.h"
#include "param_blend.h"
#include "extendVector.h"
#include "blendctrl.h"
#include "splitpart.h"
#include "queue.h"

#define CROSS_TOLERANCE 0.1

double  subtreeUpdateCost     (pastix_int_t, CostMatrix *, const EliminTree *);

typedef struct propmap_s {
    BlendCtrl         *ctrl;
    EliminTree        *etree;
    Cand              *candtab;
    CostMatrix        *costmtx;
    SymbolMatrix      *symbmtx;
    const Dof         *dofptr;
    ExtraSymbolMatrix *extrasymb;
    ExtraCostMatrix   *extracost;
    pastix_int_t       procnbr;
    int                split;
    int                nocrossproc;
    int                allcand;
} propmap_t;

static inline void
propMappSubtreeOn1P( propmap_t   *pmptr,
                     pastix_int_t rootnum,
                     pastix_int_t fcandnum,
                     pastix_int_t lcandnum,
                     pastix_int_t cluster )
{
    pastix_int_t i;
    pastix_int_t sonsnbr;

    pmptr->candtab[rootnum].fcandnum = fcandnum;
    pmptr->candtab[rootnum].lcandnum = lcandnum;
    pmptr->candtab[rootnum].cluster  = cluster;

    /* Split the treenode */
    if ( pmptr->split ) {
        splitOnProcs( pmptr->symbmtx,
                      pmptr->extrasymb,
                      pmptr->extracost,
                      pmptr->ctrl,
                      pmptr->dofptr,
                      rootnum, 1);

        /* Correct the subtree cost */
        subtreeUpdateCost(rootnum, pmptr->costmtx, pmptr->etree);
    }

    /* Recursively apply the affectation to the sons */
    sonsnbr = pmptr->etree->nodetab[rootnum].sonsnbr;
    for(i=0;i<sonsnbr;i++)
        propMappSubtreeOn1P( pmptr, eTreeSonI(pmptr->etree, rootnum, i),
                             fcandnum, lcandnum, cluster );

    return;
}

static inline void
propMappSubtree( propmap_t    *pmptr,
                 pastix_int_t  rootnum,
                 pastix_int_t  fcandnum,
                 pastix_int_t  lcandnum,
                 pastix_int_t  cluster,
                 double       *cost_remain)
{
    Queue *queue_tree;
    pastix_int_t p;
    pastix_int_t candnbr, scandnbr;
    pastix_int_t fcand = 0;
    pastix_int_t lcand = 0;
    pastix_int_t i;
    pastix_int_t sonsnbr;
    double isocost;
    double aspt_cost;
    double cumul_cost;
    double *sub_cost_remain = NULL;
    double epsilon;

    candnbr = lcandnum - fcandnum + 1;

    /* Si il n'y a qu'un candidat, tous le sous-arbre lui appartient */
    if (candnbr == 1)
    {
        memFree_null(cost_remain);
        propMappSubtreeOn1P( pmptr, rootnum, fcandnum, lcandnum, cluster );
        return;
    }

    /* Set the cand group for this tree node */
    if(pmptr->allcand)
    {
        pmptr->candtab[rootnum].fcandnum = 0;
        pmptr->candtab[rootnum].lcandnum = pmptr->procnbr - 1;
        pmptr->candtab[rootnum].cluster  = cluster;
    }
    else
    {
        pmptr->candtab[rootnum].fcandnum = fcandnum;
        pmptr->candtab[rootnum].lcandnum = lcandnum;
        pmptr->candtab[rootnum].cluster  = cluster;

    }

    /* This treenode is a leave, split it and return */
    if(pmptr->etree->nodetab[rootnum].sonsnbr == 0)
    {
        memFree_null(cost_remain);
        if ( pmptr->split ) {
            splitOnProcs( pmptr->symbmtx,
                          pmptr->extrasymb,
                          pmptr->extracost,
                          pmptr->ctrl,
                          pmptr->dofptr,
                          rootnum, candnbr );
        }
        return;
    }

    /* Work that each processor is intended to get from this treenode */
    isocost = pmptr->costmtx->cblktab[rootnum].total / candnbr;
    for(p=0;p<candnbr;p++)
        cost_remain[p] -= isocost;

    /* Split the treenode */
    if ( pmptr->split ) {
        splitOnProcs( pmptr->symbmtx,
                      pmptr->extrasymb,
                      pmptr->extracost,
                      pmptr->ctrl,
                      pmptr->dofptr,
                      rootnum, candnbr );

        /* Correct the subtree cost */
        subtreeUpdateCost(rootnum, pmptr->costmtx, pmptr->etree);
    }

    /* Get after split cost remaining in the descendance of the treenode */
    aspt_cost = pmptr->costmtx->cblktab[rootnum].subtree - pmptr->costmtx->cblktab[rootnum].total;

    /*
     * If the first and last candidate have already received more work that they
     * should, we remove them from the set
     */
    if(cost_remain[0] <= 0)
        fcand = 1;
    else
        fcand = 0;
    if (cost_remain[candnbr-1] <= 0)
        candnbr--;

    assert(fcand < candnbr);
    assert(candnbr > 1);

    /* Make sure that the sum of cost_remain in used proc is at least equals to after split cost */
    cumul_cost = 0;
    for(i=fcand; i<candnbr; i++)
        cumul_cost += cost_remain[i];

    if (cumul_cost > aspt_cost) {
        double ratio = aspt_cost / cumul_cost;
        for(i=fcand; i<candnbr; i++)
            cost_remain[i] *= ratio;
    }

    /* Compute the minimun participation rate of a candidat processor */
    epsilon = (CROSS_TOLERANCE * cumul_cost ) / (double)candnbr;

#ifdef DEBUG_BLEND_SPLIT
    {
        fprintf(stdout, "Rootnum %ld fproc %ld lproc %ld [ ", (long)rootnum, (long)fcandnum, (long)lcandnum);
        for(p=0;p<candnbr;p++)
            fprintf(stdout, "%g ", cost_remain[p]);
        fprintf(stdout, " ]\n");
        fprintf(stdout, " Sons ");
        for(i=0;i<pmptr->etree->nodetab[rootnum].sonsnbr;i++)
            fprintf(stdout, " [%ld, %g] ", (long)i, pmptr->costmtx->cblktab[eTreeSonI(pmptr->etree, rootnum, i)].subtree);
        fprintf(stdout, "\n");
    }
#endif

    /*
     * Compute the cand group for each proc
     */
    sonsnbr = pmptr->etree->nodetab[rootnum].sonsnbr;

    /* Create the list of sons sorted by descending order of cost */
    MALLOC_INTERN(queue_tree, 1, Queue);
    queueInit(queue_tree, sonsnbr);
    for(i=0; i<sonsnbr; i++)
    {
        double soncost;

        /* Cost in the current subtree to be mapped */
        cumul_cost = -pmptr->costmtx->cblktab[eTreeSonI(pmptr->etree, rootnum, i)].subtree;

        /* Cost of the root node in the subtree */
        soncost    = -pmptr->costmtx->cblktab[eTreeSonI(pmptr->etree, rootnum, i)].total;

        queueAdd2(queue_tree, i, cumul_cost, (pastix_int_t)soncost);
    }

    /* Proportionnal mapping of the subtree on remaining candidates           */
    /* The first stage deals only with nodes that require multiple candidates */
    lcand = fcand;
    while (queueSize(queue_tree) > 0)
    {
        i = queueGet2( queue_tree, &cumul_cost, NULL );
        cumul_cost = -cumul_cost;

        /*
         * If the subtree cost is less than what the first candidate can
         * process, folowing nodes will be smaller and following candidates have
         * at least the same amount of ressources available so, we skip to the
         * second stage
         */
        if (cumul_cost <= cost_remain[fcand]) {
            cost_remain[fcand] -= cumul_cost;
            propMappSubtreeOn1P( pmptr, eTreeSonI(pmptr->etree, rootnum, i),
                                 fcandnum+fcand, fcandnum+fcand, cluster);
            break;
        }

        /*
         * If a candidate cannot participate to multiple subtree, we check with
         * epsilon to avoid overflow
         */
        if ( pmptr->nocrossproc &&
             (cumul_cost <= (cost_remain[fcand]+epsilon)) )
        {
            cost_remain[fcand] -= cumul_cost;
            propMappSubtreeOn1P( pmptr, eTreeSonI(pmptr->etree, rootnum, i),
                                 fcandnum+fcand, fcandnum+fcand, cluster);
            fcand++;
            continue;
        }

        /* If the first candidate doesn't have enough ressources, we skip it */
        if( (cost_remain[lcand] <= epsilon) && (lcand < candnbr-1) )
        {
            fcand++;
        }

        /*
         * Computes how many candidate will participate to this node.  We add
         * candidate as long as we have some and they all have more than espilon
         * extra work.
         */
        lcand = fcand;
        scandnbr = 1;
        cumul_cost -= cost_remain[fcand];
        while ((cumul_cost > (scandnbr * epsilon)) &&
               (lcand < candnbr - 1))
        {
            lcand++; scandnbr++;
            assert( cost_remain[lcand] > 0 );
            cumul_cost -= cost_remain[lcand];
        }

        /*
         * Prepare the sub_cost_remain array.
         * If some cost is remaining, we distribute equally the work to each
         * candidate, otherwise we give the cost remaining to the last candidate
         * and set to zero the remaining cost of the first candidates.
         */
        MALLOC_INTERN(sub_cost_remain, lcand-fcand+1, double);
        if (cumul_cost > 0)
        {
            isocost = cumul_cost / scandnbr;
            for(p=0; p<scandnbr; p++)
            {
                sub_cost_remain[p]   = cost_remain[fcand+p] + isocost;
                cost_remain[fcand+p] = - isocost;
            }
        }
        else
        {
            for(p=0; p<scandnbr-1; p++)
            {
                sub_cost_remain[p] = cost_remain[fcand+p];
                cost_remain[fcand+p] = 0.0;
            }
            sub_cost_remain[scandnbr-1] = cost_remain[lcand];
            cost_remain[fcand+scandnbr-1] += cumul_cost;  /* cumul_cost <= 0 */
        }

        /* Go on to subtree */
        propMappSubtree( pmptr, eTreeSonI(pmptr->etree, rootnum, i),
                         fcandnum+fcand, fcandnum+lcand, cluster, sub_cost_remain);

        if ( (lcand < candnbr - 1) &&
             ( pmptr->nocrossproc ||
               (cost_remain[lcand] < epsilon) ) )
        {
            fcand = lcand+1;
        }
        else
        {
            if ( pmptr->nocrossproc )
                break;
            fcand = lcand;
        }
    }

    if (queueSize(queue_tree) > 0)
    {
        Queue *queue_proc;

        /*
         * Second stage:
         * Distribute the single candidate subtree. At each step of the algorithm we
         * associate together the candidate with maximum cost_remain and the largest
         * son.
         */

        /* Fill queue proc order by remain cost descending */
        MALLOC_INTERN(queue_proc, 1, Queue);
        queueInit(queue_proc, candnbr);
        for (i=0; i<candnbr; i++)
            queueAdd(queue_proc, i, -cost_remain[i]);

        while (queueSize(queue_tree) > 0)
        {
            /* Get the largest node */
            i = queueGet2( queue_tree, &cumul_cost, NULL );
            cumul_cost = -cumul_cost;

            /* Get the candidate with the largest cost_remain */
            fcand = queueGet(queue_proc);

            /* Map them together */
            propMappSubtreeOn1P( pmptr, eTreeSonI(pmptr->etree, rootnum, i),
                                 fcandnum+fcand, fcandnum+fcand, cluster);

            /* Update cost_remain and re-insert into the sorted queue */
            cost_remain[fcand] -= cumul_cost;
            queueAdd(queue_proc, fcand, -cost_remain[fcand]);
        }

        queueExit(queue_proc);
        memFree(queue_proc);
    }


    queueExit(queue_tree);
    memFree(queue_tree);
    memFree_null(cost_remain);
    return;
}

void
propMappTree( BlendCtrl         *ctrl,
              EliminTree        *etree,
              Cand              *candtab,
              CostMatrix        *costmtx,
              SymbolMatrix      *symbmtx,
              const Dof         *dofptr,
              ExtraSymbolMatrix *extrasymb,
              ExtraCostMatrix   *extracost,
              pastix_int_t       procnbr,
              int split, int nocrossproc, int allcand )
{
    propmap_t pmdata;
    pastix_int_t p;
    double *cost_remain = NULL;
    double isocost;

    /* Prepare the initial cost_remain array */
    MALLOC_INTERN(cost_remain, procnbr, double);
    isocost = costmtx->cblktab[ eTreeRoot(etree) ].subtree / procnbr;

    for(p=0; p<procnbr; p++)
        cost_remain[p] = isocost;

    /* Prepare the stucture */
    pmdata.ctrl        = ctrl;
    pmdata.etree       = etree;
    pmdata.candtab     = candtab;
    pmdata.costmtx     = costmtx;
    pmdata.symbmtx     = symbmtx;
    pmdata.dofptr      = dofptr;
    pmdata.extrasymb   = extrasymb;
    pmdata.extracost   = extracost;
    pmdata.procnbr     = procnbr;
    pmdata.split       = split;
    pmdata.nocrossproc = nocrossproc;
    pmdata.allcand     = allcand;

    propMappSubtree( &pmdata, eTreeRoot(etree),
                     0, procnbr-1,
                     NOCLUSTER, cost_remain);
}
