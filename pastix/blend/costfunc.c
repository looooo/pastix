#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "dof.h"
#include "cost.h"
#include "ftgt.h"
#include "symbol.h"
#include "elimin.h"
#include "perf.h"
#include "cand.h"
#include "queue.h"
#include "bulles.h"
#include "extendVector.h"
#include "blendctrl.h"
#include "updown.h"
#include "solver.h"
#include "simu.h"
#include "costfunc.h"

#define BLEND_CHOLESKY /** LLt version **/

void   subtreeSetNullCost    (pastix_int_t, const BlendCtrl * ctrl, const SymbolMatrix *, const SimuCtrl *,  pastix_int_t);
double cblkComputeCost2DLocal(pastix_int_t, const BlendCtrl * ctrl, const SymbolMatrix *, const Dof *, const SimuCtrl *);

#if 0
/*+ Summ the subtree cost local node ; do not recompute node cost +*/
double subtreeUpdateCostLocal(pastix_int_t rootnum, const BlendCtrl * ctrl, const SymbolMatrix *symbmtx,
                              const SimuCtrl *simuctrl, const Dof * dofptr,  pastix_int_t clustnum)
{

    CostMatrix *costmtx = ctrl->costmtx;
    EliminTree *etree   = ctrl->etree;
    Cand       *candtab = ctrl->candtab;
    pastix_int_t         i;

    /* Update cost of local task 1D */
    if (candtab[rootnum].distrib == D1)
    {
        /* Subtree is not local */
        if (candtab[rootnum].cluster != clustnum)
        {
            costmtx->cblktab[rootnum].total = 0.0;
        }
        /* If not, we don't touch the cost */
    }
    /* 2D */
    else
    {
        /* Update cost */
        costmtx->cblktab[rootnum].total = cblkComputeCost2DLocal(rootnum, ctrl, symbmtx, dofptr, simuctrl);
    }

    costmtx->cblktab[rootnum].subtree = costmtx->cblktab[rootnum].total;

    if ((candtab[rootnum].fccandnum <= clustnum) &&
        (candtab[rootnum].lccandnum >= clustnum))
    {
        for(i=0;i<etree->nodetab[rootnum].sonsnbr;i++)
            costmtx->cblktab[rootnum].subtree += subtreeUpdateCostLocal(eTreeSonI(etree, rootnum, i), ctrl,
                                                                        symbmtx, simuctrl, dofptr, clustnum);
    }
#ifdef DEBUG_BLEND
    else
    {
        for(i=0;i<etree->nodetab[rootnum].sonsnbr;i++)
            subtreeSetNullCost(eTreeSonI(etree, rootnum, i), ctrl,
                               symbmtx, simuctrl, clustnum);
    }
#endif

    /* Sort the sons by decreasing order */
    {
        pastix_int_t son, i, sonsnbr;
        double cumul_cost, soncost;
        Queue *queue_tree;

        sonsnbr = etree->nodetab[rootnum].sonsnbr;

        MALLOC_INTERN(queue_tree, 1, Queue);
        queueInit(queue_tree, sonsnbr);
        for(i=0;i<sonsnbr;i++)
        {
            son = eTreeSonI(etree, rootnum, i);

            /** Cost in the current subtree to be mapped **/
            cumul_cost = -ctrl->costmtx->cblktab[son].subtree;

            /* Cost of the root node in the subtree */
            soncost    = -ctrl->costmtx->cblktab[son].total;

            queueAdd2(queue_tree, son, cumul_cost, soncost);
        }

        for(i=0;i<sonsnbr;i++)
        {
            eTreeSonI(etree, rootnum, i) = queueGet(queue_tree);
        }
        queueExit(queue_tree);
        memFree(queue_tree);


        for(i=1;i<sonsnbr;i++)
        {
            assert( ctrl->costmtx->cblktab[eTreeSonI(etree, rootnum, i)].subtree
                    <= ctrl->costmtx->cblktab[eTreeSonI(etree, rootnum, i-1)].subtree );
        }
    }


    return costmtx->cblktab[rootnum].subtree;
}

void subtreeSetNullCost(pastix_int_t rootnum, const BlendCtrl * ctrl,
                        const SymbolMatrix *symbmtx, const SimuCtrl *simuctrl,
                        pastix_int_t clustnum)
{
    CostMatrix *costmtx = ctrl->costmtx;
    EliminTree *etree   = ctrl->etree;
    pastix_int_t         i;

    assert(ctrl->candtab[rootnum].cluster != clustnum);
    assert(simuctrl->bloktab[symbmtx->cblktab[rootnum].bloknum].ownerclust != clustnum );

    costmtx->cblktab[rootnum].total   = 0.0;
    costmtx->cblktab[rootnum].subtree = 0.0;
    for(i=0;i<etree->nodetab[rootnum].sonsnbr;i++)
        subtreeSetNullCost(eTreeSonI(etree, rootnum, i), ctrl, symbmtx, simuctrl, clustnum);

    return;
}
#endif

double cblkComputeCost2D(pastix_int_t cblknum, CostMatrix *costmtx, const SymbolMatrix *symbptr, const Dof * dofptr)
{
    pastix_int_t i, j;
    pastix_int_t L, h, g;
    double cost = 0.0;
    (void)costmtx;

    L    = (symbptr->cblktab[cblknum].lcolnum - symbptr->cblktab[cblknum].fcolnum + 1);
    L   *= (dofptr)->noddval;
    cost = DIAGCost(L);
    for(i=symbptr->cblktab[cblknum].bloknum+1;i<symbptr->cblktab[cblknum+1].bloknum;i++)
    {
        h = symbptr->bloktab[i].lrownum - symbptr->bloktab[i].frownum + 1;
        h *= (dofptr)->noddval;
        cost += E1Cost(L, h);
        for(j=i;j<symbptr->cblktab[cblknum+1].bloknum;j++)
        {

            g = symbptr->bloktab[j].lrownum - symbptr->bloktab[j].frownum + 1;
            g *= (dofptr)->noddval;
            cost += E2Cost(L, h, g);
#ifdef DEBUG_BLEND
            ASSERT(L > 0,MOD_BLEND);
            ASSERT(h > 0,MOD_BLEND);
            ASSERT(g > 0,MOD_BLEND);
#endif
        }
    }
#ifdef DEBUG_BLEND
    /*  ASSERT(cost >= 0,MOD_BLEND);*/
#endif
    return cost;
}

double cblkComputeCost2DLocal(pastix_int_t cblknum, const BlendCtrl * ctrl, const SymbolMatrix *symbptr,
                              const Dof * dofptr, const SimuCtrl *simuctrl)
{
    double      cost = 0.0;
    pastix_int_t         i, j;
    pastix_int_t         L, h, g;

    L  = (symbptr->cblktab[cblknum].lcolnum - symbptr->cblktab[cblknum].fcolnum + 1);
    L *= (dofptr)->noddval;

    /*  if (simuctrl->bloktab[symbptr->cblktab[cblknum].bloknum].tasknum != -1)*/
    if ( simuctrl->bloktab[symbptr->cblktab[cblknum].bloknum].ownerclust == ctrl->clustnum)
        cost = DIAGCost(L);

    for(i=symbptr->cblktab[cblknum].bloknum+1;i<symbptr->cblktab[cblknum+1].bloknum;i++)
    {
        if (simuctrl->bloktab[i].ownerclust != ctrl->clustnum)
            continue;

        h     = symbptr->bloktab[i].lrownum - symbptr->bloktab[i].frownum + 1;
        h    *= (dofptr)->noddval;
        cost += E1Cost(L, h);

        for(j=i; j<symbptr->cblktab[cblknum+1].bloknum; j++)
        {
            g     = symbptr->bloktab[j].lrownum - symbptr->bloktab[j].frownum + 1;
            g    *= (dofptr)->noddval;
            cost += E2Cost(L, h, g);
        }
    }
    return cost;
}

/*+ Compute cost of the cblk, return total cost +*/

/** Assure that cblkComputeCost and cblkCost() compute the same things !!!! **/
double cblkComputeCost(pastix_int_t cblknum, CostMatrix *costmtx, const SymbolMatrix *symbmtx, const Dof * dofptr)
{
    double contribsum;
    pastix_int_t l, h, g;
    pastix_int_t k;
    pastix_int_t noddval = symbmtx->dof;

    /** we need the height of cblk non empty lines  and the broadness
     of the cbl to compute the local compute cost **/
    l  = symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1;
    l *= noddval;

    g = 0;
    for(k = symbmtx->cblktab[cblknum].bloknum;
        k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
        g += (symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1) * noddval;
    }

    /* retrieve diag height so let g be the odb non empty lines height */
    g -= l;

    /** compute the local compute cost **/
    if( l != 0 ) {
        contribsum = computeCost(l, g);
    }
    else {
        contribsum = 0.;
    }

    costmtx->bloktab[ symbmtx->cblktab[cblknum].bloknum ].contrib = contribsum;

    /** compute for each odb its contribution compute cost and add cost **/
    for(k = symbmtx->cblktab[cblknum].bloknum+1;
        k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
#ifdef  DOF_CONSTANT
        h = (symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1)*noddval;
#endif

        /*if(l!=0 && h != 0 && g != 0)*/
        costmtx->bloktab[k].contrib     = contribCompCost(l, h, g);
        /*else
         costmtx->bloktab[k].contrib     = 0;*/
        /*if(h != 0 && g != 0)*/
        costmtx->bloktab[k].contrib    += contribAddCost(h, g);

        contribsum += costmtx->bloktab[k].contrib;
        g -= h;
    }

#ifdef DEBUG_BLEND
    {
        pastix_int_t stride=0;
        double cost2=0;

        for(k=symbmtx->cblktab[cblknum].bloknum; k<symbmtx->cblktab[cblknum+1].bloknum;k++)
            stride += symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1;

        //ASSERT(costmtx->cblktab[cblknum].total > 0,MOD_BLEND);
        cost2 = cblkCost(symbmtx->cblktab[cblknum+1].bloknum -  symbmtx->cblktab[cblknum].bloknum,
                         &(symbmtx->bloktab[symbmtx->cblktab[cblknum].bloknum]),
                         dofptr);
        /* Values should be equals but we accept the machine computational error */
        /* ASSERT(costmtx->cblktab[cblknum].total - cost2 < 10e-15, */
        /*        MOD_BLEND); */

    }
#endif

    return contribsum;
}

/** Assure that cblkComputeCost and cblkCost compute the same things !!!! **/
void cblkComputeCostLL(pastix_int_t cblknum, CostMatrix *costmtx, const SymbolMatrix *symbmtx, const Dof * dofptr)
{
    double contribsum;
    pastix_int_t fbloknum, lbloknum, fcblknum;
    pastix_int_t l, h, g;
    pastix_int_t k;
    pastix_int_t noddval = ( dofptr == NULL ) ? 1 : dofptr->noddval;

    /* Compute the height and width of the cblk */
    l =  symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1;
    l *= noddval;

    g = 0;
    fbloknum = symbmtx->cblktab[cblknum  ].bloknum;
    lbloknum = symbmtx->cblktab[cblknum+1].bloknum;
    for(k=fbloknum; k<lbloknum; k++)
    {
        g += (symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);
    }
    g *= noddval;

    /* Substract the diagonal height so let g be the off-diagonal blocks non empty lines height */
    g -= l;

    /* Compute the cost of factorizing the diagonal block and updating the panel */
    assert( l > 0 );
    contribsum = computeCost(l, g);

    costmtx->cblkcost[cblknum] += contribsum;

    /* Compute for each off-diagonal block its contribution and compute cost and add cost **/
    for(k=fbloknum+1; k<lbloknum; k++)
    {
        h  = (symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);
        h *= noddval;

        fcblknum = symbmtx->bloktab[k].cblknum;

        /* g is the number of off-diagonal blocks in the panel below k (k included) */
        costmtx->cblkcost[fcblknum]  = contribCompCost(l, h, g);
        costmtx->cblkcost[fcblknum] += contribAddCost(h, g);

        g -= h;
    }

    return;
}


/*****************************************************************************************
 *      There are the cost functions of the compute phase of factorization algorithm     *
 *****************************************************************************************/
#ifndef BLEND_CHOLESKY
double computeCost(pastix_int_t L, pastix_int_t g_total)
{
    double total = 0;
    total =(double)(L*PERF_COPY(L)+ PERF_PPF(L) + PERF_TRSM(L, g_total) + L*PERF_SCAL(g_total)
                    + L*PERF_COPY(g_total));
    return (total>0)?total:0;
}
#else
double computeCost(pastix_int_t L, pastix_int_t g_total)
{
    double total = 0;
    total =(double)(PERF_POF(L) + PERF_TRSM(L, g_total)) ;
    return (total>0)?total:0;
}
#endif


double contribCompCost(pastix_int_t L, pastix_int_t h, pastix_int_t g)
{
    double total = 0;
    assert(L> 0);
    assert(h>=0);
    assert(g>=0);

    total = (double)(PERF_GEMM(g,h,L));
    return (total>0)?total:0;
}

double contribAddCost(pastix_int_t h, pastix_int_t g)
{
    double total = 0;
#ifdef DEBUG_BLEND
    ASSERT(h>=0,MOD_BLEND);
    ASSERT(g>0,MOD_BLEND);
#endif

    total = (double)(PERF_GEAM(g, h));

    return (total>0)?total:0;
}

/**********************************************/
/*     Pour le 2D                             */
/**********************************************/
#ifndef BLEND_CHOLESKY
double DIAGCost(pastix_int_t L)
{
    return (double)(L*PERF_COPY(L)+ PERF_PPF(L));
}
double E1Cost(pastix_int_t L, pastix_int_t g)
{
    return (double)(PERF_TRSM(L, g) + L*PERF_SCAL(g)
                    + L*PERF_COPY(g));
}
double E2Cost(pastix_int_t L, pastix_int_t h, pastix_int_t g)
{
    return (double)(PERF_GEMM(g,h,L) + PERF_GEAM(g, h));
}
#else
double DIAGCost(pastix_int_t L)
{
    return (double)(PERF_POF(L));
}
double E1Cost(pastix_int_t L, pastix_int_t g)
{
    return (double)(PERF_TRSM(L, g));
}
double E2Cost(pastix_int_t L, pastix_int_t h, pastix_int_t g)
{
    return (double)(PERF_GEMM(g,h,L));
}
#endif
/*****************************************************************************************
 *                  END of cost functions                                                *
 *****************************************************************************************/

#if 0
double cblkMaxCost(pastix_int_t cblknbr, const CostMatrix *costmtx)
{
    pastix_int_t i;
    double maxcost;
    maxcost = 0;
    for(i=0;i< cblknbr;i++)
        if(costmtx->cblktab[i].total > maxcost)
            maxcost = costmtx->cblktab[i].total;
    return maxcost;
}



double totalCost(pastix_int_t cblknbr, const CostMatrix *costmtx)
{
    pastix_int_t i;
    double total=0;
    for(i=0;i<cblknbr;i++)
        total += costmtx->cblktab[i].total;
    return total;
}

#endif
double memorySpaceCost(const SolverMatrix *solvmtx)
{
    double space=0;
    space += solverSpaceCost(solvmtx);
    return space;
}


double solverSpaceCost(const SolverMatrix *solvmtx)
{
    double space=0;
    /*pastix_int_t i;*/
    space += sizeof(double)*(solvmtx->coefnbr);

    space += sizeof(SolverCblk)*(solvmtx->cblknbr);
    space += sizeof(SolverBlok)*(solvmtx->bloknbr);
    space += sizeof(FanInTarget)*(solvmtx->ftgtnbr);
    space += sizeof(pastix_int_t)*MAXINFO*(solvmtx->ftgtnbr);
    /*  for(i=0;i<solvmtx->ftgtnbr;i++)
     space += sizeof(double) * solvmtx->ftgttab[i].indtab[solvmtx->ftgttab[i].infotab[FTGT_BLOKNBR]]
     * (solvmtx->ftgttab[i].infotab[FTGT_LCOLNUM] -solvmtx->ftgttab[i].infotab[FTGT_FCOLNUM] + 1) ;
     */
    return space;
}

double symbolSpaceCost(const SymbolMatrix *symbmtx)
{
    double space=0;
    space += sizeof(SymbolCblk)*(symbmtx->cblknbr+1);
    space += sizeof(SymbolBlok)*(symbmtx->bloknbr+1);
    return space;
}


void printSolverInfo(FILE *out, const SolverMatrix * solvmtx)
{
    pastix_int_t procnum = solvmtx->clustnum;
    double totalspace = memorySpaceCost(solvmtx);
    fprintf(out,   " %ld : Number of Non Null Coeff         : %ld --> %g Octs \n",
            (long)procnum, (long)solvmtx->coefnbr, (double)solvmtx->coefnbr*sizeof(double));
    fprintf(out,   " %ld : ExtraStructure Memory Space      : %g Octs \n", (long)procnum, totalspace - (double)solvmtx->coefnbr*sizeof(double));
    fprintf(out,   " %ld : Total Memory space               : %g  Octs\n", (long)procnum, totalspace);
}



double cblk_time_fact(pastix_int_t n, pastix_int_t *ja, pastix_int_t colnbr)
{
    /*******************************************/
    /* Compute the time to compute a cblk      */
    /* according to the BLAS modelization      */
    /*******************************************/
    double cost;
    pastix_int_t i;
    pastix_int_t L, G, H;

    /** The formula are based on the costfunc.c in blend **/
    /** @@@ OIMBE: il faudra faire les DOF_CONSTANT ***/

    L = colnbr; /* Size of diagonal block        */
    G = n-L;    /* Number of extra-diagonal rows */

#define CHOLESKY
#ifndef CHOLESKY
    cost = (double)(L*PERF_COPY(L)+ PERF_PPF(L) + PERF_TRSM(L, G) + L*PERF_SCAL(G)
                    + L*PERF_COPY(G));
#else
    cost = (double)( PERF_POF(L) + PERF_TRSM(L, G) );
#endif

    /** Contributions **/
    i = colnbr;
    while(i < n)
    {
        H = 1;
        i++;
        while((i<n) && (ja[i] == ja[i-1]+1))
        {
            i++;
            H++;
        }
        cost += (double)(PERF_GEMM(G, H, L));
        G -= H;
    }
    return cost;
}


double cblk_time_solve(pastix_int_t n, pastix_int_t *ja, pastix_int_t colnbr)
{
    double cost;
    pastix_int_t L;
    (void)ja;

    L = colnbr;

    cost = (double)PERF_TRSV(L) + (double) PERF_GEMV(L, n-L);
    return cost;
}
