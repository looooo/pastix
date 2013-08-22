#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "common.h"
#include "dof.h"
#include "cost.h"
#include "symbol.h"
#include "elimin.h"
#include "extrastruct.h"
#include "extendVector.h"
#include "cand.h"
#include "queue.h"
#include "bulles.h"
#include "blendctrl.h"
#include "ftgt.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "simu.h"
#include "costfunc.h"
#include "partbuild.h"
#include "splitpart.h"

/*
  Function: splitPart

  Repartitioning of the initial symbolic factorization
  and processing of candidate processors group for
  each colum bloc

  Parameters:
    symbmtx - Symbolic matrix.
    ctrl    -
    dofptr  -
*/
void splitPart(SymbolMatrix *symbmtx,
               BlendCtrl    *ctrl,
               const Dof    *dofptr)
{
    pastix_int_t i;
    ExtraSymbolMatrix *extrasymb;
    ExtraCostMatrix   *extracost;

    /* Initialization of the extra structures */
    {
        MALLOC_INTERN(extrasymb, 1, ExtraSymbolMatrix);
        extrasymbolInit(extrasymb);

        /* initialize spt[tab] */
        MALLOC_INTERN(extrasymb->sptcblk,      symbmtx->cblknbr, pastix_int_t);
        MALLOC_INTERN(extrasymb->sptcbnb,      symbmtx->cblknbr, pastix_int_t);
        MALLOC_INTERN(extrasymb->subtreeblnbr, symbmtx->cblknbr, pastix_int_t);
        MALLOC_INTERN(extrasymb->sptblok,      symbmtx->bloknbr, pastix_int_t);
        MALLOC_INTERN(extrasymb->sptblnb,      symbmtx->bloknbr, pastix_int_t);

        /* set spt* array to -1 ,
         positive or null value means cblk/blok has been splitted
         the value is, in this case the index in the extra- cblk/blok -tab */
        for(i=0;i<symbmtx->cblknbr;i++)
        {
            extrasymb->sptcblk[i] = -1;
            extrasymb->sptcbnb[i] =  1;
        }

        for(i=0;i<symbmtx->bloknbr;i++)
            extrasymb->sptblok[i] = -1;
        bzero(extrasymb->sptblnb, symbmtx->bloknbr * sizeof(pastix_int_t));

        extrasymb->curcblk = 0;
        extrasymb->curblok = 0;

        /* We choose an arbitrary size for initial allocation of bloktab and cblktab */
        MALLOC_INTERN(extrasymb->cblktab, symbmtx->cblknbr/3 + 1, SymbolCblk);
        extrasymb->sizcblk = symbmtx->cblknbr/3 + 1;
        MALLOC_INTERN(extrasymb->bloktab, symbmtx->bloknbr/3 + 1, SymbolBlok);
        extrasymb->sizblok = symbmtx->bloknbr/3 + 1;

        MALLOC_INTERN(extracost, 1, ExtraCostMatrix);
        extracostInit(extracost);

        MALLOC_INTERN(extracost->bloktab, symbmtx->bloknbr/3 + 1, CostBlok);
    }

    /* Stupid split */
    {
        pastix_int_t cblknum;
        for(cblknum = 0; cblknum < symbmtx->cblknbr; cblknum++)
        {
            splitOnProcs( symbmtx, extrasymb, extracost, ctrl,
                          dofptr, cblknum,
                          ctrl->candtab[ cblknum ].lcandnum -
                          ctrl->candtab[ cblknum ].fcandnum + 1);
        }
    }

    /* Rebuild the symbolic matrix */
    partBuild(ctrl, dofptr, extrasymb, extracost,
              symbmtx, ctrl->costmtx, &(ctrl->candtab));

    extrasymbolExit(extrasymb);
    extracostExit(extracost);

    /*********************************/
    /* Restore the elimination tree **/
    /*********************************/
    {
        eTreeExit(ctrl->etree);
        MALLOC_INTERN(ctrl->etree, 1, EliminTree);
        eTreeInit(ctrl->etree);
        eTreeBuild(ctrl->etree, symbmtx);
    }
    candBuild( ctrl->autolevel,
               ctrl->level2D,
               ctrl->ratiolimit,
               ctrl->candtab,
               ctrl->etree,
               symbmtx,
               ctrl->costmtx );

#ifdef DEBUG_BLEND
    if(ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, " New Cost of the Matrix %g \n",
                totalCost(symbmtx->cblknbr, ctrl->costmtx));
#endif
}


/*
 Function: splitOnProcs

 Parameters:
 symbmtx    - Symbolic matrix
 extrasymb  -
 extracost  -
 ctrl       -
 dofptr     -
 cblknum    -
 procnbr    -
 */
void splitOnProcs(SymbolMatrix      *symbmtx,
                  ExtraSymbolMatrix *extrasymb,
                  ExtraCostMatrix   *extracost,
                  BlendCtrl         *ctrl,
                  const Dof         *dofptr,
                  pastix_int_t                cblknum,
                  pastix_int_t                procnbr)
{
    pastix_int_t i;
    pastix_int_t blas_min_col;
    pastix_int_t blas_max_col;
    pastix_int_t pas;
    pastix_int_t *seq;
    pastix_int_t nseq;
    pastix_int_t width;


    /* if only one proc : no need to split */
    /*if(procnbr == 1)
     return;*/


    /* Compute minimun broadness for splitting this cblk */
    if(procnbr > ctrl->ratiolimit)
    {
        blas_min_col = ctrl->blblokmin;
        blas_max_col = ctrl->blblokmax;
    }
    else
    {
        blas_min_col = ctrl->blcolmin;
        blas_max_col = ctrl->blcolmax;
    }


#ifdef DOF_CONSTANT
    blas_min_col /= dofptr->noddval;
    if(blas_min_col == 0)
        blas_min_col = 1;
#endif

    /* number of split, for instance
     we choose to split at the maximum */
#ifdef DOF_CONSTANT

    width = symbmtx->cblktab[cblknum].lcolnum
        -   symbmtx->cblktab[cblknum].fcolnum + 1;

    if(procnbr == 1)
    {
        /*** Need to split big supernode because
         the diagonal block factorization is written
         in BLAS1 (due to the pivoting in LDLt and LU) ***/
        /*
         * if the column block size is small enough there is no need to
         * split it.
         */
        if( width <= blas_max_col)
            return;

        pas  = blas_max_col;
        nseq = width / pas;
    }
    else
    {
        pastix_int_t abs = ctrl->abs;
        if(procnbr > ctrl->ratiolimit)
        {
            abs *= 2; /* Increase abs for 2D */
        }

        /***  If option adaptative block size is set then compute the size of a column block ***/
        if(abs > 0)
        {
            pas = width / (abs * procnbr);

            pas = MAX(pas, blas_min_col);
            pas = MIN(pas, blas_max_col);

            nseq = width / pas;
        }
        else
        {
            nseq = (int)ceil( (double)width / (double)blas_min_col );
        }
    }

    /*nseq = splitSeqCblk2D(cblknum, procnbr, symbmtx, extrasymb, dofptr, ctrl, P2D);*/
#endif

    /** No parallelism available above 4 splitted cblk **/
    if(nseq < 4)
        return;

#ifdef SMART_CBLK_SPLIT
    {
        int old_nseq = nseq;
        smart_cblk_split(ctrl,
                         symbmtx,
                         cblknum,
                         procnbr,
                         blas_min_col,
                         blas_max_col,
                         &nseq,
                         &seq);
        /** No parallelism available above 4 splitted cblk **/
        if(nseq < 4)
            return;
    }
#else /* SMART_CBLK_SPLIT */
    pas = (int)((symbmtx->cblktab[cblknum].lcolnum -
                 symbmtx->cblktab[cblknum].fcolnum+1.0)/nseq);
    /* extend intvec to contain seq */
    extendint_ToSize(2*nseq, ctrl->intvec);
    seq = ctrl->intvec->inttab;

    for(i=0;i<nseq-1;i++)
    {
        seq[2*i]   = symbmtx->cblktab[cblknum].fcolnum + pas*i;
        seq[2*i+1] = symbmtx->cblktab[cblknum].fcolnum + pas*(i+1)-1;
    }
    seq[2*(nseq-1)] =  symbmtx->cblktab[cblknum].fcolnum + pas*(nseq-1);
    seq[2*nseq-1]   =  symbmtx->cblktab[cblknum].lcolnum;

    ASSERTDBG(seq[2*(nseq-1)]<= seq[2*nseq-1],MOD_BLEND);

#endif /* SMART_CBLK_SPLIT */
    splitCblk(symbmtx, extrasymb, extracost, ctrl, dofptr, cblknum, nseq, seq);
}


/*
 Function: splitCblk

 Split a column bloc in nseq column blocs according to the seq array.

 Update symbolic matrix and extra symbolic matrix.

 Update cost of the node descendant of the splitted column bloc
 in the elimination tree.

 Parameters:
 symbmtx   - Symbol matrix.
 extrasymb - Extra symbol matrix.
 extracost -
 ctrl      - Blend control structure.
 dofptr    - Structure for degree of freedom,
 ! does not work anymore !.
 cblknum   - Column bloc to split.
 nseq      - Number of part of the split
 *seq      - Splitting indexes array.
 */
void  splitCblk(SymbolMatrix      *symbmtx,
                ExtraSymbolMatrix *extrasymb,
                ExtraCostMatrix   *extracost,
                BlendCtrl         *ctrl,
                const Dof         *dofptr,
                pastix_int_t         cblknum,
                pastix_int_t         nseq,
                pastix_int_t        *seq)
{
    pastix_int_t i, j, s;
    pastix_int_t bloknbr      = 0;
    pastix_int_t splitbloknbr = 0;

    /**no need to split **/
    if(nseq == 1)
        return;

    assert(symbmtx->cblktab[cblknum].fcolnum == seq[0]       );
    assert(symbmtx->cblktab[cblknum].lcolnum == seq[2*nseq-1]);

    /* number of blok in the cblk to be split (diag included)  */
    bloknbr = 0;
    for(j = symbmtx->cblktab[cblknum].bloknum;
        j < symbmtx->cblktab[cblknum+1].bloknum; j++)
    {
        if(extrasymb->sptblok[j] >= 0)
            bloknbr += extrasymb->sptblnb[j];
        else
            bloknbr++;
    }

    /*
     * XL: For Schur complement we keep the last column block.as full
     */
    if (!((ctrl->iparm[IPARM_SCHUR] == API_YES) &&
          (symbmtx->cblktab[cblknum].lcolnum == symbmtx->nodenbr -1)))
    {
        /** mark the cblk to be splitted **/
        /** NB odb of this cblk don't need to be marked
         *  because they won't be splitted any more in our
         *  top-down tree strategy **/
        extrasymb->sptcblk[cblknum] = extrasymb->curcblk;
        extrasymb->sptcbnb[cblknum] = nseq;

        /* now create the new cblk and associated blok */
        for(i=0;i<nseq;i++)
        {
            extrasymb->cblktab[extrasymb->curcblk].fcolnum = seq[2*i];
            extrasymb->cblktab[extrasymb->curcblk].lcolnum = seq[2*i+1];
            extrasymb->cblktab[extrasymb->curcblk].bloknum = extrasymb->curblok;

            /* create diag blok */
            extrasymb->bloktab[extrasymb->curblok].frownum = seq[2*i];
            extrasymb->bloktab[extrasymb->curblok].lrownum = seq[2*i+1];
            extrasymb->bloktab[extrasymb->curblok].cblknum = extrasymb->curcblk;
            extrasymb->bloktab[extrasymb->curblok].levfval = 0;

            extra_inc_blok(extrasymb, extracost);
            splitbloknbr++;

            /* create odb due to the splitting of the diag blok */
            for(j=i+1;j<nseq;j++)
            {
                extrasymb->bloktab[extrasymb->curblok].frownum = seq[2*j];
                extrasymb->bloktab[extrasymb->curblok].lrownum = seq[2*j+1];
                extrasymb->bloktab[extrasymb->curblok].cblknum = extrasymb->sptcblk[cblknum] + j;
                extrasymb->bloktab[extrasymb->curblok].levfval = 0;
                extra_inc_blok(extrasymb, extracost);
                splitbloknbr++;
            }
            /* create other odb */
            /* We have to test if some of them have been splitted before */
            for(j=symbmtx->cblktab[cblknum].bloknum+1;j<symbmtx->cblktab[cblknum+1].bloknum;j++)
            {
                /* this odb hasn't been splitted */
                if(extrasymb->sptblok[j]<0)
                {
                    extrasymb->bloktab[extrasymb->curblok].frownum = symbmtx->bloktab[j].frownum;
                    extrasymb->bloktab[extrasymb->curblok].lrownum = symbmtx->bloktab[j].lrownum;
                    extrasymb->bloktab[extrasymb->curblok].cblknum = symbmtx->bloktab[j].cblknum;
                    extrasymb->bloktab[extrasymb->curblok].levfval = 0;
                    extra_inc_blok(extrasymb, extracost);
                    splitbloknbr++;
                }
                /* this odb has been splitted before (its facing diag blok has been splitted) */
                else
                {
                    for(s=extrasymb->sptblok[j];s<extrasymb->sptblok[j]+extrasymb->sptblnb[j];s++)
                    {
                        extrasymb->bloktab[extrasymb->curblok].frownum = extrasymb->bloktab[s].frownum;
                        extrasymb->bloktab[extrasymb->curblok].lrownum = extrasymb->bloktab[s].lrownum;
                        extrasymb->bloktab[extrasymb->curblok].cblknum = extrasymb->bloktab[s].cblknum;
                        extrasymb->bloktab[extrasymb->curblok].levfval = 0;
                        extra_inc_blok(extrasymb, extracost);
                        splitbloknbr++;
                    }
                }
            }
            extra_inc_cblk(extrasymb);
        }
        /* update extracblk and extrablok */
        extrasymb->addcblk += nseq-1;
        extrasymb->addblok += splitbloknbr - bloknbr;
    }

    /** we have to add an extra cblk to extrasymb because of side effect **/
    /* NB extrasymb->cblktab have at least one extra-allocated cells
     *    so don't worry about allocated memory size */
    ASSERTDBG(extrasymb->sizcblk > extrasymb->curcblk,MOD_BLEND);
    extrasymb->cblktab[extrasymb->curcblk].bloknum = extrasymb->curblok;

    /** Now we're going to split odb that hit our diag blok **/
    /** We have to mark them as splitted bloks because they may be split again later
     if their cblk owner is splitted **/
    /** NB odb blok that hit the diag (those that we're about to split )
     can have been splitted before because of our top-down tree strategy **/
    for(i=0;i<ctrl->egraph->verttab[cblknum].innbr;i++)
    {
        pastix_int_t sptbloknbr; /* number of splitted bloks resulting */
        pastix_int_t bloknum;
        bloknum = ctrl->egraph->inbltab[ctrl->egraph->verttab[cblknum].innum+i];

        ASSERTDBG(symbmtx->bloktab[bloknum].cblknum == cblknum,MOD_BLEND);

        sptbloknbr = 0;
        /* mark this blok as splitted */
        extrasymb->sptblok[bloknum] = extrasymb->curblok;
        for(j=0;j<nseq;j++)
        {
            /* there are six possible cases of intersection
             beetween the blok and the seq we consider ,
             among these six cases, only 4 of them are
             not empty intersection */

            /* empty intersections */
            if(symbmtx->bloktab[bloknum].frownum > seq[2*j+1])
                continue;

            if(symbmtx->bloktab[bloknum].lrownum < seq[2*j])
                /* in this case there will no more splitted blok to create */
                break;

            /* not empty intersections */
            if((symbmtx->bloktab[bloknum].frownum >= seq[2*j])
               && (symbmtx->bloktab[bloknum].lrownum >= seq[2*j+1]))
            {
                extrasymb->bloktab[extrasymb->curblok].frownum = symbmtx->bloktab[bloknum].frownum;
                extrasymb->bloktab[extrasymb->curblok].lrownum = seq[2*j+1];
                goto endloop;
            }
            if((symbmtx->bloktab[bloknum].frownum <= seq[2*j])
               && (symbmtx->bloktab[bloknum].lrownum >= seq[2*j+1]))
            {
                extrasymb->bloktab[extrasymb->curblok].frownum = seq[2*j];
                extrasymb->bloktab[extrasymb->curblok].lrownum = seq[2*j+1];
                goto endloop;
            }
            if((symbmtx->bloktab[bloknum].frownum <= seq[2*j])
               && (symbmtx->bloktab[bloknum].lrownum <= seq[2*j+1]))
            {
                extrasymb->bloktab[extrasymb->curblok].frownum = seq[2*j];
                extrasymb->bloktab[extrasymb->curblok].lrownum = symbmtx->bloktab[bloknum].lrownum;
                goto endloop;
            }
            if((symbmtx->bloktab[bloknum].frownum >= seq[2*j])
               && (symbmtx->bloktab[bloknum].lrownum <= seq[2*j+1]))
            {
                extrasymb->bloktab[extrasymb->curblok].frownum = symbmtx->bloktab[bloknum].frownum;
                extrasymb->bloktab[extrasymb->curblok].lrownum = symbmtx->bloktab[bloknum].lrownum;
                goto endloop;
            }
          endloop:
            extrasymb->bloktab[extrasymb->curblok].cblknum = extrasymb->sptcblk[cblknum]+j;
            extrasymb->bloktab[extrasymb->curblok].levfval = 0;
            sptbloknbr++;
            extra_inc_blok(extrasymb, extracost);
        }
        extrasymb->sptblnb[bloknum] = sptbloknbr;
        extrasymb->addblok += sptbloknbr-1;

        /** update cost of the cblk owning the splitted bloks **/
        //blokUpdateCost(bloknum, ctrl->egraph->ownetab[bloknum], ctrl->costmtx, extracost, symbmtx, extrasymb, ctrl, dofptr);
    }
}

double maxProcCost(double *proc_cost, pastix_int_t procnbr)
{
    double maxcost = 0;
    pastix_int_t p;
    for(p=0;p<procnbr;p++)
        if(proc_cost[p]>maxcost)
            maxcost = proc_cost[p];
    return maxcost;
}

/* /\*+ Recompute cost of cblk which some odb have been splitted, return new total cost - old total cost +*\/ */
/* double blokUpdateCost(pastix_int_t bloknum, pastix_int_t cblknum, CostMatrix *costmtx, ExtraCostMatrix *extracost, const SymbolMatrix *symbmtx, const ExtraSymbolMatrix *extrasymb, BlendCtrl *ctrl, const Dof * dofptr) */
/* { */
/*     pastix_int_t L, h, g; */
/*     pastix_int_t s; */
/*     double oldcost, newcost; */
/* #ifndef DOF_CONSTANT */
/*     pastix_int_t i; */
/* #endif */

/*     ASSERTDBG(extrasymb->sptblok[bloknum] >= 0,MOD_BLEND); */

/*     /\** we need the height of the odb and the broadness */
/*      of the cbl to compute the local compute cost **\/ */

/* #ifdef DOF_CONSTANT */
/*     L = (symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1)*(dofptr)->noddval; */
/* #else */
/*     for(i=symbmtx->cblktab[cblknum].fcolnum;i<=symbmtx->cblktab[cblknum].lcolnum;i++) */
/*         L+= noddDlt(dofptr, i); */
/* #endif */
/*     /\** no need to recompute the local compute cost because odb lines number */
/*      is not changed **\/ */


/*     /\** recompute for each splitted odb its contribution compute cost and add cost **\/ */
/*     oldcost = 0; */
/*     newcost = 0; */

/*     oldcost = costmtx->bloktab[bloknum].contrib; /\* cost of blok not splitted *\/ */
/*     g = costmtx->bloktab[bloknum].linenbr; */

/*     /\* now compute the new cost of each resulting blok of the spliting *\/ */
/*     for(s = extrasymb->sptblok[bloknum];s<extrasymb->sptblok[bloknum]+extrasymb->sptblnb[bloknum];s++) */
/*     { */

/* #ifdef  DOF_CONSTANT */
/*         h = (extrasymb->bloktab[s].lrownum - extrasymb->bloktab[s].frownum + 1)*(dofptr)->noddval; */
/* #else */
/*         for(i=extrasymb->bloktab[s].frownum;i<=extrasymb->bloktab[s].lrownum;i++) */
/*             h+= noddDlt(dofptr, i); */
/* #endif */
/*         extracost->bloktab[s].linenbr = g; */
/*         extracost->bloktab[s].contrib =  contribCompCost(L, h, g) + contribAddCost(h, g); */
/*         newcost += extracost->bloktab[s].contrib; */
/*         g -= h; */
/*     } */

/*     {  */
/*         pastix_int_t bloknum = symbmtx->cblktab[ cblknum ].bloknum; */
/*         if(ctrl->candtab[cblknum].distrib == D1) */
/*             costmtx->bloktab[bloknum].contrib += newcost - oldcost; */
/*     } */
/*     return newcost - oldcost; */
/* } */


pastix_int_t countBlok(pastix_int_t cblknum, SymbolMatrix *symbptr, pastix_int_t blcolmin)
{
    pastix_int_t i;
    pastix_int_t bloknbr;
    double delta;
    double stride = 0;
    delta = (double)(symbptr->cblktab[cblknum].lcolnum - symbptr->cblktab[cblknum].fcolnum+1);
    delta = ceil(delta/blcolmin);

    for(i=symbptr->cblktab[cblknum].bloknum+1;i<symbptr->cblktab[cblknum+1].bloknum;i++)
        stride += symbptr->bloktab[i].lrownum-symbptr->bloktab[i].frownum + 1;
    stride = ceil(stride/blcolmin);
    /*fprintf(stdout, "delta %g stride %g blcolmin %ld \n", delta, stride, blcolmin); */
    bloknbr = 0;
    bloknbr += (pastix_int_t) (((delta + 1)*delta)/2);
    bloknbr += (pastix_int_t) (stride*delta);

    return bloknbr;
}

pastix_int_t setSubtreeBlokNbr(pastix_int_t rootnum, const EliminTree *etree, SymbolMatrix *symbptr, ExtraSymbolMatrix *extrasymb, pastix_int_t blcolmin)
{
    pastix_int_t i;
    extrasymb->subtreeblnbr[rootnum] =  countBlok(rootnum, symbptr, blcolmin);
    /*fprintf(stdout, "Rootnum %ld bloknbr %ld \n", rootnum, extrasymb->blnbtab[rootnum]);*/
    for(i=0;i<etree->nodetab[rootnum].sonsnbr;i++)
        extrasymb->subtreeblnbr[rootnum] += setSubtreeBlokNbr(eTreeSonI(etree, rootnum, i), etree, symbptr, extrasymb, blcolmin);
    return extrasymb->subtreeblnbr[rootnum];
}


