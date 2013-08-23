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

static inline void
splitSubtree( pastix_int_t rootnum,
              const EliminTree   *etree,
              const SymbolMatrix *symbmtx,
              ExtraSymbolMatrix  *extrasymb,
              ExtraCostMatrix    *extracost,
              const BlendCtrl    *ctrl,
              const Dof          *dofptr,
              const Cand         *candtab )
{
    pastix_int_t i, son;

    splitOnProcs( symbmtx, extrasymb, extracost, ctrl,
                  dofptr, rootnum,
                  candtab[ rootnum ].lcandnum -
                  candtab[ rootnum ].fcandnum + 1 );

    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        splitSubtree( son, etree, symbmtx, extrasymb, extracost, ctrl, dofptr, candtab );
    }
}

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

        extrasymb->curcblk = -1;
        extrasymb->curblok = -1;

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
        splitSubtree( eTreeRoot( ctrl->etree), ctrl->etree, symbmtx, extrasymb, extracost, ctrl, dofptr, ctrl->candtab );

        extrasymb->curblok++;
        extrasymb->curcblk++;
    }

    /* Rebuild the symbolic matrix */
    {
        Clock timer;
        clockStart(timer);
        symbolMerge( ctrl, dofptr,
                     symbmtx,       extrasymb,
                     ctrl->costmtx, extracost,
                     &(ctrl->candtab) );
        clockStop(timer);

        pastix_print( 0, 0, "symbolMerge perform in %e s\n",
                      clockVal(timer) );
    }

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
void splitOnProcs(const SymbolMatrix      *symbmtx,
                  ExtraSymbolMatrix *extrasymb,
                  ExtraCostMatrix   *extracost,
                  const BlendCtrl         *ctrl,
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

    splitCblk( ctrl, symbmtx, dofptr, extrasymb, extracost, cblknum, nseq, seq);
}



static inline void
extraAddBlok( pastix_int_t       frownum,
              pastix_int_t       lrownum,
              pastix_int_t       cblknum,
              ExtraSymbolMatrix *extrasymb,
              ExtraCostMatrix   *extracost )
{
    extra_inc_blok( extrasymb, extracost );

    /* Add the new block */
    extrasymb->bloktab[extrasymb->curblok].frownum = frownum;
    extrasymb->bloktab[extrasymb->curblok].lrownum = lrownum;
    extrasymb->bloktab[extrasymb->curblok].cblknum = cblknum;
    extrasymb->bloktab[extrasymb->curblok].levfval = 0;

    return;
}

static inline void
extraAddBlokCopy( const SymbolBlok  *blokptr,
                  ExtraSymbolMatrix *extrasymb,
                  ExtraCostMatrix   *extracost )
{
    extra_inc_blok( extrasymb, extracost );

    /* Add the new block */
    memcpy( &(extrasymb->bloktab[extrasymb->curblok]),
            blokptr, sizeof(SymbolBlok) );

    return;
}

static inline void
extraAddCblk( pastix_int_t fcolnum,
              pastix_int_t lcolnum,
              ExtraSymbolMatrix *extrasymb,
              ExtraCostMatrix   *extracost)
{
    extra_inc_cblk( extrasymb );

    /* Add the new diagonal block associated after having update curcblk */
    extraAddBlok( fcolnum, lcolnum, extrasymb->curcblk, extrasymb, extracost );

    /* Add new cblk after having updated curblok */
    extrasymb->cblktab[extrasymb->curcblk].fcolnum = fcolnum;
    extrasymb->cblktab[extrasymb->curcblk].lcolnum = lcolnum;
    extrasymb->cblktab[extrasymb->curcblk].bloknum = extrasymb->curblok;

    return;
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
void  splitCblk( const BlendCtrl    *ctrl,
                 const SymbolMatrix *symbmtx,
                 const Dof          *dofptr,
                 ExtraSymbolMatrix  *extrasymb,
                 ExtraCostMatrix    *extracost,
                 pastix_int_t        cblknum,
                 pastix_int_t        nseq,
                 pastix_int_t       *seq)
{
    pastix_int_t i, j, s;
    pastix_int_t bloknbr      = 0;
    pastix_int_t splitbloknbr = 0;

    /* No split to perform */
    if(nseq == 1)
        return;

    assert(symbmtx->cblktab[cblknum].fcolnum == seq[0]       );
    assert(symbmtx->cblktab[cblknum].lcolnum == seq[2*nseq-1]);

    /* Number of bloks in the cblk to be split (diag included)  */
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
     * Mark the cblk as being splitted
     * NB: off diagonal blocks of this cblk don't need to be marked because they
     * won't be splitted any more in our top-down tree strategy
     */
    extrasymb->sptcblk[cblknum] = extrasymb->curcblk + 1;
    extrasymb->sptcbnb[cblknum] = nseq;

    /* Now create the new cblk and associated blok */
    for(i=0;i<nseq;i++)
    {
        /* Add cblk AND its diagonal block */
        extraAddCblk( seq[2*i], seq[2*i+1], extrasymb, extracost );
        splitbloknbr++;

        /* Create odb due to the splitting of the diag blok */
        for(j=i+1; j<nseq; j++)
        {
            extraAddBlok( seq[2*j], seq[2*j+1], extrasymb->sptcblk[cblknum]+j,
                          extrasymb, extracost );
            splitbloknbr++;
        }

        /* Create other odb */
        /* We have to test if some of them have been splitted before */
        for(j = symbmtx->cblktab[cblknum].bloknum+1;
            j < symbmtx->cblktab[cblknum+1].bloknum; j++)
        {
            /* This odb hasn't been splitted */
            if( extrasymb->sptblok[j] < 0 )
            {
                extraAddBlokCopy( &(symbmtx->bloktab[j]), extrasymb, extracost );
                splitbloknbr++;
            }
            /* This odb has been splitted horizontally before (its facing diag blok that has been splitted) */
            else
            {
                for(s = extrasymb->sptblok[j];
                    s < extrasymb->sptblok[j] + extrasymb->sptblnb[j]; s++)
                {
                    extraAddBlokCopy( &(extrasymb->bloktab[s]), extrasymb, extracost );
                    splitbloknbr++;
                }
            }
        }
    }

    /* update extracblk and extrablok */
    extrasymb->addcblk += nseq-1;
    extrasymb->addblok += splitbloknbr - bloknbr;

    /** Now we're going to split odb that hit our diag blok **/
    /** We have to mark them as splitted bloks because they may be split again later
     if their cblk owner is splitted **/
    /** NB odb blok that hit the diag (those that we're about to split )
     can have been splitted before because of our top-down tree strategy **/
    for(i=0;i<ctrl->egraph->verttab[cblknum].innbr;i++)
    {
        pastix_int_t sptbloknbr; /* number of splitted bloks resulting */
        pastix_int_t bloknum;
        pastix_int_t frownum;
        pastix_int_t lrownum;
        bloknum = ctrl->egraph->inbltab[ctrl->egraph->verttab[cblknum].innum+i];
        frownum = symbmtx->bloktab[bloknum].frownum;
        lrownum = symbmtx->bloktab[bloknum].lrownum;

        assert( symbmtx->bloktab[bloknum].cblknum == cblknum );

        sptbloknbr = 0;
        /* mark this blok as splitted */
        extrasymb->sptblok[bloknum] = extrasymb->curblok + 1;
        for(j=0;j<nseq;j++)
        {
            pastix_int_t fcblknum = extrasymb->sptcblk[cblknum]+j;

            /* there are six possible cases of intersection
             beetween the blok and the seq we consider ,
             among these six cases, only 4 of them are
             not empty intersection */

            /* Empty intersections */
            if( frownum > seq[2*j+1] )
                continue;

            /* In this case there will no more splitted blok to create */
            if( lrownum < seq[2*j] )
                break;

            /* Not empty intersections */
            if( frownum >= seq[2*j] )
            {
                if ( lrownum <= seq[2*j+1] )
                    extraAddBlok( frownum,  lrownum,    fcblknum, extrasymb, extracost );
                else
                    extraAddBlok( frownum,  seq[2*j+1], fcblknum, extrasymb, extracost );
            }
            else
            {
                if ( lrownum <= seq[2*j+1] )
                    extraAddBlok( seq[2*j], lrownum,    fcblknum, extrasymb, extracost );
                else
                    extraAddBlok( seq[2*j], seq[2*j+1], fcblknum, extrasymb, extracost );
            }
            sptbloknbr++;
        }

        extrasymb->sptblnb[bloknum] = sptbloknbr;
        extrasymb->addblok += sptbloknbr-1;
    }
}
