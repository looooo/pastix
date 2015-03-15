/**
 *
 * @file cand.c
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

void
candInit( Cand *candtab,
          pastix_int_t cblknbr )
{
    pastix_int_t i;
    for(i=0;i<cblknbr;i++)
    {
        candtab[i].costlevel = 0.0;
        candtab[i].treelevel = 0;
        candtab[i].fcandnum  = -1;
        candtab[i].lcandnum  = -1;
        candtab[i].fccandnum = -1;
        candtab[i].lccandnum = -1;
        candtab[i].cluster   = -1;
        candtab[i].cblktype  = CBLK_1D;
    }
}

void
candSave( const Cand *candtab,
          pastix_int_t cblknbr )
{
    FILE *f = fopen("candtab.txt", "w");
    pastix_int_t i;
    fprintf(f, "%ld\n", cblknbr );
    for(i=0;i<cblknbr;i++)
    {
        fprintf(f, "%lf %ld %ld %ld %ld %ld %ld %d\n",
                candtab[i].costlevel,
                candtab[i].treelevel,
                candtab[i].fcandnum,
                candtab[i].lcandnum,
                candtab[i].fccandnum,
                candtab[i].lccandnum,
                candtab[i].cluster,
                candtab[i].cblktype );
    }
    fclose(f);
}

void
candSetSubCandidate( Cand *candtab,
                     const EliminTree *etree,
                     pastix_int_t rootnum,
                     pastix_int_t procnum )
{
    pastix_int_t i;

    candtab[rootnum].fcandnum = procnum;
    candtab[rootnum].lcandnum = procnum;

    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
        candSetSubCandidate( candtab, etree, eTreeSonI(etree, rootnum, i), procnum );
}

void
candSetClusterCand( Cand *candtab,
                    pastix_int_t  cblknbr,
                    pastix_int_t *core2clust,
                    pastix_int_t  coresnbr )
{
    pastix_int_t i;
    (void)coresnbr;

    for(i=0; i<cblknbr; i++) {
        assert( candtab[i].fcandnum >= 0 );
        assert( candtab[i].lcandnum >= 0 );
        assert( candtab[i].fcandnum < coresnbr );
        assert( candtab[i].lcandnum < coresnbr );
        candtab[i].fccandnum = core2clust[ candtab[i].fcandnum ];
        candtab[i].lccandnum = core2clust[ candtab[i].lcandnum ];
    }
}

int
candCheck( Cand *candtab, SymbolMatrix *symbmtx )
{
    pastix_int_t i, j;
    pastix_int_t facecblknum;

    for(i=0; i<symbmtx->cblknbr; i++)
    {
        for(j = symbmtx->cblktab[i].bloknum;
            j < symbmtx->cblktab[i+1].bloknum; j++)
        {
            facecblknum = symbmtx->bloktab[j].cblknum;

            if( (candtab[i].fcandnum < candtab[facecblknum].fcandnum) ||
                (candtab[i].lcandnum > candtab[facecblknum].lcandnum) )
            {
                errorPrint("bad processor candidat sets : cblk %ld candidat =[%ld %ld] father %ld candidat = [%ld %ld].",
                           (long)i, (long)candtab[i].fcandnum, (long)candtab[i].lcandnum,
                           (long)facecblknum, (long)candtab[facecblknum].fcandnum,
                           (long)candtab[facecblknum].lcandnum);
                return 0;
            }
        }
    }
    return 1;
}


static inline double
candSubTreeBuild( pastix_int_t        rootnum,
                  Cand               *candtab,
                  EliminTree         *etree,
                  const SymbolMatrix *symbmtx,
                  const CostMatrix   *costmtx )
{
    double cost;
    pastix_int_t i, son, bloknum;

    /* Compute cost of current node */
#if defined(BLEND_COST_LL)
    cost = costmtx->cblkcost[rootnum];
#else
    cost = 0.;
    for( bloknum = symbmtx->cblktab[ rootnum   ].bloknum;
         bloknum < symbmtx->cblktab[ rootnum+1 ].bloknum; bloknum++)
    {
        cost += costmtx->bloktab[ bloknum ].contrib;
    }
#endif

    etree->nodetab[ rootnum ].total   = cost;
    etree->nodetab[ rootnum ].subtree = cost;

    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        candtab[ son ].treelevel = candtab[ rootnum ].treelevel - 1;
        candtab[ son ].costlevel = candtab[ rootnum ].costlevel - cost;

        etree->nodetab[ rootnum ].subtree +=
            candSubTreeBuild( son, candtab, etree, symbmtx, costmtx );
    }

    return etree->nodetab[ rootnum ].subtree;
}

static inline void
candSubTreeDistribWithSize( pastix_int_t        rootnum,
                            pastix_int_t        cblktype,
                            pastix_int_t        ratiolimit,
                            Cand               *candtab,
                            const EliminTree   *etree,
                            const SymbolMatrix *symbmtx )
{
    pastix_int_t i, son;

    if(cblktype == CBLK_1D)
    {
        candtab[ rootnum ].cblktype = CBLK_1D;
    }
    else
    {
        pastix_int_t width = symbmtx->cblktab[ rootnum ].lcolnum - symbmtx->cblktab[ rootnum ].fcolnum + 1;

        if(width >= ratiolimit)
            candtab[ rootnum ].cblktype = CBLK_SPLIT;
        else
            candtab[ rootnum ].cblktype = CBLK_1D;
    }

    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        candSubTreeDistribWithSize( son, candtab[ rootnum ].cblktype, ratiolimit, candtab, etree, symbmtx);
    }
}

static inline void
candDistribWithDepth( pastix_int_t depth,
                      pastix_int_t cblknbr,
                      Cand        *candtab )
{
    pastix_int_t i;

    for(i=0;i<cblknbr;i++)
    {
        if( candtab[i].treelevel > depth )
            candtab[i].cblktype = CBLK_SPLIT;
        else
            candtab[i].cblktype = CBLK_1D;
    }
}

void
candBuild( pastix_int_t autolevel, pastix_int_t level2D, double ratiolimit,
           Cand               *candtab,
           EliminTree         *etree,
           const SymbolMatrix *symbmtx,
           const CostMatrix   *costmtx )
{
    pastix_int_t root = eTreeRoot(etree);

    /* Let's start with the root */
    candtab[ root ].costlevel = -1.0;
    candtab[ root ].treelevel = -1;

    candSubTreeBuild( root, candtab, etree, symbmtx, costmtx );

    /* Let's set the cblk type of each node */
    /* For now, it sets the distrib field and not the cblktype field */
    if(autolevel)
    {
        candSubTreeDistribWithSize( eTreeRoot(etree), CBLK_SPLIT, (pastix_int_t)ratiolimit,
                                    candtab, etree, symbmtx );
    }
    else
    {
        candDistribWithDepth( -level2D, symbmtx->cblknbr, candtab );
    }
}
