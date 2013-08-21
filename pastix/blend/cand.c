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
        candtab[i].distrib   = D1;
        candtab[i].cluster   = -1;
    }
}

static inline void
candSetSubTreelevel( Cand *candtab, const EliminTree *etree, pastix_int_t rootnum )
{
    pastix_int_t i, son;
    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        candtab[ son ].treelevel = candtab[rootnum].treelevel - 1;
        candSetSubTreelevel(candtab, etree, son );
    }
}

void
candSetTreelevel( Cand *candtab, const EliminTree *etree )
{
    pastix_int_t root = eTreeRoot(etree);
    candtab[root].treelevel = -1;
    candSetSubTreelevel(candtab, etree, root );
}

static inline void
candSetSubCostlevel(Cand *candtab, const EliminTree *etree, const CostMatrix *costmtx, pastix_int_t rootnum )
{
    pastix_int_t i, son;
    for(i=0; i<etree->nodetab[rootnum].sonsnbr; i++)
    {
        son = eTreeSonI(etree, rootnum, i);
        candtab[ son ].costlevel = candtab[rootnum].costlevel - costmtx->cblktab[rootnum].total;
        candSetSubCostlevel( candtab, etree, costmtx, son );
    }
}

void
candSetCostlevel(Cand *candtab, const EliminTree *etree, const CostMatrix *costmtx)
{
    pastix_int_t root = eTreeRoot(etree);
    candtab[ root ].costlevel = -1.0;
    candSetSubCostlevel( candtab, etree, costmtx, root );
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
