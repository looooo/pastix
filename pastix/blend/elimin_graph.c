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

pastix_int_t
eGraphInit( EliminGraph *egraph )
{
    egraph->baseval = 0;
    egraph->vertnbr = 0;
    egraph->verttab = NULL;
    egraph->inbltab = NULL;
    egraph->ownetab = NULL;
    return 1;
}

void
eGraphExit( EliminGraph *egraph )
{
    memFree_null(egraph->verttab);
    memFree_null(egraph->inbltab);
    memFree_null(egraph->ownetab);
    memFree_null(egraph);
}

void
eGraphBuild( EliminGraph        *egraph,
             const SymbolMatrix *symbmtx )
{
    pastix_int_t i, j;
    pastix_int_t *tmp = NULL;
    pastix_int_t cursor;

    egraph->vertnbr = symbmtx->cblknbr;
    MALLOC_INTERN(egraph->verttab, egraph->vertnbr,  EliminVertex);
    MALLOC_INTERN(egraph->ownetab, symbmtx->bloknbr, pastix_int_t);

    if((symbmtx->bloknbr - symbmtx->cblknbr) > 0 ) {
        MALLOC_INTERN(egraph->inbltab,
                      symbmtx->bloknbr - symbmtx->cblknbr,
                      pastix_int_t);
    }

    /* Initialize verttab */
    for(i=0;i < egraph->vertnbr;i++)
    {
        egraph->verttab[i].innum = -1;
        egraph->verttab[i].innbr = -1;
    }

    /* Fill ownetab */
    for(i=0;i < symbmtx->cblknbr; i++)
    {
        for(j=symbmtx->cblktab[i].bloknum;
            j<symbmtx->cblktab[i+1].bloknum; j++)
        {
            egraph->ownetab[j] = i;
        }
    }

    /*
     * Compute deg in for each Vertex.
     * Rq innbr start at -1 because diag blok cause one extra in-edge
     */
    for(i=0;i<symbmtx->bloknbr;i++)
    {
        pastix_int_t cblknum = symbmtx->bloktab[i].cblknum;
        assert( cblknum >= 0 );
        assert( cblknum < symbmtx->cblknbr );
        egraph->verttab[ cblknum ].innbr++;
    }

    /* Compute innum for each vertex */
    cursor = 0;
    for(i=0;i<egraph->vertnbr;i++)
    {
        egraph->verttab[i].innum = cursor;
        cursor += egraph->verttab[i].innbr;
    }

    assert(cursor == (symbmtx->bloknbr-symbmtx->cblknbr));

    /* Fill inbltab */
    if (cursor > 0) {
        MALLOC_INTERN(tmp, symbmtx->cblknbr, pastix_int_t);
        bzero(tmp, symbmtx->cblknbr * sizeof(pastix_int_t));
        for(i=0;i<symbmtx->bloknbr;i++)
        {
            /* If this bloc is not a diag bloc register it*/
            if(symbmtx->cblktab[symbmtx->bloktab[i].cblknum].bloknum != i)
            {
                egraph->inbltab[ egraph->verttab[symbmtx->bloktab[i].cblknum].innum
                                 + tmp[symbmtx->bloktab[i].cblknum] ] = i;
                tmp[symbmtx->bloktab[i].cblknum]++;
            }
        }
        memFree_null(tmp);
    }
}
