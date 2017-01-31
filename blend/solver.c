/**
 * solver.c -- SolverMatrix description.
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/

#include "common.h"
#include "solver.h"
#include "coeftab.h"

void
solverInit(SolverMatrix *solvmtx)
{
    memset(solvmtx, 0, sizeof (SolverMatrix));
    return;
}

void
solverExit(SolverMatrix *solvmtx)
{
    pastix_int_t i;

    coeftabExit( solvmtx );

    /* Free arrays of solvmtx */
    if(solvmtx->cblktab) {
        memFree_null(solvmtx->cblktab);
    }
    if(solvmtx->bloktab) {
        memFree_null(solvmtx->bloktab);
    }
    if(solvmtx->browtab) {
        memFree_null(solvmtx->browtab);
    }
    if(solvmtx->ftgttab) {
        memFree_null(solvmtx->ftgttab);
    }
    if(solvmtx->tasktab) {
        memFree_null(solvmtx->tasktab);
    }
    if(solvmtx->indtab) {
        memFree_null(solvmtx->indtab);
    }
    memFree_null(solvmtx->ttsknbr);
    for (i=0;i<solvmtx->bublnbr;i++)
    {
        if (solvmtx->ttsktab[i] != NULL) {
            memFree_null(solvmtx->ttsktab[i]);
        }
    }
    memFree_null(solvmtx->ttsktab);
    memFree_null(solvmtx->proc2clust);


    /*memFree_null(solvmtx);*/
#if defined(PASTIX_WITH_STARPU)
    memFree_null(solvmtx->hcblktab);
    memFree_null(solvmtx->hbloktab);
    memFree_null(solvmtx->gcblk2halo);
    if (pastix_starpu_with_fanin() == API_YES) {
        pastix_int_t clustnum;
        for (clustnum = 0; clustnum < solvmtx->clustnbr; clustnum++) {
            if (solvmtx->fcblknbr[clustnum] > 0) {
                memFree_null(solvmtx->fcblktab[clustnum]);
                memFree_null(solvmtx->fbloktab[clustnum]);
            }
        }
        memFree_null(solvmtx->fcblknbr);
        memFree_null(solvmtx->fcblktab);
        memFree_null(solvmtx->fbloktab);
    }
#endif

    /* /\* Pour l'instant uniquement si on est en 1d *\/ */
    /* if (solvmtx->updovct.cblktab) { */
    /*     for (i=0; i<solvmtx->cblknbr; i++) */
    /*     { */
    /*         if (solvmtx->updovct.cblktab[i].browcblktab) */
    /*             memFree_null(solvmtx->updovct.cblktab[i].browcblktab); */
    /*         if (solvmtx->updovct.cblktab[i].browproctab) */
    /*             memFree_null(solvmtx->updovct.cblktab[i].browproctab); */
    /*     } */
    /* } */
    /* memFree_null(solvmtx->updovct.lblk2gcblk); */
    /* memFree_null(solvmtx->updovct.listblok); */
    /* memFree_null(solvmtx->updovct.listcblk); */
    /* memFree_null(solvmtx->updovct.gcblk2list); */
    /* memFree_null(solvmtx->updovct.loc2glob); */
    /* memFree_null(solvmtx->updovct.cblktab); */
    /* memFree_null(solvmtx->updovct.listptr); */
}

/**
 *  Function: sizeofsolver
 *
 *  Computes the size in memory of the SolverMatrix.
 *
 *  Parameters:
 *    solvptr - address of the SolverMatrix
 *
 *  Returns:
 *    SolverMatrix size.
 */
static inline size_t
solver_size( const SolverMatrix *solvptr )
{
    size_t mem = sizeof(SolverMatrix);

    /* cblk and blocks arrays */
    if ( solvptr->cblktab ) {
        mem += solvptr->cblknbr * sizeof( SolverCblk );
    }
    if ( solvptr->bloktab ) {
        mem += solvptr->bloknbr * sizeof( SolverBlok );
    }
    if ( solvptr->browtab ) {
        mem += solvptr->brownbr * sizeof( pastix_int_t );
    }
#if defined(PASTIX_WITH_PARSEC)
    if ( solvptr->parsec_desc ) {
        mem += sizeof( sparse_matrix_desc_t );
    }
#endif

    /* Fanins */
    if ( solvptr->ftgttab ) {
        mem += solvptr->ftgtnbr * sizeof( FanInTarget );
    }
    if ( solvptr->indtab ) {
        mem += solvptr->indnbr  * sizeof( pastix_int_t );
    }

    /* BubbleTree */
    if ( solvptr->btree ) {
        mem += solvptr->bublnbr * sizeof( BubbleTree );
        mem += solvptr->btree->nodemax * sizeof( BubbleTreeNode );
        mem += solvptr->btree->nodemax * sizeof( int );
    }

    /* Tasks */
    if ( solvptr->tasktab ) {
        mem += solvptr->tasknbr * sizeof(Task);
    }
    if ( solvptr->ttsknbr ) {
        int i;
        mem += solvptr->thrdnbr * sizeof(pastix_int_t);
        mem += solvptr->thrdnbr * sizeof(pastix_int_t*);

        for( i=0; i<solvptr->thrdnbr; i++ ) {
            mem += solvptr->ttsknbr[i] * sizeof(pastix_int_t);
        }
    }

    /* proc2clust */
    if ( solvptr->proc2clust ) {
        mem += solvptr->procnbr * sizeof(pastix_int_t);
    }

    return mem;
}


void
solverPrintStats( const SolverMatrix *solvptr )
{
    SolverCblk *cblk;
    size_t memstruct, memcoef;
    pastix_int_t fcol2d, lcol2d;
    pastix_int_t itercblk, cblknbr;
    double avg2d;

    fcol2d = (solvptr->cblktab + solvptr->cblkmin2d )->fcolnum;
    lcol2d = (solvptr->cblktab + solvptr->cblknbr   )->fcolnum;
    avg2d  = 0.;
    if ( (lcol2d - fcol2d) > 0 ) {
        avg2d  = (double)(lcol2d - fcol2d) / (double)( solvptr->cblknbr - solvptr->cblkmin2d);
    }
    cblknbr = solvptr->cblknbr;
    cblk    = solvptr->cblktab;
    memcoef = 0;
    for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
    {
        pastix_int_t colnbr = cblk->lcolnum - cblk->fcolnum + 1;
        pastix_int_t rownbr = cblk->stride;

        memcoef += colnbr * rownbr;
    }

    memstruct = solver_size( solvptr );

    fprintf( stdout,
             "    Solver Matrix statistics:\n"
             "      Number of cblk                    %10ld\n"
             "      Number of cblks in 2D             %10ld\n"
             "      Number of blocks in 2D            %10ld\n"
             "      First cblk in 2D                  %10ld\n"
             "      First row handled in 2D           %10ld\n"
             "      Average 2D cblk size             %11.2lf\n"
             "      Structure memory space           %11.2lf %s\n"
             "      Number of coeficients stored      %10ld\n",
             solvptr->cblknbr,
             0, /* solvptr->nb2dcblk, */
             0, /* solvptr->nb2dblok, */
             solvptr->cblkmin2d,
             fcol2d,
             avg2d,
             MEMORY_WRITE( memstruct ), MEMORY_UNIT_WRITE( memstruct ),
             memcoef );
}
