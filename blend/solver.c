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

#if defined(PASTIX_WITH_PARSEC)
    {
        if ( solvmtx->parsec_desc != NULL ) {
            sparse_matrix_destroy( solvmtx->parsec_desc );
            free( solvmtx->parsec_desc );
        }
        solvmtx->parsec_desc = NULL;
    }
#endif

    /** Free arrays of solvmtx **/
    if(solvmtx->cblktab)
    {
        for (i = 0; i < solvmtx->cblknbr; i++)
        {
            if (solvmtx->cblktab[i].lcoeftab)
                memFree_null(solvmtx->cblktab[i].lcoeftab);

            if (solvmtx->cblktab[i].ucoeftab)
                memFree_null(solvmtx->cblktab[i].ucoeftab);

            if (! (solvmtx->cblktab[i].cblktype & CBLK_DENSE)) {
                SolverBlok *blok  = solvmtx->cblktab[i].fblokptr;
                SolverBlok *lblok = solvmtx->cblktab[i+1].fblokptr;

                for (; blok<lblok; blok++) {
                    core_zlrfree(blok->LRblock);
                    core_zlrfree(blok->LRblock+1);
                }
                free(solvmtx->cblktab[i].fblokptr->LRblock);
            }
        }
        memFree_null(solvmtx->cblktab);
    }
    if(solvmtx->bloktab)
        memFree_null(solvmtx->bloktab);
    if(solvmtx->browtab)
        memFree_null(solvmtx->browtab);
    if(solvmtx->ftgttab)
        memFree_null(solvmtx->ftgttab);
    if(solvmtx->tasktab)
        memFree_null(solvmtx->tasktab);
    if(solvmtx->indtab)
        memFree_null(solvmtx->indtab);
    memFree_null(solvmtx->ttsknbr);
    for (i=0;i<solvmtx->bublnbr;i++)
    {
        if (solvmtx->ttsktab[i] != NULL)
            memFree_null(solvmtx->ttsktab[i]);
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
pastix_int_t
sizeofsolver(const SolverMatrix *solvptr,
                   pastix_int_t *iparm )
{
  pastix_int_t result=sizeof(SolverMatrix);
  pastix_int_t iter;
  (void)iparm;

  /* cblk and blocks arrays */
  result += solvptr->cblknbr*sizeof(SolverCblk);
  result += solvptr->bloknbr*sizeof(SolverBlok);

  /* fanin target */
  result += solvptr->ftgtnbr*sizeof(FanInTarget);
  result += solvptr->indnbr *sizeof(pastix_int_t);

  /* TODO: Check that it is not bubbletree + bubblnbr * bubblenode */
  result += solvptr->bublnbr*sizeof(BubbleTree);

  /* task */
  result += solvptr->tasknbr*sizeof(Task);

  /* ttsktab */
  result += solvptr->thrdnbr*sizeof(pastix_int_t);
  result += solvptr->thrdnbr*sizeof(pastix_int_t*);
  for (iter=0; iter<solvptr->thrdnbr; iter++)
    {
      result += solvptr->ttsknbr[iter]*sizeof(pastix_int_t);
    }

  /* proc2clust */
  result += solvptr->procnbr*sizeof(pastix_int_t);

  /* UpDownVector */
  /* TODO: 2D    */
  /* if (iparm[IPARM_DISTRIBUTION_LEVEL] == 0) */
  /*   { */
  /*     /\* UpDownCblk *\/ */
  /*     result += solvptr->cblknbr      *sizeof(UpDownCblk); */
  /*     for (iter=0; iter<solvptr->cblknbr; iter++) */
  /*       { */
  /*         /\* browproctab / browcblktab *\/ */
  /*         result += 2*solvptr->updovct.cblktab[iter].browprocnbr*sizeof(pastix_int_t); */
  /*       } */
  /*     /\* gcblk2list *\/ */
  /*     result += solvptr->updovct.gcblk2listnbr*sizeof(pastix_int_t); */
  /*     /\* listptr *\/ */
  /*     result += solvptr->updovct.listptrnbr   *sizeof(pastix_int_t); */
  /*     /\* listcblk / listblok *\/ */
  /*     result += 2*solvptr->updovct.listnbr    *sizeof(pastix_int_t); */
  /*     /\* loc2glob *\/ */
  /*     result += solvptr->updovct.loc2globnbr  *sizeof(pastix_int_t); */
  /*     /\* lblk2gcblk *\/ */
  /*     result += solvptr->bloknbr              *sizeof(pastix_int_t); */
  /*   } */

  return result;
}
