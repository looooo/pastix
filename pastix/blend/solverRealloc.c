#include <stdlib.h>
#include <stdio.h>


#include "common.h"
#include "ftgt.h"
#include "queue.h"
#include "bulles.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
/* #include "assert.h" */
#include "solverRealloc.h"

void setBcofPtr(SolverMatrix *solvptr, const pastix_int_t * bcofind)
{
  pastix_int_t i;
  for(i=0;i<solvptr->btagnbr;i++)
    solvptr->btagtab[i].bcofptr = &(solvptr->bcoftab[bcofind[i]]);
}

void setLocalBtagPtr(SolverMatrix *solvptr)
{
  pastix_int_t i, task, taskcnt;

  for (i=0; i<solvptr->btagnbr; i++)
    {
      if (solvptr->proc2clust[solvptr->btagtab[i].infotab[BTAG_PROCDST]] == solvptr->clustnum)
        {
          task = solvptr->btagtab[i].infotab[BTAG_TASKDST];

          /*
            for (j=0;j<SOLV_TASKNBR;j++)
            if (TASK_TASKNEXT(j)==task)
            {
            fprintf(stderr,"local Fuck!!\n");
            fprintf(stderr,"%ld %ld %ld %ld %ld %ld\n",i,j,task,TASK_TASKNEXT(j),SOLV_BTAGNBR,SOLV_TASKNBR);
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
            }
          */
          ASSERT(solvptr->tasktab[task].btagptr == NULL, MOD_SOPALIN);
          if (solvptr->tasktab[task].btagptr == NULL)
            {
              taskcnt = 1;
              solvptr->tasktab[task].btagptr = &(solvptr->btagtab[i]);

              while (solvptr->tasktab[task].tasknext != solvptr->btagtab[i].infotab[BTAG_TASKDST])
                {
                  task = solvptr->tasktab[task].tasknext;
                  solvptr->tasktab[task].btagptr = &(solvptr->btagtab[i]);
                  taskcnt++;
                }
              /* BTAG_TASKCNT(i)=taskcnt; ??? */
              ASSERTDBG(taskcnt == solvptr->btagtab[i].infotab[BTAG_TASKCNT], MOD_SOPALIN);
            }
        }
    }

  /* Compute the number of block to recv */
  {
    Queue q;
    queueInit(&q, solvptr->btgsnbr * 2); /* Size arbitrary choosen */
    for(i=0; i<solvptr->tasknbr; i++)
      {
        if (((solvptr->tasktab[i].taskid == E1) ||
             (solvptr->tasktab[i].taskid == E2)) &&
            (solvptr->tasktab[i].btagptr == NULL))
          {
            /* We take the task_master to avoid to count multiple block for one communication */
            if (queuePossess(&q, solvptr->tasktab[i].taskmstr) == API_NO)
              {
                queueAdd(&q, solvptr->tasktab[i].taskmstr, (double)solvptr->tasktab[i].taskmstr);
              }
          }
      }
    solvptr->btgrnbr = queueSize(&q);
    queueExit(&q);
  }
}

/*+ Realloc the solver matrix in a contiguous way +*/
void solverRealloc(SolverMatrix *solvmtx, pastix_int_t *bcofind)
{
    SolverMatrix *tmp;
    pastix_int_t i;

    MALLOC_INTERN(tmp, 1, SolverMatrix);
    /** copy general info **/
    memcpy(tmp, solvmtx, sizeof(SolverMatrix));

    /**OIMBE il faudra faire le REALLOC pour ooc ! **/

    /** Copy tasktab **/
    MALLOC_INTERN(solvmtx->tasktab, solvmtx->tasknbr, Task);
    memcpy(solvmtx->tasktab, tmp->tasktab, solvmtx->tasknbr*sizeof(Task));
#ifdef DEBUG_BLEND
    for(i=0;i<solvmtx->tasknbr;i++)
      ASSERT((solvmtx->tasktab[i].btagptr == NULL), MOD_BLEND);
#endif

    /** Copy cblktab and bloktab **/
    MALLOC_INTERN(solvmtx->cblktab, solvmtx->cblknbr+1, SolverCblk);
    memcpy(solvmtx->cblktab, tmp->cblktab,
           (solvmtx->cblknbr+1)*sizeof(SolverCblk));

    MALLOC_INTERN(solvmtx->bloktab, solvmtx->bloknbr, SolverBlok);
    memcpy(solvmtx->bloktab, tmp->bloktab,
           solvmtx->bloknbr*sizeof(SolverBlok));

    /** Copy BlockTarget tab **/
    if(solvmtx->btagnbr>0)
      {
        MALLOC_INTERN(solvmtx->btagtab, solvmtx->btagnbr, BlockTarget);
        memcpy(solvmtx->btagtab, tmp->btagtab,
               solvmtx->btagnbr*sizeof(BlockTarget));
      }
    else
      {
        solvmtx->btagtab = NULL;
      }

    /** Copy BlockCoef tab **/
    if(solvmtx->bcofnbr>0)
      {
        MALLOC_INTERN(solvmtx->bcoftab, solvmtx->bcofnbr, BlockCoeff);
        memcpy(solvmtx->bcoftab, tmp->bcoftab,
               solvmtx->bcofnbr*sizeof(BlockCoeff));
      }
    else
      {
        solvmtx->bcoftab = NULL;
      }

    /** Copy ftgttab **/
    if (solvmtx->ftgtnbr != 0)
      {
        MALLOC_INTERN(solvmtx->ftgttab, solvmtx->ftgtnbr, FanInTarget);
        memcpy(solvmtx->ftgttab, tmp->ftgttab,
               solvmtx->ftgtnbr*sizeof(FanInTarget));
      }
    /** copy infotab of fan intarget **/
    /*for(i=0;i<tmp->ftgtnbr;i++)
      memcpy(solvmtx->ftgttab[i].infotab, tmp->ftgttab[i].infotab, MAXINFO*sizeof(pastix_int_t));*/

    /** Copy indtab **/
    MALLOC_INTERN(solvmtx->indtab, solvmtx->indnbr, pastix_int_t);
    memcpy(solvmtx->indtab, tmp->indtab, solvmtx->indnbr*sizeof(pastix_int_t));


    /** Copy ttsktab & ttsknbr **/
    if (solvmtx->bublnbr>0)
      {
        MALLOC_INTERN(solvmtx->ttsknbr, solvmtx->bublnbr, pastix_int_t);
        memcpy(solvmtx->ttsknbr, tmp->ttsknbr, solvmtx->bublnbr*sizeof(pastix_int_t));
        MALLOC_INTERN(solvmtx->ttsktab, solvmtx->bublnbr, pastix_int_t*);
        for (i=0;i<solvmtx->bublnbr;i++)
          {
            solvmtx->ttsktab[i] = NULL;
            MALLOC_INTERN(solvmtx->ttsktab[i], solvmtx->ttsknbr[i], pastix_int_t);
            memcpy(solvmtx->ttsktab[i], tmp->ttsktab[i],
                   solvmtx->ttsknbr[i]*sizeof(pastix_int_t));
          }
      }
    else
      {
        solvmtx->ttsknbr = NULL;
        solvmtx->ttsktab = NULL;
      }

    MALLOC_INTERN(solvmtx->proc2clust, solvmtx->procnbr, pastix_int_t);
    memcpy(solvmtx->proc2clust, tmp->proc2clust,
           solvmtx->procnbr * sizeof(pastix_int_t));

    /** Free the former solver matrix **/
    solverExit(tmp);
    memFree_null(tmp);

    /** Set the bcofptr **/
    setBcofPtr(solvmtx, bcofind);

    /** Set the local btagptr **/
    setLocalBtagPtr(solvmtx);

}


void solverExit(SolverMatrix *solvmtx)
{
    pastix_int_t i;

    /** Free arrays of solvmtx **/
    if(solvmtx->cblktab)
      {
        for (i = 0; i < solvmtx->cblknbr; i++)
          {
            if (solvmtx->cblktab[i].coeftab)
              memFree_null(solvmtx->cblktab[i].coeftab);

            if (solvmtx->cblktab[i].ucoeftab)
              memFree_null(solvmtx->cblktab[i].ucoeftab);
          }
        memFree_null(solvmtx->cblktab);
      }
    if(solvmtx->bloktab)
      memFree_null(solvmtx->bloktab);
    /*if(solvmtx->coeftab)
      memFree_null(solvmtx->coeftab);*/
    if(solvmtx->ftgttab)
      memFree_null(solvmtx->ftgttab);
    if(solvmtx->btagtab)
      memFree_null(solvmtx->btagtab);
    if(solvmtx->bcoftab)
      memFree_null(solvmtx->bcoftab);
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

}


void solverInit(SolverMatrix *solvmtx)
{
  solvmtx->cblktab = NULL;
  solvmtx->bloktab = NULL;
  solvmtx->coefnbr = 0;
  solvmtx->ftgtnbr = 0;
  solvmtx->btagtab = NULL;
  solvmtx->bcoftab = NULL;

  solvmtx->ftgttab = NULL;
  solvmtx->cpftmax = 0;
  solvmtx->bpftmax = 0;
  solvmtx->coefmax = 0;
  memset(solvmtx, 0, sizeof (SolverMatrix));


  solvmtx->baseval = 0;
  solvmtx->cblknbr = 0;
  solvmtx->bloknbr = 0;
  solvmtx->nodenbr = 0;

}
