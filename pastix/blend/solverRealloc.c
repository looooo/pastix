#include <stdlib.h>
#include <stdio.h>


#include "common.h"
#include "ftgt.h"
#include "queue.h"
#include "bulles.h"
#include "updown.h"
#include "solver.h"
/* #include "assert.h" */
#include "solverRealloc.h"

/*+ Realloc the solver matrix in a contiguous way +*/
void solverRealloc(SolverMatrix *solvmtx)
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

#if defined(PASTIX_WITH_STARPU)
    MALLOC_INTERN(solvmtx->hcblktab, solvmtx->hcblknbr+1, SolverCblk);
    memcpy(solvmtx->hcblktab, tmp->hcblktab,
           (solvmtx->hcblknbr+1)*sizeof(SolverCblk));
    MALLOC_INTERN(solvmtx->hbloktab, solvmtx->hcblktab[solvmtx->hcblknbr].bloknum,
                  SolverBlok);
    memcpy(solvmtx->hbloktab, tmp->hbloktab,
           solvmtx->hcblktab[solvmtx->hcblknbr].bloknum*sizeof(SolverBlok));
    MALLOC_INTERN(solvmtx->gcblk2halo, solvmtx->gcblknbr, pastix_int_t);
    memcpy(solvmtx->gcblk2halo, tmp->gcblk2halo,
           solvmtx->gcblknbr*sizeof(pastix_int_t));
#endif /* defined(PASTIX_WITH_STARPU) */

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
#endif
}


void solverInit(SolverMatrix *solvmtx)
{
  solvmtx->cblktab = NULL;
  solvmtx->bloktab = NULL;
  solvmtx->coefnbr = 0;
  solvmtx->ftgtnbr = 0;

  solvmtx->ftgttab = NULL;
  solvmtx->coefmax = 0;
  memset(solvmtx, 0, sizeof (SolverMatrix));

  solvmtx->baseval = 0;
  solvmtx->cblknbr = 0;
  solvmtx->bloknbr = 0;
  solvmtx->nodenbr = 0;

}
