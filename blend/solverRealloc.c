#include <stdlib.h>
#include <stdio.h>


#include "common.h"
#include "ftgt.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"
/* #include "assert.h" */
#include "solverRealloc.h"

/*+ Realloc the solver matrix in a contiguous way +*/
void solverRealloc(SolverMatrix *solvmtx)
{
    SolverMatrix *tmp;
    SolverBlok   *solvblok;
    SolverCblk   *solvcblk;
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

    MALLOC_INTERN(solvmtx->browtab, solvmtx->brownbr, pastix_int_t);
    memcpy(solvmtx->browtab, tmp->browtab,
           solvmtx->brownbr*sizeof(pastix_int_t));

    solvblok = solvmtx->bloktab;
    for (solvcblk = solvmtx->cblktab; solvcblk  < solvmtx->cblktab + solvmtx->cblknbr; solvcblk++) {
        pastix_int_t bloknbr = (solvcblk+1)->fblokptr - solvcblk->fblokptr;
        solvcblk->fblokptr = solvblok;
        solvblok+= bloknbr;
    }
    solvcblk->fblokptr = solvblok;

#if defined(PASTIX_WITH_STARPU)
    if ( tmp->gcblk2halo ) {
        MALLOC_INTERN(solvmtx->gcblk2halo, solvmtx->gcblknbr, pastix_int_t);
        memcpy(solvmtx->gcblk2halo, tmp->gcblk2halo,
               solvmtx->gcblknbr*sizeof(pastix_int_t));
    }
    if ( tmp->hcblktab ) {
        MALLOC_INTERN(solvmtx->hcblktab, solvmtx->hcblknbr+1, SolverCblk);
        memcpy(solvmtx->hcblktab, tmp->hcblktab,
               (solvmtx->hcblknbr+1)*sizeof(SolverCblk));
        MALLOC_INTERN(solvmtx->hbloktab, tmp->hcblktab[tmp->hcblknbr].fblokptr - tmp->hbloktab,
                      SolverBlok);
        memcpy(solvmtx->hbloktab, tmp->hbloktab,
               (tmp->hcblktab[tmp->hcblknbr].fblokptr - tmp->hbloktab)*sizeof(SolverBlok));

        solvblok = solvmtx->hbloktab;
        for (solvcblk = solvmtx->hcblktab;
             solvcblk  < solvmtx->hcblktab + solvmtx->hcblknbr;
             solvcblk++) {
            pastix_int_t bloknbr = (solvcblk+1)->fblokptr - solvcblk->fblokptr;
            solvcblk->fblokptr = solvblok;
            solvblok+= bloknbr;
        }
        solvcblk->fblokptr = solvblok;
    }

    if (pastix_starpu_with_fanin() == API_YES) {
        /* FANIN info */
        pastix_int_t clustnum;
        MALLOC_INTERN(solvmtx->fcblktab, solvmtx->clustnbr, SolverCblk*);
        MALLOC_INTERN(solvmtx->fbloktab, solvmtx->clustnbr, SolverBlok*);
        MALLOC_INTERN(solvmtx->fcblknbr, solvmtx->clustnbr, pastix_int_t);
        memset(solvmtx->fcblknbr, 0, solvmtx->clustnbr*sizeof(pastix_int_t));
        memset(solvmtx->fcblktab, 0, solvmtx->clustnbr*sizeof(SolverCblk*));
        memset(solvmtx->fbloktab, 0, solvmtx->clustnbr*sizeof(SolverBlok*));

        memcpy(solvmtx->fcblknbr, tmp->fcblknbr,
               solvmtx->clustnbr*sizeof(pastix_int_t));
        for (clustnum = 0; clustnum < solvmtx->clustnbr; clustnum++) {
            pastix_int_t bloknbr;
            MALLOC_INTERN(solvmtx->fcblktab[clustnum],
                          solvmtx->fcblknbr[clustnum]+1,
                          SolverCblk);
            if ( solvmtx->fcblknbr[clustnum] > 0 ) {
                memcpy(solvmtx->fcblktab[clustnum],
                       tmp->fcblktab[clustnum],
                       (solvmtx->fcblknbr[clustnum]+1)*sizeof(SolverCblk));
                bloknbr = tmp->fcblktab[clustnum][tmp->fcblknbr[clustnum]].fblokptr - tmp->fbloktab[clustnum];
                MALLOC_INTERN(solvmtx->fbloktab[clustnum],
                              bloknbr,
                              SolverBlok);
                memcpy(solvmtx->fbloktab[clustnum],
                       tmp->fbloktab[clustnum],
                       (bloknbr)*sizeof(SolverBlok));
                solvblok = solvmtx->fbloktab[clustnum];
                for (solvcblk = solvmtx->fcblktab[clustnum];
                     solvcblk  < solvmtx->fcblktab[clustnum] +
                         solvmtx->fcblknbr[clustnum];
                     solvcblk++) {

                    pastix_int_t bloknbr = (solvcblk+1)->fblokptr - solvcblk->fblokptr;
                    solvcblk->fblokptr = solvblok;
                    solvblok+= bloknbr;
                }
                solvcblk->fblokptr = solvblok;

            }

        }
    }

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
            if (solvmtx->cblktab[i].lcoeftab)
                memFree_null(solvmtx->cblktab[i].lcoeftab);

            if (solvmtx->cblktab[i].ucoeftab)
                memFree_null(solvmtx->cblktab[i].ucoeftab);
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


void solverInit(SolverMatrix *solvmtx)
{
    memset(solvmtx, 0, sizeof (SolverMatrix));
    return;
}
