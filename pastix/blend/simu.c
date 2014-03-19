#include <stdio.h>
#include <strings.h>
#include <math.h>

#include "common.h"
#include "ftgt.h"
#include "symbol.h"
#include "elimin.h"
#include "cost.h"
#include "cand.h"
#include "queue.h"
#include "extendVector.h"
#include "blendctrl.h"
#include "simu.h"

#define TIMEBASE 10.0

pastix_int_t
simuInit( SimuCtrl     *simuctrl,
          SymbolMatrix *symbptr,
          pastix_int_t  clustnbr,
          pastix_int_t  procnbr,
          pastix_int_t  cblknbr,
          pastix_int_t  bloknbr,
          Cand         *candtab)
{
    pastix_int_t i, j;
    pastix_int_t p;
    pastix_int_t ftgtcur;
    pastix_int_t candnbr;
    pastix_int_t step;

    simuctrl->cblknbr  = cblknbr;
    simuctrl->ftgtprio = 0;
    simuctrl->tasktab  = NULL;
    simuctrl->ftgtcnt  = 0;

    /** Processor initialisation **/
    MALLOC_INTERN(simuctrl->proctab, procnbr, SimuProc);
    for(i=0;i<procnbr;i++)
    {
        timerSet(TIMER(i), 0.0); /** for paragraph numeric tolerance **/
        simuctrl->proctab[i].prionum   = 0;
        MALLOC_INTERN(simuctrl->proctab[i].futuretask, 1, pastix_queue_t);
        MALLOC_INTERN(simuctrl->proctab[i].readytask,  1, pastix_queue_t);
        pqueueInit(simuctrl->proctab[i].futuretask, 100);
        pqueueInit(simuctrl->proctab[i].readytask,  100);

        MALLOC_INTERN(simuctrl->proctab[i].tasktab, 1, ExtendVectorINT);
        extendint_Init(simuctrl->proctab[i].tasktab, bloknbr/procnbr + 1);
    }

    /** Cluster initialization **/
    MALLOC_INTERN(simuctrl->clustab, clustnbr, SimuCluster);
    step = procnbr / clustnbr;
    for(i=0;i<clustnbr;i++)
    {
        simuctrl->clustab[i].fprocnum = i*step;
        simuctrl->clustab[i].lprocnum = simuctrl->clustab[i].fprocnum + step - 1;
        MALLOC_INTERN(simuctrl->clustab[i].ftgtsend, clustnbr, ExtendVectorINT);
        simuctrl->clustab[i].prionum  = 0;
        for(p=0;p<clustnbr;p++)
            extendint_Init(&(simuctrl->clustab[i].ftgtsend[p]), cblknbr/(2*clustnbr)+1);
    }
    simuctrl->clustab[clustnbr-1].lprocnum = procnbr-1;

    MALLOC_INTERN(simuctrl->ownetab, cblknbr, pastix_int_t);
    MALLOC_INTERN(simuctrl->blprtab, bloknbr, pastix_int_t);

    /* affect a negative value to cblk not mapped */
    for(i=0;i<cblknbr;i++)
        simuctrl->ownetab[i] = -1;
    for(i=0;i<bloknbr;i++)
        simuctrl->blprtab[i] = -1;

    MALLOC_INTERN(simuctrl->cblktab, cblknbr+1, SimuCblk);
    MALLOC_INTERN(simuctrl->bloktab, bloknbr+1, SimuBlok);
    ftgtcur = 0;

    for(i=0;i<cblknbr;i++)
    {
        candnbr = candtab[i].lccandnum - candtab[i].fccandnum + 1;
        simuctrl->cblktab[i].ctrbcnt = 0;

        for(j=symbptr->cblktab[i].bloknum;j<symbptr->cblktab[i+1].bloknum;j++)
        {
            simuctrl->bloktab[j].ftgtnum   = ftgtcur;
            simuctrl->bloktab[j].tasknum   = -1;
            simuctrl->bloktab[j].fccandnum = candtab[i].fccandnum;
            simuctrl->bloktab[j].ctrbcnt   = 0;
            /*if(candnbr > 1)*/
            ftgtcur += candnbr;
        }
    }
    /* one extracblk for avoiding side effect */
    simuctrl->bloktab[bloknbr].ftgtnum = ftgtcur;
    simuctrl->ftgtnbr = ftgtcur;

    if(simuctrl->ftgtnbr > 0)
    {
        /** Allocate and Initialize the timer for the reception of each ftgt on a candidate cluster **/
        MALLOC_INTERN(simuctrl->ftgttimetab, simuctrl->ftgtnbr, SimuTimer);
        for(i=0;i<simuctrl->ftgtnbr;i++)
            timerSet(&(simuctrl->ftgttimetab[i]), 0.0);

        MALLOC_INTERN(simuctrl->ftgttab, ftgtcur, SimuFtgt);
        for(i=0;i<simuctrl->ftgtnbr;i++)
        {
            simuctrl->ftgttab[i].clustnum = -1;
            timerSet(&(simuctrl->ftgttab[i].timerecv), 0.0);
            simuctrl->ftgttab[i].costsend = 0.0;
            simuctrl->ftgttab[i].costadd  = 0.0;
            bzero(simuctrl->ftgttab[i].ftgt.infotab,MAXINFO*sizeof(pastix_int_t));
            simuctrl->ftgttab[i].ftgt.infotab[FTGT_FCOLNUM] = INTVALMAX;
            simuctrl->ftgttab[i].ftgt.infotab[FTGT_FROWNUM] = INTVALMAX;
            simuctrl->ftgttab[i].ftgt.infotab[FTGT_CTRBNBR] = 0;
            simuctrl->ftgttab[i].ftgt.infotab[FTGT_CTRBCNT] = 0;
        }
    }
    else
    {
        simuctrl->ftgttab     = NULL;
        simuctrl->tasktimetab = NULL;
        simuctrl->ftgttimetab = NULL;
    }

    return 1;
}


pastix_int_t simuRealloc(SimuCtrl *simuctrl, pastix_int_t procnbr, pastix_int_t local_nbthrds)
{
    pastix_int_t i;

    /* Free processor structure */
    for(i=0;i<procnbr;i++)
    {
        pqueueExit    (simuctrl->proctab[i].readytask);
        memFree_null  (simuctrl->proctab[i].readytask);
        pqueueExit    (simuctrl->proctab[i].futuretask);
        memFree_null  (simuctrl->proctab[i].futuretask);
        extendint_Exit(simuctrl->proctab[i].tasktab);
        memFree_null  (simuctrl->proctab[i].tasktab);
    }
    memFree_null(simuctrl->proctab);

    /** Initialisation for local thread **/
    MALLOC_INTERN(simuctrl->proctab, local_nbthrds, SimuProc);
    for(i=0;i<local_nbthrds;i++)
    {
        MALLOC_INTERN(simuctrl->proctab[i].tasktab, 1, ExtendVectorINT);
        /* On initialise pas les vecteur d'entier, car il nous faudrait la structure d'arbre de bulles */
    }

    return 1;
}

void simuExit(SimuCtrl *simuctrl, pastix_int_t clustnbr, pastix_int_t procnbr, pastix_int_t local_nbthrds)
{
    pastix_int_t i,j;
    (void)local_nbthrds; (void)procnbr;

#ifndef PASTIX_DYNSCHED
    for(i=0;i<procnbr;i++)
    {
        pqueueExit    (simuctrl->proctab[i].readytask);
        memFree_null  (simuctrl->proctab[i].readytask);
        pqueueExit    (simuctrl->proctab[i].futuretask);
        memFree_null  (simuctrl->proctab[i].futuretask);
        extendint_Exit(simuctrl->proctab[i].tasktab);
        memFree_null  (simuctrl->proctab[i].tasktab);
    }
#else
    for(i=0;i<local_nbthrds;i++)
    {
        extendint_Exit(simuctrl->proctab[i].tasktab);
        memFree_null(simuctrl->proctab[i].tasktab);
    }
#endif

    for(i=0;i<clustnbr;i++)
    {
        for(j=0;j<clustnbr;j++)
            extendint_Exit(&(simuctrl->clustab[i].ftgtsend[j]));
        memFree_null(simuctrl->clustab[i].ftgtsend);
    }

    if(simuctrl->ftgttab != NULL)
    {
        memFree_null(simuctrl->ftgttab);
        memFree_null(simuctrl->ftgttimetab);
    }
    /* memFree_null(simuctrl->tasktimetab); */
    memFree_null(simuctrl->tasktab);
    memFree_null(simuctrl->proctab);
    memFree_null(simuctrl->clustab);
    memFree_null(simuctrl->ownetab);
    memFree_null(simuctrl->blprtab);
    memFree_null(simuctrl->cblktab);
    memFree_null(simuctrl->bloktab);
    memFree_null(simuctrl);
}


pastix_int_t compTimer(SimuTimer *t1, SimuTimer *t2)
{
    /** Return 1 if t1 < t2 **/
    /** 0 in other cases **/
    if(t1->s < t2->s)
        return 1;
    /*if(t1->s == t2->s && t1->ms < t2->ms)
     return 1;*/
    return 0;
}

void timerAdd(SimuTimer *timer, double t)
{
    timer->s += t;

    /*timer->s += floor(t);
     timer->ms += (t-floor(t))*TIMEBASE;
     if(timer->ms >= TIMEBASE)
     {
     timer->s += 1;
     timer->ms = timer->ms - TIMEBASE;
     }*/
}

double timerVal(const SimuTimer *t)
{
    /*#ifdef DEBUG_BLEND
     ASSERT(t->ms < TIMEBASE,MOD_BLEND);
     #endif*/
    return (t->s /* + t->ms/TIMEBASE */);
}

void timerSet(SimuTimer *timer, double t)
{
    timer->s = t;
    /*timer->s = floor(t);
     timer->ms = (t-floor(t))*TIMEBASE;*/
}

void timerSetMax(SimuTimer *timer, double t)
{
    if ( t > timer->s ) {
        timer->s = t;
    }
}


/** OIMBE pour optimisation faire un timerMax pour distrib **/
/** OIMBE pour optimisation faire un timerAffect(timer) pour distrib **/
