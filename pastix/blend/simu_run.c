#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "common.h"
#include "symbol.h"
#include "extendVector.h"
#include "queue.h"
#include "ftgt.h"
#include "elimin.h"
#include "cost.h"
#include "cand.h"
#include "bulles.h"
#include "blendctrl.h"
#include "simu.h"
#include "dof.h"
#include "task.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "costfunc.h"

/** OIMBE a t on  tjs besoins de egraph->ownetab ???? **/

static inline void
computeBlockCtrbNbr(const BlendCtrl    *ctrl,
                    const SymbolMatrix *symbptr,
                          SimuCtrl     *simuctrl )
{
    pastix_int_t i, j, k;
    pastix_int_t facebloknum;

    /* Setup the block ctrbcnt only if 2D is enabled */
    {
        SymbolCblk *curcblk;
        SimuBlok   *curblok;

        curcblk = symbptr->cblktab;
        curblok = simuctrl->bloktab;
        for(i=0; i<symbptr->cblknbr; i++, curcblk++)
        {
            pastix_int_t fbloknum = curcblk[0].bloknum;
            pastix_int_t lbloknum = curcblk[1].bloknum;

            /* Skip diagonal block */
            curblok++;

            /* 1D cblk computed */
            for(j=fbloknum+1; j<lbloknum; j++)
            {
                facebloknum = 0;
                /* Add contribution due to E2 */
                for(k=j; k<lbloknum; k++)
                {
                    /* We don't care if facing task is 1D or 2D */
                    facebloknum = getFaceBlockE2(facebloknum, j, k, symbptr, ctrl->ricar);

                    if(facebloknum >= 0)
                        simuctrl->bloktab[facebloknum].ctrbcnt++;
                }
            }
        }
    }

    /* Set up the task ctrbcnt and cblkcnt */
    {
        SimuTask *task = simuctrl->tasktab;

        for(i=0;i<simuctrl->tasknbr;i++)
        {
            pastix_int_t fbloknum = symbptr->cblktab[task->cblknum  ].bloknum;
            pastix_int_t lbloknum = symbptr->cblktab[task->cblknum+1].bloknum;

            task->ctrbcnt = 0;
            for(j=fbloknum; j<lbloknum; j++)
                task->ctrbcnt += simuctrl->bloktab[j].ctrbcnt;

            simuctrl->cblktab[task->cblknum].ctrbcnt = task->ctrbcnt;
            task++;
        }
    }
}

static inline void
simu_putInAllReadyQueues(const BlendCtrl *ctrl,
                         SimuCtrl        *simuctrl,
                         pastix_int_t     tasknum )
{
    /*---------------------------------------------------------/
     / This function according to the ready date of a task      /
     / put this task on the ready queue of a processor          /
     / NOTE: when the ready date of a task is inferior to the   /
     / proc timer then he task is ordered according to its      /
     / priorities in the elimination tree                       /
     /---------------------------------------------------------*/
    const SimuTask *task     = simuctrl->tasktab + tasknum;
    const Cand     *cblkcand = ctrl->candtab + task->cblknum;
    double ready_date = 0.0;
    pastix_int_t procnum;
    pastix_int_t bloknum = task->bloknum;
    pastix_int_t treelevel = cblkcand->treelevel;

    /* Get the ready date of the task on the processor passed in parameter */
    if( cblkcand->fccandnum == cblkcand->lccandnum )
    {
        ready_date = timerVal( &(task->time) );

        for(procnum =  cblkcand->fcandnum;
            procnum <= cblkcand->lcandnum; procnum++)
        {
            if(ready_date > timerVal(TIMER(procnum)))
                queueAdd2(simuctrl->proctab[procnum].taskheap, tasknum, ready_date, treelevel);
            else
                queueAdd2(simuctrl->proctab[procnum].taskheap2, tasknum, (double)treelevel, bloknum);
        }
    }
    else
    {
        for(procnum =  cblkcand->fcandnum;
            procnum <= cblkcand->lcandnum; procnum++)
        {
            ready_date = timerVal( simuctrl->ftgttimetab + CLUST2INDEX(bloknum, ctrl->core2clust[procnum]) );

            if( ready_date > timerVal(TIMER(procnum)) )
                queueAdd2(simuctrl->proctab[procnum].taskheap, tasknum, ready_date, treelevel);
            else
                queueAdd2(simuctrl->proctab[procnum].taskheap2, tasknum, (double)treelevel, bloknum);
        }
    }
}

static inline pastix_int_t
getNextTaskNextProc(SimuCtrl *simuctrl, BlendCtrl *ctrl, pastix_int_t *procnumptr)
{
    /*----------------------------------------------------------------------------------------------------/
     /  Get the next task and the next proc in order that they are the first that can compute something    /
     /  On return : earlier task index                                                                     /
     /----------------------------------------------------------------------------------------------------*/

    pastix_int_t p;
    pastix_int_t procnum = -1;
    pastix_int_t tasknum;
    double earlytimeready = INTVALMAX;
    double earlyproctimer = INTVALMAX;
    double timeready;
    pastix_int_t earlytask = -1;

    /** Find the earlier task in the processor heaps **/
    for(p=0;p<ctrl->total_nbcores;p++)
    {
        tasknum = -1;
        /** First we search the earlier task in the set of task whose ready date is < proc timer **/
        while(queueSize(simuctrl->proctab[p].taskheap2)>0)
        {
            tasknum = queueRead(simuctrl->proctab[p].taskheap2);
            if( simuctrl->blprtab[simuctrl->tasktab[tasknum].bloknum]>=0 )
            {
                /** This task have to be remove from the heap (already mapped) **/
                queueGet(simuctrl->proctab[p].taskheap2);
                tasknum = -1;
            }
            else
                break;
        }
        /** We found no task which ready date is < proc timer so we search one that minimizes ready date - proc-timer **/
        if(tasknum == -1)
        {
            while(queueSize(simuctrl->proctab[p].taskheap)>0)
            {
                tasknum = queueRead(simuctrl->proctab[p].taskheap);
                if( simuctrl->blprtab[simuctrl->tasktab[tasknum].bloknum]>=0 )
                {
                    /** This task have to be remove from the heap (already mapped) **/
                    queueGet(simuctrl->proctab[p].taskheap);
                    tasknum = -1;
                }
                else
                    break;
            }
        }

        if(tasknum != -1)
        {
            timeready = MAX(timerVal(TIMER(p)), timerVal(&(simuctrl->ftgttimetab[CLUST2INDEX(simuctrl->tasktab[tasknum].bloknum, ctrl->core2clust[p])])));

            /** We prevent to distribute on the same processor set when all time are equals **/
            if((timeready == earlytimeready) && (timerVal(TIMER(p)) < earlyproctimer))
            {
                procnum = p;
                earlyproctimer = timerVal(TIMER(p));
                earlytask = tasknum;
                earlytimeready = timeready;
            }

            if(timeready < earlytimeready)
            {
                procnum  = p;
                earlytask = tasknum;
                earlytimeready = timeready;
            }
        }
    }
    if(procnum != -1)
    {
        if(queueSize(simuctrl->proctab[procnum].taskheap2)>0)
        {assert(earlytask == queueGet(simuctrl->proctab[procnum].taskheap2));}
        else
            assert(earlytask == queueGet(simuctrl->proctab[procnum].taskheap));
    }
    *procnumptr = procnum;
    return earlytask;
}

static inline void
computeTaskReceiveTime(const pastix_int_t tasknum, SymbolMatrix *symbptr, SimuCtrl *simuctrl, BlendCtrl *ctrl, const Dof * dofptr)
{
    /*-------------------------------------------------------------------------/
     / Compute the time the cblk would have RECEIVED and ADDED                  /
     / all its contributions if it was mapped on a given cand CLUSTER           /
     / !! These times not include add time for fan in target !!                 /
     /--------------------------------------------------------------------------*/

    pastix_int_t i, j;
    double lftgttime = 0;
    double sftgttime = 0;
    pastix_int_t   lftgtnum  = -1;
    pastix_int_t   cblknum;
    pastix_int_t   bloknum;
    pastix_int_t   clustdst;

    bloknum = simuctrl->tasktab[tasknum].bloknum;
    cblknum = simuctrl->tasktab[tasknum].cblknum;

    /* no fan_in_target-> no need treatment */
    if(ctrl->candtab[cblknum].fccandnum == ctrl->candtab[cblknum].lccandnum)
        return;

    /*------------------------------------------------------------------------------------------------/
     / Compute the cblk on proc timer that is time the cblk would have received                        /
     / all its contributions if it was mapped on a given cand processor                                /
     / These times INCLUDE add time for fan in target !!                                               /
     /------------------------------------------------------------------------------------------------*/

    /** Compute receive time (time at which a non-local processor should received the target **/
    /* find the latest ftgt receive time and the second latest*/
    for(i=simuctrl->bloktab[bloknum].ftgtnum; i<simuctrl->bloktab[bloknum+1].ftgtnum;i++)
    {
        /* Source of this ftgt */
        clustdst = INDEX2CLUST(i,bloknum);

        /** Task COMP_1D with several cand proc **/
        /** The information about ftgt costs are in the ftgt of the diagonal block;
         this loop sums the cost of all the ftgt received by the blocks in this column block **/
        if(simuctrl->ftgttab[i].ftgt.infotab[FTGT_CTRBNBR]>0)
            for(j=bloknum;j<symbptr->cblktab[cblknum+1].bloknum;j++)
            {
                if(simuctrl->ftgttab[simuctrl->bloktab[j].ftgtnum + i-simuctrl->bloktab[bloknum].ftgtnum].ftgt.infotab[FTGT_CTRBNBR]>0)
                {
                    simuctrl->ftgttab[i].costadd +=
                        costFtgtAdd(&(simuctrl->ftgttab[CLUST2INDEX(j, clustdst)].ftgt), dofptr);

                    simuctrl->ftgttab[i].costsend +=
                        costFtgtSend(clustdst, ctrl->candtab[cblknum].lccandnum-ctrl->candtab[cblknum].fccandnum+1,
                                     &(simuctrl->ftgttab[CLUST2INDEX(j, clustdst)].ftgt), ctrl, dofptr);
                }
            }

#ifdef DEBUG_BLEND
        if(!(simuctrl->ftgttab[i].costsend >= 0.0))
            errorPrint("ftgt %ld costsend %f", (long)i, simuctrl->ftgttab[i].costsend);
        if(!(simuctrl->ftgttab[i].costadd >= 0.0))
            errorPrint("ftgt %ld costadd %f", (long)i, simuctrl->ftgttab[i].costadd);

        assert(simuctrl->ftgttab[i].costsend >= 0.0);
        assert(simuctrl->ftgttab[i].costadd >= 0.0);
#endif

        /** ftgttab[].timerecv is the time this ftgt will be receive **/
        timerSet(&(simuctrl->ftgttab[i].timerecv), timerVal(&(simuctrl->ftgttimetab[i])) + simuctrl->ftgttab[i].costsend + simuctrl->ftgttab[i].costadd);

        /** If this ftgt the last reveived or the second last received ?? **/
        if(timerVal(&(simuctrl->ftgttab[i].timerecv)) > lftgttime)
        {
            lftgttime = timerVal(&(simuctrl->ftgttab[i].timerecv));
            lftgtnum  = i;
        }
        else
            if(timerVal(&(simuctrl->ftgttab[i].timerecv)) > sftgttime)
                sftgttime = timerVal(&(simuctrl->ftgttab[i].timerecv));
    }


    /*------------------------------------------------------/
     / Put in ftgttimetab[] the date at which the cluster    /
     / would have received and add all the ftgt if the task  /
     /  was mapped on it                                     /
     /------------------------------------------------------*/
    for(i=simuctrl->bloktab[bloknum].ftgtnum; i<simuctrl->bloktab[bloknum+1].ftgtnum;i++)
    {
        if(i != lftgtnum)
            timerSet(&(simuctrl->ftgttimetab[i]), lftgttime);
        else
            timerSet(&(simuctrl->ftgttimetab[i]), MAX(timerVal(&(simuctrl->ftgttimetab[i])), sftgttime));
    }
}

static inline void
updateFtgtStruct(pastix_int_t bloknum, pastix_int_t bloknum2, pastix_int_t ftgtnum, SymbolMatrix *symbptr, SimuCtrl *simuctrl, BlendCtrl *ctrl)
{
    SymbolBlok * blokptr;
    SymbolBlok * blokptr2;
    (void)ctrl;

    blokptr  = &(symbptr->bloktab[bloknum]);
    blokptr2 = &(symbptr->bloktab[bloknum2]);
    simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_CTRBNBR]++;
    if(blokptr2->frownum < simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_FCOLNUM])
        simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_FCOLNUM] = blokptr2->frownum;
    if(blokptr2->lrownum > simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_LCOLNUM])
        simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_LCOLNUM] = blokptr2->lrownum;
    if(blokptr->frownum < simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_FROWNUM])
        simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_FROWNUM] = blokptr->frownum;
    if(blokptr->lrownum > simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_LROWNUM])
        simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_LROWNUM] = blokptr->lrownum;

    assert(simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_LCOLNUM] - simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_FCOLNUM]+1 > 0);
    assert(simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_LROWNUM] - simuctrl->ftgttab[ftgtnum].ftgt.infotab[FTGT_FROWNUM]+1 > 0);

}

static inline void
taskExec_COMP1D(pastix_int_t tasknum, SymbolMatrix *symbptr, SimuCtrl *simuctrl, BlendCtrl *ctrl, const Dof * dofptr)
{
    pastix_int_t          i, j;
    pastix_int_t          cblknum;
    pastix_int_t          facebloknum;
    pastix_int_t          facecblknum;
    pastix_int_t          local;
    pastix_int_t          facetasknum;
    pastix_int_t          ftgtnum;
    pastix_int_t          procnum;
    /*pastix_int_t          clustnum;*/
    SimuProc     *proc;
    CostMatrix   *costmtx;

    cblknum = simuctrl->tasktab[tasknum].cblknum; /* in case of COMP1D bloknum in a SimuTask struct means cblknum */
    procnum = simuctrl->ownetab[cblknum];
    proc    = &(simuctrl->proctab[procnum]);
    costmtx = ctrl->costmtx;

#ifdef DEBUG_BLEND
    if (procnum < ctrl->candtab[cblknum].fcandnum || procnum > ctrl->candtab[cblknum].lcandnum)
        fprintf(stderr, "procnum : %ld, fcandnum : %ld, lcandnum : %ld\n",
                (long)procnum, (long)ctrl->candtab[cblknum].fcandnum, (long)ctrl->candtab[cblknum].lcandnum);
    assert(procnum >= ctrl->candtab[cblknum].fcandnum && procnum <= ctrl->candtab[cblknum].lcandnum);
#endif

    /** Add time for factorizatoin of the diag blok and repercution on the off diag bloks **/
    i = symbptr->cblktab[cblknum].bloknum;
    timerAdd(&(proc->timer), costmtx->bloktab[i].contrib);

    for(i++;i<symbptr->cblktab[cblknum+1].bloknum;i++)
    {
        /** Add time for compute of the contrib due to this odb **/
        /** OIMBE: pour l'instant je considere que les contribs sont calculees
         en un bloc **/
        timerAdd(&(proc->timer), costmtx->bloktab[i].contrib);

        facecblknum = symbptr->bloktab[i].cblknum;

        if(ctrl->candtab[facecblknum].fccandnum ==  ctrl->candtab[facecblknum].lccandnum)
            local = 1; /** Facing task is COMP1D and is on the local cluster subtree **/
        else
            local = 0;

        if(!local || ctrl->candtab[facecblknum].cblktype != CBLK_1D )
        {
            facebloknum = 0;

            for(j=i;j<symbptr->cblktab[cblknum+1].bloknum;j++)
            {
                /* OIMBE trop couteux !! */
                facebloknum = getFaceBlockE2(facebloknum, i, j, symbptr, ctrl->ricar);

#ifdef DEBUG_M
                if(ctrl->ricar == 0)
                    assert(facebloknum >= 0);
#endif
                if(facebloknum >= 0)
                {
                    if(ctrl->candtab[facecblknum].cblktype == CBLK_1D)
                        simuctrl->cblktab[facecblknum].ctrbcnt--;
                    else
                        simuctrl->bloktab[facebloknum].ctrbcnt--;

                    if(!local)
                    {
                        ftgtnum = CLUST2INDEX(facebloknum, ctrl->core2clust[procnum]);
                        updateFtgtStruct(j, i, ftgtnum, symbptr, simuctrl, ctrl);

                        /** Update timer ready for receiver of the ftgt **/
                        ftgtnum = CLUST2INDEX(symbptr->cblktab[facecblknum].bloknum, ctrl->core2clust[procnum]);
                        timerSet(&(simuctrl->ftgttimetab[ftgtnum]) , MAX( timerVal(&(simuctrl->ftgttimetab[ftgtnum])) ,
                                                                          timerVal(&(proc->timer))));
                    }
                    else
                    {
                        /*** LOCAL task ***/
                        facetasknum = simuctrl->bloktab[symbptr->cblktab[facecblknum].bloknum].tasknum;
                        timerSet(&(simuctrl->tasktab[facetasknum].time), MAX(timerVal(&(proc->timer)),
                                                                             timerVal(&(simuctrl->tasktab[facetasknum].time))));
                    }

                    if( (ctrl->candtab[facecblknum].cblktype == CBLK_1D && simuctrl->cblktab[facecblknum].ctrbcnt == 0) )
                    {
                        facetasknum = simuctrl->bloktab[facebloknum].tasknum;

                        if(!local)
                            computeTaskReceiveTime(facetasknum, symbptr, simuctrl, ctrl, dofptr);

                        assert( facecblknum == simuctrl->tasktab[facetasknum].cblknum );
                        simu_putInAllReadyQueues( ctrl, simuctrl, facetasknum );
                        assert(facetasknum < simuctrl->tasknbr);
                    }
                }
            }
        }
        else
        {
            /** The facing task is local COMP_1D**/
            if(ctrl->ricar == 1)
            {
                facebloknum = 0;
                for(j=i;j<symbptr->cblktab[cblknum+1].bloknum;j++)
                {
                    /* OIMBE trop couteux ON PEUT FAIRE MIEUX EN PARCOURANT EN DESCENDANT!! */
                    facebloknum = getFaceBlockE2(facebloknum, i, j, symbptr, ctrl->ricar);
                    if(facebloknum>=0)
                        simuctrl->cblktab[facecblknum].ctrbcnt--;
                }
            }
            else
                simuctrl->cblktab[facecblknum].ctrbcnt -= symbptr->cblktab[cblknum+1].bloknum - i;            /** A REVOIR POUR NAPA **/

            assert(simuctrl->cblktab[facecblknum].ctrbcnt >= 0);

            if(simuctrl->cblktab[facecblknum].ctrbcnt == 0)
            {
                facetasknum = simuctrl->bloktab[symbptr->cblktab[facecblknum].bloknum].tasknum;
                /** Update timer of the task (owned by the diag block**/
                assert(ctrl->candtab[facecblknum].fccandnum == ctrl->candtab[facecblknum].lccandnum);

                /* NB: The facing cblk is local !! */
                timerSet(&(simuctrl->tasktab[facetasknum].time), MAX(timerVal(&(proc->timer)),
                                                                     timerVal(&(simuctrl->tasktab[facetasknum].time))) );
#ifdef DEBUG_BLEND
                assert(simuctrl->tasktab[facetasknum].taskid == COMP_1D);
                if(ctrl->ricar == 1)
                    assert(ctrl->candtab[facecblknum].lccandnum == ctrl->candtab[facecblknum].fccandnum);
                if(ctrl->core2clust[procnum] != ctrl->candtab[facecblknum].fccandnum)
                {
                    fprintf(stderr, "clustnum %ld  face proc cand %ld \n", (long)ctrl->core2clust[procnum], (long)ctrl->candtab[facecblknum].fccandnum);
                    fprintf(stderr, "%ld candidat [%ld %ld] => %ld candidat [%ld %ld ]\n", (long)cblknum, (long)ctrl->candtab[cblknum].fccandnum, (long)ctrl->candtab[cblknum].lccandnum,
                            (long)facecblknum, (long)ctrl->candtab[facecblknum].fccandnum, (long)ctrl->candtab[facecblknum].lccandnum);
                }
                assert(ctrl->core2clust[procnum] == ctrl->candtab[facecblknum].fccandnum);
                assert(ctrl->core2clust[procnum] == ctrl->candtab[facecblknum].lccandnum);
                assert(facetasknum<simuctrl->tasknbr);
#endif

                /** Put the task in the ready heap of its local candidat processor **/
                assert( facecblknum == simuctrl->tasktab[facetasknum].cblknum );
                simu_putInAllReadyQueues( ctrl, simuctrl, facetasknum );
            }
        }
    }
}

pastix_int_t comp_int(const pastix_int_t * a, const pastix_int_t * b)
{
    if(a[0]>b[0])
        return 1;
    if(a[0]<b[0])
        return -1;

    /*a == b*/
    return 0;
}

/* OIMBE utiliser un tas pour getNextProc --> cout log(P) au lieu de P */
pastix_int_t getNextProc(SimuProc *proctab, pastix_int_t procnbr)
{
    double min;
    pastix_int_t pr;
    pastix_int_t procnum = -1;

    min = (double)INTVALMAX;
    for(pr=0;pr<procnbr;pr++)
        if((timerVal(&(proctab[pr].timer)) < min) && ( (queueSize(proctab[pr].taskheap)>0) || (queueSize(proctab[pr].taskheap2)>0) ))
        {
            min = timerVal(&(proctab[pr].timer));
            procnum = pr;
        }
    return procnum;
}

pastix_int_t getTaskUnmapped(Queue *q1, Queue *q2, SimuCtrl *simuctrl)
{
    pastix_int_t next = -1;
    while(queueSize(q2)>0)
    {
        next = queueGet(q2);
        if(simuctrl->blprtab[simuctrl->tasktab[next].bloknum]<0)
        {
            assert(simuctrl->ownetab[simuctrl->tasktab[next].cblknum]<0);
            goto end;
        }
    }
    while(queueSize(q1)>0)
    {
        next = queueGet(q1);
        if(simuctrl->blprtab[simuctrl->tasktab[next].bloknum]<0)
        {
            assert(simuctrl->ownetab[simuctrl->tasktab[next].cblknum]<0);
            goto end;
        }
    }

    /** no unmapped task found **/
    return -1;

  end:
    return next;
}



static inline void
queueReorder(const BlendCtrl *ctrl, SimuCtrl *simuctrl, pastix_int_t t )
{
    pastix_int_t tasknum;
    pastix_int_t cblknum;
    pastix_int_t procnum;

    {
        procnum = simuctrl->blprtab[simuctrl->tasktab[t].bloknum];
        while(queueSize(simuctrl->proctab[procnum].taskheap)>0)
        {
            tasknum = queueRead(simuctrl->proctab[procnum].taskheap);
            cblknum = simuctrl->tasktab[tasknum].cblknum;

            if(!compTimer(TIMER(procnum), &(simuctrl->ftgttimetab[CLUST2INDEX(simuctrl->tasktab[tasknum].bloknum, ctrl->core2clust[procnum])])))
            {
                tasknum = queueGet(simuctrl->proctab[procnum].taskheap);
                cblknum = simuctrl->tasktab[tasknum].cblknum;
                queueAdd2(simuctrl->proctab[procnum].taskheap2, tasknum,(double)ctrl->candtab[cblknum].treelevel, simuctrl->tasktab[tasknum].bloknum );
            }
            else
                break;
        }
    }
}


void simuRun( SymbolMatrix *symbptr,
              SimuCtrl *simuctrl,
              BlendCtrl *ctrl,
              const Dof * dofptr )
{

    pastix_int_t             i, j, b;
    pastix_int_t             cblknum, bloknum;
    /*pastix_int_t             c;*/
    pastix_int_t             pr;


    for(i=0;i<symbptr->cblknbr;i++)
        simuctrl->cblktab[i].ctrbcnt = ctrl->egraph->verttab[i].innbr;

    /* OIMBE attention les ctrbcnt des cblk en COMP1D sont recalculee dans computeBlockCtrbNbr */
    /** Compute number of contributions for blocks **/
    computeBlockCtrbNbr(ctrl, symbptr, simuctrl);

    /*
     * All ready tasks are put in the task heaps of their respective candidates
     */
    for(i=0;i<symbptr->cblknbr;i++)
    {
        pastix_int_t tasknum;
        if(simuctrl->cblktab[i].ctrbcnt == 0)
        {
            tasknum = simuctrl->bloktab[symbptr->cblktab[i].bloknum].tasknum;
            assert(ctrl->candtab[i].treelevel < 0);

            if(ctrl->costlevel)
                assert(ctrl->candtab[i].costlevel < 0);
            assert( simuctrl->tasktab[tasknum].taskid == COMP_1D );

            assert(simuctrl->tasktab[tasknum].cblknum == i);
            assert(ctrl->candtab[i].cblktype == CBLK_1D);

            simu_putInAllReadyQueues( ctrl, simuctrl, tasknum );
        }
    }

    /*
     * Run simulation and mapp the task onto a single candidate
     */
    while(1)
    {
        /* Get the next earlier task index and the processor on which it is mapped */
        i = getNextTaskNextProc(simuctrl, ctrl, &pr);

        /* No more tasks */
        if( i == -1 )
            break;

        bloknum = simuctrl->tasktab[i].bloknum;
        cblknum = simuctrl->tasktab[i].cblknum;

        /*
         * Compute the time at which each proc cand will have added its ftgt and
         * received block target if the task is mapped on
         */
        assert( simuctrl->ownetab[cblknum]<0 );

        /* Set processor owners */
        simuctrl->ownetab[cblknum] = pr;
        for(j=symbptr->cblktab[cblknum].bloknum;
            j<symbptr->cblktab[cblknum+1].bloknum;j++)
        {
            simuctrl->blprtab[j] = pr;
        }

        simuctrl->tasktab[i].prionum = simuctrl->clustab[ctrl->core2clust[pr]].prionum;
        simuctrl->clustab[ctrl->core2clust[pr]].prionum++;

        /* Ajout de la tache a la file du proc pour version standard */
        extendint_Add(simuctrl->proctab[pr].tasktab, i);

        /* Sauvegarde du processus MPI devant executer la tache pour version MARCEL */
        assert(simuctrl->tasktab[i].cblknum < symbptr->cblknbr);
        ctrl->candtab[simuctrl->tasktab[i].cblknum].cluster = ctrl->core2clust[pr];

        /*-------------------------------------------------------------/
         /   UPDATE TIMER OF THE PROC ON WHICH IS MAPPED THE TASK       /
         /   TO THE DATE THE PROC WILL BEGIN THE TASK INNER COMPUTATION /
         /-------------------------------------------------------------*/
        if( ctrl->candtab[simuctrl->tasktab[i].cblknum].fccandnum == ctrl->candtab[simuctrl->tasktab[i].cblknum].lccandnum )
            /** Time do not depend on the reception of a ftgt **/
            timerSet(TIMER(pr), MAX(timerVal(TIMER(pr)), timerVal(&(simuctrl->tasktab[i].time))));
        else
            /** Time depends on the reception of a ftgt **/
            timerSet(TIMER(pr),
                     MAX(timerVal(TIMER(pr)),
                         timerVal(&(simuctrl->ftgttimetab[CLUST2INDEX(bloknum, ctrl->core2clust[pr])]))));

        /*------------------------------------------------------------------------/
         /  Fill some fanintarget info (task of type E2 does not have any ftgt)    /
         /------------------------------------------------------------------------*/
        if(simuctrl->bloktab[bloknum].ftgtnum< simuctrl->bloktab[bloknum+1].ftgtnum)
        {
            /** Task COMP_1D with several cand cluster **/
            for(b=bloknum;b<symbptr->cblktab[cblknum+1].bloknum;b++)
            {
                for(j=simuctrl->bloktab[b].ftgtnum;j<simuctrl->bloktab[b+1].ftgtnum;j++)
                {
                    if((simuctrl->ftgttab[j].ftgt.infotab[FTGT_CTRBNBR] >0)
                       && (j != CLUST2INDEX(b, ctrl->core2clust[pr])))
                    {
                        simuctrl->ftgttab[j].clustnum = INDEX2CLUST(j, b);
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_PRIONUM] = simuctrl->tasktab[i].prionum;
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_PROCDST] = pr;
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_BLOKDST] = b;
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_TASKDST] = simuctrl->bloktab[bloknum].tasknum;
#ifdef OOC_FTGT
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_GCBKDST] = simuctrl->tasktab[simuctrl->bloktab[bloknum].tasknum].cblknum;
#endif
                        extendint_Add(&(simuctrl->clustab[INDEX2CLUST(j,b)].ftgtsend[ctrl->core2clust[pr]]), j);

                        simuctrl->tasktab[simuctrl->bloktab[bloknum].tasknum].ftgtcnt++;

                        if (ctrl->core2clust[pr] == ctrl->clustnum)
                            simuctrl->ftgtcnt++;
                    }
                }

            }
            simuctrl->ftgtprio++;

        }
        else
            assert(ctrl->candtab[cblknum].fccandnum == ctrl->candtab[cblknum].lccandnum);

        /* Simule the computing of the task */
        taskExec_COMP1D(i, symbptr, simuctrl, ctrl, dofptr);

        queueReorder(ctrl, simuctrl, i);
    }

    double maxtime = 0;
    for(pr=0;pr<ctrl->total_nbcores;pr++)
    {
        if(timerVal(TIMER(pr)) > maxtime)
            maxtime = timerVal(TIMER(pr));
    }
    set_dparm(ctrl->dparm, DPARM_PRED_FACT_TIME, maxtime);

#ifdef DEBUG_BLEND

    for(i=0;i<simuctrl->cblknbr;i++)
        if(ctrl->candtab[i].cblktype == CBLK_1D)
            if(simuctrl->ownetab[i] < 0) /** Check valid for 1D distribution only **/
                fprintf(stderr, "CBLK %ld has no processor \n", (long)i);

    for(i=0;i<symbptr->bloknbr;i++)
        if(!(simuctrl->blprtab[i]>=0))
        {
            fprintf(stderr, "BLOCK %ld has no processor \n", (long)i);
            fprintf(stdout, "blprtab [ %ld ] = %ld type %ld \n", (long)i,
                    (long)simuctrl->blprtab[i],
                    (long)simuctrl->tasktab[simuctrl->bloktab[i].tasknum].taskid);
            EXIT(MOD_BLEND,INTERNAL_ERR);
        }
    for(i=0;i<symbptr->bloknbr;i++)
        assert(simuctrl->blprtab[i]>=0);
#endif
}
