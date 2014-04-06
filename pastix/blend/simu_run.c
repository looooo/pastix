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
#include "updown.h"
#include "solver.h"
#include "costfunc.h"

static inline void
simu_computeBlockCtrbNbr(const SymbolMatrix *symbptr,
                         SimuCtrl     *simuctrl,
                         pastix_int_t  ricar )
{
    pastix_int_t i, j, k;
    pastix_int_t facebloknum, firstbloknum;

    /* Compute the number of contribution per block to each block */
    /* Might be optimized if we computed the input graph before */
    {
        SymbolCblk *curcblk;

        curcblk = symbptr->cblktab;
        for(i=0; i<symbptr->cblknbr; i++, curcblk++)
        {
            pastix_int_t fbloknum = curcblk[0].bloknum + 1;
            pastix_int_t lbloknum = curcblk[1].bloknum;

            /* 1D cblk computed */
            for(j=fbloknum; j<lbloknum; j++)
            {
                firstbloknum = 0;

                /* Add contribution due to E2 */
                for(k=j; k<lbloknum; k++)
                {
                    facebloknum = symbolGetFacingBloknum( symbptr, j, k, firstbloknum, ricar );
                    if(facebloknum >= 0) {
                        simuctrl->bloktab[facebloknum].ctrbcnt++;
                        firstbloknum = facebloknum;
                    }
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
    SimuProc  *sproc;
    double ready_date = 0.0;
    pastix_int_t procnum;
    pastix_int_t bloknum = task->bloknum;
    pastix_int_t treelevel = cblkcand->treelevel;

    /* Get the ready date of the task on the processor passed in parameter */
    if( cblkcand->fccandnum == cblkcand->lccandnum )
    {
        ready_date = timerVal( &(task->time) );
        sproc = &(simuctrl->proctab[cblkcand->fcandnum]);

        for(procnum =  cblkcand->fcandnum;
            procnum <= cblkcand->lcandnum; procnum++, sproc++)
        {
            if(ready_date > timerVal(TIMER(procnum)))
                pqueuePush2( sproc->futuretask, tasknum, ready_date, treelevel);
            else
                pqueuePush2( sproc->readytask, tasknum, treelevel, bloknum);
        }
    }
    else
    {
        sproc = &(simuctrl->proctab[cblkcand->fcandnum]);

        for(procnum =  cblkcand->fcandnum;
            procnum <= cblkcand->lcandnum; procnum++, sproc++)
        {
            ready_date = timerVal( simuctrl->ftgttimetab + CLUST2INDEX(bloknum, ctrl->core2clust[procnum]) );

            if(ready_date > timerVal(TIMER(procnum)))
                pqueuePush2( sproc->futuretask, tasknum, ready_date, treelevel);
            else
                pqueuePush2( sproc->readytask, tasknum, treelevel, bloknum);
        }
    }
}

static inline pastix_int_t
simu_getNextTaskNextProc( const BlendCtrl *ctrl,
                          SimuCtrl        *simuctrl,
                          pastix_int_t    *procnumptr )
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
        while(pqueueSize(simuctrl->proctab[p].readytask)>0)
        {
            tasknum = pqueueRead(simuctrl->proctab[p].readytask);
            if( simuctrl->bloktab[simuctrl->tasktab[tasknum].bloknum].ownerclust >= 0 )
            {
                /** This task have to be remove from the heap (already mapped) **/
                pqueuePop(simuctrl->proctab[p].readytask);
                tasknum = -1;
            }
            else
                break;
        }
        /** We found no task which ready date is < proc timer so we search one that minimizes ready date - proc-timer **/
        if(tasknum == -1)
        {
            while(pqueueSize(simuctrl->proctab[p].futuretask)>0)
            {
                tasknum = pqueueRead(simuctrl->proctab[p].futuretask);
                if( simuctrl->bloktab[simuctrl->tasktab[tasknum].bloknum].ownerclust >= 0 )
                {
                    /** This task have to be remove from the heap (already mapped) **/
                    pqueuePop(simuctrl->proctab[p].futuretask);
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
        if(pqueueSize(simuctrl->proctab[procnum].readytask)>0)
        {
            assert(earlytask == pqueuePop(simuctrl->proctab[procnum].readytask));
        }
        else
        {
            assert(earlytask == pqueuePop(simuctrl->proctab[procnum].futuretask));
        }
    }
    *procnumptr = procnum;
    return earlytask;
}

static inline void
simu_computeTaskReceiveTime( const BlendCtrl    *ctrl,
                             const SymbolMatrix *symbptr,
                             const Dof          *dofptr,
                             SimuCtrl           *simuctrl,
                             pastix_int_t        tasknum )
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

    /* If the task is local, all sons sending contributions are local => no treatment */
    if( ctrl->candtab[cblknum].fccandnum == ctrl->candtab[cblknum].lccandnum )
        return;

    /*------------------------------------------------------------------------------------------------/
     / Compute the cblk on proc timer that is time the cblk would have received                        /
     / all its contributions if it was mapped on a given cand processor                                /
     / These times INCLUDE add time for fan in target !!                                               /
     /------------------------------------------------------------------------------------------------*/

    /** Compute receive time (time at which a non-local processor should received the target **/
    /* find the latest ftgt receive time and the second latest*/
    for(i=simuctrl->bloktab[bloknum].ftgtnum; i<simuctrl->bloktab[bloknum+1].ftgtnum; i++)
    {
        /* Source of this ftgt */
        clustdst = INDEX2CLUST(i, bloknum);

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
                        costFtgtSend( ctrl, dofptr,
                                      &(simuctrl->ftgttab[CLUST2INDEX(j, clustdst)].ftgt),
                                      clustdst, ctrl->candtab[cblknum].lccandnum-ctrl->candtab[cblknum].fccandnum+1 );
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
            timerSetMax( &(simuctrl->ftgttimetab[i]), sftgttime );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_simulation
 *
 * simu_updateFtgt - Update the Fan In target structure by incrementing the
 * contribution counter and integrating to the ftgt area the new contribution.
 *
 *******************************************************************************
 *
 * @param[in] symbptr
 *          The pointer to the symbolic matrix structure.
 *
 * @param[in,out] simuctrl
 *          The pointer to the simulation structure. On exit, data regarding the
 *          computational unit pr are updated.
 *
 * @param[in] ftgtnum
 *          Index of the fanin target to update.
 *
 * @param[in] bloknum
 *          Index of the first of diagonal block generating a contribution to
 *          the ftgtnum Fan In.
 *
 * @param[in] fbloknum
 *          Index of the facing blok of bloknum that will receive the final
 *          contribution.
 *
 *******************************************************************************/
static inline void
simu_updateFtgt( const SymbolMatrix *symbptr,
                 SimuCtrl           *simuctrl,
                 pastix_int_t        ftgtnum,
                 pastix_int_t        bloknum,
                 pastix_int_t        fbloknum )
{
    FanInTarget  *ftgt     = &(simuctrl->ftgttab[ftgtnum].ftgt);
    pastix_int_t *infotab  = ftgt->infotab;
    SymbolBlok   *blokptr  = (symbptr->bloktab) + bloknum;
    SymbolBlok   *fblokptr = (symbptr->bloktab) + fbloknum;

    infotab[FTGT_CTRBNBR]++;

    /* Update ftgt dimensions to the maximum area covering all contributions */
    if( blokptr->frownum < infotab[FTGT_FCOLNUM] )
        infotab[FTGT_FCOLNUM] = blokptr->frownum;

    if( blokptr->lrownum > infotab[FTGT_LCOLNUM] )
        infotab[FTGT_LCOLNUM] = blokptr->lrownum;

    if( fblokptr->frownum < infotab[FTGT_FROWNUM] )
        infotab[FTGT_FROWNUM] = fblokptr->frownum;

    if( fblokptr->lrownum > infotab[FTGT_LROWNUM] )
        infotab[FTGT_LROWNUM] = fblokptr->lrownum;

    assert( (infotab[FTGT_LCOLNUM] - infotab[FTGT_FCOLNUM] + 1) > 0 );
    assert( (infotab[FTGT_LROWNUM] - infotab[FTGT_FROWNUM] + 1) > 0 );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_simulation
 *
 * simu_computetTask - This routine simulates the task execution by updating the
 * timers of selected processor, and cblk, as well as cblk targetted by the
 * updates.
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The pointer to the global blend control structure.
 *
 * @param[in] symbptr
 *          The pointer to the symbolic matrix structure.
 *
 * @param[in,out] simuctrl
 *          The pointer to the simulation structure. On exit, data regarding the
 *          computational unit pr are updated.
 *
 * @param[in] tasknum
 *          The task index of the one, we want to simulate the execution.
 *
 *******************************************************************************
 *
 * Remark: In this function, we use the standard [f|l]blocknum for first and
 * last bloknum, and facingcblk, facingblok for the facing block and column
 * block.
 *
 *******************************************************************************/
static inline void
simu_computeTask( const BlendCtrl    *ctrl,
                  const SymbolMatrix *symbptr,
                  const Dof          *dofptr,
                  SimuCtrl           *simuctrl,
                  pastix_int_t        tasknum )
{
    pastix_int_t  i, j;
    pastix_int_t  cblknum;
    pastix_int_t  fbloknum;
    pastix_int_t  lbloknum;
    pastix_int_t  firstfacingblok;
    pastix_int_t  facingblok;
    pastix_int_t  facingcblk;
    pastix_int_t  local;
    pastix_int_t  ftgtnum;
    pastix_int_t  procnum;
    pastix_int_t  clustnum;
    SimuProc     *sproc;
    CostMatrix   *costmtx;

    cblknum  = simuctrl->tasktab[tasknum].cblknum;
    procnum  = simuctrl->ownetab[cblknum];
    clustnum = ctrl->core2clust[procnum];
    sproc    = &(simuctrl->proctab[procnum]);
    costmtx  = ctrl->costmtx;

    fbloknum = symbptr->cblktab[cblknum  ].bloknum;
    lbloknum = symbptr->cblktab[cblknum+1].bloknum;

    assert( (procnum >= ctrl->candtab[cblknum].fcandnum) &&
            (procnum <= ctrl->candtab[cblknum].lcandnum) );

    /* Add factorization time to the diagonal blok */
    timerAdd(&(sproc->timer), costmtx->bloktab[fbloknum].contrib);

    for(i=fbloknum+1; i<lbloknum; i++)
    {
        /* Add trsm time of this off-diagonal block */
        timerAdd(&(sproc->timer), costmtx->bloktab[i].contrib);

        facingcblk = symbptr->bloktab[i].cblknum;

        /*
         * If only one candidate cluster, we can consider the facingcblk as
         * local because it is an ancestor of the current cblk in the
         * elimination tree.
         */
        local = ( ctrl->candtab[facingcblk].fccandnum == ctrl->candtab[facingcblk].lccandnum ) ? 1 : 0;

        firstfacingblok = symbptr->cblktab[facingcblk].bloknum;

        for(j=i; j<lbloknum; j++)
        {
            /* TODO: symbolGetFacingBloknum is too expensive !! */
            facingblok = symbolGetFacingBloknum(symbptr, i, j, firstfacingblok, ctrl->ricar);

            /* If the couple (i, j) generates a contribution, applies it */
            if( facingblok >= 0 ) {
                pastix_int_t facingdiagblok;
                pastix_int_t facingtask;

                /* Decrease contributions on block and column block */
                simuctrl->cblktab[facingcblk].ctrbcnt--;
                simuctrl->bloktab[facingblok].ctrbcnt--;

                /* Checks */
                assert(simuctrl->cblktab[facingcblk].ctrbcnt >= 0);
                assert(simuctrl->bloktab[facingblok].ctrbcnt >= 0);

                /* Update to start next search from the last facing block */
                firstfacingblok = facingblok;

                facingdiagblok = symbptr->cblktab[facingcblk].bloknum;
                facingtask     = simuctrl->bloktab[facingdiagblok].tasknum;

                assert( facingcblk == simuctrl->tasktab[facingtask].cblknum );
                assert( facingtask < simuctrl->tasknbr );

                if(!local)
                {
                    ftgtnum = CLUST2INDEX(facingblok, clustnum);
                    simu_updateFtgt( symbptr, simuctrl, ftgtnum, i, j );

                    /* Update timer ready for receiver of the ftgt */
                    ftgtnum = CLUST2INDEX( facingdiagblok, clustnum );
                    timerSetMax( &(simuctrl->ftgttimetab[ftgtnum]),
                                 timerVal(&(sproc->timer)) );

                }
                else {

                    /* Update timer of the task (associated to the diagonal block) */
                    timerSetMax( &(simuctrl->tasktab[facingtask].time),
                                 timerVal(&(sproc->timer)) );
                }

                if( simuctrl->cblktab[facingcblk].ctrbcnt == 0 ) {
                    if (!local)
                        simu_computeTaskReceiveTime(ctrl, symbptr, dofptr, simuctrl, facingtask );

                    /* Put the task in the ready heap of its local candidat processor */
                    simu_putInAllReadyQueues( ctrl, simuctrl, facingtask );
                }
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_simulation
 *
 * simu_pushToReadyHeap - This routine pushes all future task from the future
 * task heap to the ready one, if the time at which the task will be ready is
 * already passed by the computation unit.
 *
 *******************************************************************************
 *
 * @param[in] ctrl
 *          The pointer to the global blend control structure.
 *
 * @param[in,out] simuctrl
 *          The pointer to the simulation structure. On exit, data regarding the
 *          computational unit pr are updated.
 *
 * @param[in] pr
 *          The computational unit index for which, the data need to be transfer
 *          from the future task heap to ready task heap if the computatuional
 *          unit timer is more advanced than the ready time of the tasks.
 *
 *******************************************************************************/
static inline void
simu_pushToReadyHeap(const BlendCtrl *ctrl,
                     SimuCtrl        *simuctrl,
                     pastix_int_t     procnum )
{
    SimuProc    *sproc;
    pastix_int_t tasknum;
    pastix_int_t cblknum;
    pastix_int_t clustnum;

    clustnum = ctrl->core2clust[procnum];
    sproc    = &(simuctrl->proctab[procnum]);

    /*
     * Move each task from future task heap to ready heap if the timer is
     * further in the future than the ready time
     */
    while( pqueueSize(sproc->futuretask) > 0 )
    {
        tasknum = pqueueRead(sproc->futuretask);

        if(! timerComp( &(sproc->timer),
                        &(simuctrl->ftgttimetab[CLUST2INDEX(simuctrl->tasktab[tasknum].bloknum, clustnum )])) )
        {
            tasknum = pqueuePop(sproc->futuretask);
            cblknum = simuctrl->tasktab[tasknum].cblknum;

            pqueuePush2(sproc->readytask, tasknum,
                        ctrl->candtab[cblknum].treelevel,
                        simuctrl->tasktab[tasknum].bloknum );
        }
        else
            break;
    }
}


void
simuRun( SimuCtrl           *simuctrl,
         const BlendCtrl    *ctrl,
         const SymbolMatrix *symbptr,
         const Dof          *dofptr )
{

    pastix_int_t             i, j, b;
    pastix_int_t             cblknum, bloknum;
    /*pastix_int_t             c;*/
    pastix_int_t             pr;

    /* Compute number of contributions per blocks, cblks, tasks */
    simu_computeBlockCtrbNbr( symbptr, simuctrl, ctrl->ricar );

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
        SimuTask    *task;
        pastix_int_t clustnum;

        /* Get the next earlier task index and the processor on which it is mapped */
        i = simu_getNextTaskNextProc(ctrl, simuctrl, &pr);

        /* No more tasks */
        if( i == -1 )
            break;

        task    = &(simuctrl->tasktab[i]);
        bloknum = task->bloknum;
        cblknum = task->cblknum;
        clustnum= ctrl->core2clust[pr];

        assert(cblknum < symbptr->cblknbr);
        assert(bloknum < symbptr->bloknbr);

        /* Make sure the cblk is not already attibuted to someone and give it to the selected proc */
        assert( simuctrl->ownetab[cblknum] < 0 );
        simuctrl->ownetab[cblknum] = pr;
        for(j=symbptr->cblktab[cblknum].bloknum;
            j<symbptr->cblktab[cblknum+1].bloknum;j++)
        {
            simuctrl->bloktab[j].ownerclust = clustnum;
        }
        task->prionum = simuctrl->clustab[clustnum].prionum;
        simuctrl->clustab[clustnum].prionum++;

        /* Add task to the selected processor list */
        extendint_Add(simuctrl->proctab[pr].tasktab, i);

        /* Backup which cluster will get the data for the second run of proportionnal mapping */
        ctrl->candtab[cblknum].cluster = clustnum;

        /*
         * Compute the time at which each proc cand will have added its ftgt and
         * received block target if the task is mapped on
         */
        if( ctrl->candtab[cblknum].fccandnum == ctrl->candtab[cblknum].lccandnum ) {
            /*
             * All contributions come from the same node
             * Time do not depend on the reception of a ftgt
             */
            timerSetMax( TIMER(pr), timerVal(&(task->time)));
        }
        else {
            /*
             * Contributions might come from different nodes
             * Time depends on the reception of a ftgt
             */
            timerSetMax( TIMER(pr),
                         timerVal(&(simuctrl->ftgttimetab[CLUST2INDEX(bloknum, clustnum)])));
        }

        /*
         * Fill some fanintarget info (task of type E2 does not have any ftgt)
         */
        if(simuctrl->bloktab[bloknum].ftgtnum < simuctrl->bloktab[bloknum+1].ftgtnum)
        {
            /* Task COMP_1D with several cand cluster */
            for(b=bloknum; b<symbptr->cblktab[cblknum+1].bloknum; b++)
            {
                for(j=simuctrl->bloktab[b].ftgtnum; j<simuctrl->bloktab[b+1].ftgtnum; j++)
                {
                    if( (simuctrl->ftgttab[j].ftgt.infotab[FTGT_CTRBNBR] >0)
                        && (j != CLUST2INDEX(b, clustnum)))
                    {
                        simuctrl->ftgttab[j].clustnum = INDEX2CLUST(j, b);
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_PRIONUM] = task->prionum;
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_PROCDST] = pr;
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_BLOKDST] = b;
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_TASKDST] = simuctrl->bloktab[bloknum].tasknum;
#ifdef OOC_FTGT
                        simuctrl->ftgttab[j].ftgt.infotab[FTGT_GCBKDST] = simuctrl->tasktab[simuctrl->bloktab[bloknum].tasknum].cblknum;
#endif
                        extendint_Add(&(simuctrl->clustab[INDEX2CLUST(j,b)].ftgtsend[clustnum]), j);

                        simuctrl->tasktab[simuctrl->bloktab[bloknum].tasknum].ftgtcnt++;

                        if (clustnum == ctrl->clustnum)
                            simuctrl->ftgtcnt++;
                    }
                }
            }
            simuctrl->ftgtprio++;
        }
        else {
            assert(ctrl->candtab[cblknum].fccandnum == ctrl->candtab[cblknum].lccandnum);
        }

        /* Simule the computing of the task */
        simu_computeTask( ctrl, symbptr, dofptr, simuctrl, i );

        simu_pushToReadyHeap(ctrl, simuctrl, pr);
    }

    /* Compute maximum time */
    {
        double maxtime = 0;
        for(pr=0; pr<ctrl->total_nbcores; pr++)
        {
            if(timerVal(TIMER(pr)) > maxtime)
                maxtime = timerVal(TIMER(pr));
        }
        set_dparm(ctrl->dparm, DPARM_PRED_FACT_TIME, maxtime);
    }

#ifdef DEBUG_BLEND
    for(i=0;i<simuctrl->cblknbr;i++)
        if(ctrl->candtab[i].cblktype == CBLK_1D)
            if(simuctrl->ownetab[i] < 0) /** Check valid for 1D distribution only **/
                fprintf(stderr, "CBLK %ld has no processor \n", (long)i);

    for(i=0;i<symbptr->bloknbr;i++)
        if(!(simuctrl->bloktab[i].ownerclust>=0))
        {
            fprintf(stderr, "BLOCK %ld has no processor \n", (long)i);
            fprintf(stdout, "ownerclust [ %ld ] = %ld type %ld \n", (long)i,
                    (long)simuctrl->bloktab[i].ownerclust,
                    (long)simuctrl->tasktab[simuctrl->bloktab[i].tasknum].taskid);
            EXIT(MOD_BLEND,INTERNAL_ERR);
        }
    for(i=0;i<symbptr->bloknbr;i++)
        assert(simuctrl->bloktab[i].ownerclust>=0);
#endif
}
