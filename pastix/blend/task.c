#include <stdio.h>
#include <math.h>

#include "common.h"
#include "symbol.h"
#include "ftgt.h"
#include "queue.h"
#include "extendVector.h"
#include "elimin.h"
#include "cost.h"
#include "cand.h"
#include "dof.h"
#include "bulles.h"
#include "blendctrl.h"
#include "simu.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "costfunc.h"
/* #include "extrastruct.h" */
/* #include "param_comm.h" */
/* #include "perf.h" */
/* #include "assert.h" */
#include "task.h"

void taskBuild(SimuCtrl *simuctrl, SymbolMatrix *symbptr, Cand *candtab,
               const Dof * dofptr, BlendCtrl *ctrl)
{
    pastix_int_t i, j;
    pastix_int_t tasknbr = 0;
    SimuTask *task = NULL;

    /* Count number of task */
    for(i=0;i<symbptr->cblknbr;i++)
    {
        switch(candtab[i].cblktype)
        {
        case CBLK_1D:
        case CBLK_SPLIT:
            tasknbr++;
            break;
        default:
            fprintf(stderr, "Task No %ld has wrong type \n", (long)i);
            EXIT(MOD_BLEND,INTERNAL_ERR);
        }
    }

    simuctrl->tasknbr = tasknbr;

    MALLOC_INTERN(simuctrl->tasktab, tasknbr, SimuTask);
#ifdef DEBUG_BLEND
    ASSERT(simuctrl->tasktab != NULL,MOD_BLEND);
#endif
    task = simuctrl->tasktab;
    tasknbr = 0;
    for(i=0;i<symbptr->cblknbr;i++)
    {
        switch(candtab[i].cblktype)
        {
        case CBLK_1D:
        case CBLK_SPLIT:
            task->taskid   = COMP_1D;
            task->prionum  = -1;
            task->cblknum  = i;
            task->bloknum  = symbptr->cblktab[i].bloknum;
            task->bloknum2 = -1;
            task->ctrbcnt  = 0;
            task->ftgtcnt  = 0;
            task->facebloknum = -1;
            task->cost     = -1;
            timerSet(&(task->time), 0.0);
            task->mesglen  = 0.0;
            task->tasknext = -1;
            for(j = symbptr->cblktab[i].bloknum;
                j < symbptr->cblktab[i+1].bloknum; j++ )
            {
                simuctrl->bloktab[j].tasknum = tasknbr;
            }
            tasknbr++;
            task++;
            break;

        default:
            fprintf(stderr, "Task No %ld has wrong type \n", (long)i);
            EXIT(MOD_BLEND,INTERNAL_ERR);
        }
    }

#ifdef DEBUG_BLEND
    ASSERT(simuctrl->tasknbr == tasknbr,MOD_BLEND);
    for(i=0;i<simuctrl->tasknbr;i++)
        if(simuctrl->tasktab[i].tasknext != -1)
            if(simuctrl->tasktab[simuctrl->tasktab[i].tasknext].taskid != -1)
                ASSERT(simuctrl->tasktab[simuctrl->tasktab[i].tasknext].taskid == simuctrl->tasktab[i].taskid,MOD_BLEND);
#endif
    ;
}

double taskSendCost(SimuTask *taskptr, const pastix_int_t clustsrc, const pastix_int_t clustdst, BlendCtrl *ctrl)
{
    double startup, bandwidth;

    getCommunicationCosts( ctrl, clustsrc, clustdst,
                           ctrl->candtab[taskptr->cblknum].lccandnum -
                           ctrl->candtab[taskptr->cblknum].fccandnum + 1,
                           &startup, &bandwidth );

    assert( taskptr->taskid != COMP_1D );

    return (startup + bandwidth * taskptr->mesglen);
}
