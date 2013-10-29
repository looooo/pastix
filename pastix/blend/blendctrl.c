#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>

#include "common.h"
#include "elimin.h"
#include "cost.h"
#include "extendVector.h"
#include "cand.h"
#include "queue.h"
#include "bulles.h"
#include "blendctrl.h"
#include "perf.h"

void getCommunicationCosts( BlendCtrl   *ctrl,
                            pastix_int_t clustsrc,
                            pastix_int_t clustdst,
                            pastix_int_t sync_comm_nbr,
                            double      *startup,
                            double      *bandwidth)
{
    assert((clustsrc >= 0) && (clustsrc < ctrl->clustnbr));
    assert((clustdst >= 0) && (clustdst < ctrl->clustnbr));
    assert((sync_comm_nbr > 0) && (sync_comm_nbr <= ctrl->clustnbr));

    if(clustsrc == clustdst)
    {
        *startup   = 0.;
        *bandwidth = 0.;
        return;
    }

    /* Shared Memory */
    if( ctrl->clust2smp[clustsrc] == ctrl->clust2smp[clustdst] )
    {
        switch (sync_comm_nbr)
        {
        case 1:
        case 2:
            *startup   = SHARED_STARTUP_1;
            *bandwidth = SHARED_BANDWIDTH_1;
            return;
        case 3:
        case 4:
            *startup   = SHARED_STARTUP_2;
            *bandwidth = SHARED_BANDWIDTH_2;
            return;
        case 5:
        case 6:
        case 7:
        case 8:
            *startup   = SHARED_STARTUP_4;
            *bandwidth = SHARED_BANDWIDTH_4;
            return;
        default:
            *startup   = SHARED_STARTUP_8;
            *bandwidth = SHARED_BANDWIDTH_8;
            return;
        }
    }
    else
    {
        switch (sync_comm_nbr)
        {
        case 1:
        case 2:
            *startup   = CLUSTER_STARTUP_1;
            *bandwidth = CLUSTER_BANDWIDTH_1;
            return;
        case 3:
        case 4:
            *startup   = CLUSTER_STARTUP_2;
            *bandwidth = CLUSTER_BANDWIDTH_2;
            return;
        case 5:
        case 6:
        case 7:
        case 8:
            *startup   = CLUSTER_STARTUP_4;
            *bandwidth = CLUSTER_BANDWIDTH_4;
            return;
        default:
            *startup   = CLUSTER_STARTUP_8;
            *bandwidth = CLUSTER_BANDWIDTH_8;
            return;
        }
    }
}

int
blendCtrlInit(BlendCtrl    *ctrl,
              pastix_int_t  procnum,
              pastix_int_t  procnbr,
              pastix_int_t  local_coresnbr,
              pastix_int_t  local_thrdsnbr,
              pastix_int_t *iparm)
{
    pastix_int_t i;

    /* Check parameters */
    if( ctrl == NULL )
    {
        errorPrint("blendCtrlInit: Illegal ctrl parameter\n");
        return PASTIX_ERR_BADPARAMETER;
    }
    if( procnum < 0 )
    {
        errorPrint("blendCtrlInit: Illegal procnum parameter\n");
        return PASTIX_ERR_BADPARAMETER;
    }
    if( procnbr < 1 )
    {
        errorPrint("blendCtrlInit: Illegal procnbr parameter\n");
        return PASTIX_ERR_BADPARAMETER;
    }
    if( local_coresnbr < 1 )
    {
        errorPrint("blendCtrlInit: Illegal local_coresnbr parameter\n");
        return PASTIX_ERR_BADPARAMETER;
    }
    if( local_thrdsnbr < 1 )
    {
        errorPrint("blendCtrlInit: Illegal local_thrdsnbr parameter\n");
        return PASTIX_ERR_BADPARAMETER;
    }
    if( procnum >= procnbr )
    {
        errorPrint("blendCtrlInit: Incompatible values of procnum(%d) and procnbr (%d)\n",
                   (int) procnum, (int) procnbr);
        return PASTIX_ERR_BADPARAMETER;
    }
    if( ctrl == NULL )
    {
        errorPrint("blendCtrlInit: Illegal ctrl parameter\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Initialize options */
    ctrl->count_ops = 1;
#if !defined(PASTIX_DEBUG_BLEND)
    ctrl->debug     = 1;
#else
    ctrl->debug     = 0;
#endif
    ctrl->timer     = 1;
    ctrl->ooc       = 0;
    ctrl->ricar     = iparm[IPARM_INCOMPLETE];
    ctrl->leader           = 0;

    /* Proportional Mapping options */
    ctrl->allcand     = 0;
    ctrl->nocrossproc = 0;
    ctrl->costlevel   = 1;

    /* Spliting options */
    ctrl->blcolmin = 60;
    ctrl->blcolmax = 120;
    ctrl->blcolmin = iparm[IPARM_MIN_BLOCKSIZE];
    ctrl->blcolmax = iparm[IPARM_MAX_BLOCKSIZE];
    ctrl->abs      = iparm[IPARM_ABS];
    ctrl->updatecandtab = 0;
    if(ctrl->blcolmin > ctrl->blcolmax)
    {
        errorPrint("Parameter error : blocksize max < blocksize min (cf. iparm.txt).");
        assert(ctrl->blcolmin <= ctrl->blcolmax);
    }

    /* 2D options */
    ctrl->autolevel  = 0;
    ctrl->level2D    = 100000000;
    ctrl->level2D    = iparm[IPARM_DISTRIBUTION_LEVEL];
    ctrl->ratiolimit = 0.0;
    ctrl->ratiolimit = (double)(iparm[IPARM_DISTRIBUTION_LEVEL]);
    ctrl->blblokmin  = 90;
    ctrl->blblokmax  = 140;
    ctrl->blblokmin  = iparm[IPARM_MIN_BLOCKSIZE];
    ctrl->blblokmax  = iparm[IPARM_MAX_BLOCKSIZE];
    /* OOC works only with 1D structures */
    if(ctrl->ooc)
    {
        pastix_print( procnum, 0, "Force 1D distribution because of OOC \n" );
        ctrl->ratiolimit = INTVALMAX;
    }

    if (ctrl->autolevel)
        printf("ratiolimit=%lf\n", ctrl->ratiolimit );
    else
        printf("level2D=%ld\n", (long) ctrl->level2D );

    /* Save iparm for other options */
    ctrl->iparm = iparm;

    /*
     * Initialize architecture description
     */

    /* Id and number of MPI processes */
    ctrl->clustnum = procnum;
    ctrl->clustnbr = procnbr;

    /* Local informations */
    ctrl->local_nbcores = local_coresnbr;
    ctrl->local_nbthrds = local_thrdsnbr;
    ctrl->local_nbctxts = ctrl->local_nbthrds;

    /* Total information (should require a MPI_Reduce if different informations on each node) */
    ctrl->total_nbcores = ctrl->local_nbcores * procnbr;
    ctrl->total_nbthrds = ctrl->local_nbthrds * procnbr;

    /* Create the array of associativity bewteen MPI process ids and SMP node ids */
    /* Rq: We could use a MPI reduction for irregular pattern                     */
    /* TODO: insert back the number of MPI processes per node                     */
    MALLOC_INTERN(ctrl->clust2smp, ctrl->clustnbr, pastix_int_t);
    for(i=0; i < ctrl->clustnbr; i++)
        ctrl->clust2smp[i] = i;

    /* Create the array of associativity bewteen core ids and MPI process ids */
    /* Rq: We could use a MPI reduction for irregular pattern             */
    MALLOC_INTERN(ctrl->core2clust, ctrl->total_nbcores, pastix_int_t);
    for(i=0; i < ctrl->total_nbcores; i++)
        ctrl->core2clust[i] = i / ctrl->local_nbcores;

    ctrl->egraph  = NULL;
    ctrl->etree   = NULL;
    ctrl->costmtx = NULL;
    ctrl->candtab = NULL;

    MALLOC_INTERN(ctrl->lheap, 1, Queue);
    queueInit(ctrl->lheap, 1000);

    MALLOC_INTERN(ctrl->intvec,  1, ExtendVectorINT);
    MALLOC_INTERN(ctrl->intvec2, 1, ExtendVectorINT);
    extendint_Init(ctrl->intvec,  10);
    extendint_Init(ctrl->intvec2, 10);

#ifdef PASTIX_DYNSCHED
    MALLOC_INTERN(ctrl->btree, 1, BubbleTree);
#endif

    return PASTIX_SUCCESS;
}


void blendCtrlExit(BlendCtrl *ctrl)
{
    queueExit(ctrl->lheap);
    memFree_null(ctrl->lheap);

    extendint_Exit(ctrl->intvec);
    memFree_null(ctrl->intvec);

    extendint_Exit(ctrl->intvec2);
    memFree_null(ctrl->intvec2);

    if(ctrl->clust2smp)
        memFree_null(ctrl->clust2smp);
    if(ctrl->core2clust)
        memFree_null(ctrl->core2clust);
    if(ctrl->candtab)
        memFree_null(ctrl->candtab);
}
