#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "common.h"
#include "param_blend.h"

/* we set default option */
pastix_int_t blendParamInit(BlendParam  *param,
                            pastix_int_t procnum,
                            pastix_int_t *iparm )
{
    param->count_ops        = 1;
    param->debug            = 0;
    param->timer            = 1;
    param->leader           = 0;
    param->allcand          = 0;
    param->nocrossproc      = 0;
    param->level2D          = 100000000;
    param->costlevel        = 1;
    param->autolevel        = 0;
    /* On suppose que tous les noeuds smp utilisés ont la même configuration */
    param->procnbr          = sysconf(_SC_NPROCESSORS_ONLN);
    param->ratiolimit       = 0.0;
    param->ooc              = 0;

    param->smpnbr = iparm[IPARM_NB_SMP_NODE_USED];
    param->ricar  = iparm[IPARM_INCOMPLETE];
    param->abs    = iparm[IPARM_ABS];

    /* TODO: differentiate the iparm */
    param->blcolmin  = 60;
    param->blcolmax  = 120;
    param->blblokmin = 90;
    param->blblokmax = 140;
    param->blcolmin  = iparm[IPARM_MIN_BLOCKSIZE];
    param->blcolmax  = iparm[IPARM_MAX_BLOCKSIZE];
    param->blblokmin = iparm[IPARM_MIN_BLOCKSIZE];
    param->blblokmax = iparm[IPARM_MAX_BLOCKSIZE];

    if(param->blcolmin > param->blcolmax)
    {
        errorPrint("Parameter error : blocksize max < blocksize min (cf. iparm.txt).");
        ASSERT(param->blcolmin <=  param->blcolmax, MOD_SOPALIN);
    }

    param->level2D    = iparm[IPARM_DISTRIBUTION_LEVEL];
    param->ratiolimit = (double)(iparm[IPARM_DISTRIBUTION_LEVEL]);

    if (param->autolevel)
        printf("ratiolimit=%lf\n",param->ratiolimit);
    else
        printf("level2D=%ld\n", (long) param->level2D);

    if(param->ooc)
    {
        /* OOC works only with 1D structures */
        pastix_print( procnum, 0, "Force 1D distribution because of OOC \n" );
        param->ratiolimit = INTVALMAX;
    }

    param->iparm = iparm;
    /* param->dparm = dparm; */

    return 1;
}

void blendParamExit(BlendParam *param)
{
    memFree_null(param);
}
