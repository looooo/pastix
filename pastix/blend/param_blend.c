#include <stdio.h>
#include <stdlib.h>

#include "common_pastix.h"
#include "param_blend.h"

/* we set default option */
PASTIX_INT blendParamInit(BlendParam *param)
{
  param->hpf_filename     = NULL;
  param->trace_filename   = "traceBlend.trf";
  param->ps_filename      = "matrix.ps";
  param->hpf              = 0;
  param->tracegen         = 0;
  param->ps               = 0;
  param->assembly         = 0;
  param->solvmtx_filename = "solvmtx.";
  param->sequentiel       = 0;
  param->count_ops        = 1;
  param->debug            = 0;
  param->timer            = 1;
  param->recover          = 1;
  param->blcolmin         = 60;
  param->blcolmax         = 120;
  param->blblokmin        = 90;
  param->blblokmax        = 140;
  param->leader           = 0;
  param->allcand          = 0;
  param->nocrossproc      = 0;
  param->forceE2          = 0;
  param->level2D          = 100000000;
  param->candcorrect      = 0;
  param->clusterprop      = 0;
  param->costlevel        = 1;
  param->autolevel        = 0;
  param->malt_limit       = -1;
  param->smpnbr           = 0;
  /* On suppose que tous les noeuds smp utilisés ont la même configuration */
  param->procnbr          = sysconf(_SC_NPROCESSORS_ONLN);
  param->ratiolimit       = 0.0;
  param->dense_endblock   = 0;
  param->ooc              = 0;
  param->oocmemlimit      = 4e7;
  param->abs              = 4;
  param->ricar            = 0;
  
  return 1;
}

void blendParamExit(BlendParam *param)
{
  memFree_null(param->hpf_filename);
  memFree_null(param->trace_filename);
  memFree_null(param->ps_filename);
  memFree_null(param->solvmtx_filename);
  memFree_null(param);
}


