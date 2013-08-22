#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "cost.h"



pastix_int_t costInit(CostMatrix *costmtx)
{
    costmtx->bloktab = NULL;
    return 1;
}

void costExit(CostMatrix *costmtx)
{
    if(costmtx->bloktab != NULL)
	memFree_null(costmtx->bloktab);
    memFree_null(costmtx);
}
