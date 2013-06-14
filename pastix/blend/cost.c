#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "cost.h"



PASTIX_INT costInit(CostMatrix *costmtx)
{
    costmtx->cblktab = NULL;
    costmtx->bloktab = NULL;
    return 1;
}
	
void costExit(CostMatrix *costmtx)
{
    if(costmtx->cblktab != NULL)
	memFree_null(costmtx->cblktab);
    if(costmtx->bloktab != NULL)
	memFree_null(costmtx->bloktab);
    memFree_null(costmtx);
}
