#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "cost.h"

pastix_int_t
costMatrixInit( CostMatrix *costmtx )
{
    costmtx->bloktab = NULL;
    costmtx->cblkcost = NULL;
    return 1;
}

void
costMatrixExit( CostMatrix *costmtx )
{
    if(costmtx->bloktab != NULL)
	memFree_null(costmtx->bloktab);

    if(costmtx->cblkcost != NULL)
	memFree_null(costmtx->cblkcost);
}


CostMatrix *
costMatrixBuild( const SymbolMatrix * symbmtx,
                 const Dof * dofptr)
{
    CostMatrix *costmtx = NULL;
    pastix_int_t i;

    MALLOC_INTERN(costmtx, 1, CostMatrix);
    costMatrixInit(costmtx);

    MALLOC_INTERN(costmtx->bloktab,  symbmtx->bloknbr, CostBlok);
    MALLOC_INTERN(costmtx->cblkcost, symbmtx->cblknbr, double);

    memset( costmtx->cblkcost, 0, symbmtx->cblknbr * sizeof(double) );

    /* Init block cost */
    for(i=0;i<symbmtx->cblknbr;i++)
        cblkComputeCost(i, costmtx, symbmtx, dofptr);

    /* Init cblk cost */
    for(i=0;i<symbmtx->cblknbr;i++)
        cblkComputeCostLL(i, costmtx, symbmtx, dofptr);

    return costmtx;
}
