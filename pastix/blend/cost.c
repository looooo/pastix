#include <stdio.h>
#include <stdlib.h>

#include "common.h"
#include "symbol.h"
#include "cost.h"

pastix_int_t
costMatrixInit( CostMatrix *costmtx )
{
    costmtx->blokcost = NULL;
    costmtx->cblkcost = NULL;
    return 1;
}

void
costMatrixExit( CostMatrix *costmtx )
{
    if(costmtx->blokcost != NULL)
	memFree_null(costmtx->blokcost);

    if(costmtx->cblkcost != NULL)
	memFree_null(costmtx->cblkcost);
}


CostMatrix *
costMatrixBuild( const SymbolMatrix *symbmtx,
                 pastix_coeftype_t   flttype,
                 pastix_factotype_t  factotype )
{
    CostMatrix *costmtx = NULL;

    MALLOC_INTERN(costmtx, 1, CostMatrix);
    costMatrixInit(costmtx);

    MALLOC_INTERN( costmtx->cblkcost, symbmtx->cblknbr, double );
    MALLOC_INTERN( costmtx->blokcost, symbmtx->bloknbr, double );

    symbolGetTimes( symbmtx, flttype, factotype,
                    costmtx->cblkcost, costmtx->blokcost );

    return costmtx;
}
