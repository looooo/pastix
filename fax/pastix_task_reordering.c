#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <pastix.h>
#include "common.h"
#include "order.h"

int
pastix_task_reordering(pastix_data_t *pastix_data)
{
    Clock timer;

    /* Compute the Reordering complexity */
    symbolReorderingPrintComplexity( pastix_data->symbmtx );

    clockStart(timer);

    symbolReordering( pastix_data->symbmtx, pastix_data->ordemesh,
                      pastix_data->iparm[IPARM_REORDERING_SPLIT],
                      pastix_data->iparm[IPARM_REORDERING_STOP] );

    /* !!! We have to save the new ordering !!! */
    /* int retval = orderSave( pastix_data->ordemesh, NULL ); */

    symbolExit(pastix_data->symbmtx);
    memFree_null(pastix_data->symbmtx);
    pastix_data->symbmtx = NULL;

    pastix_task_symbfact( pastix_data, NULL, NULL );

    clockStop(timer);
    printf("Total time for reordering (with extra symbolic factorization) %.3g s\n", (double)clockVal(timer));

    if (pastix_data->graph != NULL)
    {
        graphExit( pastix_data->graph );
        memFree_null( pastix_data->graph );
    }

    return PASTIX_SUCCESS;
}
