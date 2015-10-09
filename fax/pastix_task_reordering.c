#include "common.h"
#include "order.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_reordering
 * @ingroup pastix
 *
 * pastix_task_reordering - Reorders unknowns given a Scotch parition and the
 * corresponding symbolc factoriszation
 *
 * The function is an algorithm to reorder the problem. It takes as input the
 * ordemesh structure, as well as the symbolic factorization and returns both a
 * new ordemesh and a new symbolic factorization.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_IO_STRATEGY, IPARM_REORDERING_SPLIT,
 *   IPARM_REORDERING_STOP
 *
 *******************************************************************************
 *
 * @param[in,out] pastix_data
 *          The pastix_data structure that describes the solver instance.  On
 *          exit, the field symbmtx is updated with the new symbol matrix, and
 *          the field ordemesh is updated with the new ordering.
 *          - IPARM_IO_STRATEGY will enable to load/store the result to files.
 *          If set to API_IO_SAVE, the symbmtx and the generated ordemesh is
 *          dump to file.
 *          If set to APÏ_IO_LOAD, the symbmtx (only) is loaded from the files.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on successful exit
 *          \retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/

int
pastix_task_reordering(pastix_data_t *pastix_data)
{
    Clock         timer;
    pastix_int_t *iparm;
    Order        *ordemesh;
    pastix_int_t  procnum;

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_task_symbfact: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    iparm    = pastix_data->iparm;
    procnum  = pastix_data->procnum;
    ordemesh = pastix_data->ordemesh;

    /* Compute the Reordering complexity */
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        symbolReorderingPrintComplexity( pastix_data->symbmtx );


    clockStart(timer);

    pastix_print(procnum, 0, OUT_REORDERING_PARAMS,
                 iparm[IPARM_REORDERING_SPLIT],
                 iparm[IPARM_REORDERING_STOP]);

    symbolReordering( pastix_data->symbmtx, ordemesh,
                      iparm[IPARM_REORDERING_SPLIT],
                      iparm[IPARM_REORDERING_STOP] );

    if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE))
    {
        if (procnum == 0) {
            orderSave( ordemesh, NULL );
        }
    }

    symbolExit(pastix_data->symbmtx);
    memFree_null(pastix_data->symbmtx);
    pastix_data->symbmtx = NULL;

    pastix_task_symbfact( pastix_data, NULL, NULL );

    clockStop(timer);
    pastix_print(procnum, 0, OUT_REORDERING_TIME,
                 (double)clockVal(timer));

    if (pastix_data->graph != NULL)
    {
        graphExit( pastix_data->graph );
        memFree_null( pastix_data->graph );
    }

    return PASTIX_SUCCESS;
}
