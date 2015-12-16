#include "common.h"
#include "order.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_reordering
 * @ingroup pastix
 *
 * pastix_task_reordering - Apply the reordering step, to reorder unknowns
 * within each supernode. It takes as input the ordering provided by
 * pastix_task_order(), and the symbolic factorization computed by
 * pastix_task_symbfact(). It returns boths structures Order and Symbol updated.
 *
 * This function reorders the unknowns of the problem based on traveller
 * salesman problem to gather together the contibutions facing each supernodes.
 *
 * See reordering paper (TODO: Put link to the paper when published)
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
 *          If set to API_IO_SAVE, the symbmtx and the generated ordemesh are
 *          dump to file.
 *          If set to API_IO_LOAD, the symbmtx (only) is loaded from the files.
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
    //EliminTree   *etree;
    Clock         timer;
    pastix_int_t *iparm;
    Order        *ordemesh;
    pastix_int_t  procnum;

    /**
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_task_reordering: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    iparm    = pastix_data->iparm;
    procnum  = pastix_data->procnum;
    ordemesh = pastix_data->ordemesh;

    assert(ordemesh->rangtab);
    assert(ordemesh->treetab);

    /* Start the step */
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO) {
        pastix_print(procnum, 0, OUT_STEP_REORDER,
                     iparm[IPARM_REORDERING_SPLIT],
                     iparm[IPARM_REORDERING_STOP]);
    }

    /* Print the reordering complexity */
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
        symbolReorderingPrintComplexity( pastix_data->symbmtx );

    clockStart(timer);

    /**
     * Build the elimination tree from the symbolic partition to reorder the
     * supernodes
     */
    //etree = eTreeBuild(symbmtx);

    /**
     * Reorder the rows of each supernode in order to compact coupling blocks
     */
    symbolReordering( pastix_data->symbmtx, ordemesh,
                      iparm[IPARM_REORDERING_SPLIT],
                      iparm[IPARM_REORDERING_STOP] );

    /* Backup the new ordering */
    if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE))
    {
        if (procnum == 0) {
            orderSave( ordemesh, NULL );
        }
    }

    symbolExit(pastix_data->symbmtx);
    memFree_null(pastix_data->symbmtx);
    pastix_data->symbmtx = NULL;

    /* Re-build the symbolic structure */
    /* TODO: Create a function to update the symbolic factorization, instead of computing it from scratch,
     this can be easily made in // per column, as opposed to the symbolic factorization itself */
    pastix_task_symbfact( pastix_data, NULL, NULL );

    clockStop(timer);

#if !defined(NDEBUG)
    if ( orderCheck( ordemesh ) != 0) {
        errorPrint("pastix_task_reordering: orderCheck on final ordering failed !!!");
        assert(0);
    }
    if( symbolCheck(pastix_data->symbmtx) != 0 ) {
        errorPrint("pastix_task_reordering: symbolCheck on final symbol matrix failed !!!");
        assert(0);
    }
#endif

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO) {
        pastix_print(procnum, 0, OUT_REORDERING_TIME,
                     (double)clockVal(timer));
    }
    return PASTIX_SUCCESS;
}
