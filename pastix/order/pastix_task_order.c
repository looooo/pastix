/*
 *  File: pastix_task_order.c
 *
 *  Function that handle the ordering step of PaStiX. It is a wrapper
 *  to the different ordering library that can be used by PaStiX.
 *
 *  Authors:
 *    Mathieu  Faverge    - faverge@labri.fr
 *    Xavier   LACOSTE    - lacoste@labri.fr
 *    Pierre   RAMET      - ramet@labri.fr
 *
 *  Dates:
 *    Version 0.0 - from 08 may 1998
 *                  to   08 jan 2001
 *    Version 1.0 - from 06 jun 2002
 *                  to   06 jun 2002
 */
#include "common.h"
#include "csc_utils.h"
#include "cscd_utils_intern.h"

/*
  Function: pastix_task_order

  Execute ordering task, with a centralised graph.

  Free *col2*  and *row2* entries of pastix_data if <pastix_task_order>
  has already been called.

  Set *col2*, *row2* and *n2* to a copy of user's CSC.

  Symmetrize this CSC.

  Remove diagonal elements from it.

  Clean last oredering if it exists.
  Depending on *IPARM_ORDERING* :
  - Calls Scotch ordering,
  - Calls Metis ordering,
  - Uses user oredering,
  - Loads oredering stored on disk in a Scotch format.

  Can save computed ordering on disk.

  returns compuited ordering into user arrays.

  Parameters:
  pastix_data - PaStiX data structure.
  perm        - permutation tabular
  invp        - reverse permutation tabular
*/
int pastix_task_order(pastix_data_t *pastix_data,
                      pastix_int_t   n,
                      pastix_int_t  *colptr,
                      pastix_int_t  *rows,
                      pastix_int_t  *loc2glob,
                      pastix_int_t  *perm,
                      pastix_int_t  *invp)
{
    pastix_int_t    schur_n;
    pastix_int_t   *schur_colptr;
    pastix_int_t   *schur_rows;
    pastix_int_t   *schur_perm;
    pastix_int_t    zeros_n;
    pastix_int_t   *zeros_colptr;
    pastix_int_t   *zeros_rows;
    pastix_int_t   *zeros_perm;
    pastix_int_t   *iparm = pastix_data->iparm;
    pastix_graph_t  subgraph;
    pastix_graph_t *graph;
    Order          *ordemesh;
    Clock           timer;
    int             procnum;
    int             retval = PASTIX_SUCCESS;
    int             retval_rcv;
    int             do_schur = 1;
    int             do_zeros = 1;

    /*
     * Check parameters
     */
    if ((iparm[IPARM_SCHUR] == API_YES) &&
        (pastix_data->schur_n > 0) )
    {
        /*
         * If ordering is set to API_ORDER_PERSONAL or API_ORDER_LOAD, we
         * consider that the schur complement is already isolated at the end of
         * permutation array
         */
        if ((iparm[IPARM_ORDERING] == API_ORDER_PERSONAL) ||
            (iparm[IPARM_ORDERING] == API_ORDER_LOAD) ) {
            do_schur = 0;
        }
    } else {
        do_schur = 0;
    }
    if ((iparm[IPARM_ISOLATE_ZEROS] == API_YES) &&
        (pastix_data->zeros_n > 0) )
    {
        /*
         * If ordering is set to API_ORDER_PERSONAL or API_ORDER_LOAD, we
         * consider that the schur complement is already isolated at the end of
         * permutation array
         */
        if ((iparm[IPARM_ORDERING] == API_ORDER_PERSONAL) ||
            (iparm[IPARM_ORDERING] == API_ORDER_LOAD) ) {
            do_zeros = 0;
        }
    } else {
        do_zeros = 0;
    }

    /*
     * Clean ordering if it exists
     */
    if (pastix_data->ordemesh != NULL) {
        orderExit(pastix_data->ordemesh);
    } else {
        MALLOC_INTERN( pastix_data->ordemesh, 1, Order );
    }

    ordemesh = pastix_data->ordemesh;
    procnum  = pastix_data->procnum;
    orderInit( ordemesh, 0, 0 );

    print_debug(DBG_STEP, "-> pastix_task_order\n");
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        pastix_print(procnum, 0, "%s", OUT_STEP_ORDER);

    /*
     * Prepare a copy of user's CSC
     * Copy the given csc in pastix_data structure and performs basic required
     * operations as symmetrizing the graph and removing the diagonal
     * coefficients
     */
    graphPrepare( pastix_data, n, colptr, rows, loc2glob, &(pastix_data->csc) );
    graph = pastix_data->csc;
    pastix_data->n2 = graph->n;
    pastix_data->gN = graph->gN;

    /*
     * Isolate Shur elements
     */
    if ( do_schur )
    {
        assert( graph->loc2glob == NULL );
        assert( pastix_data->schur_list != NULL );
        graphIsolate(graph->n,
                     graph->colptr,
                     graph->rows,
                     pastix_data->schur_n,
                     pastix_data->schur_list,
                     &schur_colptr,
                     &schur_rows,
                     &schur_perm,
                     NULL);

        schur_n = n - pastix_data->schur_n;
    } else {
        schur_n      = n;
        schur_colptr = graph->colptr;
        schur_rows   = graph->rows;
    }

    /*
     * Isolate diagonal elements close to 0.
     */
    if ( do_zeros )
    {
        assert( graph->loc2glob == NULL );
        assert( pastix_data->zeros_list != NULL );
        graphIsolate(schur_n,
                     schur_colptr,
                     schur_rows,
                     pastix_data->zeros_n,
                     pastix_data->zeros_list,
                     &zeros_colptr,
                     &zeros_rows,
                     &zeros_perm,
                     NULL);

        zeros_n = schur_n - pastix_data->zeros_n;
    } else {
        zeros_n      = schur_n;
        zeros_colptr = schur_colptr;
        zeros_rows   = schur_rows;
    }

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
        pastix_print(procnum, 0, "%s", OUT_ORDERINIT);

    clockStart(timer);

    subgraph.gN       = graph->gN;
    subgraph.n        = zeros_n;
    subgraph.colptr   = zeros_colptr;
    subgraph.rows     = zeros_rows;
    subgraph.loc2glob = graph->loc2glob;

    /* Select the ordering method chosen by the user */
    switch (iparm[IPARM_ORDERING]) {
        /*
         * Scotch Ordering
         */
    case API_ORDER_SCOTCH:
#if defined(HAVE_SCOTCH)
        retval = orderComputeScotch( pastix_data, &subgraph );
#else
        errorPrint("PASTIX Ordering: Ordering with Scotch requires to enable -DPASTIX_ORDERING_SCOTCH option");
        retval = PASTIX_ERR_BADPARAMETER;
#endif
        break;

        /*
         * PT-Scotch Ordering
         */
    case API_ORDER_PTSCOTCH:
#if defined(HAVE_PTSCOTCH)
        retval = orderComputePTScotch( pastix_data, &subgraph );
#else
        errorPrint("PASTIX Ordering: Ordering with PT-Scotch requires to enable -DPASTIX_ORDERING_PTSCOTCH option");
        retval = PASTIX_ERR_BADPARAMETER;
#endif
        break;

        /*
         *  METIS ordering
         */
    case API_ORDER_METIS:
#if defined(HAVE_METIS)
        retval = orderComputeMetis( pastix_data, &subgraph );
        assert( ordemesh->rangtab == NULL );
#else
        errorPrint("PASTIX Ordering: Ordering with Metis requires -DHAVE_METIS flag at compile time");
        retval = PASTIX_ERR_BADPARAMETER;
#endif
        break;

        /*
         * Personal Ordering
         */
    case API_ORDER_PERSONAL:
        {
            /* orderInit(ordemesh, n, 0); */
            /* memcpy(ordemesh->permtab, perm, n*sizeof(pastix_int_t)); */
            /* memcpy(ordemesh->peritab, invp, n*sizeof(pastix_int_t)); */

            orderLoadFiles( pastix_data );
            memFree_null( ordemesh->rangtab );
            ordemesh->cblknbr = 0;
        }
        break;

        /*
         * Load ordering
         */
    case API_ORDER_LOAD:
        retval = orderLoadFiles( pastix_data );
        break;

    default:
        errorPrint("Ordering not available");
        retval = PASTIX_ERR_BADPARAMETER;
        break;
    }

    if (retval != PASTIX_SUCCESS )
        return retval;

    /*
     * Add the isolated elements to the ordering structure
     */
    if ( do_zeros )
    {
        orderAddIsolate( ordemesh, schur_n, zeros_perm );

        if ( zeros_colptr != schur_colptr ) { memFree_null( zeros_colptr ); }
        if ( zeros_rows   != schur_rows   ) { memFree_null( zeros_rows   ); }
        if ( zeros_perm   != NULL         ) { memFree_null( zeros_perm   ); }
    }

    /*
     * Add the isolated elements to the ordering structure
     */
    if ( do_schur )
    {
        orderAddIsolate( ordemesh, n, schur_perm );

        if ( schur_colptr != graph->colptr ) { memFree_null( schur_colptr ); }
        if ( schur_rows   != graph->rows   ) { memFree_null( schur_rows   ); }
        if ( schur_perm   != NULL          ) { memFree_null( schur_perm   ); }
    }

    /* Reduce the error code */
    MPI_Allreduce(&retval, &retval_rcv, 1, MPI_INT, MPI_MAX,
                  pastix_data->pastix_comm);
    if (retval_rcv != PASTIX_SUCCESS)
        return retval_rcv;

    /* Rebase the ordering to 0 */
    orderBase(ordemesh, 0);

    clockStop(timer);
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        pastix_print(procnum, 0, TIME_COMPUTE_ORDERING, clockVal(timer));

    /* Save i/o strategy */
    if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE)) {
        retval = orderSaveFiles( pastix_data );
        if (retval != PASTIX_SUCCESS)
            return retval;
    }

    /*
     * Return the ordering to user if perm/invp are not NULL
     * Remark: No need to copy back for personal
     */
    if (iparm[IPARM_ORDERING] != API_ORDER_PERSONAL) {
        if (perm != NULL) memcpy(perm, ordemesh->permtab, n*sizeof(pastix_int_t));
        if (invp != NULL) memcpy(invp, ordemesh->peritab, n*sizeof(pastix_int_t));
    }

    iparm[IPARM_START_TASK]++;
    return PASTIX_SUCCESS;
}
