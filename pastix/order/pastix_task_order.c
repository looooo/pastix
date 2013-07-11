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

int orderAddIsolate( Order        *ordemesh,
                     pastix_int_t  new_n,
                     pastix_int_t *perm )
{
    Order ordesave;
    pastix_int_t i, ip;
    pastix_int_t n       = ordemesh->vertnbr;
    pastix_int_t cblknbr = ordemesh->cblknbr;
    int          baseval = ordemesh->baseval;

    assert( n <= new_n );

    /* Quick return */
    if ( n == new_n )
        return PASTIX_SUCCESS;

    memcpy( &ordesave, ordemesh, sizeof(Order) );
    orderInit( ordemesh, new_n, cblknbr + 1 );
    ordemesh->baseval = baseval;

    for(i=0; i< new_n; i++) {
        ip = perm[i];
        if (ip < n-baseval)
            ordemesh->permtab[i] = ordesave.permtab[ ip ];
        else
            ordemesh->permtab[i] = ip+baseval;
    }
    for(i=0; i<new_n; i++) {
        ip = ordemesh->permtab[i] - baseval;
        assert( (ip > -1) && (ip < new_n) );
        ordemesh->peritab[ip] = i + baseval;
    }

    /* Copy the cblknbr+1 first element of old rangtab and add last element */
    memcpy( &(ordemesh->rangtab), &(ordesave.rangtab), ordemesh->cblknbr * sizeof(pastix_int_t) );
    ordemesh->rangtab[ ordemesh->cblknbr ] = new_n + baseval;

    orderExit( &ordesave );
    return PASTIX_SUCCESS;
}

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
    pastix_int_t  schur_n;
    pastix_int_t *schur_colptr;
    pastix_int_t *schur_rows;
    pastix_int_t *schur_perm;
    pastix_int_t  zeros_n;
    pastix_int_t *zeros_colptr;
    pastix_int_t *zeros_rows;
    pastix_int_t *zeros_perm;
    pastix_int_t *iparm = pastix_data->iparm;
    pastix_csc_t  csc;
    Order        *ordemesh;
    Clock         timer;
    int           procnum;
    int           retval = PASTIX_SUCCESS;
    int           retval_rcv;

    /* Clean ordering if it exists */
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
     *
     */
    if (!(PASTIX_MASK_ISTRUE(iparm[IPARM_ORDERING], API_ORDER_LOAD))) {
        orderPrepareCSC( pastix_data, n, colptr, rows, loc2glob );
        pastix_data->bmalcolrow = 1;
    }

    /*
     * Isolate Shur elements
     */
    if ((iparm[IPARM_SCHUR] == API_YES) &&
        (pastix_data->schur_n > 0) )
    {
        assert( pastix_data->loc2glob2 == NULL );
        assert( pastix_data->schur_list != NULL );
        cscIsolate(pastix_data->n2,
                   pastix_data->col2,
                   pastix_data->row2,
                   pastix_data->schur_n,
                   pastix_data->schur_list,
                   &schur_colptr,
                   &schur_rows,
                   &schur_perm,
                   NULL);

        schur_n   = n - pastix_data->schur_n;
    } else {
        schur_n      = n;
        schur_colptr = pastix_data->col2;
        schur_rows   = pastix_data->row2;
    }

    /*
     * Isolate diagonal elements close to 0.
     */
    if ((iparm[IPARM_ISOLATE_ZEROS] == API_YES) &&
        (pastix_data->zeros_n > 0) )
    {
        assert( pastix_data->loc2glob2 == NULL );
        assert( pastix_data->zeros_list != NULL );
        cscIsolate(schur_n,
                   schur_colptr,
                   schur_rows,
                   pastix_data->zeros_n,
                   pastix_data->zeros_list,
                   &zeros_colptr,
                   &zeros_rows,
                   &zeros_perm,
                   NULL);

        zeros_n   = schur_n - pastix_data->zeros_n;
    } else {
        zeros_n      = schur_n;
        zeros_colptr = schur_colptr;
        zeros_rows   = schur_rows;
    }

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
        pastix_print(procnum, 0, "%s", OUT_ORDERINIT);

    clockStart(timer);

    csc.gN       = pastix_data->gN;
    csc.n        = zeros_n;
    csc.colptr   = zeros_colptr;
    csc.rows     = zeros_rows;
    csc.loc2glob = pastix_data->loc2glob2;

    /* Select the ordering method chosen by the user */
    switch (iparm[IPARM_ORDERING]) {
        /*
         * Scotch Ordering
         */
    case API_ORDER_SCOTCH:
#if defined(HAVE_SCOTCH)
        retval = orderComputeScotch( pastix_data, &csc );
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
        retval = orderComputePTScotch( pastix_data, &csc );
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
        retval = orderComputeMetis( pastix_data, &csc );
        /* orderFindSupernodes( pastix_data->n2, */
        /*                      pastix_data->col2, */
        /*                      pastix_data->row2, */
        /*                      ordemesh, NULL); */
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
            orderInit(ordemesh, n, 0);
            memcpy(ordemesh->permtab, perm, n*sizeof(pastix_int_t));
            memcpy(ordemesh->peritab, invp, n*sizeof(pastix_int_t));

            /* orderLoadFiles( pastix_data ); */
            /* orderFindSupernodes( pastix_data->n2, */
            /*                      pastix_data->col2, */
            /*                      pastix_data->row2, */
            /*                      ordemesh, NULL); */
        }
        break;

        /*
         * Load ordering with Scotch Format
         */
    case API_ORDER_LOAD:
        retval = orderLoadFiles( pastix_data, &csc );
        pastix_data->n2        = csc.n;
        pastix_data->col2      = csc.colptr;
        pastix_data->row2      = csc.rows;
        pastix_data->loc2glob2 = csc.loc2glob;
        break;

    default:
        errorPrint("Ordering not available");
        retval = PASTIX_ERR_BADPARAMETER;
        break;
    }

    /*
     * Add the isolated elements to the ordering structure
     */
    if ((iparm[IPARM_ISOLATE_ZEROS] == API_YES) &&
        (pastix_data->zeros_n > 0) )
    {
        orderAddIsolate( ordemesh, schur_n, zeros_perm );

        if ( zeros_colptr != schur_colptr ) { memFree_null( zeros_colptr ); }
        if ( zeros_rows   != schur_rows   ) { memFree_null( zeros_rows   ); }
        if ( zeros_perm   != NULL         ) { memFree_null( zeros_perm   ); }
    }

    /*
     * Add the isolated elements to the ordering structure
     */
    if ((iparm[IPARM_SCHUR] == API_YES) &&
        (pastix_data->schur_n > 0) )
    {
        orderAddIsolate( ordemesh, n, schur_perm );

        if ( schur_colptr != pastix_data->col2 ) { memFree_null( schur_colptr ); }
        if ( schur_rows   != pastix_data->row2 ) { memFree_null( schur_rows   ); }
        if ( schur_perm   != NULL              ) { memFree_null( schur_perm   ); }
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
