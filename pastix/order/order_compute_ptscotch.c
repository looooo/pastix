/*
 *  File: order_compute_scotch.c
 *
 *  Wrapper to compute the ordering with Scotch Library.
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
#include <scotch.h>
#include <ptscotch.h>
#include "csc_utils.h"
#include "scotch_strats.h"
#include "cscd_utils_intern.h"

int orderComputePTScotch( pastix_data_t *pastix_data, const pastix_csc_t *csc )
{
    SCOTCH_Dordering ordedat;
    SCOTCH_Ordering  ordering;
    SCOTCH_Dgraph    dgraph;

    Order        *ordemesh = &(pastix_data->ordemesh);
    SCOTCH_Graph *grafmesh = &(ordemesh->grafmesh);
    SCOTCH_Strat  stratdat;
    char          strat[1024];
    pastix_int_t *colptr, *colptr_schur;
    pastix_int_t *rows, *rows_schur;
    pastix_int_t *perm_schur, *invp_schur;
    pastix_int_t *iparm = pastix_data->iparm;
    pastix_int_t  procnum;
    pastix_int_t n;
    pastix_int_t nnz;
    int ret;

    procnum   = pastix_data->procnum;

    /* Check integer compatibility */
    if (sizeof(pastix_int_t) != sizeof(SCOTCH_Num)) {
        errorPrint("Inconsistent integer type\n");
        return INTEGER_TYPE_ERR;
    }

    gN     = csc->gN;
    n      = csc->n;
    colptr = csc->colptr;
    rows   = csc->rows;
    nnz    = colptr[n] - 1;

    /* Build distributed graph */
    print_debug(DBG_SCOTCH, "> SCOTCH_dgraphInit <\n");
    SCOTCH_dgraphInit(dgraph, pastix_comm);

    print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_graphBuild <\n");
    if ( SCOTCH_dgraphBuild (&dgraph,
                             colptr[0],    /* baseval */
                             n,            /* number of local vertices */
                             n,            /* Maximum number of local vertices     */
                             colptr,
                             NULL,
                             NULL,         /* Local vertex load array (if any)     */
                             NULL,         /* Local vertex label array (if any)    */
                             nnz,
                             nnz,
                             rows,         /* Local edge array                     */
                             NULL,         /* Ghost edge array (if any); not const */
                             NULL))
        {
            errorPrint("SCOTCH_dgraphBuild");
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }

    print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_dgraphCheck <\n");
    if (SCOTCH_dgraphCheck(&dgraf)) {
        errorPrint("pastix: SCOTCH_dgraphCheck");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    print_debug(DBG_SCOTCH, "> SCOTCH_stratInit <\n");
    if (SCOTCH_stratInit(&stratdat))
    {
        errorPrint("pastix : SCOTCH_stratInit");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    /* TODO : Add default strategies for PT-Scotch */
    if (iparm[IPARM_DEFAULT_ORDERING] == API_YES)
    {
        if (iparm[IPARM_INCOMPLETE] == API_NO)
            sprintf(strat, PTSCOTCH_STRAT_DIRECT);
        else
            sprintf(strat, PTSCOTCH_STRAT_INCOMP);
    }
    else /* Personal strategy */
    {
        sprintf(strat, PTSCOTCH_STRAT_PERSO,
                (long) iparm[IPARM_ORDERING_SWITCH_LEVEL],
                (long) iparm[IPARM_ORDERING_CMIN],
                (long) iparm[IPARM_ORDERING_CMAX],
                ((float)iparm[IPARM_ORDERING_FRAT])/100.,
                (long) iparm[IPARM_ORDERING_SWITCH_LEVEL],
                (long) iparm[IPARM_ORDERING_CMIN],
                (long) iparm[IPARM_ORDERING_CMAX],
                ((float)iparm[IPARM_ORDERING_FRAT])/100.);

        if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
            print_onempi("PT-Scotch Strategy |%s|\n", strat);
    }

    clockStart(timer1);

    /*    print_debug(DBG_SCOTCH, "> SCOTCH_stratDgraphOrder <\n"); */
    /*    if (SCOTCH_stratDgraphOrder(&stratdat, strat)) */
    /*      { */
    /*        errorPrint("pastix : SCOTCH_stratDgraphOrder"); */
    /*        EXIT(MOD_SOPALIN,INTERNAL_ERR); */
    /*      } */
    if (procnum == 0 && iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        errorPrintW("PaStiX works only with PT-Scotch default strategy");

    print_debug(DBG_SCOTCH, "> SCOTCH_dgraphOrderInit <\n");
    if (0 != SCOTCH_dgraphOrderInit(dgraph, ordedat))
    {
        errorPrint("pastix : SCOTCH_dgraphOrderInit");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    print_debug(DBG_SCOTCH, "> SCOTCH_dgraphOrderCompute <\n");
    if (0 != SCOTCH_dgraphOrderCompute(dgraph, ordedat, &stratdat))
    {
        errorPrint("pastix : SCOTCH_dgraphOrderCompute");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    print_debug(DBG_SCOTCH, "> SCOTCH_stratExit <\n");
    SCOTCH_stratExit(&stratdat);

    /* print_debug(DBG_SCOTCH, "> SCOTCH_dgraphOrderPerm <\n"); */
    /*       if (0 != SCOTCH_dgraphOrderPerm(dgraph, ordedat, perm)) */
    /*  { */
    /*    errorPrint("pastix : SCOTCH_dgraphOrderPerm"); */
    /*    EXIT(MOD_SOPALIN,INTERNAL_ERR); */
    /*  } */

    clockStop(timer1);
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        print_onempi(TIME_COMPUTE_ORDERING, clockVal(timer1));

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
        print_onempi("%s", OUT_ORDERINIT);

    orderInit(ordemesh, gN, gN);
    memset( ordemesh->rangtab, 0, (gN+1)*sizeof(pastix_int_t));

    SCOTCH_dgraphCorderInit (&dgraph,
                             &ordering,
                             (SCOTCH_Num *)(ordemesh->permtab),
                             (SCOTCH_Num *)(ordemesh->peritab),
                             &ordemesh->cblknbr,
                             ordemesh->rangtab,
                             NULL);

    if (procnum == 0) {
        SCOTCH_dgraphOrderGather (&dgraph, &ordedat, &ordering);
    }
    else {
        SCOTCH_dgraphOrderGather (&dgraph, &ordedat, NULL);
    }

    MPI_Bcast(&ordemesh->cblknbr, 1,                     PASTIX_MPI_INT, 0, pastix_comm);
    MPI_Bcast( ordemesh->rangtab, (ordemesh->cblknbr+1), PASTIX_MPI_INT, 0, pastix_comm);
    MPI_Bcast( ordemesh->permtab, gN,                    PASTIX_MPI_INT, 0, pastix_comm);
    MPI_Bcast( ordemesh->peritab, gN,                    PASTIX_MPI_INT, 0, pastix_comm);

    global2localperm(n, perm, ((*pastix_data)->ordemesh).permtab, loc2glob);

    /* Gathering graph */
    print_debug(DBG_SCOTCH, "> SCOTCH_dgraphGather <\n");
    SCOTCH_dgraphGather( &dgraph, &(ordemesh->grafmesh) );
    SCOTCH_dgraphCorderExit( &dgraph, &ordering );
    SCOTCH_dgraphOrderExit( &dgraph, ordedat );
    SCOTCH_dgraphExit( &dgraph );

    SCOTCH_graphBase(&(ordemesh->grafmesh), 0);
    orderBase(ordemesh, 0);

#if defined(FORGET_PARTITION)
    ordemesh->cblknbr = 0;
    if (ordemesh->rangtab != NULL) memFree_null(ordemesh->rangtab);
#endif

    return PASTIX_SUCCESS;
}
