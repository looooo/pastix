/**
 *
 * @file order_compute_ptscotch.c
 *
 *  PaStiX order routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains functions to perform ordering with PT-Scotch library.
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "graph.h"
#include "order.h"
#if (defined PASTIX_ORDERING_SCOTCH) || (defined PASTIX_ORDERING_PTSCOTCH)
#  ifdef    PASTIX_ORDERING_PTSCOTCH
#    include <ptscotch.h>
#  else  /* PASTIX_ORDERING_PTSCOTCH */
#    include <scotch.h>
#  endif /* PASTIX_ORDERING_PTSCOTCH */
#endif /* PASTIX_ORDERING_PTSCOTCH || PASTIX_ORDERING_SCOTCH */
#include "order_scotch_strats.h"
#include "d_cscd_utils_intern.h"

/* TODO: take care of this */
void global2localperm(pastix_int_t  lN,
                      pastix_int_t *lperm,
                      pastix_int_t *gperm,
                      pastix_int_t *loc2glob);
/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering
 *
 * orderComputePTScotch - Compute the ordering of the graph given as parameter
 * with PT-Scotch library.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_ORDERING_DEFAULT, IPARM_SCOTCH_SWITCH_LEVEL,
 *   IPARM_SCOTCH_CMIN, IPARM_SCOTCH_CMAX, IPARM_SCOTCH_FRAT
 *
 *******************************************************************************
 *
 * @param[in,out] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the field oerdemesh is initialize with the result of the
 *          ordering realized by Scotch.
 *
 * @param[in] graph
 *          The graph prepared by graphPrepare function on which wwe want to
 *          perform the ordering.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on successful exit
 *          \retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *          \retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *          \retval PASTIX_ERR_INTEGER_TYPE if Scotch integer type is not the
 *                  same size as PaStiX ones.
 *          \retval PASTIX_ERR_INTERNAL if an error occurs internally to Scotch.
 *
 *******************************************************************************/
int
orderComputePTScotch(       d_pastix_data_t  *pastix_data,
                      const pastix_graph_t *graph )
{
    SCOTCH_Dordering ordedat;
    SCOTCH_Ordering  ordering;
    SCOTCH_Dgraph    dgraph;
    SCOTCH_Strat     stratdat;
#if defined(PERSONAL_PTSCOTCH_STRATEGY)
    char             strat[1024];
#endif
    MPI_Comm      pastix_comm;
    Order        *ordemesh = pastix_data->ordemesh;
    Clock         timer;
    pastix_int_t *colptr;
    pastix_int_t *rows;
    pastix_int_t *iparm = pastix_data->iparm;
    pastix_int_t  procnum;
    pastix_int_t  n, gN, nnz;

    procnum     = pastix_data->procnum;
    pastix_comm = pastix_data->pastix_comm;

    /* Check integer compatibility */
    if (sizeof(pastix_int_t) != sizeof(SCOTCH_Num)) {
        errorPrint("orderComputePTScotch: Inconsistent integer type between Pastix and PT-Scotch\n");
        return PASTIX_ERR_INTEGER_TYPE;
    }

    gN     = graph->gN;
    n      = graph->n;
    colptr = graph->colptr;
    rows   = graph->rows;
    nnz    = colptr[n] - 1;

    /* Build distributed graph */
    print_debug(DBG_SCOTCH, "> SCOTCH_dgraphInit <\n");
    SCOTCH_dgraphInit(&dgraph, pastix_comm);

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
    if (SCOTCH_dgraphCheck(&dgraph)) {
        errorPrint("pastix: SCOTCH_dgraphCheck");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    print_debug(DBG_SCOTCH, "> SCOTCH_stratInit <\n");
    if (SCOTCH_stratInit(&stratdat))
    {
        errorPrint("pastix : SCOTCH_stratInit");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    /*
     * Create Strategy string for Scotch
     */
    /* TODO : Add default strategies for PT-Scotch */
#if defined(PERSONAL_PTSCOTCH_STRATEGY)
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
            pastix_print(procnum, 0, "PT-Scotch Strategy |%s|\n", strat);
    }

    print_debug(DBG_SCOTCH, "> SCOTCH_stratDgraphOrder <\n");
    if (SCOTCH_stratDgraphOrder(&stratdat, strat))
    {
        errorPrint("pastix : SCOTCH_stratDgraphOrder");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }
#else
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        pastix_print(procnum, 0, "PaStiX works only with PT-Scotch default strategy %s", "");
#endif

    clockStart(timer);
    print_debug(DBG_SCOTCH, "> SCOTCH_dgraphOrderInit <\n");
    if (0 != SCOTCH_dgraphOrderInit(&dgraph, &ordedat))
    {
        errorPrint("pastix : SCOTCH_dgraphOrderInit");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    print_debug(DBG_SCOTCH, "> SCOTCH_dgraphOrderCompute <\n");
    if (0 != SCOTCH_dgraphOrderCompute(&dgraph, &ordedat, &stratdat))
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

    SCOTCH_dgraphCorderExit( &dgraph, &ordering );
    SCOTCH_dgraphOrderExit( &dgraph, &ordedat );
    SCOTCH_dgraphExit( &dgraph );

    orderBase(ordemesh, 0);

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        pastix_print(procnum, 0, TIME_COMPUTE_ORDERING, clockVal(timer));

    clockStop(timer);
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
        pastix_print(procnum, 0, "%s", OUT_ORDERINIT);

    return PASTIX_SUCCESS;
}
