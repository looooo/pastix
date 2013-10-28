/**
 *
 * @file order_compute_scotch.c
 *
 *  PaStiX order routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains functions to perform ordering with Scotch library.
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
#include <scotch.h>
#include "order_scotch_strats.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering
 *
 * orderComputeScotch - Compute the ordering of the graph given as parameter
 * with Scotch library.
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
orderComputeScotch(       pastix_data_t  *pastix_data,
                    const pastix_graph_t *graph )
{
    Order        *ordemesh = pastix_data->ordemesh;
    SCOTCH_Graph  scotchgraph;
    SCOTCH_Strat  stratdat;
    char          strat[1024];
    pastix_int_t *colptr;
    pastix_int_t *rows;
    pastix_int_t *iparm = pastix_data->iparm;
    pastix_int_t  procnum;
    pastix_int_t n;
    pastix_int_t nnz;
    int ret;

    procnum = pastix_data->procnum;

    /* Check integer compatibility */
    if (sizeof(pastix_int_t) != sizeof(SCOTCH_Num)) {
        errorPrint("Inconsistent integer type\n");
        return PASTIX_ERR_INTEGER_TYPE;
    }

    n      = graph->n;
    colptr = graph->colptr;
    rows   = graph->rows;
    nnz    = colptr[n] - 1;

    print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_graphInit <\n");
    orderInit(ordemesh, n, n);
    SCOTCH_graphInit( &scotchgraph );

    print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_graphBuild <\n");
    if (SCOTCH_graphBuild(&scotchgraph,   /* Graph to build     */
                          1,              /* baseval            */
                          n,              /* Number of vertices */
                          colptr,         /* Vertex array       */
                          NULL,
                          NULL,           /* Array of vertex weights (DOFs) */
                          NULL,
                          nnz,            /* Number of arcs     */
                          rows,           /* Edge array         */
                          NULL))
        {
            errorPrint("pastix : graphBuildGraph");
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }

#if defined(PASTIX_DEBUG_ORDERING)
    {
        Clock timer;
        clockStart(timer);
        print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_graphCheck <\n");
        if (SCOTCH_graphCheck(&scotchgraph)) {
            errorPrint("pastix: graphCheck");
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }
        clockStop(timer);
        pastix_print( procnum, 0, "SCOTCH_graphCheck done in %lf second\n", clockVal(timer) );
    }
#endif
    SCOTCH_graphBase(&scotchgraph, 0);

    /* The graph is build, let's compute the ordering */
    SCOTCH_stratInit(&stratdat);

    /*
     * Create Strategy string for Scotch
     */
    /* default ordering */
    if (iparm[IPARM_DEFAULT_ORDERING] == API_YES) {
        if (iparm[IPARM_INCOMPLETE] == API_NO) {
            if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
                pastix_print(procnum, 0, "%s", "Scotch direct strategy\n");
            sprintf(strat, SCOTCH_STRAT_DIRECT);
        }
        else {
            if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
                pastix_print(procnum, 0, "%s", "Scotch incomplete strategy\n");
            sprintf(strat, SCOTCH_STRAT_INCOMP);
        }
    }
    /* personal ordering */
    else {
        sprintf(strat, SCOTCH_STRAT_PERSO,
                (long)  iparm[IPARM_ORDERING_SWITCH_LEVEL],
                (long)  iparm[IPARM_ORDERING_CMIN],
                (long)  iparm[IPARM_ORDERING_CMAX],
                ((float)iparm[IPARM_ORDERING_FRAT])/100,
                (long)  iparm[IPARM_ORDERING_SWITCH_LEVEL],
                (long)  iparm[IPARM_ORDERING_CMIN],
                (long)  iparm[IPARM_ORDERING_CMAX],
                ((float)iparm[IPARM_ORDERING_FRAT])/100);
        if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
            pastix_print(procnum, 0, "Scotch personal strategy |%s|\n", strat);
    }

    ret = SCOTCH_stratGraphOrder (&stratdat, strat);
    if (ret == 0) {
        /* Compute graph ordering */
        ret = SCOTCH_graphOrderList(&scotchgraph,
                                    (SCOTCH_Num)   n,
                                    (SCOTCH_Num *) NULL,
                                    &stratdat,
                                    (SCOTCH_Num *) ordemesh->permtab,
                                    (SCOTCH_Num *) ordemesh->peritab,
                                    (SCOTCH_Num *)&ordemesh->cblknbr,
                                    (SCOTCH_Num *) ordemesh->rangtab,
                                    NULL);

    }

    SCOTCH_stratExit (&stratdat);
    SCOTCH_graphExit( &scotchgraph );

    if (ret != 0) {           /* If something failed in Scotch */
        orderExit (ordemesh);    /* Free ordering arrays          */
        return PASTIX_ERR_INTERNAL;
    }

#if defined(FORGET_PARTITION)
    ordemesh->cblknbr = 0;
    if (ordemesh->rangtab != NULL) memFree_null(ordemesh->rangtab);
#else
    /* Redimensionnement de rangtab a cblknbr */
    ordemesh->rangtab =
        (pastix_int_t *) memRealloc (ordemesh->rangtab,
                                     (ordemesh->cblknbr + 1)*sizeof (pastix_int_t));
#endif

    return PASTIX_SUCCESS;
}
