/**
 *
 * @file order_compute_scotch.c
 *
 * PaStiX order driver to perform ordering with Scotch library.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "graph.h"
#include "order.h"
#if defined(PASTIX_ORDERING_PTSCOTCH)
#include <ptscotch.h>
#elif defined(PASTIX_ORDERING_SCOTCH)
#include <scotch.h>
#endif /* defined(PASTIX_ORDERING_PTSCOTCH) */
#include "order_scotch_strats.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Compute the ordering of the graph given as parameter
 * with Scotch library.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_ORDERING_DEFAULT, IPARM_SCOTCH_SWITCH_LEVEL,
 *   IPARM_SCOTCH_CMIN, IPARM_SCOTCH_CMAX, IPARM_SCOTCH_FRAT
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the field oerdemesh is initialize with the result of the
 *          ordering realized by Scotch.
 *
 * @param[inout] graph
 *          The graph prepared by graphPrepare function on which wwe want to
 *          perform the ordering. On exit, the graph might be rebased.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed,
 * @retval PASTIX_ERR_INTEGER_TYPE if Scotch integer type is not the
 *         same size as PaStiX ones,
 * @retval PASTIX_ERR_INTERNAL if an error occurs internally to Scotch.
 *
 *******************************************************************************/
int
pastixOrderComputeScotch( pastix_data_t  *pastix_data,
                          pastix_graph_t *graph )
{
    pastix_order_t *ordemesh = pastix_data->ordemesh;
    SCOTCH_Graph    scotchgraph;
    SCOTCH_Strat    stratdat;
    char            strat[1024];
    pastix_int_t   *colptr;
    pastix_int_t   *rows;
    pastix_int_t   *iparm = pastix_data->iparm;
    pastix_int_t    procnum;
    pastix_int_t    n, nnz, baseval;
    int ret;

    procnum = pastix_data->procnum;

    /* Check integer compatibility */
    if (sizeof(pastix_int_t) != sizeof(SCOTCH_Num)) {
        errorPrint("Inconsistent integer type\n");
        return PASTIX_ERR_INTEGER_TYPE;
    }

     /* Enable this define to fix the SCOTCH random generator */
#if defined(SCOTCH_FIX_SEED)
    _SCOTCHintRandInit();
#endif

    /* Centralized */
#if 0
    if (iparm[IPARM_GRAPHDIST] == 0) {
#endif
        n      = graph->n;
        colptr = graph->colptr;
        rows   = graph->rows;
        baseval= colptr[0];
        nnz    = colptr[n] - baseval;
#if 0
    }
    /* Distributed */
    else {
        cscd2csc_int( graph->n,
                      graph->colptr,
                      graph->rows,
                      NULL, NULL, NULL, NULL,
                      &n, &colptr, &rows,
                      NULL, NULL, NULL, NULL,
                      graph->loc2glob,
                      pastix_data->pastix_comm,
                      0, /* DoF to 0 as we have no values */
                      1);
    }
#endif
    print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_graphInit <\n");
    pastixOrderAlloc(ordemesh, n, n);
    SCOTCH_graphInit( &scotchgraph );

    print_debug(DBG_ORDER_SCOTCH, "> SCOTCH_graphBuild <\n");
    if (SCOTCH_graphBuild(&scotchgraph,   /* Graph to build     */
                          baseval,        /* baseval            */
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
        if (iparm[IPARM_VERBOSE] > PastixVerboseYes)
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
    if (iparm[IPARM_ORDERING_DEFAULT] == 1) {
        if (iparm[IPARM_INCOMPLETE] == 0) {
            if (iparm[IPARM_VERBOSE] > PastixVerboseNo)
                pastix_print(procnum, 0, "      Scotch direct strategy\n");
            sprintf(strat, SCOTCH_STRAT_DIRECT);
        }
        else {
            if (iparm[IPARM_VERBOSE] > PastixVerboseNo)
                pastix_print(procnum, 0, "      Scotch incomplete strategy\n");
            sprintf(strat, SCOTCH_STRAT_INCOMP);
        }
    }
    /* personal ordering */
    else {
        sprintf(strat, SCOTCH_STRAT_PERSO,
                (long)  iparm[IPARM_SCOTCH_SWITCH_LEVEL],
                (long)  iparm[IPARM_SCOTCH_CMIN],
                (long)  iparm[IPARM_SCOTCH_CMAX],
                ((float)iparm[IPARM_SCOTCH_FRAT])/100,
                (long)  iparm[IPARM_SCOTCH_SWITCH_LEVEL],
                (long)  iparm[IPARM_SCOTCH_CMIN],
                (long)  iparm[IPARM_SCOTCH_CMAX],
                ((float)iparm[IPARM_SCOTCH_FRAT])/100);
        if (iparm[IPARM_VERBOSE] > PastixVerboseNo)
            pastix_print(procnum, 0, "Scotch personal strategy |%s|\n", strat);
    }

    {
        static volatile pastix_atomic_lock_t strat_lock = PASTIX_ATOMIC_UNLOCKED;
        pastix_atomic_lock( &strat_lock );
        ret = SCOTCH_stratGraphOrder (&stratdat, strat);
        pastix_atomic_unlock( &strat_lock );
    }
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
                                    (SCOTCH_Num *) ordemesh->treetab);
    }

    SCOTCH_stratExit (&stratdat);
    SCOTCH_graphExit( &scotchgraph );
#if 0
    if (iparm[IPARM_GRAPHDIST] == 1) {
        memFree_null(colptr);
        memFree_null(rows);
    }
#endif

    if (ret != 0) {           /* If something failed in Scotch */
        pastixOrderExit (ordemesh);    /* Free ordering arrays          */
        return PASTIX_ERR_INTERNAL;
    }

#if defined(FORGET_PARTITION)
    ordemesh->cblknbr = 0;
    if (ordemesh->rangtab != NULL) memFree_null(ordemesh->rangtab);
    if (ordemesh->treetab != NULL) memFree_null(ordemesh->treetab);
#else
    /**
     * Adapt size of rangtab and treetab to the new cblknbr
     * WARNING: If no nodes in the graph, nothing has been initialized.
     */
    ordemesh->rangtab =
        (pastix_int_t *) memRealloc (ordemesh->rangtab,
                                     (ordemesh->cblknbr + 1)*sizeof (pastix_int_t));
    ordemesh->treetab =
        (pastix_int_t *) memRealloc (ordemesh->treetab,
                                     (ordemesh->cblknbr)*sizeof (pastix_int_t));
    if (ordemesh->cblknbr == 0) {
        ordemesh->rangtab[0] = 0;
    }
#endif

    return PASTIX_SUCCESS;
}
