/**
 *
 * @file graph_base.c
 *
 *  PaStiX graph routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
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

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * graphClean - Free the graph structure given in parameter.
 *
 *******************************************************************************
 *
 * @param[in,out] graph
 *          The graph structure to free.
 *
 *******************************************************************************/
void graphClean( pastix_graph_t *graph )
{
    graph->gN = 0;
    graph->n  = 0;

    /* Parameter checks */
    if ( graph == NULL ) {
        errorPrint("graphClean: graph pointer is NULL");
        return;
    }
    if ( (graph->colptr == NULL) ||
         (graph->rows   == NULL) )
    {
        errorPrint("graphClean: graph pointer is not correctly initialized");
        return;
    }

    memFree_null( graph->colptr);
    memFree_null( graph->rows);
    if (graph->loc2glob != NULL) { memFree_null( graph->loc2glob);}

    memFree_null( graph );
    return;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * graphBase - Rebase the graph to the given value.
 *
 *******************************************************************************
 *
 * @param[in,out] graph
 *          The graph to rebase.
 *
 * @param[in] baseval
 *          The base value to use in the graph (0 or 1).
 *
 *******************************************************************************/
void graphBase( pastix_graph_t *graph,
                int             baseval )
{
    pastix_int_t baseadj;
    pastix_int_t i, n, nnz;

    /* Parameter checks */
    if ( graph == NULL ) {
        errorPrint("graphBase: graph pointer is NULL");
        return;
    }
    if ( (graph->colptr == NULL) ||
         (graph->rows   == NULL) )
    {
        errorPrint("graphBase: graph pointer is not correctly initialized");
        return;
    }
    if ( (baseval != 0) &&
         (baseval != 1) )
    {
        errorPrint("graphBase: baseval is incorrect, must be 0 or 1");
        return;
    }

    baseadj = baseval - graph->colptr[0];
    if (baseadj == 0)
	return;

    n   = graph->n;
    nnz = graph->colptr[n] - graph->colptr[0];

    for (i = 0; i <= n; i++) {
        graph->colptr[i]   += baseadj;
    }
    for (i = 0; i < nnz; i++) {
        graph->rows[i] += baseadj;
    }

    if (graph->loc2glob != NULL) {
        for (i = 0; i < n; i++) {
            graph->loc2glob[i] += baseadj;
        }
    }
    return;
}
