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
 * graphBase - Rebase the graph to the given value.
 *
 *******************************************************************************
 *
 * @param[in,out] graph
 *          The graph to rebase.
 *
 * @param[in] baseval
 *          The base value to use in the graph.
 *
 *******************************************************************************/
void graphBase( pastix_graph_t *graph,
                int             baseval )
{
    pastix_int_t baseadj;
    pastix_int_t i, n, nnz;

    baseadj = baseval - graph->colptr[0];
    if (baseadj == 0)
	return;

    n   = graph->n;
    nnz = graph->colptr[n] - graph->colptr[0];

    if (graph->colptr != NULL) {
        for (i = 0; i <= n; i++) {
            graph->colptr[i]   += baseadj;
        }
    }
    if (graph->rows != NULL) {
        for (i = 0; i < nnz; i++)
            graph->rows[i] += baseadj;
    }

    if (graph->loc2glob != NULL) {
        for (i = 0; i < n; i++) {
            graph->loc2glob[i] += baseadj;
        }
    }
    return;
}

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

    if (graph->colptr   != NULL) { memFree_null( graph->colptr);  }
    if (graph->rows     != NULL) { memFree_null( graph->rows);    }
    if (graph->loc2glob != NULL) { memFree_null( graph->loc2glob);}

    memFree_null( graph );
    return;
}
