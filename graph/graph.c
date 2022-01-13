/**
 *
 * @file graph.c
 *
 * PaStiX graph structure routines
 *
 * @copyright 2004-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2021-01-25
 *
 * @addtogroup pastix_graph
 * @{
 *
 **/
#include "common.h"
#include "graph/graph.h"

#define assert_graph( _graph_ )                         \
    do {                                                \
        assert( (_graph_)->fmttype == SpmCSC     );     \
        assert( (_graph_)->flttype == SpmPattern );     \
        assert( (_graph_)->values  == NULL       );     \
    } while (0)

/**
 *******************************************************************************
 *
 * @brief Initialize an empty graph
 *
 *******************************************************************************
 *
 * @param[inout] graph
 *          The empty graph to init.
 *
 * @param[in] comm
 *          The MPI communicator used for the graph.
 *
 *******************************************************************************/
void
graphInitEmpty( pastix_graph_t *graph,
                PASTIX_Comm     comm )
{
    spmInitDist( graph, comm );

    graph->flttype = SpmPattern;

    graph->colptr = malloc( sizeof( pastix_int_t ) );
    graph->colptr[0] = 0;

    graphUpdateComputedFields( graph );
}

/**
 *******************************************************************************
 *
 * @brief Free the content of the graph structure.
 *
 *******************************************************************************
 *
 * @param[inout] graph
 *          The pointer graph structure to free.
 *
 *******************************************************************************/
void
graphExit( pastix_graph_t *graph )
{
    /* Parameter checks */
    if ( graph == NULL ) {
        pastix_print_error( "graphExit: graph pointer is NULL" );
        return;
    }
    assert_graph( graph );

    spmExit( graph );

    return;
}

/**
 *******************************************************************************
 *
 * @brief Rebase the graph to the given value.
 *
 *******************************************************************************
 *
 * @param[inout] graph
 *          The graph to rebase.
 *
 * @param[in] baseval
 *          The base value to use in the graph (0 or 1).
 *
 *******************************************************************************/
void
graphBase( pastix_graph_t *graph,
           pastix_int_t    baseval )
{
    /* Parameter checks */
    if ( graph == NULL ) {
        pastix_print_error( "graphBase: graph pointer is NULL" );
        return;
    }
    if ( (baseval != 0) &&
         (baseval != 1) )
    {
        pastix_print_error( "graphBase: baseval is incorrect, must be 0 or 1" );
        return;
    }

    assert_graph( graph );

    spmBase( graph, baseval );

    return;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief This routine copy a given ordering in a new one.
 *
 * This function copies a graph structure into another one. If all subpointers
 * are NULL, then they are all allocated and contains the original graphsrc
 * values on exit. If one or more array pointers are not NULL, then, only those
 * are copied to the graphdst structure.
 *
 *******************************************************************************
 *
 * @param[output] graphdst
 *          The destination graph. Must be allocated on entry.
 *
 * @param[in] graphsrc
 *          The source graph
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
graphCopy( pastix_graph_t       *graphdst,
           const pastix_graph_t *graphsrc )
{
    /* Parameter checks */
    if ( graphdst == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( graphsrc == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( graphsrc == graphdst ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    assert_graph( graphsrc );

    /* Copy the source graph */
    spmCopy( graphsrc, graphdst );

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief This routine scatter a graph from node root to the other nodes
 *
 *******************************************************************************
 *
 * @param[inout] graph
 *          On entry, the graph to scatter.
 *          On exit, the scattered graph
 *
 * @param[in] comm
 *          MPI communicator.
 *
 *******************************************************************************
 *
 * @retval 1 if the graph has been scattered, 0 if untouched
 *
 *******************************************************************************/
int
graphScatterInPlace( pastix_graph_t *graph,
                     PASTIX_Comm     comm )
{
    pastix_graph_t newgraph;
    int            rc;

    assert_graph( graph );
    if ( graph->loc2glob != NULL ) {
        return 0;
    }

    /* Scatter the graph */
    rc = spmScatter( &newgraph, -1, graph, -1, NULL, 1, comm );

    if ( rc == SPM_SUCCESS ) {
        graphExit( graph );
        memcpy( graph, &newgraph, sizeof( pastix_graph_t ) );

        assert_graph( graph );
        return 1;
    }
    else {
        return 0;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief This routine gather a distributed graph on each note in place.
 *
 *******************************************************************************
 *
 * @param[in] graph
 *          On entry, the distributed graph.
 *          On exit, the gathered graph
 *
 * @param[in] root
 *          The root where we want to gather the graph
 *
 *******************************************************************************
 *
 * @retval 1 if the graph has been gathered, 0 if untouched
 *
 ********************************************************************************/
int
graphGatherInPlace( pastix_graph_t *graph )
{
    pastix_graph_t newgraph;
    int            rc;

    assert_graph( graph );
    if ( graph->loc2glob == NULL ) {
        return 0;
    }

    rc = spmGather( graph, -1, &newgraph );

    if ( rc == SPM_SUCCESS ) {
        graphExit( graph );
        memcpy( graph, &newgraph, sizeof( pastix_graph_t ) );

        assert_graph( graph );
        return 1;
    }
    else {
        return 0;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief This routine sortes the subarray of edges of each vertex.
 *
 * WARNING: The sort is always performed, do not call this routine
 * when it is not required.
 *
 *******************************************************************************
 *
 * @param[inout] graph
 *          On entry, the pointer to the graph structure.
 *          On exit, the same graph with subarrays of edges sorted by ascending
 *          order.
 *
 *******************************************************************************/
void
graphSort( pastix_graph_t *graph )
{
    assert_graph( graph );
    spmSort( graph );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief Symmetrize a given graph
 *
 *******************************************************************************
 *
 * @param[inout] graph
 *          The initialized graph structure that will be symmetrized.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on success,
 * @retval PASTIX_ERR_BADPARAMETER if incorrect parameters are given.
 *
 *******************************************************************************/
int
graphSymmetrize( pastix_graph_t *graph )
{
    /* Parameter checks */
    if ( graph == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    assert_graph( graph );

    spmSymmetrize( graph );
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief Update dofs, nnz, nnzexp, gnnz, n, nexp, gN of a given graph
 *
 *******************************************************************************
 *
 * @param[inout] graph
 *          The initialized graph structure that will be updated.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on success,
 * @retval PASTIX_ERR_BADPARAMETER if incorrect parameters are given.
 *
 *******************************************************************************/
int
graphUpdateComputedFields( pastix_graph_t *graph )
{
    /* Parameter checks */
    if ( graph == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    assert_graph( graph );

    spmUpdateComputedFields( graph );
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief This routine build a graph thanks to an spm;
 *
 * This function copies an spm structure into a graph one. If all subpointers
 * are NULL, then they are all allocated and contains the original spm
 * values on exit. If one or more array pointers are not NULL, then, only those
 * are copied to the graphdst structure.
 * We will take care that our graph does not contain coefficients, therefore has
 * SpmPattern floating type and is a SpmCSC format type.
 *
 *******************************************************************************
 *
 * @param[inout] graph
 *          The destination graph.
 *
 * @param[in] graphsrc
 *          The source Sparse Matrix.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
graphSpm2Graph( pastix_graph_t   *graph,
                const spmatrix_t *spm )
{
    /* Parameter checks */
    if ( graph == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( spm == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    /*
     * Clear the prexisting graph
     * Might be uninitialized, so we call spmExit instead of graphExit.
     */
    spmExit( graph );

    /* Copy existing datas */
    spmCopy( spm, graph );

    /* Enforce Pattern type to the graph */
    graph->flttype = SpmPattern;
    if ( graph->values ) {
        free( graph->values );
        graph->values = NULL;
    }

    /* Make sure the graph is in CSC format */
    spmConvert( SpmCSC, graph );

    return PASTIX_SUCCESS;
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief Build the vertex weight array out of the dof array.
 *
 *******************************************************************************
 *
 * @param[in] graph
 *          Pointer to the graph structure.
 *
 *******************************************************************************
 *
 * @retval The vertex weight array if graph->dof != 1, NULL otherwise.
 *
 *******************************************************************************/
pastix_int_t *
graphGetWeights( const pastix_graph_t *graph )
{
    pastix_int_t  i, n;
    pastix_int_t *weights, *wptr;

    if ( graph->dof == 1 ) {
        return NULL;
    }

    n = graph->n;
    MALLOC_INTERN( weights, n, pastix_int_t );

    wptr = weights;
    /* Constant dof */
    if ( graph->dof > 1 ) {
        for ( i=0; i<n; i++, wptr++ ) {
            *wptr = graph->dof;
        }
    }
    /* Variadic dof */
    else {
        pastix_int_t *dofptr = graph->dofs;
        if ( graph->loc2glob == NULL ) {
            for ( i=0; i<n; i++, wptr++, dofptr++ ) {
                *wptr = dofptr[1] - dofptr[0];
            }
        }
        else {
            pastix_int_t *dofptr   = graph->dofs - graph->baseval;
            pastix_int_t *loc2glob = graph->loc2glob;

            for ( i=0; i<n; i++, wptr++, loc2glob++ ) {
                *wptr = dofptr[ *loc2glob + 1 ] - dofptr[ *loc2glob ];
            }
        }
    }

    return weights;
}

/**
 * @}
 */
