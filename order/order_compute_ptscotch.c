/**
 *
 * @file order_compute_ptscotch.c
 *
 * PaStiX order driver to perform ordering with PT-Scotch library.
 *
 * @copyright 2004-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2021-01-25
 *
 **/
#include "common.h"
#include "graph/graph.h"
#include "pastix/order.h"
#if defined(PASTIX_ORDERING_PTSCOTCH)
#include <ptscotch.h>
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Check the PT-Scotch distributed graph.
 *
 *******************************************************************************
 *
 * @param[in] graph
 *          Pointer to the SCOTCH_Dgraph structure.
 *
 * @param[in] procnum
 *          Procnum of the process. Output purpose.
 *
 *******************************************************************************/
static inline void
ocpts_graph_check( const SCOTCH_Dgraph *graph,
                   int                  procnum )
{
#if defined(PASTIX_DEBUG_ORDERING)
    Clock timer;
    clockStart(timer);
    if ( SCOTCH_dgraphCheck( graph ) ) {
        errorPrint( "pastix: dgraphCheck" );
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }
    clockStop(timer);
    pastix_print( procnum, 0, "SCOTCH_dgraphCheck done in %lf second\n", clockVal(timer) );
#else
    (void)graph;
    (void)procnum;
#endif /* defined(PASTIX_DEBUG_ORDERING) */
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Build the PT-Scotch distributed graph.
 *
 *******************************************************************************
 *
 * @param[inout] scotchgraph
 *          The SCOTCH_Dgraph structure that will be build.
 *
 * @param[inout] graph
 *          The graph prepared by graphPrepare function on which we want to
 *          perform the ordering.
 *
 * @param[in] comm
 *          PaStiX communicator.
 *
 *******************************************************************************/
static inline void
ocpts_graph_init( SCOTCH_Dgraph  *scotchgraph,
                  pastix_graph_t *graph,
                  PASTIX_Comm     comm )
{
    pastix_int_t *colptr, *rowptr;
    pastix_int_t *weights;
    pastix_int_t  n, nnz, baseval;

    /* Make sure graph is 0-based */
    graphBase( graph, 0 );

    n       = graph->n;
    nnz     = graph->nnz;
    colptr  = graph->colptr;
    rowptr  = graph->rowptr;
    baseval = graph->baseval;
    weights = NULL;

    assert( baseval == colptr[0] );
    assert( nnz == colptr[n] - colptr[0] );

    SCOTCH_dgraphInit( scotchgraph, comm );

    /*
     * Generate the vertex load array if dof != 1
     */
    weights = graphGetWeights( graph );

    if ( SCOTCH_dgraphBuild( scotchgraph,
                             baseval,      /* baseval */
                             n,            /* Number of local vertices             */
                             n,            /* Maximum number of local vertices     */
                             colptr,
                             NULL,
                             weights,      /* Local vertex load array (if any)     */
                             NULL,         /* Local vertex label array (if any)    */
                             nnz,
                             nnz,
                             rowptr,       /* Local edge array                     */
                             NULL,         /* Ghost edge array (if any); not const */
                             NULL ) )
    {
        errorPrint( "pastix : SCOTCH_dgraphBuild\n" );
        EXIT(MOD_SOPALIN,PASTIX_ERR_INTERNAL);
    }

    /* Check the generated Scotch graph structure */
    ocpts_graph_check( scotchgraph, graph->clustnum );

    return;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Cleanup the PT-Scotch graph structure and free the associated data.
 *
 *******************************************************************************
 *
 * @param[inout] scotchgraph
 *          The SCOTCH_Dgraph structure that will be clean.
 *
 * @param[inout] comm
 *          The PaStiX communicator.
 *
 *******************************************************************************/
static inline void
ocpts_graph_exit( SCOTCH_Dgraph *scotchgraph,
                  PASTIX_Comm    comm )
{
    pastix_int_t *colptr  = NULL;
    pastix_int_t *rowptr  = NULL;
    pastix_int_t *weights = NULL;

    SCOTCH_dgraphData (
        scotchgraph,
        NULL,     /* Base value                            */
        NULL,     /* Global number of vertices             */
        NULL,     /* Local number of vertices              */
        NULL,     /* Max number of local vertices          */
        NULL,     /* Number of local and ghost vertices    */
        &colptr,  /* Vertex array [vertnbr+1]              */
        NULL,     /* Vertex array [vertnbr]                */
        &weights, /* Vertex load array                     */
        NULL,     /* Vertex label array                    */
        NULL,     /* Global number of edges (arcs)         */
        NULL,     /* Declared size of the local edge array */
        NULL,     /* Max number of local edges             */
        &rowptr,  /* Edge array [edgenbr]                  */
        NULL,     /* Ghost adjency array                   */
        NULL,     /* Edge load array                       */
        &comm );

    /* Free the vertex load array */
    if ( weights != NULL ) {
        memFree_null( weights );
    }

    return;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Compute the graph ordering
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          Pointer to the pastix_data instance.
 *
 * @param[inout] scotchgraph
 *          The SCOTCH_Dgraph structure that will be ordered.
 *
 ********************************************************************************
 *
 * @retval The return value of SCOTCH_graphDOrderList.
 *
 *******************************************************************************/
static inline int
ocpts_compute_graph_ordering( pastix_data_t  *pastix_data,
                              SCOTCH_Dgraph  *scotchgraph )
{
    pastix_order_t  *ordemesh = pastix_data->ordemesh;
    SCOTCH_Strat     stratdat;
    SCOTCH_Dordering ordedat;
    SCOTCH_Ordering  ordering;

    /* Create Strategy string for Scotch */
    SCOTCH_stratInit( &stratdat );

    /* Init distributed ordering */
    if ( SCOTCH_dgraphOrderInit(scotchgraph, &ordedat) )
    {
        pastix_print_error("pastix : SCOTCH_dgraphOrderInit\n");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    /* Compute distributed ordering */
    if ( SCOTCH_dgraphOrderCompute( scotchgraph, &ordedat, &stratdat ) )
    {
        pastix_print_error( "pastix : SCOTCH_dgraphOrderCompute" );
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    SCOTCH_stratExit( &stratdat );

    /* Init centralized ordering */
    if ( SCOTCH_dgraphCorderInit( scotchgraph,
                                  &ordering,
                                  (SCOTCH_Num *) ordemesh->permtab,
                                  (SCOTCH_Num *) ordemesh->peritab,
                                  (SCOTCH_Num *)&ordemesh->cblknbr,
                                  (SCOTCH_Num *) ordemesh->rangtab,
                                  (SCOTCH_Num *) ordemesh->treetab) )
    {
        pastix_print_error( "pastix : SCOTCH_dgraphCorderInit" );
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
    }

    /* Gather distributed ordering on node 0 */
    if ( pastix_data->procnum == 0 ) {
        SCOTCH_dgraphOrderGather( scotchgraph, &ordedat, &ordering );
    }
    else {
        SCOTCH_dgraphOrderGather( scotchgraph, &ordedat, NULL );
    }

    /* Broadcast node 0 datas on all nodes */
    {
        int rangnbr, permnbr;
        MPI_Bcast( &ordemesh->cblknbr, 1, PASTIX_MPI_INT, 0, pastix_data->pastix_comm );

        rangnbr = ordemesh->cblknbr;
        MPI_Bcast( ordemesh->rangtab, rangnbr+1, PASTIX_MPI_INT, 0, pastix_data->pastix_comm );
        MPI_Bcast( ordemesh->treetab, rangnbr,   PASTIX_MPI_INT, 0, pastix_data->pastix_comm );

        permnbr = pastix_data->graph->gN;
        MPI_Bcast( ordemesh->permtab, permnbr, PASTIX_MPI_INT, 0, pastix_data->pastix_comm );
        MPI_Bcast( ordemesh->peritab, permnbr, PASTIX_MPI_INT, 0, pastix_data->pastix_comm );
    }

    /* Exit PT-Scotch ordering structures */
    SCOTCH_dgraphCorderExit( scotchgraph, &ordering );
    SCOTCH_dgraphOrderExit( scotchgraph, &ordedat );

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Compute the ordering of the graph given as parameter
 * with PT-Scotch library.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_ORDERING_DEFAULT, IPARM_SCOTCH_SWITCH_LEVEL,
 *   IPARM_SCOTCH_CMIN, IPARM_SCOTCH_CMAX, IPARM_SCOTCH_FRAT
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the field ordemesh is initialized with the result of the
 *          ordering realized by PT-Scotch.
 *
 * @param[inout] graph
 *          The graph prepared by graphPrepare function on which we want to
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
pastixOrderComputePTScotch( pastix_data_t  *pastix_data,
                            pastix_graph_t *graph )
{
    SCOTCH_Dgraph   scotchdgraph;
    int             ret      = PASTIX_SUCCESS;
    pastix_order_t *ordemesh = pastix_data->ordemesh;
    pastix_int_t    baseval  = graph->baseval;

    /* Check integer compatibility */
    if ( sizeof(pastix_int_t) != sizeof(SCOTCH_Num) ) {
        errorPrint("pastixOrderComputePTScotch: Inconsistent integer type between Pastix and PT-Scotch\n");
        return PASTIX_ERR_INTEGER_TYPE;
    }

    /* Enable this define to fix the SCOTCH random generator */
#if defined(PASTIX_ORDERING_FIX_SEED) && defined(PASTIX_ORDERING_SCOTCH_MT)
    SCOTCH_randomSeed( (SCOTCH_Num)(pastix_data->id) );
#endif

    /* Build The Scotch Dgraph */
    ocpts_graph_init( &scotchdgraph, graph, pastix_data->pastix_comm );

    /* Allocate the ordering structure */
    pastixOrderAlloc( ordemesh, graph->gN, graph->gN );

    /* Compute the ordering */
    ret = ocpts_compute_graph_ordering( pastix_data, &scotchdgraph );

    /*
     * SCOTCH_graphBase may have had side effects on our graph.
     * Call this routine again with our saved baseval.
     */
    graphBase( graph, baseval );

    /* Free the Scotch graph structure */
    ocpts_graph_exit( &scotchdgraph, pastix_data->pastix_comm );

    /* If something has failed in PT-Scotch */
    if ( ret != PASTIX_SUCCESS ) {
        pastixOrderExit( ordemesh );
        return PASTIX_ERR_INTERNAL;
    }
    order_scotch_reallocate_ordemesh( ordemesh );

    return ret;
}
