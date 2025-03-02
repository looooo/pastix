/**
 *
 * @file pastix_subtask_order.c
 *
 * PaStiX ordering task.
 * Contains wrappers to build a good ordering for sparse direct solvers.
 * Affected by the compilation time options:
 *    - PASTIX_ORDERING_SCOTCH: Enable Scotch graph partitioning library.
 *    - PASTIX_ORDERING_PTSCOTCH: Enable PT-Scotch graph partitioning library.
 *    - PASTIX_ORDERING_METIS: Enable Metis graph partitioning library.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Gregoire Pichon
 * @author Tony Delarue
 * @date 2024-07-05
 *
 **/
#include "common.h"
#include <spm.h>
#include "graph/graph.h"
#include "blend/elimintree.h"
#include "order_internal.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_analyze
 *
 * @brief Computes the ordering of the given graph in parameters.
 *
 * The graph given by the user is used to generate a graph that can be used by
 * ordering tools and symbolic factorization. This graph is stored in the
 * pastix_data->graph to be pass over to the symbolic factorization. If it exists
 * before to call this routine, then the current structure is cleaned and a new
 * one is created. From this structure, an ordering is computed by the ordering
 * tool chosen by IPARM_ORDERING and stored in pastix_data->ordemesh. At the end
 * the full ordering stucture: pemutation, inverse permutation, partition, and
 * partion tree is generated such that it can be used by any symbolic
 * factorization algorithm.
 * The user can get back the permutation generated by providing allocated ordering
 * structure where the results is stored after computation.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_ORDERING, IPARM_IO_STRATEGY
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the field ordemesh is initialize with the result of the
 *          ordering.
 *          - IPARM_ORDERING will determine which ordering tool is used.
 *          - IPARM_IO_STRATEGY:
 *             - If set to PastixIOSave, the results will be written on files on
 *               exit.
 *             - If set to PastixIOLoad and IPARM_ORDERING is set to personal,
 *               then the ordering is loaded from files and no ordering is
 *               called.
 *          - If the function pastix_setSchurUnknownList() function has been
 *          previously called to set the list of vertices to isolate in the
 *          schur complement, those vertices are isolated at the end of the
 *          matrix in a dedicated supernode..
 *          - If the function pastix_setZerosUnknownList() function has been
 *          previously called to set the list of diagonal elements that may
 *          cause problem during the factorization, those vertices are isolated
 *          at the end of the matrix in a dedicated supernode..
 *
 * @param[in] spm
 *          The sparse matrix given by the user on which the ordering will be
 *          computed.
 *
 *
 * @param[inout] myorder
 *          On entry, the permutation provide by the user if IPARM_ORDERING
 *          parameter set to PastixOrderPersonal. Not read otherwise.
 *          On exit, if the structure attributs != NULL and IPARM_ORDERING parameter
 *          is not set to PastixOrderPersonal, contains the permutation generated.
 *          Otherwise, it is not referenced.
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
pastix_subtask_order(       pastix_data_t  *pastix_data,
                      const spmatrix_t     *spm,
                            pastix_order_t *myorder )
{
    pastix_int_t    schur_gN;
    pastix_int_t   *schur_perm = NULL;
    pastix_int_t   *zeros_perm = NULL;
    pastix_int_t   *iparm;
    pastix_graph_t *subgraph;
    pastix_graph_t  schurgraph;
    pastix_graph_t  zerosgraph;
    pastix_graph_t *graph;
    pastix_order_t *ordemesh;
    Clock           timer;
    int             procnum;
    int             retval = PASTIX_SUCCESS;
    int             retval_rcv;
    int             do_schur = 1;
    int             do_zeros = 1;
    int             subgraph_is_a_copy = 0;
    int             spmbase;

    /*
     * Check parameters
     */
    if ( pastix_data == NULL ) {
        pastix_print_error( "pastix_subtask_order: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( spm == NULL ) {
        pastix_print_error( "pastix_subtask_order: wrong spm parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( !(pastix_data->steps & STEP_INIT) ) {
        pastix_print_error( "pastix_subtask_order: pastixInit() has to be called before calling this function" );
        return PASTIX_ERR_BADPARAMETER;
    }

    iparm = pastix_data->iparm;
    /* Backup flttype from the spm into iparm[IPARM_FLOAT] for later use */
    iparm[IPARM_FLOAT] = spm->flttype;

    if ( pastix_data->schur_n > 0 )
    {
        /*
         * If ordering is set to PastixOrderPersonal, we consider that the schur
         * complement is already isolated at the end of the permutation array.
         */
        if ( iparm[IPARM_ORDERING] == PastixOrderPersonal ) {
            do_schur = 0;
        }
    }
    else {
        do_schur = 0;
    }

    if ( pastix_data->zeros_n > 0 )
    {
        /*
         * If ordering is set to PastixOrderPersonal, we consider that the zeros
         * on diagonal are already isolated at the end of the permutation array.
         */
        if ( iparm[IPARM_ORDERING] == PastixOrderPersonal ) {
            do_zeros = 0;
        }
    }
    else {
        do_zeros = 0;
    }

    /*
     * Clean ordering if it exists
     */
    if ( pastix_data->ordemesh != NULL ) {
        pastixOrderExit(pastix_data->ordemesh);
    }
    else {
        MALLOC_INTERN( pastix_data->ordemesh, 1, pastix_order_t );
    }

    ordemesh = pastix_data->ordemesh;
    procnum  = pastix_data->procnum;
    pastixOrderAlloc( ordemesh, 0, 0 );

    if ( iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
        pastix_print( procnum, 0, "%s", OUT_STEP_ORDER );
    }

    /*
     * Prepare a copy of user's SPM
     * Copy the given spm in pastix_data structure and performs basic required
     * operations as symmetrizing the graph and removing the diagonal
     * coefficients
     */
    if ( pastix_data->graph != NULL ) {
        graphExit( pastix_data->graph );
        memFree_null( pastix_data->graph );
    }
    graphPrepare( pastix_data, spm, &(pastix_data->graph) );
    graph   = pastix_data->graph;
    spmbase = spmFindBase( spm );

    graphBase( graph, 0 );
    subgraph = graph;

    /*
     * Isolate Shur elements
     */
    if ( do_schur )
    {
        assert( pastix_data->schur_list != NULL );

        if ( spmbase != 0 ) {
            /* We need to rebase the schur unknown list */
            pastix_int_t i;

            for( i=0; i<pastix_data->schur_n; i++ ) {
                pastix_data->schur_list[i] -= spmbase;
            }
        }

        graphIsolate( graph, &schurgraph,
                      pastix_data->schur_n,
                      pastix_data->schur_list,
                      &schur_perm,
                      NULL );

        subgraph_is_a_copy = 1;
        subgraph = &schurgraph;
    }

    /* Backup only the needed part */
    schur_gN = subgraph->gN;

    /*
     * Isolate diagonal elements close to 0.
     */
    if ( do_zeros )
    {
        assert( pastix_data->zeros_list != NULL );

        if ( spmbase != 0 ) {
            /* We need to rebase the zeros unknown list */
            pastix_int_t i;

            for( i=0; i<pastix_data->zeros_n; i++ ) {
                pastix_data->zeros_list[i] -= spmbase;
            }
        }

        graphIsolate( subgraph, &zerosgraph,
                      pastix_data->zeros_n,
                      pastix_data->zeros_list,
                      &zeros_perm,
                      NULL );

        if ( subgraph_is_a_copy ) {
            /* We do not longer require arrays from the schur subgraph */
            graphExit( subgraph );
        }
        subgraph_is_a_copy = 1;
        subgraph = &zerosgraph;
    }

    if ( iparm[IPARM_VERBOSE] > PastixVerboseYes ) {
        pastix_print( procnum, 0, "%s", OUT_ORDER_INIT );
    }

    clockStart(timer);

    /* Select the ordering method chosen by the user */
    switch ( iparm[IPARM_ORDERING] ) {
        /*
         * Scotch Ordering
         */
    case PastixOrderScotch:
        if ( iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
            pastix_print( procnum, 0, OUT_ORDER_METHOD, "Scotch" );
        }
#if defined(PASTIX_ORDERING_SCOTCH)
        graphGatherInPlace( subgraph );
        retval = orderComputeScotch( pastix_data, subgraph );
#else
        pastix_print_error( "pastix_subtask_order: Ordering with Scotch requires to enable -DPASTIX_ORDERING_SCOTCH option" );
        retval = PASTIX_ERR_BADPARAMETER;
#endif
        break;

        /*
         * PT-Scotch Ordering
         */
    case PastixOrderPtScotch:
        if (iparm[IPARM_VERBOSE] > PastixVerboseNot) {
            pastix_print(procnum, 0, OUT_ORDER_METHOD, "PT-Scotch" );
        }
#if defined(PASTIX_ORDERING_PTSCOTCH)
        graphScatterInPlace( subgraph, pastix_data->pastix_comm );
        retval = orderComputePTScotch( pastix_data, subgraph );
#else
        pastix_print_error( "pastix_subtask_order: Ordering with PT-Scotch requires to enable -DPASTIX_ORDERING_PTSCOTCH option" );
        retval = PASTIX_ERR_BADPARAMETER;
#endif
        break;

        /*
         *  METIS ordering
         */
    case PastixOrderMetis:
        if (iparm[IPARM_VERBOSE] > PastixVerboseNot) {
            pastix_print(procnum, 0, OUT_ORDER_METHOD, "Metis" );
        }
#if defined(PASTIX_ORDERING_METIS)
        graphGatherInPlace( subgraph );
        retval = orderComputeMetis( pastix_data, subgraph );
        assert( ordemesh->rangtab == NULL );
#else
        pastix_print_error( "pastix_subtask_order: Ordering with Metis requires -DPASTIX_ORDERING_METIS flag at compile time" );
        retval = PASTIX_ERR_BADPARAMETER;
#endif
        break;

        /*
         * Personal Ordering
         */
    case PastixOrderPersonal:
        if ( iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
            pastix_print( procnum, 0, OUT_ORDER_METHOD, "Personal" );
        }
        retval = orderComputePersonal( pastix_data, subgraph, myorder );
        break;

    default:
        pastix_print_error( "pastix_subtask_order: Ordering not available (iparm[IPARM_ORDERING]=%d)\n",
                    (int)iparm[IPARM_ORDERING] );
        retval = PASTIX_ERR_BADPARAMETER;
        break;
    }

    /*
     * Reduce the error code
     */
    MPI_Allreduce( &retval, &retval_rcv, 1, MPI_INT, MPI_MAX,
                   pastix_data->pastix_comm );
    if (retval_rcv != PASTIX_SUCCESS) {

        /* Cleanup memory */
        if ( do_zeros ) {
            graphExit( &zerosgraph );
            memFree_null( zeros_perm );
        }
        if ( do_schur ) {
            graphExit( &schurgraph );
            memFree_null( schur_perm );
        }

        return retval_rcv;
    }

    /* Rebase the ordering to 0 (for orderFindSupernodes) */
    pastixOrderBase( ordemesh, 0 );

    /*
     * If the rangtab or the treetab are not initialized, let's find it ourself
     */
    if (( ordemesh->rangtab == NULL ) ||
        ( ordemesh->treetab == NULL ) )
    {
        /* The ordering step may be distributed, but supernodes routines are not */
        graphGatherInPlace( subgraph );

        /* TODO: if rangtab is provided, treetab could be easily calculated */
        orderFindSupernodes( subgraph, ordemesh );

#if !defined(NDEBUG) && defined(PASTIX_DEBUG_ORDERING)
        assert( pastixOrderCheck( ordemesh ) == PASTIX_SUCCESS );
#endif

        orderAmalgamate( iparm[IPARM_VERBOSE],
                               iparm[IPARM_INCOMPLETE],
                               iparm[IPARM_LEVEL_OF_FILL],
                               iparm[IPARM_AMALGAMATION_LVLCBLK],
                               iparm[IPARM_AMALGAMATION_LVLBLAS],
                               subgraph,
                               ordemesh,
                               pastix_data->pastix_comm );
    }

#if !defined(NDEBUG) && defined(PASTIX_DEBUG_ORDERING)
    assert( pastixOrderCheck( ordemesh ) == PASTIX_SUCCESS );
#endif

    /*
     * The subgraph is no longer needed, let's free it if it was a copy
     */
    if ( subgraph_is_a_copy ) {
        graphExit( subgraph );
    }

    /*
     * Reorder supernodes by level to get a better order for runtime systems,
     * and for the reordering algorithm
     */
    orderApplyLevelOrder( ordemesh,
                          iparm[IPARM_TASKS2D_LEVEL],
                          iparm[IPARM_TASKS2D_WIDTH] );

#if !defined(NDEBUG) && defined(PASTIX_DEBUG_ORDERING)
    assert( pastixOrderCheck( ordemesh ) == PASTIX_SUCCESS );
#endif

    /*
     * Add the isolated elements to the ordering structure
     */
    if ( do_zeros )
    {
        orderAddIsolate( ordemesh, schur_gN, zeros_perm );

        if ( zeros_perm != NULL ) { memFree_null( zeros_perm ); }
    }

    /*
     * Add the isolated elements to the ordering structure
     */
    if ( do_schur )
    {
        orderAddIsolate( ordemesh, graph->gN, schur_perm );

        if ( schur_perm != NULL ) { memFree_null( schur_perm ); }
    }

    /*
     * Backup of the original supernodes
     */
    ordemesh->sndenbr = ordemesh->cblknbr;
    ordemesh->sndetab = malloc( (ordemesh->sndenbr+1) * sizeof(pastix_int_t) );
    memcpy( ordemesh->sndetab, ordemesh->rangtab, (ordemesh->sndenbr+1) * sizeof(pastix_int_t) );

    /*
     * Block Low-Rank clustering
     */
    if ( ( iparm[IPARM_COMPRESS_WHEN] != PastixCompressNever ) &&
         ( iparm[IPARM_SPLITTING_STRATEGY] != PastixSplitNot ) )
    {
#if !defined(PASTIX_ORDERING_SCOTCH)
        pastix_print_warning( "Clustering is not available yet when Scotch is disabled" );
#else
        EliminTree  *etree;
        pastix_int_t min_cblk = ordemesh->rangtab[ordemesh->cblknbr-1];
        pastix_int_t ret;

        graphBase( graph, 0 );

        etree = orderBuildEtree( ordemesh );

        ret = orderSupernodes( graph, ordemesh,
                               etree, iparm, do_schur );

        eTreeExit( etree );

        (void)ret;
        (void)min_cblk;
#endif /* !defined(PASTIX_ORDERING_SCOTCH) */
    }

#if defined(PASTIX_ORDER_DRAW_LASTSEP)
    /*
     * Draw last separator graph and xyz
     */
    orderDraw( pastix_data, NULL, ordemesh->sndenbr-1,
               orderDrawGraph | orderDrawCoordinates );
#endif

    /* Reduce the error code */
    MPI_Allreduce( &retval, &retval_rcv, 1, MPI_INT, MPI_MAX,
                   pastix_data->pastix_comm );
    if ( retval_rcv != PASTIX_SUCCESS ) {
        return retval_rcv;
    }

    clockStop(timer);
    pastix_data->dparm[DPARM_ORDER_TIME] = clockVal(timer);
    if ( iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
        pastix_print(procnum, 0, OUT_ORDER_TIME, clockVal(timer));
    }

    /*
     * Save i/o strategy
     */
    if ( iparm[IPARM_IO_STRATEGY] & PastixIOSave ) {
        if ( procnum == 0 ) {
            retval = pastixOrderSave( pastix_data, ordemesh );
        }

        MPI_Allreduce( &retval, &retval_rcv, 1, MPI_INT, MPI_MAX,
                       pastix_data->pastix_comm );
        if ( retval_rcv != PASTIX_SUCCESS ) {
            return retval_rcv;
        }
    }

    /*
     * Return the ordering to user if structure is not NULL
     * Remark: No need to copy back for personal
     */
    if ( ( iparm[IPARM_ORDERING] != PastixOrderPersonal ) &&
         ( myorder != NULL ) )
    {
        if ( graph->loc2glob == NULL )
        {
            retval = pastixOrderCopy( myorder, ordemesh );
            MPI_Allreduce( &retval, &retval_rcv, 1, MPI_INT, MPI_MAX,
                           pastix_data->pastix_comm );
            if ( retval_rcv != PASTIX_SUCCESS ) {
                return retval_rcv;
            }
        }
        else
        {
            pastix_int_t *loc2glob;
            pastix_int_t  baseval = graph->baseval;
            pastix_int_t  i, n;

            n = graph->n;
            if ( myorder->permtab != NULL ) {
                pastix_int_t *permtab = ordemesh->permtab - baseval;

                loc2glob = graph->loc2glob;
                for( i = 0; i < n; i++, loc2glob++ ) {
                    myorder->permtab[i] = permtab[*loc2glob];
                }
            }
            if (myorder->peritab != NULL) {
                pastix_int_t *peritab = ordemesh->peritab - baseval;
                pastix_int_t i;

                loc2glob = graph->loc2glob;
                for( i = 0; i < n; i++, loc2glob++ ) {
                    myorder->peritab[i] = peritab[*loc2glob];
                }
            }
            /* TODO: Copy also rangtab and treetab ? */
        }
    }

    /*
     * For now, what the rank 0 has will overwrite what the others have, even if
     * they all computed something
     */
    pastixOrderBcast( pastix_data->ordemesh, 0, pastix_data->pastix_comm );

#if !defined(NDEBUG)
    assert( pastixOrderCheck( pastix_data->ordemesh ) == PASTIX_SUCCESS );
#endif

    /* Backup the spm pointer for further information */
    pastix_data->csc = spm;

    /* Invalidate following steps, and add order step to the ones performed */
    pastix_data->steps &= ~( STEP_SYMBFACT  |
                             STEP_ANALYSE   |
                             STEP_CSC2BCSC  |
                             STEP_BCSC2CTAB |
                             STEP_NUMFACT   |
                             STEP_SOLVE     |
                             STEP_REFINE    );
    pastix_data->steps |= STEP_ORDERING;

    return PASTIX_SUCCESS;
}
