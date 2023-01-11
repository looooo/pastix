/**
 *
 * @file graph_isolate.c
 *
 * PaStiX graph isolate routine
 *
 * @copyright 2004-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2021-12-21
 *
 **/
#include "common.h"
#include "pastix/order.h"
#include "blend/extendVector.h"
#include "graph/graph.h"

/**
 *******************************************************************************
 *
 * @brief Assign the new pointer to the temporary one
 *
 *******************************************************************************
 *
 * @param[in] newptr
 *          TODO
 *
 * @param[in] tmpptr
 *          TODO
 *
 *******************************************************************************/
static inline void
graph_isolate_assign_ptr( pastix_int_t **newptr,
                          pastix_int_t  *tmpptr )
{
    if ( newptr != NULL ) {
        *newptr = tmpptr;
    } else {
        memFree_null( tmpptr );
    }
}

/**
 *******************************************************************************
 * @brief If isolate_n == n, everything needs to be isolated.
 *        That means that we just need an id permuation.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          TODO
 *
 * @param[in] newperm
 *          TODO
 *
 * @param[in] newinvp
 *          TODO
 *
 *******************************************************************************/
static inline void
graph_isolate_everything( pastix_int_t   n,
                          pastix_int_t **newperm,
                          pastix_int_t **newinvp )
{
    pastix_int_t i;

    if( newperm != NULL ) {
        pastix_int_t *perm;

        MALLOC_INTERN( *newperm, n, pastix_int_t );

        perm = *newperm;
        for (i = 0; i < n; i++, perm++) {
            *perm = i;
        }
    }
    if( newinvp != NULL ) {
        pastix_int_t *invp;

        MALLOC_INTERN( *newinvp, n, pastix_int_t );

        invp = *newinvp;
        for (i = 0; i < n; i++, invp++) {
            *invp = i;
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Init and fill the inverse permutation array.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          TODO
 *
 * @param[in] isolate_n
 *          TODO
 *
 * @param[in] isolate_list
 *          TODO
 *
 * @param[in] permtab
 *          TODO
 *
 * @param[in] invptab
 *          TODO
 *
 * @param[in] baseval
 *          TODO
 *
 *******************************************************************************/
static inline void
graph_isolate_permutations( pastix_int_t        n,
                            pastix_int_t        isolate_n,
                            const pastix_int_t *isolate_list,
                            pastix_int_t       *permtab,
                            pastix_int_t       *invptab,
                            pastix_int_t        baseval )
{
    pastix_int_t *invp_isolate, *invp_nonisol;
    pastix_int_t  i;

    /* First, fill inverse permutations array */
    invp_nonisol = invptab;
    invp_isolate = invptab + (n - isolate_n);
    for (i=0; i<n; i++)
    {
        if ( ((invp_isolate - invptab) <  n) &&
             (i == (*isolate_list - baseval)) )
        {
            *invp_isolate = i;
            invp_isolate++;
            isolate_list++;
        }
        else {
            *invp_nonisol = i;
            invp_nonisol++;
        }
    }
    assert( (invp_nonisol - invptab) == (n - isolate_n) );
    assert( (invp_isolate - invptab) ==  n );

    /* Second, fill permutation array */
    invp_nonisol = invptab;
    for( i = 0; i < n; i++, invp_nonisol++ )
    {
        permtab[ *invp_nonisol ] = i;
    }

#if defined(PASTIX_DEBUG_GRAPH)
    for(i = 0; i < n; i++)
    {
        assert(permtab[i] < n );
        assert(permtab[i] > -1);
    }
#endif
}

/**
 *******************************************************************************
 *
 * @brief Compress the local copy of the graph with only the kept unknowns.
 *
 *******************************************************************************
 *
 * @param[in] oldgraph
 *          TODO
 *
 * @param[in] newgraph
 *          TODO
 *
 * @param[in] new_gn
 *          TODO
 *
 * @param[in] permtab
 *          TODO
 *
 *******************************************************************************/
static inline void
graph_isolate_compress( const pastix_graph_t *oldgraph,
                        pastix_graph_t       *newgraph,
                        pastix_int_t          new_gn,
                        const pastix_int_t   *permtab )
{
    pastix_int_t *newcol  = newgraph->colptr;
    pastix_int_t *oldcol  = oldgraph->colptr;
    pastix_int_t *newrow  = newgraph->rowptr;
    pastix_int_t *oldrow  = oldgraph->rowptr;
    pastix_int_t *newdof  = newgraph->dofs;
    pastix_int_t *olddof  = oldgraph->dofs;
    pastix_int_t *newl2g  = newgraph->loc2glob;
    pastix_int_t *oldl2g  = oldgraph->loc2glob;
    pastix_int_t  baseval = oldgraph->baseval;

    pastix_int_t  i, j, k, ip, jp;
    pastix_int_t  nbrows;
    pastix_int_t  n = oldgraph->n;

    /* Make sure the glob2loc is destroyed */
    if ( newgraph->glob2loc ) {
        free( newgraph->glob2loc );
        newgraph->glob2loc = NULL;
    }

    /* Initial values */
    *newcol = baseval;
    if ( oldgraph->dofs ) {
        *newdof = baseval;
    }

    for ( i=0; i<n; i++, oldcol++, olddof++, oldl2g++ )
    {
        /* Get the new global index after permutation */
        k  = oldgraph->loc2glob ? *oldl2g : i;
        ip = permtab[k];

        /* The vertex is isolated, we remove the full column */
        if( ip >= new_gn ) {
            oldrow += (oldcol[1] - oldcol[0]);
            continue;
        }

        nbrows = 0;
        for( k=oldcol[0]; k<oldcol[1]; k++, oldrow++ )
        {
            j  = *oldrow - baseval;
            jp = permtab[j];

            /* Copy only edges that are kept */
            if( jp < new_gn ) {
                *newrow = jp + baseval;
                newrow++;
                nbrows++;
            }
        }

        /* Update the new colptr */
        newcol[1] = newcol[0] + nbrows;
        newcol++;

        /* Update the new loc2glob if any */
        if ( oldgraph->loc2glob ) {
            *newl2g = ip;
            newl2g++;
        }

        /* Update the new dofs if any */
        if ( oldgraph->dofs ) {
            newdof[1] = newdof[0] + (olddof[1] - olddof[0]);
            newdof++;
        }
    }

    /* Update sizes */
    {
        pastix_int_t new_n, new_nnz;
        new_n   = newcol - newgraph->colptr;
        new_nnz = *newcol - baseval;

        assert( new_n <= new_gn );
        assert( new_nnz <= oldgraph->nnz );

        newgraph->n   = new_n;
        newgraph->nnz = new_nnz;

        graphUpdateComputedFields( newgraph );
    }

    /* Realloc what needs to be */
    newgraph->colptr = realloc( newgraph->colptr, (newgraph->n+1) * sizeof(pastix_int_t) );
    newgraph->rowptr = realloc( newgraph->rowptr, newgraph->nnz * sizeof(pastix_int_t) );
    if ( newgraph->loc2glob ) {
        newgraph->loc2glob = realloc( newgraph->loc2glob, newgraph->n * sizeof(pastix_int_t) );
    }
    if ( newgraph->dofs ) {
        newgraph->dofs = realloc( newgraph->dofs, (newgraph->gN+1) * sizeof(pastix_int_t) );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief Isolate a subset of vertices from a given graph.
 *
 * Return a new graph cleaned from those vertices.
 *
 *******************************************************************************
 *
 * @param[in] ingraph
 *          TODO
 *
 * @param[in] outgraph
 *          TODO
 *
 * @param[in] isolate_n
 *          The number of columns to isolate from the original graph.
 *
 * @param[inout] isolate_list
 *          Array of size isolate_n.
 *          List of columns to isolate. On exit, the list is sorted by ascending
 *          indexes. Must be based as the graph.
 *
 * @param[out] new_perm
 *          Array of size n-isolate_n.
 *          Contains permutation generated to isolate the columns at the end of
 *          the graph that is 0-based.
 *          If new_perm == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 * @param[out] new_invp
 *          Array of size n-isolate_n.
 *          Contains the inverse permutation generated to isolate the columns
 *          at the end of the graph that is 0-based.
 *          If new_invp == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on success,
 * @retval PASTIX_ERR_ALLOC if allocation went wrong,
 * @retval PASTIX_ERR_BADPARAMETER if incorrect parameters are given.
 *
 *******************************************************************************/
int graphIsolate( const pastix_graph_t *ingraph,
                  pastix_graph_t       *outgraph,
                  pastix_int_t          isolate_n,
                  pastix_int_t         *isolate_list,
                  pastix_int_t        **new_perm,
                  pastix_int_t        **new_invp )
{
    pastix_int_t *tmpperm = NULL;
    pastix_int_t *tmpinvp = NULL;
    pastix_int_t  baseval = ingraph->baseval;
    pastix_int_t  gN      = ingraph->gN;
    pastix_int_t  new_n   = gN - isolate_n;

    if ( (isolate_n > gN)  || (isolate_n < 0) ) {
        pastix_print_warning( "Number of columns to isolate greater than the columns in the graph matrix\n" );
        return PASTIX_ERR_BADPARAMETER;
    }

    /* For the moment, make sure the graph is 0 based */
    assert( baseval == 0 );

    /* We isolate the whole graph, so the new graph is empty */
    if ( isolate_n == gN )
    {
        graphInitEmpty( outgraph, ingraph->comm );
        graph_isolate_everything( gN, new_perm, new_invp );
        return PASTIX_SUCCESS;
    }

    graphCopy( outgraph, ingraph );

    /* Quick Return */
    if ( isolate_n == 0 ) {
        pastix_print_warning( "graphIsolate: the graph is beeing duplicated to isolate no unknowns, are you sure you wanted to do that ?\n" );
        return PASTIX_SUCCESS;
    }

    /* Sort the list of vertices */
    intSort1asc1( isolate_list, isolate_n );

    /* Init invp and perm array */
    MALLOC_INTERN( tmpinvp, gN, pastix_int_t );
    MALLOC_INTERN( tmpperm, gN, pastix_int_t );
    graph_isolate_permutations( gN, isolate_n, isolate_list,
                                tmpperm, tmpinvp, baseval );

    /* Create the new graph */
    graph_isolate_compress( ingraph, outgraph, new_n, tmpperm );

    /* Backup the permutation is asked by the caller */
    graph_isolate_assign_ptr( new_perm, tmpperm );
    graph_isolate_assign_ptr( new_invp, tmpinvp );

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @brief Fill the isolated colptr and rowptr.
 *
 *******************************************************************************
 *
 * @param[in] graph
 *          The original graph associated from which vertices and edges must be
 *          extracted.
 *
 * @param[in] order
 *          The ordering structure associated to the graph.
 *
 * @param[in] vec
 *          TODO
 *
 * @param[inout] out_colptr
 *          TODO
 *
 * @param[inout] out_rowptr
 *          TODO
 *
 * @param[in] fnode
 *          The index of the first node to extract in the inverse permutation.
 *
 * @param[in] lnode
 *          The index (+1) of the last node to extract in the inverse permutation.
 *
 * @param[in] distance
 *          Distance considered in number of edges to create an edge in isolated
 *          graph.
 *
 *******************************************************************************/
static inline void
graph_iRange_fill_outptr( const pastix_graph_t *graph,
                          const pastix_order_t *order,
                          ExtendVectorINT      *vec,
                          pastix_int_t         *out_colptr,
                          pastix_int_t         *out_rowptr,
                          pastix_int_t          fnode,
                          pastix_int_t          lnode,
                          pastix_int_t          distance )
{
    pastix_int_t *out_connected;
    pastix_int_t *colptr ,*rowptr, *rowtmp;
    pastix_int_t *perm   = order->permtab;
    pastix_int_t *invp   = order->peritab;
    pastix_int_t  out_n  = lnode - fnode;
    pastix_int_t  ip, jp, k, jj, i, j;
    pastix_int_t  d, sze, baseval;

    MALLOC_INTERN( out_connected, out_n, pastix_int_t );

    /* The first loop counts the number of edges */
    baseval = graph->baseval;
    colptr  = graph->colptr;
    rowptr  = graph->rowptr - baseval;
    rowtmp  = out_rowptr;
    for( ip=fnode; ip<lnode; ip++ )
    {
        /* Clear the previous computations */
        extendint_Clear( vec );
        memset(out_connected, 0, out_n * sizeof(pastix_int_t));
        out_connected[ip-fnode] = 1;

        /* i^th vertex in the initial numbering */
        extendint_Add( vec, invp[ip] );
        sze =  1;
        d   = -1;
        k   =  0;

        while( d < distance ) {
            for(; k<sze; k++) {
                i = extendint_Read( vec, k );

                for (jj = colptr[i]; jj < colptr[i+1]; jj++)
                {
                    j  = rowptr[jj] - baseval;
                    jp = perm[j];

                    /* Count edges in each column of the new graph */
                    if ( ( jp >= fnode ) && ( jp < lnode ) ) {
                        if (out_connected[jp-fnode] == 0){
                            out_connected[jp-fnode] = 1;

                            /* Fill the out_colptr */
                            out_colptr[ip-fnode+1]++;

                            /* Fill the out_rowptr */
                            *rowtmp = jp-fnode;
                            rowtmp++;
                        }
                    }
                    else {
                        extendint_Add( vec, j );
                    }
                }
            }
            d++;
            sze = extendint_Size( vec );
        }
    }
    free(out_connected);
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * @brief Isolate the subgraph associated to a range of unknowns in the permuted
 * graph.
 *
 * This routine isolates a continuous subset of vertices from a given graph, and
 * returns a new graph made of those vertices and internal connexions. Extra
 * edges are created between vertices if they are connected through a halo at a
 * distance d given in parameter.
 *
 *******************************************************************************
 *
 * @param[in] graph
 *          The original graph associated from which vertices and edges must be
 *          extracted.
 *
 * @param[in] order
 *          The ordering structure associated to the graph.
 *
 * @param[inout] out_graph
 *          The extracted graph. If the graph is allocated, it is freed before usage.
 *          On exit, contains the subgraph of the vertices invp[fnode] to invp[lnode-1].
 *
 * @param[in] fnode
 *          The index of the first node to extract in the inverse permutation.
 *
 * @param[in] lnode
 *          The index (+1) of the last node to extract in the inverse permutation.
 *
 * @param[in] distance
 *          Distance considered in number of edges to create an edge in isolated
 *          graph.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on success.
 * @retval PASTIX_ERR_ALLOC if allocation went wrong.
 * @retval PASTIX_ERR_BADPARAMETER if incorrect parameters are given.
 *
 *******************************************************************************/
int
graphIsolateRange( const pastix_graph_t *graph,
                   const pastix_order_t *order,
                   pastix_graph_t       *out_graph,
                   pastix_int_t          fnode,
                   pastix_int_t          lnode,
                   pastix_int_t          distance )
{
    ExtendVectorINT vec;
    pastix_int_t   *out_rowptr;
    pastix_int_t   *out_colptr;
    pastix_int_t    baseval = graph->baseval;
    pastix_int_t    n       = graph->n;
    pastix_int_t    out_n   = lnode - fnode;
    pastix_int_t    i;

    if ( out_graph == NULL ) {
        pastix_print_warning( "graphIsolateSupernode: Incorrect pointer for the output graph\n" );
        return PASTIX_ERR_BADPARAMETER;
    }

    n             = graph->n;
    out_graph->n  = out_n;

    /* Copy dofs datas, as they are global */
    out_graph->dof  = graph->dof;
    if( out_graph->dof < 0 ) {
        MALLOC_INTERN( out_graph->dofs, graph->gN, pastix_int_t );
        memcpy( out_graph->dofs, graph->dofs, graph->gN * sizeof(pastix_int_t) );
    }

    if ( out_n == 0 ) {
        pastix_print_warning( "graphIsolateSupernode: Empty supernode\n" );
        return PASTIX_ERR_BADPARAMETER;
    }

    /* For the moment, make sure the graph is 0 based */
    assert( baseval == 0 );

    /* Quick Return */
    if ( out_n == n ) {
        /* Only one supernode */
        assert( order->cblknbr == 1 );
        graphCopy( out_graph, graph );
        return PASTIX_SUCCESS;
    }

    /* Create the new_colptr array */
    MALLOC_INTERN( out_graph->colptr, out_n + 1, pastix_int_t );
    memset( out_graph->colptr, 0, (out_n + 1) * sizeof(pastix_int_t) );
    out_colptr = out_graph->colptr;

    /* Create the new_rowptr array, will be reallocated later */
    MALLOC_INTERN( out_graph->rowptr, graph->nnz, pastix_int_t );
    out_rowptr = out_graph->rowptr;

    /*
     * (i,j) in permutated ordering
     * (ip,jp) in initial ordering
     */
    out_graph->baseval = baseval;
    out_colptr[0]      = baseval;

    extendint_Init( &vec, 100 );

    /* Fill the out_colptr and the out_rowptr */
    graph_iRange_fill_outptr( graph, order, &vec,
                              out_colptr, out_rowptr,
                              fnode, lnode, distance );

    /* Update the colptr */
    for (i = 0; i < out_n; i++){
        out_colptr[i+1] += out_colptr[i];
    }

    out_graph->nnz = out_colptr[out_n] - out_colptr[0];

    /* Allocation will fail if matrix is diagonal and no off-diagonal elements are found */
    if ( out_graph->nnz == 0 ){
        fprintf( stderr, "Diagonal matrix cannot be correcly managed here!\n" );
        return EXIT_FAILURE;
    }

    /* Create the new rowptr array */
    out_graph->rowptr =
        (pastix_int_t *) memRealloc( out_graph->rowptr, out_graph->nnz * sizeof(pastix_int_t) );

    extendint_Exit( &vec );
    graphBase( out_graph, 0 );

    /* Update mainly gN and gnnz */
    graphUpdateComputedFields( out_graph );

    return PASTIX_SUCCESS;
}
