/**
 *
 * @file fax_csr_iluk.c
 *
 * This file contains croutines to create the graph of the incomplete
 * factorization of a graph A.
 *
 * @copyright 2004-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2024-07-05
 *
 **/
#include "common.h"
#include "pastix/order.h"
#include "fax_csr.h"

/**
 *******************************************************************************
 *
 * @ingroup symbol_dev_csr
 *
 * @brief Compute the non zero pattern of the levelized
 * incomplete factor for a sparse lower triangular matrix in CSC format. This
 * pattern is exact iff the matrix has a SYMMETRIC non zero structure.
 *
 * @remark This algorithm has been implemented according to the paper of David
 *     Hysom and Alex Pothen : Level-based Incomplete LU factorization: Graph
 *     Model and Algorithm
 *
 *******************************************************************************
 *
 * @param[in] graphA
 *          The graph structure of the original matrix A to factorize.
 *
 * @param[in] order
 *          The order structure holding the number of cblk and the partition
 *
 * @param[in] level
 *          It is the level desired for the ilu(k) factorization. (level >= 0)
 *
 * @param[inout] graphL
 *          The graph the structure of the non zero pattern of the factorized matrix.
 *          On entry, a pointer to a graph structure. No need for initialization.
 *          On exit, the structure contains the computed graph.
 *
 *******************************************************************************
 *
 * @retval >=0, the number of non zero entries in the generated graph.
 * @retval -i, if the i^th parameter is incorrect
 *
 *******************************************************************************/
pastix_int_t
faxCSRFactILUk( const fax_csr_t      *graphA,
                const pastix_order_t *order,
                pastix_int_t          level,
                fax_csr_t            *graphL )
{
    pastix_int_t *visited = NULL;
    pastix_int_t *length  = NULL;
    pastix_int_t *stack   = NULL;
    pastix_int_t *adj     = NULL;
    pastix_int_t *ja      = NULL;
    pastix_int_t  used;
    pastix_int_t  h, i, j, k, t;
    long          nnz;
    fax_csr_t     tmpgrah;

    /* Check parameters */
    if ( graphA == NULL ) {
        return -1;
    }
    if ( level < 0 ) {
        return -2;
    }
    if ( graphL == NULL ) {
        return -3;
    }

    /* Quick return */
    if ( graphA->n == 0 ) {
        return 0;
    }

    /** Allocated the working array **/
    MALLOC_INTERN( visited, graphA->n, pastix_int_t );
    MALLOC_INTERN( length,  graphA->n, pastix_int_t );
    MALLOC_INTERN( stack,   graphA->n, pastix_int_t );
    MALLOC_INTERN( ja,      graphA->n, pastix_int_t );
    nnz = 0;

    /** Initialized visited ***/
    for ( j = 0; j < graphA->n; j++ ) {
        visited[j] = -1;
        length[j]  = 0;
    }

    /** Apply GS_Urow for each row **/
    faxCSRInit( graphA->n, graphL );
    for ( i = 0; i < graphA->n; i++ ) {
        /** Reset the stack number of elements **/
        stack[0] = i;
        used     = 1;

        length[i]  = 0;
        visited[i] = i;

        ja[0] = i; /** Put the diagonal term **/
        k     = 1;

        /** BFS phase **/
        while ( used > 0 ) {
            used--;
            h   = stack[used];
            adj = graphA->rows[h];
            for ( j = 0; j < graphA->nnz[h]; j++ ) {
                t = adj[j];
                if ( visited[t] != i ) {
                    visited[t] = i;
                    if ( ( t < i ) && ( length[h] < level ) ) {
                        stack[used] = t;
                        used++;
                        length[t] = length[h] + 1;
                    }
                    if ( t > i ) {
                        ja[k++] = t;
                    }
                }
            }
        }

        assert( k > 0 );

        graphL->nnz[i] = k;
        MALLOC_INTERN( graphL->rows[i], k, pastix_int_t );
        memcpy( graphL->rows[i], ja, k * sizeof( pastix_int_t ) );

        intSort1asc1( graphL->rows[i], graphL->nnz[i] );

        nnz += k;
    }

    /* Compress with the elementary partition */
    tmpgrah = *graphL;
    faxCSRCblkCompress( &tmpgrah, order, graphL, visited );
    faxCSRClean( &tmpgrah );

    memFree_null( ja );
    memFree_null( visited );
    memFree_null( length );
    memFree_null( stack );

    graphL->total_nnz = nnz;
    return nnz;
}
