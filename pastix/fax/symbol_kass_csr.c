/**
 *
 * @file kass_genPA.c
 *
 *  PaStiX symbolic factorization routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "graph.h"
#include "kass.h"

void
kass_csrInit( pastix_int_t n, kass_csr_t *csr )
{
    csr->n = n;
    MALLOC_INTERN( csr->nnz,  n, pastix_int_t  );
    MALLOC_INTERN( csr->rows, n, pastix_int_t* );

    memset( csr->nnz,  0, n*sizeof(pastix_int_t ) );
    memset( csr->rows, 0, n*sizeof(pastix_int_t*) );
}


void
kass_csrClean( kass_csr_t *csr )
{
    pastix_int_t i;
    for(i=0; i< csr->n; i++) {
        if ( csr->nnz[i] != 0 ) {
            memFree_null( csr->rows[i] );
        }
    }
    memFree_null( csr->rows );
    memFree_null( csr->nnz  );
}

pastix_int_t
kass_csrGetNNZ( kass_csr_t *csr )
{
    pastix_int_t i, nnz;
    nnz  = 0;
    for(i=0; i< csr->n; i++) {
        nnz += csr->nnz[i];
    }
    return nnz;
}

void
kass_csrCompact( kass_csr_t *csr )
{
    pastix_int_t n = csr->n;
    pastix_int_t i, j;

    /* Look for first deleted node */
    for(i=0,j=0; i<n; i++,j++) {
        if (csr->nnz[i] == 0)
            break;
    }

    /* Compact the nodes */
    for(; i<n; i++) {
        if ( csr->nnz[i] > 0 ) {
            assert( j < i );
            csr->nnz[j]  = csr->nnz[i];
            csr->rows[j] = csr->rows[i];
            csr->nnz[i]  = 0;
            csr->rows[i] = NULL;
            j++;
        }
    }

    csr->n = j;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * graphApplyPerm - Generate the graph of P*A from the graph of A and the
 * permutation vector.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of vertex of the original graph.
 *
 * @param[in] ia
 *          Array of size n+1
 *          Index of first edge for each vertex in ja array.
 *
 * @param[in] ja
 *          Array of size nnz = ia[n] - ia[0].
 *          Edges for each vertex.
 *
 * @param[in] loc2glob
 *          Array of size n
 *          Global numbering of each local vertex.
 *
 * @param[in,out] newgraph
 *          The initialized graph structure where the symmetrized graph will be
 *          stored.
 *          The allocated data must be freed with graphClean.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on success.
 *          \retval PASTIX_ERR_ALLOC if allocation went wrong.
 *          \retval PASTIX_ERR_BADPARAMETER if incorrect parameters are given.
 *
 *******************************************************************************/
int
kass_csrGenPA( const pastix_graph_t *graphA,
               const pastix_int_t   *perm,
                     kass_csr_t     *graphPA )
{
    pastix_int_t *rowsPA, *rowsA;
    pastix_int_t i, j, ip;
    pastix_int_t baseval;
    pastix_int_t n = graphPA->n = graphA->n;

    baseval = graphA->colptr[0];

    /* Compute the number of nnz per vertex */
    for(i=0; i<n; i++)
    {
        ip = perm[i];
        /* Add diagonal (could be removed fro direct) */
        graphPA->nnz[ip] = graphA->colptr[i+1] - graphA->colptr[i] + 1;
    }

    /* Create the row vector */
    for(i=0; i<n; i++)
    {
        ip = perm[i] - baseval;

        MALLOC_INTERN( graphPA->rows[ip],
                       graphPA->nnz[ip],
                       pastix_int_t );

        rowsPA = graphPA->rows[ip];
        rowsA  = graphA->rows + graphA->colptr[i] - baseval;

        /* Add diagonal */
        *rowsPA = ip;
        rowsPA++;

        for(j=1; j<graphPA->nnz[ip]; j++, rowsPA++) {
            *rowsPA = perm[ *rowsA ];
            rowsA++;
        }

        intSort1asc1( graphPA->rows[ip],
                      graphPA->nnz[ip]);
    }
    return EXIT_SUCCESS;
}
