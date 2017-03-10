/**
 *
 * @file graph_apply_perm.c
 *
 *  PaStiX graph apply permutation routine
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Pascal Henon
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
 * @brief Generate the graph of P*A from the graph of A and the
 * permutation vector P.
 *
 *******************************************************************************
 *
 * @param[in] graphA
 *          The original graph.
 *
 * @param[in] perm
 *          Array of size n.
 *          The permutation to apply.
 *
 * @param[inout] graphPA
 *          The initialized graph structure where the permuted graph will be
 *          stored.
 *          The allocated data must be freed with graphClean.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on success,
 * @retval PASTIX_ERR_ALLOC if allocation went wrong,
 * @retval PASTIX_ERR_BADPARAMETER if incorrect parameters are given.
 *
 *******************************************************************************/
int
graphApplyPerm( const pastix_graph_t *graphA,
                const pastix_int_t   *perm,
                      pastix_graph_t *graphPA )
{
    pastix_int_t *rowsPA, *rowsA;
    pastix_int_t i, j, ip;
    pastix_int_t nnz;
    pastix_int_t baseval;
    pastix_int_t n = graphPA->n = graphA->n;
    graphPA->gN  = graphA->gN;
    graphPA->dof = graphA->dof;
    baseval = graphA->colptr[0];

    MALLOC_INTERN( graphPA->colptr, n+1, pastix_int_t );
    MALLOC_INTERN( graphPA->nnz,    n,   pastix_int_t );

    /* Compute the number of nnz per vertex */
    for(i=0; i<n; i++)
    {
        ip = perm[i];
        graphPA->nnz[ip] = graphA->colptr[i+1] - graphA->colptr[i] + 1;
    }

    /* Create the colptr */
    graphPA->colptr[0] = baseval;
    for(i=0; i<n; i++)
    {
        graphPA->colptr[i+1] = graphPA->colptr[i] + graphPA->nnz[i];
    }

    /* Create the row vector */
    nnz = graphPA->colptr[n] - baseval;
    MALLOC_INTERN( (graphPA->rows), nnz, pastix_int_t );
    for(i=0; i<n; i++)
    {
        ip = perm[i];
        rowsPA = graphPA->rows + graphPA->colptr[ip] - baseval;
        rowsA  = graphA->rows  + graphA->colptr[i]   - baseval;

        /* Add the missing diagonal */
        *rowsPA = ip;
        rowsPA++;
        for(j=1; j<graphPA->nnz[ip]; j++, rowsPA++) {
            *rowsPA = perm[ *rowsA ];
            rowsA++;
        }

        intSort1asc1( graphPA->rows + graphPA->colptr[ip] - baseval,
                      graphPA->nnz[ip]);
    }
    return EXIT_SUCCESS;
}
