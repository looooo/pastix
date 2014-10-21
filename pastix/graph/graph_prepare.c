/**
 *
 * @file graph_prepare.c
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
#if defined(PASTIX_DISTRIBUTED)
#include "cscd_utils_intern.h"
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * graphNoDiag - This routine removes the diagonal edges from a centralized
 * graph.
 *
 *******************************************************************************
 *
 * @param[in,out] graph
 *          On entry, the pointer to the graph structure with possible diagonal
 *          edges (i,i).
 *          On exit, all entries of type (i,i) are removed from the graph.
 *
 *******************************************************************************/
void
graphNoDiag( pastix_graph_t *graph )
{
    pastix_int_t  i, j, indj;
    pastix_int_t  n   = graph->n;
    pastix_int_t *ia  = graph->colptr;
    pastix_int_t *ja  = graph->rows;
    pastix_int_t *ja2 = graph->rows;
    int baseval = ia[0];

    indj = ia[0];
    for(i=0; i<n; i++, ia++)
    {
        for (j = ia[0]; j < ia[1]; j++, ja++ )
        {
            /* If diagonal element, we skip it */
            if ( (ja[0]-baseval) == i ) {
                continue;
            }
            /* Otherwise we save it */
            else {
                *ja2 = *ja;
                ja2++;
            }
        }
        /* Save the new ia[i] */
        ia[0] = indj;

        /* Store the new ia[i+1] */
        indj = (ja2 - graph->rows) + baseval;
    }
    ia[0] = indj;

    graph->rows =
        (pastix_int_t *) memRealloc ( graph->rows,
                                      (ia[0]-baseval)*sizeof (pastix_int_t) );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * graphSort - This routine sortes the subarray of edges of each vertex in a
 * centralized graph. WARNING: the sort is always performed, so be carefull to
 * not call this routine when it is not required.
 *
 *******************************************************************************
 *
 * @param[in,out] graph
 *          On entry, the pointer to the graph structure.
 *          On exit, the same graph with subarrays of edges sorted by ascending
 *          order.
 *
 *******************************************************************************/
void
graphSort( pastix_graph_t *graph )
{
    pastix_int_t *ia = graph->colptr;
    pastix_int_t *ja = graph->rows;
    pastix_int_t  n = graph->n;
    pastix_int_t  itercol;
    int baseval = ia[0];

    /* Sort in place each subset */
    for (itercol=0;itercol<n;itercol++)
    {
        pastix_int_t frow = ia[itercol]   - baseval;
        pastix_int_t lrow = ia[itercol+1] - baseval;

        intSort1asc1( (ja+frow), (lrow-frow) );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * graphPrepare - This routine initializes the graph for future call to ordering
 * and symbol matrix generation tools: symmetrize the graph, remove duplicates,
 * ...
 *
 *******************************************************************************
 *
 * @param[in,out] pastix_data
 *          The pointer to the solver instance. On exit, the fields n, cols,
 *          rows and loc2glob are initialized for future steps of the solver.
 *
 * @param[in] n
 *          The number of vertices.
 *
 * @param[in] colptr
 *          Array of size n+1
 *          The array of indirection to the rows array for each vertex.
 *          rows[ colptr[i] ] to rows[ colptr[i+1] are the edges of the
 *          ith vertex.
 *          Can be equal to NULL if graph load is asked.
 *
 * @param[in] rows
 *          Array of size nnz = colptr[n] - colptr[0]. The array of edges.
 *          rows[ colptr[i]   - colptr[0] ] to
 *          rows[ colptr[i+1] - colptr[0] ] are the edges of the ith vertex.
 *          Can be equal to NULL if graph load is asked.
 *
 * @param[in] loc2glob
 *          Array of size n
 *          Global numbering of each local vertex.
 *          Can be equal to NULL if graph load is asked.
 *
 * @param[out] graph
 *          On exit, the pointer to the allocated graph structure is returned.
 *          The graph can then be used with ordering and symbol factorization
 *          tools.
 *          The graph is symmetrized without diagonal elements and rows array is
 *          sorted.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on success.
 *          \retval !0 on failure.
 *
 *******************************************************************************/
int
graphPrepare(      pastix_data_t   *pastix_data,
                   pastix_int_t     n,
             const pastix_int_t    *colptr,
             const pastix_int_t    *rows,
             const pastix_int_t    *loc2glob,
                   pastix_graph_t **graph )
{
    pastix_graph_t *tmpgraph  = NULL;
    pastix_int_t *iparm   = pastix_data->iparm;
    pastix_int_t  procnum = pastix_data->procnum;
    int io_strategy = iparm[IPARM_IO_STRATEGY];

    MALLOC_INTERN( tmpgraph, 1, pastix_graph_t );
    memset( tmpgraph, 0, sizeof(pastix_graph_t) );

    if (PASTIX_MASK_ISTRUE(io_strategy, API_IO_LOAD_GRAPH))
    {
        graphLoad( pastix_data, tmpgraph );
    }
    else
    {
        /* Check that we use Fortran ordering */
        assert( colptr[0] == 1 );

        /*
         * Centralized graph
         */
        //if (iparm[IPARM_GRAPHDIST] == API_NO)
        if (loc2glob == NULL)
        {
            tmpgraph->gN = n;

            /*
             * TODO: change test for requirement from the user to correct his
             * mistakes
             */
            if ( (iparm[IPARM_SYM] == API_SYM_YES) ||
                 (iparm[IPARM_SYM] == API_SYM_HER) )
            {
                graphSymmetrize( n, colptr, rows, loc2glob, tmpgraph );
                assert( n == tmpgraph->n );
                if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO) {
                    fprintf(stdout, "SY - N=%ld, NNZ=%ld\n", n, tmpgraph->colptr[n] - tmpgraph->colptr[0]);
                }
            }
            else
            {
                pastix_int_t nnz = colptr[n]-colptr[0];
                tmpgraph->n = n;
                MALLOC_INTERN(tmpgraph->colptr, (n+1), pastix_int_t);
                MALLOC_INTERN(tmpgraph->rows,   nnz,   pastix_int_t);
                memcpy(tmpgraph->colptr, colptr, (n+1)*sizeof(pastix_int_t));
                memcpy(tmpgraph->rows,   rows,     nnz*sizeof(pastix_int_t));
                if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
                    fprintf(stdout, "GE - N=%ld, NNZ=%ld\n", n, nnz);

                graphSort( tmpgraph );
            }

            if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
                pastix_print(procnum, 0, "%s", OUT_NODIAG);

            graphNoDiag( tmpgraph );
        }
#if defined(PASTIX_DISTRIBUTED)
        /*
         * Distributed graph
         */
        else
        {
            MPI_Comm     pastix_comm = pastix_data->pastix_comm;
            pastix_int_t gN = 0;
            int copy_l2g = 1;

            MPI_Allreduce(&n, &gN, 1, PASTIX_MPI_INT, MPI_SUM, pastix_comm);
            if (iparm[IPARM_SYM]==API_SYM_YES || iparm[IPARM_SYM] == API_SYM_HER) {
                cscd_symgraph_int(n, colptr, rows, NULL,
                                  &(tmpgraph->n),
                                  &(tmpgraph->colptr),
                                  &(tmpgraph->rows), NULL,
                                  loc2glob,
                                  pastix_comm, API_YES );
                assert( n == tmpgraph->n );
            } else {
                pastix_int_t nnz = colptr[n]-colptr[0];
                tmpgraph->n = n;
                MALLOC_INTERN(tmpgraph->colptr, (n+1), pastix_int_t);
                MALLOC_INTERN(tmpgraph->rows,   nnz,   pastix_int_t);
                memcpy(tmpgraph->colptr, colptr, (n+1)*sizeof(pastix_int_t));
                memcpy(tmpgraph->rows,   rows,     nnz*sizeof(pastix_int_t));
            }
            MALLOC_INTERN(tmpgraph->loc2glob,   n,   pastix_int_t);
            memcpy(tmpgraph->loc2glob, loc2glob, n*sizeof(pastix_int_t));

            cscd_noDiag(tmpgraph->n,
                        tmpgraph->colptr,
                        tmpgraph->rows,
                        NULL,
                        loc2glob);

            /* Create contiguous partitions for ordering tools */
            {
                pastix_int_t i;
                int ok  = 0;
                int gok = 0;

                /* Check if matrix is allready partitionned in contiguous blocks */
                for (i = 0; i < n-1; i++)
                    if (loc2glob[i] != (loc2glob[i+1] - 1) )
                        ok = 1;

                MPI_Allreduce(&ok, &gok, 1, MPI_INT, MPI_SUM, pastix_comm);

                /*
                 * If the partition is incorrect, we create a permutation to linearize the sets
                 */
                if ( !gok ) {
                    pastix_int_t  ldisp;
                    pastix_int_t *all_n;
                    pastix_int_t *displs;

                    /* Gather the locals n */
                    MALLOC_INTERN(all_n,  pastix_data->procnbr, pastix_int_t);
                    MALLOC_INTERN(displs, pastix_data->procnbr, pastix_int_t);

                    MPI_Allgather(&n,    1, PASTIX_MPI_INT,
                                  all_n, 1, PASTIX_MPI_INT,
                                  pastix_comm);

                    displs[0] = 0;
                    for (i = 1; i < pastix_data->procnbr; i++)
                        displs[i] = displs[i-1] + all_n[i-1];
                    ldisp = displs[ pastix_data->procnum ] + 1;

                    /* Collect the locals loc2glob */
                    MALLOC_INTERN(pastix_data->PTS_peritab, gN, pastix_int_t);
                    MPI_Allgatherv((void*)loc2glob, n, PASTIX_MPI_INT,
                                   pastix_data->PTS_peritab, all_n, displs, PASTIX_MPI_INT,
                                   pastix_comm);

                    memFree_null(displs);
                    memFree_null(all_n);

                    MALLOC_INTERN(pastix_data->PTS_permtab, gN, pastix_int_t);
                    for (i = 0; i < gN; i++)
                        pastix_data->PTS_permtab[pastix_data->PTS_peritab[i]-1] = i+1;

                    /* Apply the new permutation to the local graph */
                    for (i = 0; i < (tmpgraph->colptr)[n] - 1; i++)
                        tmpgraph->rows[i] = pastix_data->PTS_permtab[(tmpgraph->rows)[i]-1];

                    /* Initialize loc2glob */
                    copy_l2g = 0;
                    MALLOC_INTERN(tmpgraph->loc2glob, n, pastix_int_t);
                    for (i = 0; i < n; i++,ldisp++)
                        tmpgraph->loc2glob[i] = ldisp;
                }

                tmpgraph->gN = gN;
            }

            if (copy_l2g)
            {
                MALLOC_INTERN(tmpgraph->loc2glob, n, pastix_int_t);
                memcpy(tmpgraph->loc2glob, loc2glob, n*sizeof(pastix_int_t));
            }
        }
#else
        assert(loc2glob == NULL );
#endif
    }

    graphBase( tmpgraph, 1 );

    *graph = tmpgraph;
    return PASTIX_SUCCESS;
}
