/**
 *
 * @file graph_isolate.c
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

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * graphIsolate - This routine isolate a subset of vertices from a given graph,
 * and return a new GRAPH cleaned from those vertices.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of columns of the original GRAPH matrix.
 *
 * @param[in] colptr
 *          Array of size n+1
 *          Index of first element of each column in rows array.
 *
 * @param[in] rows
 *          Array of size nnz = colptr[n] - colptr[0].
 *          Rows of each non zero entries.
 *
 * @param[in] isolate_n
 *          The number of columns to isolate from the original GRAPH matrix.
 *
 * @param[in,out] isolate_list
 *          Array of size isolate_n.
 *          List of columns to isolate. On exit, the list is sorted by ascending
 *          indexes.
 *
 * @param[out] new_colptr
 *          Array of size n-isolate_n+1
 *          Index of first element of each column in rows array for the new GRAPH
 *          matrix.
 *          If new_colptr == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 * @param[out] new_rows
 *          Array of size new_nnz = (*new_colptr)[n] - (*new_colptr)[0].
 *          Rows of each non zero entries for the new GRAPH matrix.
 *          If new_rows == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 * @param[out] new_perm
 *          Array of size n-isolate_n.
 *          Contains permutation generated to isolate the columns at the end of
 *          the GRAPH.
 *          If new_perm == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 * @param[out] new_invp
 *          Array of size n-isolate_n.
 *          Contains the inverse permutation generated to isolate the columns
 *          at the end of the GRAPH.
 *          If new_invp == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on success.
 *          \retval PASTIX_ERR_ALLOC if allocation went wrong.
 *          \retval PASTIX_ERR_BADPARAMETER if incorrect parameters are given.
 *
 *******************************************************************************/
int graphIsolate(       pastix_int_t   n,
                  const pastix_int_t  *colptr,
                  const pastix_int_t  *rows,
                        pastix_int_t   isolate_n,
                        pastix_int_t  *isolate_list,
                        pastix_int_t **new_colptr,
                        pastix_int_t **new_rows,
                        pastix_int_t **new_perm,
                        pastix_int_t **new_invp )
{
    pastix_int_t *tmpcolptr = NULL;
    pastix_int_t *tmprows   = NULL;
    pastix_int_t *tmpperm   = NULL;
    pastix_int_t *tmpinvp   = NULL;
    pastix_int_t  baseval = colptr[0];
    pastix_int_t  nnz = colptr[n] - baseval;
    pastix_int_t  new_n = n - isolate_n;
    pastix_int_t  new_nnz;
    pastix_int_t  i, j, ip, k;
    pastix_int_t  iter_isolate = 0;
    pastix_int_t  iter_non_isolate  = 0;

    if (isolate_n > n) {
        errorPrintW( "Number of columns to isolate greater than the columns in the GRAPH matrix\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Quick Return */
    if (isolate_n == 0) {
        if (new_colptr != NULL) *new_colptr = (pastix_int_t*)colptr;
        if (new_rows   != NULL) *new_rows   = (pastix_int_t*)rows;
        return PASTIX_SUCCESS;
    }

    if (isolate_n == n) {
        if (new_colptr != NULL) {
            MALLOC_INTERN(*new_colptr, n, pastix_int_t);
            memcpy( *new_colptr, colptr, n*sizeof(pastix_int_t) );
        }
        if (new_rows != NULL) {
            MALLOC_INTERN(*new_rows, nnz, pastix_int_t);
            memcpy( *new_rows, rows, nnz*sizeof(pastix_int_t) );
        }
        return PASTIX_SUCCESS;
    }

    /* Sort the lost of vertices */
    intSort1asc1(isolate_list, isolate_n);

    /* Init invp array */
    MALLOC_INTERN(tmpinvp, n, pastix_int_t);
    for (i = 0; i <n; i++) {
        if (i == isolate_list[iter_isolate]-baseval)
        {
            tmpinvp[new_n+iter_isolate] = i;
            iter_isolate++;
        }
        else
        {
            tmpinvp[iter_non_isolate] = i;
            iter_non_isolate++;
        }
    }

    assert(iter_non_isolate == new_n    );
    assert(iter_isolate     == isolate_n);

    /* Init perm array */
    MALLOC_INTERN(tmpperm, n, pastix_int_t);
    for(i = 0; i < n; i++)
        tmpperm[tmpinvp[i]] = i;

#if defined(PASTIX_DEBUG_GRAPH)
    for(i = 0; i < n; i++)
    {
        assert(tmpperm[i] < n );
        assert(tmpperm[i] > -1);
    }
#endif

    /* Create the new_colptr array */
    MALLOC_INTERN(tmpcolptr, new_n + 1, pastix_int_t);
    memset(tmpcolptr, 0, (new_n + 1)*sizeof(pastix_int_t));

    tmpcolptr[0] = baseval;
    for (i=0; i<n; i++)
    {
        ip = tmpperm[i];
        if (ip < new_n)
        {
            for (j = colptr[i]-baseval; j < colptr[i+1]-baseval; j++)
            {
                /* Count edges in each column of the new graph */
                if (tmpperm[rows[j]-baseval] < new_n)
                {
                    tmpcolptr[ip+1]++;
                }
            }
        }
    }

    for (i = 0; i < new_n; i++)
        tmpcolptr[i+1] += tmpcolptr[i];

    new_nnz = tmpcolptr[new_n] - tmpcolptr[0];
    /* TODO: be careful, allocation will fail if matrix is diagonal and no off-diagonal elements are found */

    /* Create the new rows array */
    MALLOC_INTERN(tmprows, new_nnz, pastix_int_t);
    for (i = 0; i <n; i++)
    {
        ip = tmpperm[i];
        if (ip < new_n)
        {
            k = tmpcolptr[ip]-baseval;
            for (j = colptr[i]-baseval; j < colptr[i+1]-baseval; j ++)
            {
                /* Count edges in each column of the new graph */
                if (tmpperm[rows[j]-baseval] < new_n)
                {
                    tmprows[k] = tmpperm[rows[j]-baseval] + baseval;
                    k++;
                }
            }
            assert( k == tmpcolptr[ip+1]-baseval );
        }
    }

    if (new_colptr != NULL) {
        *new_colptr = tmpcolptr;
    } else {
        memFree_null( tmpcolptr );
    }
    if (new_rows != NULL) {
        *new_rows = tmprows;
    } else {
        memFree_null( tmprows );
    }
    if (new_perm != NULL) {
        *new_perm = tmpperm;
    } else {
        memFree_null( tmpperm );
    }
    if (new_invp != NULL) {
        *new_invp = tmpinvp;
    } else {
        memFree_null( tmpinvp );
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * graphIsolateSupernode - This routine isolates a subset of vertices belonging
 * to a supernode from a given graph, and returns a new graph made up of those
 * vertices and internal connexions.  Extra edges are created between vertices
 * if they are connected through a halo at 1.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of columns of the original matrix, ie the number of
 *          vertices of the graph.
 *
 * @param[in] colptr
 *          Array of size n+1
 *          Index of first element of each column in rows array, ie the first
 *          element of each vertex in the edges array.
 *
 * @param[in] rows
 *          Array of size nnz = colptr[n] - colptr[0].
 *          Rows of each non zero entries, ie the list of neighbours of each entry.
 *
 * @param[in] perm
 *          Array of size n.
 *          Contains the permutation generated by the ordering of the initial matrix.
 *
 * @param[in] invp
 *          Array of size n.
 *          Contains the inverse permutation generated by the ordering of the initial matrix.
 *
 * @param[in] isolate_n
 *          The number of vertices to isolate from the original graph.
 *
 * @param[in,out] isolate_list
 *          Array of size isolate_n.
 *          List of columns/vertices to isolate. On exit, the list is sorted by
 *          ascending indexes.
 *
 * @param[out] new_colptr
 *          Array of size isolate_n+1
 *          Index of first element of each column in rows array for the new GRAPH
 *          matrix.
 *          If new_colptr == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 * @param[out] new_rows
 *          Array of size new_nnz = (*new_colptr)[n] - (*new_colptr)[0].
 *          Rows of each non zero entries for the new GRAPH matrix.
 *          If new_rows == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on success.
 *          \retval PASTIX_ERR_ALLOC if allocation went wrong.
 *          \retval PASTIX_ERR_BADPARAMETER if incorrect parameters are given.
 *
 *******************************************************************************/
int graphIsolateSupernode(       pastix_int_t   n,
			   const pastix_int_t  *colptr,
			   const pastix_int_t  *rows,
			   const pastix_int_t  *perm,
			   const pastix_int_t  *invp,
                                 pastix_int_t   fnode,
                                 pastix_int_t   lnode,
				 pastix_int_t   isolate_n,
			   const pastix_int_t  *isolate_list,
				 pastix_int_t **new_colptr,
				 pastix_int_t **new_rows )
{
    pastix_int_t *tmpcolptr = NULL;
    pastix_int_t *tmprows   = NULL;
    pastix_int_t *tmpperm   = NULL;
    pastix_int_t *tmpinvp   = NULL;
    pastix_int_t  baseval = colptr[0];
    pastix_int_t  nnz = colptr[n] - baseval;
    pastix_int_t  other_n = n - isolate_n;
    pastix_int_t  new_nnz;
    pastix_int_t  i, j, ip, k;
    pastix_int_t  iter_isolate = 0;
    pastix_int_t  iter_non_isolate  = 0;

    if (isolate_n > n) {
        errorPrintW( "Number of columns to isolate greater than the columns in the GRAPH matrix\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Quick Return */
    if (isolate_n == 0) {
        if (new_colptr != NULL) *new_colptr = (pastix_int_t*)colptr;
        if (new_rows   != NULL) *new_rows   = (pastix_int_t*)rows;
        return PASTIX_SUCCESS;
    }

    if (isolate_n == n) {
        if (new_colptr != NULL) {
            MALLOC_INTERN(*new_colptr, n, pastix_int_t);
            memcpy( *new_colptr, colptr, n*sizeof(pastix_int_t) );
        }
        if (new_rows != NULL) {
            MALLOC_INTERN(*new_rows, nnz, pastix_int_t);
            memcpy( *new_rows, rows, nnz*sizeof(pastix_int_t) );
        }
        return PASTIX_SUCCESS;
    }

    /* Create the new_colptr array */
    MALLOC_INTERN(tmpcolptr, isolate_n + 1, pastix_int_t);
    memset(tmpcolptr, 0, (isolate_n + 1)*sizeof(pastix_int_t));

    pastix_int_t *sn_connected;
    /* temporary array of connections to avoid double counting when extending */
    MALLOC_INTERN(sn_connected, isolate_n, pastix_int_t);

    /* NOTE: seems repetitive to go over the data twice. Could construct individual rows for each column
             in temporary arrays, then combine them all in the end when we know exactly how large rows
             should be */

    /* (i,j) in permutated ordering */
    /* (ip,jp) in initial ordering */
    tmpcolptr[0] = baseval;
    for (ip=fnode; ip<=lnode; ip++)
    {
        memset(sn_connected, 0, (isolate_n)*sizeof(pastix_int_t));
        sn_connected[ip-fnode] = 1;

        /* i^th vertex in the initial numbering */
        i = invp[ip];
        for (j = colptr[i]-baseval; j < colptr[i+1]-baseval; j++)
        {
            pastix_int_t jp = perm[ rows[j]-baseval ];

            /* Count edges in each column of the new graph */
            if ( jp >= fnode && jp <= lnode )
            {
                assert(sn_connected[jp-fnode] == 0);
                sn_connected[jp-fnode] = 1;
                tmpcolptr[ip-fnode+1]++;
            }
            else
            { /* Look for connection at distance 1 */
                for(k = colptr[rows[j]-baseval]-baseval; k < colptr[rows[j]+1-baseval]-baseval; k++)
                {
                    pastix_int_t kp = perm[ rows[k]-baseval ];
                    if ( kp >= fnode && kp <= lnode )
                    {
                        if(!sn_connected[kp-fnode])
                        {
                            sn_connected[kp-fnode] = 1;
                            tmpcolptr[ip-fnode+1]++;
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i < isolate_n; i++)
        tmpcolptr[i+1] += tmpcolptr[i];

    new_nnz = tmpcolptr[isolate_n] - tmpcolptr[0];
    /* TODO: be careful, allocation will fail if matrix is diagonal and no off-diagonal elements are found */

    /* Create the new rows array */
    MALLOC_INTERN(tmprows, new_nnz, pastix_int_t);
    pastix_int_t row_counter = 0;

    for (ip=fnode; ip<=lnode; ip++)
    {
        memset(sn_connected, 0, (isolate_n)*sizeof(pastix_int_t));
        sn_connected[ip-fnode] = 1;

        /* i^th vertex in the initial numbering */
        i = invp[ip];
        for (j = colptr[i]-baseval; j < colptr[i+1]-baseval; j++)
        {
            pastix_int_t jp = perm[ rows[j]-baseval ];

            /* Count edges in each column of the new graph */
            if ( jp >= fnode && jp <= lnode )
            {
                assert(sn_connected[jp-fnode] == 0);
                sn_connected[jp-fnode] = 1;
                tmprows[row_counter] = jp-fnode;
                row_counter++;
            }
            else
            { /* Look for connection at distance 1 */
                for(k = colptr[rows[j]-baseval]-baseval; k < colptr[rows[j]+1-baseval]-baseval; k++)
                {
                    pastix_int_t kp = perm[ rows[k]-baseval ];
                    if ( kp >= fnode && kp <= lnode )
                    {
                        if(!sn_connected[kp-fnode])
                        {
                            sn_connected[kp-fnode] = 1;
                            tmprows[row_counter] = kp-fnode;
                            row_counter++;
                        }
                    }
                }
            }
        }
    }

    if (new_colptr != NULL) {
        *new_colptr = tmpcolptr;
    } else {
        memFree_null( tmpcolptr );
    }
    if (new_rows != NULL) {
        *new_rows = tmprows;
    } else {
        memFree_null( tmprows );
    }

    return PASTIX_SUCCESS;
}

