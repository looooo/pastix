/**
 *
 * @file order_prepare_csc.h
 *
 *  PaStiX order routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "csc_utils.h"
#if defined(PASTIX_DISTRIBUTED)
#include "cscd_utils_intern.h"
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * orderPrepareCSC - This routine prepares the csc for future call to ordering
 * tools: symmetrize the graph, remove duplicates, ...
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
 *
 * @param[in] rows
 *          Array of size nnz = colptr[n] - colptr[0]. The array of edges.
 *          rows[ colptr[i]   - colptr[0] ] to
 *          rows[ colptr[i+1] - colptr[0] ] are the edges of the ith vertex.
 *
 * @param[in] loc2glob
 *          Array of size n
 *          Global numbering of each local vertex.
 *
 *******************************************************************************
 *
 * @return
 *          \retval 0 on success.
 *          \retval !0 on failure.
 *
 *******************************************************************************/
int orderPrepareCSC(pastix_data_t *pastix_data,
                    pastix_int_t   n,
                    const pastix_int_t *colptr,
                    const pastix_int_t *rows,
                    const pastix_int_t *loc2glob)
{
    pastix_int_t *iparm   = pastix_data->iparm;
    pastix_int_t  procnum = pastix_data->procnum;

    /* Check that we use Fortran ordering */
    assert( colptr[0] == 1 );

    /* Clean space already used */
    if (pastix_data->col2      != NULL) { memFree_null(pastix_data->col2); }
    if (pastix_data->row2      != NULL) { memFree_null(pastix_data->row2); }
    if (pastix_data->loc2glob2 != NULL) { memFree_null(pastix_data->loc2glob2); }
    pastix_data->bmalcolrow = 0;

    /*
     * Centralized csc
     */
    if (loc2glob == NULL) {
        if ((iparm[IPARM_SYM] == API_SYM_YES) ||
            (iparm[IPARM_SYM] == API_SYM_HER) )
        {
            csc_symgraph_int(n, colptr, rows, NULL,
                             &(pastix_data->n2),
                             &(pastix_data->col2),
                             &(pastix_data->row2), NULL, API_YES );
            assert( n == pastix_data->n2 );
        }
        else
        {
            pastix_int_t nnz = colptr[n]-colptr[0];
            pastix_data->n2 = n;
            MALLOC_INTERN(pastix_data->col2, (n+1), pastix_int_t);
            MALLOC_INTERN(pastix_data->row2, nnz,   pastix_int_t);
            memcpy(pastix_data->col2, colptr, (n+1)*sizeof(pastix_int_t));
            memcpy(pastix_data->row2, rows,     nnz*sizeof(pastix_int_t));
        }

        if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
            pastix_print(procnum, 0, "%s", OUT_NODIAG);

        csc_noDiag(pastix_data->col2[0], n,
                   pastix_data->col2,
                   pastix_data->row2, NULL);

        pastix_data->gN = n;
    }
#if defined(PASTIX_DISTRIBUTED)
    /*
     * Distributed csc
     */
    else
    {
        MPI_Comm      pastix_comm = pastix_data->pastix_comm;
        pastix_int_t  gN = 0;
        int copy_l2g = 1;

        MPI_Allreduce(&n, &gN, 1, PASTIX_MPI_INT, MPI_SUM, pastix_comm);
        if (iparm[IPARM_SYM]==API_SYM_YES || iparm[IPARM_SYM] == API_SYM_HER) {
            cscd_symgraph_int(n, colptr, rows, NULL,
                              &(pastix_data->n2),
                              &(pastix_data->col2),
                              &(pastix_data->row2), NULL,
                              loc2glob,
                              pastix_comm, API_YES );
            assert( n == pastix_data->n2 );
        }

        cscd_noDiag(pastix_data->n2,
                    pastix_data->col2,
                    pastix_data->row2,
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
                pastix_int_t ldisp;
                int *all_n;
                int *displs;

                /* Gather the locals n */
                MALLOC_INTERN(all_n,  pastix_data->procnbr, int);
                MALLOC_INTERN(displs, pastix_data->procnbr, int);

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

                /* Apply the new permutation to the local csc */
                for (i = 0; i < (pastix_data->col2)[n] - 1; i++)
                    pastix_data->row2[i] = pastix_data->PTS_permtab[(pastix_data->row2)[i]-1];

                /* Initialize loc2glob */
                copy_l2g = 0;
                MALLOC_INTERN(pastix_data->loc2glob2, n, pastix_int_t);
                for (i = 0; i < n; i++,ldisp++)
                    pastix_data->loc2glob2[i] = ldisp;
            }
        }

        if (copy_l2g)
        {
            MALLOC_INTERN(pastix_data->loc2glob2, n, pastix_int_t);
            memcpy((pastix_data->loc2glob2), loc2glob, n*sizeof(pastix_int_t));
        }
    }
#else
    assert(loc2glob == NULL );
#endif

    pastix_data->bmalcolrow = 1;
    return PASTIX_SUCCESS;
}
