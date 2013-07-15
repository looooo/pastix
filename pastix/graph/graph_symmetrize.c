/**
 *
 * @file graph_symmetrize.c
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
#include "csc_utils.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_graph
 *
 * graphSymmetrize - Symmetrize a given graph
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pointer to the solver instance to get options as rank,
 *          communicators, ...
 *
 * @param[in,out] graph
 *          The graph structure to store the loaded graph.
 *
 *******************************************************************************/
/*
  Function: csc_symgraph_int


  Modify the CSC to a symetric graph one.
  Don't use it on a lower symetric CSC
  it would give you all the CSC upper + lower.

  Parameters:
    n           - Number of columns/vertices
    ia          - Starting index of each column in *ja* and *a*
    ja          - Row index of each element
    a           - Value of each element,can be NULL
    newn        - New number of column
    newia       - Starting index of each column in *ja* and *a*
    newja       - Row index of each element
    newa        - Value of each element,can be NULL
    malloc_flag - flag to indicate if function call is intern to pastix or extern.
 */
int graphSymmetrize_int( pastix_int_t n,
                         const pastix_int_t *ia,
                         const pastix_int_t *ja,
                         const pastix_int_t *loc2glob,
                         pastix_graph_t *newgraph,
                         int malloc_flag )
{
    pastix_int_t *nbrEltCol = NULL; /* nbrEltCol[i] = Number of elt to add in column i */
    pastix_int_t  itercol, iterrow, iterrow2; /* iterators */
    pastix_int_t *newia;
    pastix_int_t *newja;
    pastix_int_t  nnz;
    int baseval = ia[0];
    (void)loc2glob;

    MALLOC_INTERN(nbrEltCol, n, pastix_int_t);
    /* !! Need check for malloc */

    /* Init nbrEltCol */
    for (itercol=0; itercol<n; itercol++)
    {
        nbrEltCol[itercol] = 0;
    }

    /*
     * Compute number of elements by col to add for correcting the CSC
     */
    for (itercol=0; itercol<n; itercol++)
    {
        pastix_int_t frow = ia[itercol]   - baseval;
        pastix_int_t lrow = ia[itercol+1] - baseval;
        for (iterrow=frow; iterrow<lrow; iterrow++)
        {
            pastix_int_t rowidx = ja[iterrow]  - baseval;
            if ( rowidx != itercol )
            {
                /* It is not a diagonal element, so we have (i,j) and we look for (j,i) element */
                /* i = itercol+1, j=ja[iterrow] */
                pastix_int_t frow2 = ia[rowidx]   - baseval;
                pastix_int_t lrow2 = ia[rowidx+1] - baseval;
                int flag = 0;

                for (iterrow2=frow2; iterrow2<lrow2; iterrow2++)
                {
                    if (ja[iterrow2] == itercol+baseval)
                    {
                        /* We found (j,i) so let's stop this madness */
                        flag = 1;
                        break;
                    }
                }

                if (flag == 0)
                {
                    /* We never found (j,i) so we increase nbrEltCol[j] */
                    (nbrEltCol[rowidx])++;
                }
            }
        }
    }

    /* Let's compute the new ia in C numbering */
    if (malloc_flag == API_YES) {
        MALLOC_INTERN(newia, n+1, pastix_int_t);
    }
    else {
        MALLOC_EXTERN(newia, n+1, pastix_int_t);
    }

    newia[0] = ia[0];
    for (itercol=0;itercol<n;itercol++)
    {
        newia[itercol+1] = newia[itercol] + ia[itercol+1] - ia[itercol] + nbrEltCol[itercol];
    }

    assert( newia[n] >= ia[n] );
    nnz = newia[n] - baseval;

    /* Let's build the new ja */
    if (malloc_flag == API_YES) {
        MALLOC_INTERN(newja, nnz, pastix_int_t);
    }
    else {
        MALLOC_EXTERN(newja, nnz, pastix_int_t);
    }

    if ( newia[n] > ia[n])
    {
        for (itercol=0;itercol<nnz;itercol++)
        {
            newja[itercol] = -1;
        }

        for (itercol=0;itercol<n;itercol++)
        {
            pastix_int_t frow  = ia[itercol]   - baseval;
            pastix_int_t lrow  = ia[itercol+1] - baseval;
            pastix_int_t nbelt = ia[itercol+1] - ia[itercol];

            assert( newia[itercol] >= ia[itercol] );
            assert( newia[itercol] < newia[itercol+1] );

            /* Let's copy what we already have at the end of the space reserved,
             * we will add the new elements in front of it */
            memcpy( newja + (newia[itercol+1] - baseval) - nbelt,
                    ja    + frow,
                    nbelt * sizeof(pastix_int_t) );

            /* We add (j, i) edges missing */
            for (iterrow=frow; iterrow<lrow; iterrow++)
            {
                pastix_int_t rowidx = ja[iterrow] - baseval;
                if ( rowidx != itercol )
                {
                    /* It is not a diagonal element, so we have (i,j) and we look for (j,i) element */
                    /* i = itercol+1, j=ja[iterrow] */
                    pastix_int_t frow2 = ia[rowidx]   - baseval;
                    pastix_int_t lrow2 = ia[rowidx+1] - baseval;
                    int flag = 0;

                    for (iterrow2=frow2; iterrow2<lrow2; iterrow2++)
                    {
                        if (ja[iterrow2] == itercol+baseval)
                        {
                            /* We found (j,i) so let's stop this madness */
                            flag = 1;
                            break;
                        }
                    }

                    if (flag == 0)
                    {
                        /* We never found (j,i) so we increase nbrEltCol[j] */
                        nbrEltCol[rowidx]--;
                        assert( nbrEltCol[rowidx] >= 0 );
                        newja[ newia[rowidx] - baseval + nbrEltCol[rowidx] ] = itercol + baseval;
                    }
                }
            }
        }

        for (itercol=0;itercol<nnz;itercol++)
        {
            assert( newja[itercol] != -1);
        }

        /* Sort in place each subset */
        for (itercol=0;itercol<n;itercol++)
        {
            pastix_int_t frow = newia[itercol]   - baseval;
            pastix_int_t lrow = newia[itercol+1] - baseval;

            intSort1asc1( (newja+frow), (lrow-frow));
            assert( nbrEltCol[itercol] == 0);
        }
    }
    else
    {
        memcpy( newja, ja, newia[n] * sizeof(pastix_int_t) );
    }

    memFree_null( nbrEltCol );

    newgraph->gN = n;
    newgraph->n  = n;
    newgraph->colptr   = newia;
    newgraph->rows     = newja;
    newgraph->loc2glob = NULL;

    graphBase( newgraph, 1 );

    return EXIT_SUCCESS;
}
