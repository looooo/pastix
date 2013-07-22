
/************************************************************/
/**                                                        **/
/**   NAME       : amalgamate.c                            **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 15/08/2006      **/
/**                                                        **/
/**                                                        **/
/************************************************************/
#include "common.h"
#include "kass.h"

pastix_int_t
SF_level( const kass_csr_t   *graphA,
                pastix_int_t  level,
                kass_csr_t   *graphL)
{
    /***********************************************************************************************/
    /* This function computes the non zero pattern of the levelized incomplete factor              */
    /* for a sparse lower triangular                                                               */
    /* matrix in CSC format.  This pattern is exact iff the matrix has a SYMMETRIC non zero        */
    /* structure.                                                                                  */
    /* On entry:                                                                                   */
    /* job : = 0 alloc the matrix ; =1 fill the indices of the sparse pattern  =2 alloc (not cont) */
    /*       and fill the coefficients                                                             */
    /* A : pointer to the matrix                                                                   */
    /* level : level desired for the ilu(k) factorization                                          */
    /* P   : an empty csr matrix (initialized with dimension n)                                    */
    /* On return:                                                                                  */
    /*     P : a csr matrix containing the non zero pattern of the factorized matrix               */
    /*         the memory for the numerical values is allocated and initialized to zero            */
    /*     The total number of nnz in the lower part is returned                                   */
    /*                                                                                             */
    /* NOTE:                                                                                       */
    /*   1) This algorithm has been implemented according to the paper of David Hysom and          */
    /*     Alex Pothen : Level-based Incomplete LU factorization: Graph Model and Algorithm        */
    /***********************************************************************************************/
    pastix_int_t *visited = NULL;
    pastix_int_t *length  = NULL;
    pastix_int_t *stack   = NULL;
    pastix_int_t *adj     = NULL;
    pastix_int_t *ja      = NULL;
    pastix_int_t used;
    pastix_int_t h, i,j,k, t;
    long nnz;

    if(graphA->n == 0)
        return 0;

    /** Allocated the working array **/
    MALLOC_INTERN(visited, graphA->n, pastix_int_t);
    MALLOC_INTERN(length,  graphA->n, pastix_int_t);
    MALLOC_INTERN(stack,   graphA->n, pastix_int_t);
    MALLOC_INTERN(ja,      graphA->n, pastix_int_t);
    nnz = 0;

    /** Initialized visited ***/
    for(j=0;j<graphA->n;j++)
    {
        visited[j] = -1;
        length[j]  = 0;
    }

    /** Apply GS_Urow for each row **/
    kass_csrInit( graphA->n, graphL );
    for(i=0;i<graphA->n;i++)
    {
        /** Reset the stack number of elements **/
        stack[0] = i;
        used = 1;

        length[i]  = 0;
        visited[i] = i;

        ja[0] = i; /** Put the diagonal term **/
        k = 1;
        /** BFS phase **/
        while(used > 0)
        {
            used--;
            h   = stack[used];
            adj = graphA->rows[h];
            for(j=0;j<graphA->nnz[h];j++)
            {
                t = adj[j];
                if(visited[t] != i)
                {
                    visited[t] = i;
                    if( (t < i) && (length[h] < level) )
                    {
                        stack[used] = t;
                        used++;
                        length[t] = length[h]+1;
                    }
                    if( t > i )
                    {
                        ja[k++] = t;
                    }
                }
            }
        }

        assert( k > 0 );

        graphL->nnz[i] = k;
        MALLOC_INTERN(graphL->rows[i], k, pastix_int_t);
        memcpy(graphL->rows[i], ja, k * sizeof(pastix_int_t));

        intSort1asc1( graphL->rows[i],
                      graphL->nnz[i]);

        nnz += k;
    }

    memFree_null(ja);
    memFree_null(visited);
    memFree_null(length);
    memFree_null(stack);

    return nnz;
}
