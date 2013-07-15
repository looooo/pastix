/**
 *
 * @file pastix_task_symbfact.c
 *
 *  PaStiX symbolic factorizations routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains wrapper to the symolic factorization step.
 * Affetcted by the compilation time options:
 *    - PASTIX_SYMBOL_FORCELOAD: Force to load the symbol matrix from file
 *    - PASTIX_SYMBOL_DUMP_SYMBMTX: Dump the symbol matrix in a postscript file.
 *    - COMPACT_SMX: Optimization for solve step (TODO: check if not obsolete)
 *    - FORGET_PARTITION: Force to forget the precomputed partition
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> s d c
 *
 **/

/*
    Fileunction: pastix_task_fax

 Symbolic factorisation.

 Parameters:
 pastix_data - PaStiX data structure
 pastix_comm - PaStiX MPI communicator
 n           - Size of the matrix
 perm        - permutation tabular
 invp        - reverse permutation tabular
 flagWinvp   - flag to indicate if we have to print warning concerning perm and invp modification.

 */
#include "common.h"
#ifdef WITH_SCOTCH
#  ifdef    PASTIX_DISTRIBUTED
#    include <ptscotch.h>
#  else
#    include <scotch.h>
#  endif /* PASTIX_DISTRIBUTED */
#endif /* WITH_SCOTCH */

#include "kass.h"
#include "csc_utils.h"
#include "cscd_utils_intern.h"
#include "fax.h"

void global2localperm(pastix_int_t  lN,
                      pastix_int_t *lperm,
                      pastix_int_t *gperm,
                      pastix_int_t *loc2glob);


#ifdef PASTIX_DISTRIBUTED
/*
 Function: dpastix_task_fax

 Symbolic factorisation.

 Parameters:
 pastix_data - PaStiX data structure
 pastix_comm - PaStiX MPI communicator
 n           - Size of the local matrix
 perm        - local permutation tabular
 loc2glob    - Global number of local columns (NULL if not ditributed)
 flagWinvp   - flag to indicate if we have to print warning concerning perm and invp modification.

 */
void dpastix_task_fax(pastix_data_t *pastix_data, MPI_Comm pastix_comm, pastix_int_t n,
                      pastix_int_t * perm, pastix_int_t * loc2glob, int flagWinvp)
{
    pastix_int_t   * gperm   = NULL;
    pastix_int_t   * ginvp   = NULL;
    pastix_int_t     gN;
    pastix_int_t     my_n;
    pastix_int_t   * my_perm = NULL;
    (void)pastix_comm;

    /* Note: for AUTOSPLIT_COMM
     *
     * perm is given by user, of size n,
     * we have to allocate it so that fax can write into it,
     * only written, data can be anything...
     *
     * loc2glob is also given by user and we need
     * a gathered one to build gperm from perm returned by fax.
     *
     * Anyway perm can't be used by user as it does not correspond to
     * the columns it gaves us. We just need to allocate data to write trash.
     */
    MPI_Allreduce(&n, &my_n, 1, PASTIX_MPI_INT, MPI_SUM,
                  pastix_data->intra_node_comm);
    if (pastix_data->intra_node_procnum == 0)
    {
        if (my_n != n)
        {
            MALLOC_INTERN(my_perm, my_n, pastix_int_t);
            perm = my_perm;
            if (pastix_data->procnum == 0)
                errorPrintW("User's perm array is invalid and will be unusable with IPARM_AUTOSPLIT_COMM");
        }

        if (!(PASTIX_MASK_ISTRUE(pastix_data->iparm[IPARM_IO_STRATEGY], API_IO_LOAD)))
        {
            gN = 0;
            MPI_Allreduce(&my_n, &gN, 1, PASTIX_MPI_INT, MPI_SUM, pastix_data->inter_node_comm);
            MALLOC_INTERN(gperm, gN, pastix_int_t);
            MALLOC_INTERN(ginvp, gN, pastix_int_t);
        }

        {
            MPI_Comm tmpcomm = pastix_data->pastix_comm;
            pastix_data->pastix_comm = pastix_data->inter_node_comm;
            pastix_task_symbfact( pastix_data,
                                  gperm, ginvp, flagWinvp);
            pastix_data->pastix_comm = tmpcomm;
        }

        if (my_n == n)
        {
            if (!(PASTIX_MASK_ISTRUE(pastix_data->iparm[IPARM_IO_STRATEGY], API_IO_LOAD)))
            {
                /* permtab may have been changed */
                global2localperm(my_n, perm, gperm, loc2glob);

                memFree_null(gperm);
                memFree_null(ginvp);
            }
        }
        else
        {
            memFree_null(my_perm);
        }
    }
    else
    {
        if ( ( pastix_data->iparm[IPARM_INCOMPLETE] == API_NO &&
               ( pastix_data->iparm[IPARM_ORDERING] == API_ORDER_PERSONAL ||
                 pastix_data->iparm[IPARM_ORDERING] == API_ORDER_METIS ||
                 pastix_data->iparm[IPARM_LEVEL_OF_FILL] == -1) &&
               pastix_data->iparm[IPARM_GRAPHDIST] == API_YES ) ||
             ( pastix_data->iparm[IPARM_INCOMPLETE] == API_YES &&
               pastix_data->iparm[IPARM_GRAPHDIST] == API_YES ) )
        {
            pastix_int_t  nkass;
            pastix_int_t *colptrkass;
            pastix_int_t *rowkass;

            CSC_sort(pastix_data->n2, pastix_data->col2, pastix_data->row2, NULL);

            cscd2csc_int(pastix_data->n2, pastix_data->col2, pastix_data->row2,
                         NULL,
                         NULL, NULL, NULL,
                         &nkass, &colptrkass, &rowkass, NULL,
                         NULL, NULL, NULL,
                         pastix_data->loc2glob2,
                         pastix_data->pastix_comm, pastix_data->iparm[IPARM_DOF_NBR], API_YES);
            memFree_null(colptrkass);
            memFree_null(rowkass);
        }

        if (pastix_data->bmalcolrow == 1)
        {
            if (pastix_data->col2      != NULL) memFree_null(pastix_data->col2);
            if (pastix_data->row2      != NULL) memFree_null(pastix_data->row2);
            if (pastix_data->loc2glob2 != NULL) memFree_null(pastix_data->loc2glob2);
            pastix_data->bmalcolrow = 0;
        }

    }
}

#endif
/*
 Function: pastix_task_symbfact

 Symbolic factorisation.

 Parameters:
 pastix_data - PaStiX data structure
 perm        - permutation tabular
 invp        - reverse permutation tabular
 flagWinvp   - flag to indicate if we have to print warning concerning perm and
               invp modification.

 */
void pastix_task_symbfact(pastix_data_t *pastix_data,
                          pastix_int_t  *perm,
                          pastix_int_t  *invp,
                          int flagWinvp)
{
    pastix_int_t   *iparm = pastix_data->iparm;
    pastix_graph_t *graph = pastix_data->csc;
    Order          *ordemesh;
    pastix_int_t    n;
    Clock           timer;
    int             procnum;
    int             retval = PASTIX_SUCCESS;
    int             retval_rcv;

#ifdef PASTIX_DISTRIBUTED
    PASTIX_INT           * PTS_perm     = pastix_data->PTS_permtab;
    PASTIX_INT           * PTS_rev_perm = pastix_data->PTS_peritab;
    PASTIX_INT           * tmpperm      = NULL;
    PASTIX_INT           * tmpperi      = NULL;
    PASTIX_INT             gN;
    PASTIX_INT             i;
#endif

    n        = graph->n;
    ordemesh = pastix_data->ordemesh;
    procnum  = pastix_data->procnum;

    print_debug(DBG_STEP, "-> pastix_task_symbfact\n");
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        pastix_print(procnum, 0, "%s", OUT_STEP_FAX);

    /* Force Load of symbmtx */
#if defined(PASTIX_SYMBOL_FORCELOAD)
    iparm[IPARM_IO_STRATEGY] = API_IO_LOAD;
#endif

    if (pastix_data->symbmtx == NULL) {
        MALLOC_INTERN( pastix_data->symbmtx, 1, SymbolMatrix );
    }
    else {
        errorPrint("PASTIX SymbFact: Symbol Matrix already allocated !!!");
    }

    /*
     * Load i/o strategy
     */
    if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD))
    {
        FILE *stream;

        /* Load graph if not already defined */
        /* if (ordemesh == NULL) { */
        /*     pastix_csc_t csc; */
        /*     orderLoadFiles( pastix_data ); */
        /* } */

        /* Load ordering if not already defined */
        if (ordemesh == NULL) {
            orderLoadFiles( pastix_data );
        }

        /* Load symbol */
        PASTIX_FOPEN(stream, "symbname", "r" );
        symbolLoad( pastix_data->symbmtx, stream );
        fclose(stream);

        /* Rebase to 0 */
        symbolBase( pastix_data->symbmtx, 0 );

    }
    /* not API_IO_LOAD */
    else
    {
        pastix_int_t  nfax;
        pastix_int_t *colptrfax;
        pastix_int_t *rowfax;

        /* Check correctness of parameters */
        if (iparm[IPARM_INCOMPLETE] == API_NO)
        {
#ifdef COMPACT_SMX
            if (procnum == 0)
                errorPrintW("COMPACT_SMX only works with incomplete factorisation, forcing incomplete factorisation.");
            iparm[IPARM_INCOMPLETE] = API_YES;
#endif /* COMPACT_SMX */
        }
        /* End of parameters check */

        /*
         * Fax works with centralized interface, we convert the cscd to csc if required
         */
        if (iparm[IPARM_GRAPHDIST] == API_YES)
        {
            CSC_sort( graph->n,
                      graph->colptr,
                      graph->rows,
                      NULL);

            cscd2csc_int( graph->n,
                          graph->colptr,
                          graph->rows,
                          NULL, NULL, NULL, NULL,
                          &nfax, &colptrfax, &rowfax,
                          NULL, NULL, NULL, NULL,
                          graph->loc2glob,
                          pastix_data->pastix_comm,
                          iparm[IPARM_DOF_NBR], API_YES);
        }
        else
        {
            nfax      = graph->n;
            colptrfax = graph->colptr;
            rowfax    = graph->rows;
        }

        symbolInit(pastix_data->symbmtx);

        /*
         * The amalgamate supernodes partition has been found with (PT-)Scotch,
         * we use it to generate the symbol matrix structure.
         */
        if ( (iparm[IPARM_INCOMPLETE]    == API_NO) &&
             (iparm[IPARM_LEVEL_OF_FILL] != -1    ) &&
             (ordemesh->rangtab != NULL) )
        {
            symbolFaxGraph(pastix_data->symbmtx, /* Symbol Matrix   */
                           nfax,                 /* Number of nodes */
                           colptrfax,            /* Nodes list      */
                           rowfax,               /* Edges list      */
                           ordemesh);
        }
        /*
         * The amalgamate supernodes partition doesn't not exist. (PT-)Scotch
         * has not been used, ILU(k) factorization is performed or we dropped
         * the partition found by Scotch.
         * We use Kass to generate both the amalgamate supernode partition and
         * the symbol matrix stucture.
         */
        else
        {
            fprintf(stderr, "HELLO I'm GOING THROUGH NEW KASS\n");
            {
                pastix_graph_t tmpgraph;
                tmpgraph.n      = nfax;
                tmpgraph.colptr = colptrfax;
                tmpgraph.rows   = rowfax;
                tmpgraph.loc2glob = NULL;

                kass2(iparm[IPARM_INCOMPLETE],
                      iparm[IPARM_LEVEL_OF_FILL],
                      iparm[IPARM_AMALGAMATION_LEVEL],
                      pastix_data->symbmtx,
                      &tmpgraph,
                      ordemesh,
                      pastix_data->pastix_comm);
            }
        }

        if (iparm[IPARM_GRAPHDIST] == API_YES)
        {
            memFree_null(colptrfax);
            memFree_null(rowfax);
        }

        symbolBase(pastix_data->symbmtx,0);

#ifdef PASTIX_DISTRIBUTED
        if (PTS_perm != NULL)
        {
            gN = n;

            MALLOC_INTERN(tmpperm, gN, PASTIX_INT);
            MALLOC_INTERN(tmpperi, gN, PASTIX_INT);
            for (i = 0; i < gN; i++)
                tmpperm[i] = ordemesh->permtab[PTS_perm[i]-1];

            memFree_null(ordemesh->permtab);
            ordemesh->permtab = tmpperm;

            for (i = 0; i < gN; i++)
                tmpperi[i] = PTS_rev_perm[ordemesh->peritab[i]]-1;
            memFree_null(ordemesh->peritab);
            ordemesh->peritab = tmpperi;

            memFree_null(PTS_perm);
            memFree_null(PTS_rev_perm);
        }
#endif /* PASTIX_DISTRIBUTED */

        /* WARNING : perm and invp can now be modified during symbolic factorization ??? */
        if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
            if (flagWinvp)
                if (procnum == 0)
                    errorPrintW("perm and invp can be modified during symbolic factorization.");

        memcpy(perm, ordemesh->permtab, n*sizeof(PASTIX_INT));
        memcpy(invp, ordemesh->peritab, n*sizeof(PASTIX_INT));
    } /* not API_IO_LOAD */

    if (pastix_data->csc != NULL)
    {
        graphClean( pastix_data->csc );
        pastix_data->csc = NULL;
    }

    /*
     * Save the symbolic factorization
     */
    if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE))
    {
        FILE *stream;
        PASTIX_FOPEN(stream, "symbgen", "w");
        symbolSave(pastix_data->symbmtx, stream);
        fclose(stream);
    }

#if defined(PASTIX_SYMBOL_DUMP_SYMBMTX)
    if (pastix_data->procnum == 0)
    {
        FILE *stream;
        PASTIX_FOPEN(stream, "symbol.eps", "w");
        symbolDraw(pastix_data->symbmtx,
                   stream);
        fclose(stream);
    }
#endif

    iparm[IPARM_START_TASK]++;
}
