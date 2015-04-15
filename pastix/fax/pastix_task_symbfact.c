/**
 *
 * @file pastix_task_symbfact.c
 *
 *  PaStiX symbolic factorizations routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains wrappers to the symbolic factorization step.
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
 **/
#include "common.h"
#include "csc.h"
#include "order.h"
#include "fax.h"
#include "kass.h"
#if defined(PASTIX_DISTRIBUTED)
#include "csc_utils.h"
#include "cscd_utils_intern.h"
#endif /* defined(PASTIX_DISTRIBUTED) */

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbfact
 * @ingroup pastix
 *
 * pastix_task_symbfact - Computes the the symbolic factorization of the matrix
 * and if required the amalgamated supernode partition.
 *
 * The function is a *centralized* algorithm to generate the symbol matrix
 * structure associated to the problem. It takes as input the ordemesh structure
 * (permutaion array, inverse permutation array, and optionnal supernodes
 * array) and returns the modified ordemesh structure if changed, and the
 * symbolic structure.
 *  - If (PT-)Scotch has been used, it generates the structure with
 * symbolFaxGraph() thanks to the supernode partition given by Scotch.
 *  - If ILU(k) factorization will be performed or if the ordering tools didn't
 * provide the supernode partition, symbolKass() is used to generate both
 * supernode partition and associated symbol matrix structure.
 *
 * Both algorithms are working with a centralized version of the graph and are
 * on every nodes. If a distributed graph has been used, it is gather on each
 * node to compute the symbol matrix.
 * If symbolKass() is used, the perm and invp vector will be modified and
 * returned to the user. BE CAREFULL if you give your own ordering and wants to
 * keep it because it will be overwritten by the updated one.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_INCOMPLETE, IPARM_LEVEL_OF_FILL,
 *   IPARM_AMALGAMATION_LVLCBLK, IPARM_AMALGAMATION_LVLBLAS,
 *   IPARM_IO_STRATEGY, IPARM_FLOAT, IPARM_FACTORIZATION
 *
 * On exit, the following parameters are set:
 *   IPARM_NNZEROS, DPARM_FACT_THFLOPS, DPARM_FACT_RLFLOPS
 *
 *******************************************************************************
 *
 * @param[in,out] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the field symbmtx is initialized with the symbol matrix,
 *          and the field ordemesh is updated if the supernode partition is
 *          computed.
 *          - IPARM_INCOMPLETE switches the factorization mode from direct to ILU(k).
 *          - IPARM_LEVEL_OF_FILL defines the level of incomplete factorization
 *            if IPARM_INCOMPLETE == API_YES. If IPARM_LEVEL_OF_FILL < 0, the
 *            full pattern is generated as for direct factorization.
 *          - IPARM_AMALGAMATION_LVLCBLK is the ratio of amalgamation allowed
 *            based on reducing the number of supernodes only.
 *          - IPARM_AMALGAMATION_LVLBLAS is the ratio of amalgamation allowed
 *            based on reducing the computational cost (solve for ILU(k), or
 *            factorization for direct factorization).
 *          - IPARM_IO_STRATEGY will enable to load/store the result to files.
 *          If set to API_IO_SAVE, the symbmtx and the generated ordemesh is
 *          dump to file.
 *          If set to APÏ_IO_LOAD, the symbmtx (only) is loaded from the files.
 *
 * @param[in,out] perm
 *          Array of size n.
 *          On entry, unused.
 *          On exit, if perm != NULL, contains the permutation array generated.
 *
 * @param[in,out] invp
 *          Array of size n.
 *          On entry, unused.
 *          On exit, if invp != NULL, contains the inverse permutation array generated.
 *
 *******************************************************************************
 *
 * @return
 *          \retval PASTIX_SUCCESS on successful exit
 *          \retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *          \retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *          \retval PASTIX_ERR_INTEGER_TYPE if Scotch integer type is not the
 *                  same size as PaStiX ones.
 *          \retval PASTIX_ERR_INTERNAL if an error occurs internally to Scotch.
 *
 *******************************************************************************/
int
pastix_task_symbfact(pastix_data_t *pastix_data,
                     pastix_int_t  *perm,
                     pastix_int_t  *invp )
{
    pastix_int_t   *iparm;
    double         *dparm;
    pastix_graph_t *graph;
    Order          *ordemesh;
    pastix_int_t    n;
    int             procnum;

#if defined(PASTIX_DISTRIBUTED)
    pastix_int_t           * PTS_perm     = pastix_data->PTS_permtab;
    pastix_int_t           * PTS_rev_perm = pastix_data->PTS_peritab;
    pastix_int_t           * tmpperm      = NULL;
    pastix_int_t           * tmpperi      = NULL;
    pastix_int_t             gN;
    pastix_int_t             i;
#endif

    /*
     * Check parameters
     */
    if (pastix_data == NULL) {
        errorPrint("pastix_task_symbfact: wrong pastix_data parameter");
        return PASTIX_ERR_BADPARAMETER;
    }
    iparm = pastix_data->iparm;
    dparm = pastix_data->dparm;

    if ( !(pastix_data->steps & STEP_ORDERING) ) {
        errorPrint("pastix_task_symbfact: pastix_task_order() has to be called before calling this function");
        return PASTIX_ERR_BADPARAMETER;
    }

    procnum  = pastix_data->procnum;
    graph    = pastix_data->graph;
    ordemesh = pastix_data->ordemesh;

    if (graph == NULL) {
        errorPrint("pastix_task_symbfact: the pastix_data->graph field has not been initialized, pastix_task_order should be called first");
        return PASTIX_ERR_BADPARAMETER;
    }
    if (ordemesh == NULL) {
        errorPrint("pastix_task_symbfact: the pastix_data->ordemesh field has not been initialized, pastix_task_order should be called first");
        return PASTIX_ERR_BADPARAMETER;
    }
    n = ordemesh->vertnbr;

    print_debug(DBG_STEP, "-> pastix_task_symbfact\n");
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        pastix_print(procnum, 0, OUT_STEP_FAX );

    /* Allocate the symbol matrix structure */
    if (pastix_data->symbmtx == NULL) {
        MALLOC_INTERN( pastix_data->symbmtx, 1, SymbolMatrix );
    }
    else {
        errorPrint("pastix_task_symbfact: Symbol Matrix already allocated !!!");
    }

    /* Force Load of symbmtx */
#if defined(PASTIX_SYMBOL_FORCELOAD)
    iparm[IPARM_IO_STRATEGY] = API_IO_LOAD;
#endif

    /*Symbol matrix loaded from file */
    if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD))
    {
        FILE *stream;
        PASTIX_FOPEN(stream, "symbname", "r" );
        symbolLoad( pastix_data->symbmtx, stream );
        fclose(stream);
    }
    /* Symbol matrix computed through Fax or Kass */
    else
    {
        pastix_int_t  nfax;
        pastix_int_t *colptrfax;
        pastix_int_t *rowfax;

        /* Check correctness of parameters */
        if (iparm[IPARM_INCOMPLETE] == API_NO)
        {
#if defined(COMPACT_SMX)
            if (procnum == 0)
                errorPrintW("COMPACT_SMX only works with incomplete factorization, force ILU(%d) factorization.",
                            iparm[IPARM_LEVEL_OF_FILL]);
            iparm[IPARM_INCOMPLETE] = API_YES;
#endif
        }
        /* End of parameters check */

        /*
         * Fax works with centralized interface, we convert the cscd to csc if required
         */
#if defined(PASTIX_DISTRIBUTED)
        if (graph->loc2glob != NULL)
        {
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
#endif
        {
            nfax      = graph->n;
            colptrfax = graph->colptr;
            rowfax    = graph->rows;
        }

        symbolInit(pastix_data->symbmtx);
        pastix_data->symbmtx->dof = graph->dof;

        /*
         * The amalgamate supernodes partition has been found with (PT-)Scotch,
         * we use it to generate the symbol matrix structure.
         * This works only if direct factorization will be performed.
         */
        if ( (iparm[IPARM_INCOMPLETE]    == API_NO) &&
             (iparm[IPARM_LEVEL_OF_FILL] != -1    ) &&
             (ordemesh->rangtab != NULL) )
        {
            pastix_print(procnum, 0, OUT_FAX_METHOD, "Fax " );
            symbolFaxGraph(pastix_data->symbmtx, /* Symbol Matrix   */
                           nfax,                 /* Number of nodes */
                           colptrfax,            /* Nodes list      */
                           rowfax,               /* Edges list      */
                           ordemesh);
        }
        /*
         * The amalgamate supernodes partition doesn't exist. (PT-)Scotch has
         * not been used, or ILU(k) factorization is performed and then, we
         * dropped the partition found by Scotch.
         * Kass is used to generate both the amalgamate supernode partition and
         * the symbol matrix stucture in this case.
         */
        else
        {
            pastix_graph_t tmpgraph;
            tmpgraph.gN     = nfax;
            tmpgraph.n      = nfax;
            tmpgraph.colptr = colptrfax;
            tmpgraph.rows   = rowfax;
            tmpgraph.loc2glob = NULL;

            pastix_print(procnum, 0, OUT_FAX_METHOD, "Kass" );
            symbolKass(iparm[IPARM_INCOMPLETE],
                       iparm[IPARM_LEVEL_OF_FILL],
                       iparm[IPARM_AMALGAMATION_LEVEL],
                       iparm[IPARM_AMALGAMATION_LEVEL],
                       pastix_data->symbmtx,
                       &tmpgraph,
                       ordemesh,
                       pastix_data->pastix_comm);
        }

        if ( graph->loc2glob != NULL )
        {
            memFree_null(colptrfax);
            memFree_null(rowfax);
        }

#if defined(PASTIX_DISTRIBUTED)
        if (PTS_perm != NULL)
        {
            gN = n;

            MALLOC_INTERN(tmpperm, gN, pastix_int_t);
            MALLOC_INTERN(tmpperi, gN, pastix_int_t);
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
#endif /* defined(PASTIX_DISTRIBUTED) */

        /*
         * Save the new ordering structure
         */
        if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE))
        {
            if (procnum == 0) {
                orderSave( ordemesh, NULL );
            }
        }

        if (perm != NULL) memcpy(perm, ordemesh->permtab, n*sizeof(pastix_int_t));
        if (invp != NULL) memcpy(invp, ordemesh->peritab, n*sizeof(pastix_int_t));
    } /* not API_IO_LOAD */

    /*
     * The graph is not useful anymore, we clean it
     */

    /* TODO: FREE SOMEWHERE */
    /* The graph is required to compute twice the symbolic factorization */
    /* if (pastix_data->graph != NULL) */
    /* { */
    /*     graphExit( pastix_data->graph ); */
    /*     memFree_null( pastix_data->graph ); */
    /* } */

    /* Rebase to 0 */
    symbolBase( pastix_data->symbmtx, 0 );

    symbolPrintStats( pastix_data->symbmtx );
    symbolCheckProperties( pastix_data->symbmtx );

    /* Rustine to be sure we have a tree
     * TODO: check difference with kassSymbolPatch */
#define RUSTINE
#ifdef RUSTINE
    symbolRustine( pastix_data->symbmtx,
                   pastix_data->symbmtx );
#endif

    /* Build the browtabs and Realign data structure */
    symbolBuildRowtab( pastix_data->symbmtx );
    symbolRealloc( pastix_data->symbmtx );

    /*
     * Save the symbolic factorization
     */
    if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE))
    {
        if (procnum == 0) {
            FILE *stream;
            PASTIX_FOPEN(stream, "symbgen", "w");
            symbolSave(pastix_data->symbmtx, stream);
            fclose(stream);
        }
    }

    /*
     * Dump an eps file of the symbolic factorization
     */
#if defined(PASTIX_SYMBOL_DUMP_SYMBMTX)
    if (procnum == 0)
    {
        FILE *stream;
        PASTIX_FOPEN(stream, "symbol.eps", "w");
        symbolDraw(pastix_data->symbmtx,
                   stream);
        fclose(stream);
    }
#endif

    /*
     * Computes statistics and print informations
     */
    iparm[IPARM_NNZEROS] = symbolGetNNZ( pastix_data->symbmtx );
    symbolGetFlops( pastix_data->symbmtx,
                    iparm[IPARM_FLOAT], iparm[IPARM_FACTORIZATION],
                    &(dparm[DPARM_FACT_THFLOPS]),
                    &(dparm[DPARM_FACT_RLFLOPS]) );

    if ( procnum == 0 ) {
        if ( iparm[IPARM_VERBOSE] > API_VERBOSE_YES ) {
            double fillin = (double)iparm[ IPARM_NNZEROS ]
                / (double)( (pastix_data->csc)->gnnz );

            fprintf(stdout, OUT_GLOBAL_NNZL,   iparm[ IPARM_NNZEROS ] );
            fprintf(stdout, OUT_GLOBAL_FILLIN, fillin );
            if ( iparm[IPARM_VERBOSE] > API_VERBOSE_CHATTERBOX ) {
                fprintf(stdout, OUT_GLOBAL_FILLIN, fillin );
                fprintf(stdout, OUT_GLOBAL_FILLIN, fillin );
            }
        }
    }

    /* Invalidate following steps, and add order step to the ones performed */
    pastix_data->steps &= ~( STEP_ANALYSE |
                             STEP_NUMFACT |
                             STEP_SOLVE   |
                             STEP_REFINE  );
    pastix_data->steps |= STEP_SYMBFACT;

    iparm[IPARM_START_TASK]++;

    return PASTIX_SUCCESS;
}
