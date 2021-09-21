/**
 *
 * @file order_compute_common.c
 *
 * PaStiX order common routine between order_compute_scotch.c and order_compute_ptscotch.c.
 *
 * @copyright 2004-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Vincent Bridonneau
 * @date 2021-06-16
 *
 **/
#include "common.h"
#include "graph/graph.h"
#include "pastix/order.h"
#include "order_scotch_strats.h"

#define STRAT_STR_MAX 1024

#define DIRECT(_isPTScotch_) (_isPTScotch_) ? PTSCOTCH_STRAT_DIRECT : SCOTCH_STRAT_DIRECT
#define INCOMP(_isPTScotch_) (_isPTScotch_) ? PTSCOTCH_STRAT_INCOMP : SCOTCH_STRAT_INCOMP
#define PERSON(_isPTScotch_) (_isPTScotch_) ? PTSCOTCH_STRAT_PERSO : SCOTCH_STRAT_PERSO

#define DIRECT_STR(_isPTScotch_) (_isPTScotch_) ? "      PT-Scotch direct strategy\n"       : "      Scotch direct strategy\n"
#define INCOMP_STR(_isPTScotch_) (_isPTScotch_) ? "      PT-Scotch incomplete strategy\n"   : "      Scotch incomplete strategy\n"
#define PERSON_STR(_isPTScotch_) (_isPTScotch_) ? "      PT-Scotch personal strategy |%s|\n": "      Scotch personal strategy |%s|\n"

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Build the vertex weight array out of the dof array.
 *
 *******************************************************************************
 *
 * @param[in] graph
 *          Pointer to the graph structure.
 *
 *******************************************************************************
 *
 * @retval The vertex weight array if graph->dof != 1, NULL otherwise.
 *
 *******************************************************************************/
pastix_int_t *
order_compute_build_weights( const pastix_graph_t *graph )
{
    pastix_int_t  i, n;
    pastix_int_t *weights, *wptr;

    if ( graph->dof == 1 ) {
        return NULL;
    }

    n = graph->n;
    MALLOC_INTERN( weights, n, pastix_int_t );

    wptr = weights;
    /* Constant dof */
    if ( graph->dof > 1 ) {
        for (i = 0; i < n; i++, wptr++ ) {
            *wptr = graph->dof;
        }
    }
    /* Variadic dof */
    else {
        pastix_int_t *dofptr = graph->dofs;
        if ( graph->loc2glob == NULL ) {
            for (i = 0; i < n; i++, wptr++, dofptr++) {
                *wptr = dofptr[1] - dofptr[0];
            }
        }
        else {
            pastix_int_t *dofptr   = graph->dofs - graph->baseval;
            pastix_int_t *loc2glob = graph->loc2glob;

            for (i = 0; i < n; i++, wptr++, loc2glob++) {
                *wptr = dofptr[ *loc2glob + 1 ] - dofptr[ *loc2glob ];
            }
        }
    }

    return weights;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Generate the ordering strategy string based on the input parameters.
 *
 *******************************************************************************
 *
 * @param[inout] strat
 *          The preallocated ordering strategy string to initialize.
 *
 * @param[in] iparm
 *          Pointer to the iparm array.
 *
 * @param[in] procnum
 *          Procnum of the process. Output purpose.
 *
 *******************************************************************************/
void
order_compute_build_strategy( char               *strat,
                              const pastix_int_t *iparm,
                              pastix_int_t        procnum,
                              int                 isPTscotch )
{
    int rc;

    /* Default ordering */
    if (iparm[IPARM_ORDERING_DEFAULT] == 1) {
        if (iparm[IPARM_INCOMPLETE] == 0) {
            if (iparm[IPARM_VERBOSE] > PastixVerboseNo) {
                pastix_print( procnum, 0, DIRECT_STR(isPTscotch) );
            }
            snprintf( strat, STRAT_STR_MAX, DIRECT(isPTscotch) );
        }
        else {
            if (iparm[IPARM_VERBOSE] > PastixVerboseNo) {
                pastix_print( procnum, 0, INCOMP_STR(isPTscotch) );
            }
            snprintf(strat, STRAT_STR_MAX, INCOMP(isPTscotch) );
        }
    }
    /* Personal ordering */
    else {
        rc = snprintf( strat, STRAT_STR_MAX, PERSON(isPTscotch),
                       (long)  iparm[IPARM_SCOTCH_SWITCH_LEVEL],
                       (long)  iparm[IPARM_SCOTCH_CMIN],
                       (long)  iparm[IPARM_SCOTCH_CMAX],
                       ((float)iparm[IPARM_SCOTCH_FRAT])/100.,
                       (long)  iparm[IPARM_SCOTCH_SWITCH_LEVEL],
                       (long)  iparm[IPARM_SCOTCH_CMIN],
                       (long)  iparm[IPARM_SCOTCH_CMAX],
                       ((float)iparm[IPARM_SCOTCH_FRAT])/100. );
        if ( rc > STRAT_STR_MAX ) {
            pastix_print_error( "Order_compute_build_strategy: Strategy string too long\n" );
            exit(-1);
        }

        if (iparm[IPARM_VERBOSE] > PastixVerboseNo) {
            pastix_print( procnum, 0, PERSON_STR(isPTscotch), strat );
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Reallocate the ordering structure.
 *
 * If we decide to drop the Scothc partition to recompute it later, then
 * partition information is freed, otherwise its memory space is compressed.
 *
 *******************************************************************************
 *
 * @param[inout] ordemesh
 *          Pointer to the ordemesh structure to reallocate.
 *
 *******************************************************************************/
void
order_compute_reallocate_ordemesh( pastix_order_t *ordemesh )
{
#if defined(FORGET_PARTITION)
    ordemesh->cblknbr = 0;
    if (ordemesh->rangtab != NULL) {
        memFree_null(ordemesh->rangtab);
    }
    if (ordemesh->treetab != NULL) {
        memFree_null(ordemesh->treetab);
    }
#else
    /**
     * Adapt size of rangtab and treetab to the new cblknbr
     * WARNING: If no nodes in the graph, nothing has been initialized.
     */
    ordemesh->rangtab =
        (pastix_int_t *) memRealloc( ordemesh->rangtab,
                                    (ordemesh->cblknbr + 1)*sizeof(pastix_int_t) );
    ordemesh->treetab =
        (pastix_int_t *) memRealloc( ordemesh->treetab,
                                    (ordemesh->cblknbr)*sizeof(pastix_int_t) );
    if (ordemesh->cblknbr == 0) {
        ordemesh->rangtab[0] = ordemesh->baseval;
    }
#endif
}
