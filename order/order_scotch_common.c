/**
 *
 * @file order_scotch_common.c
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
#include "order/order_internal.h"
#include "order_scotch_strats.h"

#define STRAT_STR_MAX 1024

#define STRAT_DIRECT(_isPTScotch_) ( (_isPTScotch_) ? PTSCOTCH_STRAT_DIRECT : SCOTCH_STRAT_DIRECT )
#define STRAT_INCOMP(_isPTScotch_) ( (_isPTScotch_) ? PTSCOTCH_STRAT_INCOMP : SCOTCH_STRAT_INCOMP )
#define STRAT_PERSON(_isPTScotch_) ( (_isPTScotch_) ? PTSCOTCH_STRAT_PERSO  : SCOTCH_STRAT_PERSO )

#define OUTPUT_DIRECT(_isPTScotch_) ( (_isPTScotch_) ? "      PT-Scotch direct strategy\n"       : "      Scotch direct strategy\n" )
#define OUTPUT_INCOMP(_isPTScotch_) ( (_isPTScotch_) ? "      PT-Scotch incomplete strategy\n"   : "      Scotch incomplete strategy\n" )
#define OUTPUT_PERSON(_isPTScotch_) ( (_isPTScotch_) ? "      PT-Scotch personal strategy |%s|\n": "      Scotch personal strategy |%s|\n" )

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Generate the ordering strategy string based on the input parameters.
 *
 *******************************************************************************
 *
 * @param[in] iparm
 *          Pointer to the iparm array.
 *
 * @param[in] procnum
 *          Procnum of the process. Output purpose.
 *
 * @param[in] isPTscotch
 *          Boolean that indicates if we use Scotch or PT-Scotch for the order
 *          step.
 *
 ********************************************************************************
 *
 * @retval The strategy string of the order step.
 *
 *******************************************************************************/
char *
order_scotch_build_strategy( const pastix_int_t *iparm,
                             pastix_int_t        procnum,
                             int                 isPTscotch )
{
    char *strat;
    MALLOC_INTERN( strat, STRAT_STR_MAX, char );

    /* Default ordering */
    if ( iparm[IPARM_ORDERING_DEFAULT] == 1 ) {
        if ( iparm[IPARM_INCOMPLETE] == 0 ) {
            if ( iparm[IPARM_VERBOSE] > PastixVerboseNo ) {
                pastix_print( procnum, 0, OUTPUT_DIRECT(isPTscotch) );
            }
            snprintf( strat, STRAT_STR_MAX, STRAT_DIRECT(isPTscotch) );
        }
        else {
            if ( iparm[IPARM_VERBOSE] > PastixVerboseNo ) {
                pastix_print( procnum, 0, OUTPUT_INCOMP(isPTscotch) );
            }
            snprintf( strat, STRAT_STR_MAX, STRAT_INCOMP(isPTscotch) );
        }
    }
    /* Personal ordering */
    else {
        int rc;
        rc = snprintf( strat, STRAT_STR_MAX, STRAT_PERSON(isPTscotch),
                       (long)  iparm[IPARM_SCOTCH_SWITCH_LEVEL],
                       (long)  iparm[IPARM_SCOTCH_CMIN],
                       (long)  iparm[IPARM_SCOTCH_CMAX],
                       ((float)iparm[IPARM_SCOTCH_FRAT])/100.,
                       (long)  iparm[IPARM_SCOTCH_SWITCH_LEVEL],
                       (long)  iparm[IPARM_SCOTCH_CMIN],
                       (long)  iparm[IPARM_SCOTCH_CMAX],
                       ((float)iparm[IPARM_SCOTCH_FRAT])/100. );
        if ( rc > STRAT_STR_MAX ) {
            pastix_print_error( "Order_scotch_build_strategy: Strategy string too long\n" );
            exit(-1);
        }

        if ( iparm[IPARM_VERBOSE] > PastixVerboseNo ) {
            pastix_print( procnum, 0, OUTPUT_PERSON(isPTscotch), strat );
        }
    }

    return strat;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Reallocate the ordering structure.
 *
 * If we decide to drop the Scotch partition to recompute it later, then
 * partition information is freed, otherwise its memory space is compressed.
 *
 *******************************************************************************
 *
 * @param[inout] ordemesh
 *          Pointer to the ordemesh structure to reallocate.
 *
 *******************************************************************************/
void
order_scotch_reallocate_ordemesh( pastix_order_t *ordemesh )
{
#if defined(FORGET_PARTITION)
    ordemesh->cblknbr = 0;
    if ( ordemesh->rangtab != NULL ) {
        memFree_null( ordemesh->rangtab );
    }
    if ( ordemesh->treetab != NULL ) {
        memFree_null( ordemesh->treetab );
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
    if ( ordemesh->cblknbr == 0 ) {
        ordemesh->rangtab[0] = ordemesh->baseval;
    }
#endif
}
