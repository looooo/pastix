/**
 *
 * @file order.h
 *
 * PaStiX order structure routines
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 *
 * @addtogroup pastix_order
 * @{
 *   @brief Functions to generate and manipulate the order structure.
 *
 *   This module provides the set of function to prepare the order structure
 *   associated to a given sparse matrix. It is possible to call Scotch,
 *   PT-Scotch, Metis and ParMetis to build a new ordering that minimize the
 *   fill-in and maximize the level of parallelism.
 *
 **/
#ifndef _PASTIX_ORDER_H_
#define _PASTIX_ORDER_H_

#include "graph.h"

/**
 * @brief Order structure.
 *
 * This structure stores the permutation (and inverse permutation) associated to the ordering.
 * It also stores the partitioning tree and the set of supernodes.
 */
typedef struct pastix_order_s {
    pastix_int_t  baseval;   /**< base value used for numbering       */
    pastix_int_t  vertnbr;   /**< Number of vertices                  */
    pastix_int_t  cblknbr;   /**< Number of column blocks             */
    pastix_int_t *permtab;   /**< Permutation array [based]           */
    pastix_int_t *peritab;   /**< Inverse permutation array [based]   */
    pastix_int_t *rangtab;   /**< Supernode array [based,+1]          */
    pastix_int_t *treetab;   /**< Partitioning tree [based]           */
} pastix_order_t;

/**
 * @name Order basic subroutines
 * @{
 */
int  pastixOrderInit (      pastix_order_t * const ordeptr, pastix_int_t baseval, pastix_int_t cblknbr, pastix_int_t vertnbr,
                            pastix_int_t *perm, pastix_int_t *peri, pastix_int_t *rang, pastix_int_t *tree );
int  pastixOrderAlloc(      pastix_order_t * const ordeptr, pastix_int_t cblknbr, pastix_int_t vertnbr);
void pastixOrderExit (      pastix_order_t * const ordeptr);
void pastixOrderBase (      pastix_order_t * const ordeptr, pastix_int_t baseval);
int  pastixOrderCheck(const pastix_order_t * const ordeptr);
int  pastixOrderCopy (      pastix_order_t * const ordedst, const pastix_order_t * const ordesrc);

/**
 * @}
 * @name Order IO subroutines
 * @{
 */
int  pastixOrderLoad( const pastix_data_t *pastix_data,       pastix_order_t *ordeptr );
int  pastixOrderSave(       pastix_data_t *pastix_data, const pastix_order_t *ordeptr );

/**
 * @}
 * @name Order compute subroutines
 * @{
 */
int  pastixOrderComputeScotch(   pastix_data_t *pastix_data, pastix_graph_t *graph );
int  pastixOrderComputePTScotch( pastix_data_t *pastix_data, pastix_graph_t *graph );
int  pastixOrderComputeMetis(    pastix_data_t *pastix_data, pastix_graph_t *graph );
int  pastixOrderComputeParMetis( pastix_data_t *pastix_data, pastix_graph_t *graph );

/**
 * @}
 * @name Order manipulation subroutines
 * @{
 */
void pastixOrderFindSupernodes( const pastix_graph_t *graph,
                                pastix_order_t * const ordeptr );

int  pastixOrderApplyLevelOrder( pastix_order_t *ordeptr,
                                 pastix_int_t    distribution_level );

int  pastixOrderAddIsolate( pastix_order_t     *ordeptr,
                            pastix_int_t        new_n,
                            const pastix_int_t *perm );

/**
 * @}
 */

#endif /* _PASTIX_ORDER_H_ */

/**
 * @}
 */
