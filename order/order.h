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
#ifndef _ORDER_H_
#define _ORDER_H_

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
} Order;

/**
 * @name Order basic subroutines
 * @{
 */
int  orderInit (      Order * const ordeptr, pastix_int_t baseval, pastix_int_t cblknbr, pastix_int_t vertnbr,
                      pastix_int_t *perm, pastix_int_t *peri, pastix_int_t *rang, pastix_int_t *tree );
int  orderAlloc(      Order * const ordeptr, pastix_int_t cblknbr, pastix_int_t vertnbr);
void orderExit (      Order * const ordeptr);
void orderBase (      Order * const ordeptr, pastix_int_t baseval);
int  orderCheck(const Order * const ordeptr);
int  orderCopy (      Order * const ordedst, const Order * const ordesrc);

/**
 * @}
 * @name Order IO subroutines
 * @{
 */
int  orderLoad( pastix_data_t *pastix_data,       Order * const ordeptr );
int  orderSave( pastix_data_t *pastix_data, const Order * const ordeptr );

/**
 * @}
 * @name Order compute subroutines
 * @{
 */
int  orderComputeScotch(   pastix_data_t *pastix_data, pastix_graph_t *graph );
int  orderComputePTScotch( pastix_data_t *pastix_data, pastix_graph_t *graph );
int  orderComputeMetis(    pastix_data_t *pastix_data, pastix_graph_t *graph );
int  orderComputeParMetis( pastix_data_t *pastix_data, pastix_graph_t *graph );

/**
 * @}
 * @name Order manipulation subroutines
 * @{
 */
void orderFindSupernodes( const pastix_graph_t *graph,
                          Order * const ordeptr );

int  orderApplyLevelOrder( Order *ordeptr, pastix_int_t distribution_level );

int  orderAddIsolate( Order              *ordeptr,
                      pastix_int_t        new_n,
                      const pastix_int_t *perm );

/**
 * @}
 */

#endif /* _ORDER_H_ */

/**
 * @}
 */
