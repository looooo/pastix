/**
 *
 * @file order.h
 *
 * PaStiX order structure routines
 *
 * @copyright 2004-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.1
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @author Gregoire Pichon
 * @author Pierre Ramet
 * @date 2023-07-21
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
#ifndef _pastix_order_h_
#define _pastix_order_h_

#include "pastix/config.h"
#include "pastix/datatypes.h"

BEGIN_C_DECLS

#ifndef DOXYGEN_SHOULD_SKIP_THIS
struct etree_s;
typedef struct etree_s EliminTree;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @brief Order structure.
 *
 * This structure stores the permutation (and inverse permutation) associated to the ordering.
 * It also stores the partitioning tree and the set of supernodes.
 */
typedef struct pastix_order_s {
    pastix_int_t  baseval;     /**< base value used for numbering                             */
    pastix_int_t  vertnbr;     /**< Number of vertices                                        */
    pastix_int_t  cblknbr;     /**< Number of column blocks                                   */
    pastix_int_t *permtab;     /**< Permutation array of size vertnbr [based]                 */
    pastix_int_t *peritab;     /**< Inverse permutation array of size vertnbr [based]         */
    pastix_int_t *rangtab;     /**< Supernode array of size cblknbr+1 [based,+1]              */
    pastix_int_t *treetab;     /**< Partitioning tree of size cblknbr+1 [based]               */
    int8_t       *selevtx;     /**< Selected vertices for low-rank clustering of size cblknbr */
    pastix_int_t  sndenbr;     /**< The number of original supernodes before clustering       */
    pastix_int_t *sndetab;     /**< Original supernode array of size sndenbr [based,+1]       */
    pastix_int_t *peritab_exp; /**< Computed field that should not be modified by the user    */
} pastix_order_t;

/**
 * @name Order basic subroutines
 * @{
 */
int pastixOrderInit( pastix_order_t *ordeptr,
                     pastix_int_t    baseval,
                     pastix_int_t    vertnbr,
                     pastix_int_t    cblknbr,
                     pastix_int_t   *perm,
                     pastix_int_t   *invp,
                     pastix_int_t   *rang,
                     pastix_int_t   *tree );

int pastixOrderAlloc( pastix_order_t *ordeptr,
                      pastix_int_t    vertnbr,
                      pastix_int_t    cblknbr );

int pastixOrderAllocId( pastix_order_t *ordeptr,
                        pastix_int_t    vertnbr );

void pastixOrderExit( pastix_order_t *ordeptr );

void pastixOrderBase( pastix_order_t *ordeptr,
                      pastix_int_t    baseval );

int pastixOrderCheck( const pastix_order_t *ordeptr );

void pastixOrderExpand( pastix_order_t   *ordeptr,
                        const spmatrix_t *spm );

int pastixOrderCopy( pastix_order_t       *ordedst,
                     const pastix_order_t *ordesrc );

pastix_order_t *pastixOrderGet( const pastix_data_t *pastix_data );

void pastixOrderBcast( pastix_order_t *ordemesh,
                       int             root,
                       PASTIX_Comm     pastix_comm );

int pastixOrderGrid( pastix_order_t **myorder,
                     pastix_int_t     nx,
                     pastix_int_t     ny,
                     pastix_int_t     nz );

/**
 * @}
 * @name Order IO subroutines
 * @{
 */
int pastixOrderLoad( const pastix_data_t *pastix_data,
                     pastix_order_t      *ordeptr );

int pastixOrderSave( pastix_data_t        *pastix_data,
                     const pastix_order_t *ordeptr );

/**
 * @}
 */

END_C_DECLS

#endif /* _pastix_order_h_ */

/**
 * @}
 */
