/**
 *
 * @file order.h
 *
 *  PaStiX order routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef _ORDER_H_
#define _ORDER_H_

#include "graph.h"

/**
 * @ingroup pastix_ordering
 * @struct Order - Ordering structure.
 */
struct Order_ {
    pastix_int_t  baseval;   /*< base value used for numbering       */
    pastix_int_t  vertnbr;   /*< Number of vertices                  */
    pastix_int_t  cblknbr;   /*< Number of column blocks             */
    pastix_int_t *permtab;   /*< Permutation array [based]           */
    pastix_int_t *peritab;   /*< Inverse permutation array [based]   */
    pastix_int_t *rangtab;   /*< Column block range array [based,+1] */
    pastix_int_t *treetab;   /*< Partitioning tree [based]           */
};

/*
 * The function prototypes.
 */
int  orderInit (      Order * const ordeptr, pastix_int_t cblknbr, pastix_int_t vertnbr);
void orderExit (      Order * const ordeptr);
void orderBase (      Order * const ordeptr, pastix_int_t baseval);
int  orderCheck(const Order * const ordeptr);

int  orderComputeScotch(   pastix_data_t *pastix_data, const pastix_graph_t *graph );
int  orderComputePTScotch( pastix_data_t *pastix_data, const pastix_graph_t *graph );
int  orderComputeMetis(    pastix_data_t *pastix_data, const pastix_graph_t *graph );
int  orderComputeParMetis( pastix_data_t *pastix_data, const pastix_graph_t *graph );
int  orderComputeOptimal(  pastix_data_t *pastix_data, pastix_int_t n );

int  orderLoad(       Order * const ordeptr, char *filename );
int  orderSave( const Order * const ordeptr, char *filename );

int  orderPrepareCSC(pastix_data_t *pastix_data,
                     pastix_int_t   n,
                     const pastix_int_t  *colptr,
                     const pastix_int_t  *rows,
                     const pastix_int_t  *loc2glob);

void orderFindSupernodes( const pastix_graph_t *graph,
                          Order * const ordeptr );

int  orderApplyLevelOrder( Order *ordeptr, pastix_int_t distribution_level );

int  orderAddIsolate( Order              *ordeptr,
                      pastix_int_t        new_n,
                      const pastix_int_t *perm );

void order_grid2D_wide(pastix_int_t *rangtab,
                       pastix_int_t *peritab,
                       pastix_int_t *cblknbr,
                       pastix_int_t x0,
                       pastix_int_t xn,
                       pastix_int_t y0,
                       pastix_int_t yn,
                       pastix_int_t max_number,
                       pastix_int_t lda,
                       pastix_int_t *current_rangtab);

void order_grid3D_wide(pastix_int_t *rangtab,
                       pastix_int_t *peritab,
                       pastix_int_t *cblknbr,
                       pastix_int_t x0,
                       pastix_int_t xn,
                       pastix_int_t y0,
                       pastix_int_t yn,
                       pastix_int_t z0,
                       pastix_int_t zn,
                       pastix_int_t *max_number,
                       pastix_int_t lda,
                       pastix_int_t *current_rangtab,
                       pastix_int_t *treetab,
                       pastix_int_t current_treetab);

void order_grid3D_classic(pastix_int_t *rangtab,
                          pastix_int_t *peritab,
                          pastix_int_t *cblknbr,
                          pastix_int_t x0,
                          pastix_int_t xn,
                          pastix_int_t y0,
                          pastix_int_t yn,
                          pastix_int_t z0,
                          pastix_int_t zn,
                          pastix_int_t *max_number,
                          pastix_int_t lda,
                          pastix_int_t *current_rangtab,
                          pastix_int_t *treetab,
                          pastix_int_t current_treetab);


#endif /* _ORDER_H_ */

