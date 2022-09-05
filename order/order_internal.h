/**
 *
 * @file order_internal.h
 *
 * PaStiX order internal routines
 *
 * @copyright 2004-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @date 2021-06-28
 *
 **/
#ifndef _order_internal_h_
#define _order_internal_h_

#include "pastix/order.h"

BEGIN_C_DECLS

#define orderDrawGraph       ( 0x1 << 0 )
#define orderDrawCoordinates ( 0x1 << 1 )
#define orderDrawMapping     ( 0x1 << 2 )

/**
 * @}
 * @name Order compute subroutines
 * @{
 */
int orderComputeScotch(   pastix_data_t *pastix_data, pastix_graph_t *graph );
int orderComputePTScotch( pastix_data_t *pastix_data, pastix_graph_t *graph );
int orderComputeMetis(    pastix_data_t *pastix_data, pastix_graph_t *graph );
int orderComputeParMetis( pastix_data_t *pastix_data, pastix_graph_t *graph );
int orderComputePersonal( pastix_data_t *pastix_data, pastix_graph_t *graph, pastix_order_t *myorder );

/**
 * @}
 * @name Order manipulation subroutines
 * @{
 */
void orderFindSupernodes( const pastix_graph_t *graph,
                                pastix_order_t * const ordeptr );

int  orderAmalgamate( int             verbose,
                      int             ilu,
                      int             levelk,
                      int             rat_cblk,
                      int             rat_blas,
                      pastix_graph_t *graph,
                      pastix_order_t *orderptr,
                      PASTIX_Comm     pastix_comm );

int  orderApplyLevelOrder( pastix_order_t *ordeptr,
                           pastix_int_t    level_tasks2d,
                           pastix_int_t    width_tasks2d );

int  orderAddIsolate( pastix_order_t     *ordeptr,
                      pastix_int_t        new_n,
                      const pastix_int_t *perm );

pastix_int_t orderSupernodes( const pastix_graph_t *graph,
                              pastix_order_t       *order,
                              EliminTree           *etree,
                              pastix_int_t         *iparm,
                              int                   do_schur );

pastix_int_t *orderGetExpandedPeritab( pastix_order_t   *ordeptr,
                                       const spmatrix_t *spm );

/**
 * @}
 * @name Order manipulation subroutines
 * @{
 */
void
orderDraw( pastix_data_t *pastix_data,
           const char    *filename,
           pastix_int_t   sndeidx,
           int            dump );

EliminTree *orderBuildEtree( const pastix_order_t *order );

char *order_scotch_build_strategy( const pastix_int_t *iparm,
                                   pastix_int_t        procnum,
                                   int                 isPTscotch );
void  order_scotch_reallocate_ordemesh( pastix_order_t *ordemesh );

END_C_DECLS

#endif /* _order_internal_h_ */

/**
 * @}
 */
