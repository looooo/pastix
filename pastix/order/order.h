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

#if defined(HAVE_SCOTCH)
#include <scotch.h>
#endif
#if defined(HAVE_PTSCOTCH)
#include <ptscotch.h>
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_ordering
 * @struct Order - Ordering structure.
 *
 *******************************************************************************/
typedef struct Order_ {
    pastix_int_t  baseval;   /*< base value used for numbering       +*/
    pastix_int_t  vertnbr;   /*< Number of vertices                  +*/
    pastix_int_t  cblknbr;   /*< Number of column blocks             +*/
    pastix_int_t *permtab;   /*< Permutation array [based]           +*/
    pastix_int_t *peritab;   /*< Inverse permutation array [based]   +*/
    pastix_int_t *rangtab;   /*< Column block range array [based,+1] +*/
#if defined(HAVE_SCOTCH)
    SCOTCH_Graph  grafmesh;  /*< Graph                                             +*/
    int           malgrf;    /*< boolean indicating if grafmesh has been allocated +*/
#endif
} Order;


/*
 * The function prototypes.
 */
int  orderInit (      Order * const ordeptr, pastix_int_t cblknbr, pastix_int_t vertnbr);
void orderExit (      Order * const ordeptr);
int  orderLoad (      Order * const ordeptr, FILE * const stream);
int  orderSave (const Order * const ordeptr, FILE * const stream);
void orderBase (      Order * const ordeptr, pastix_int_t baseval);
int  orderCheck(const Order * const ordeptr);

int orderComputeScotch(pastix_data_t *pastix_data);
int orderComputeMetis( pastix_data_t *pastix_data);
int orderLoadFiles(    pastix_data_t *pastix_data);
int orderSaveFiles(    pastix_data_t *pastix_data);

int orderPrepareCSC(pastix_data_t *pastix_data,
                    pastix_int_t   n,
                    const pastix_int_t  *colptr,
                    const pastix_int_t  *rows,
                    const pastix_int_t  *loc2glob);

void orderFindSupernodes( pastix_int_t  n,
                          pastix_int_t *ia,
                          pastix_int_t *ja,
                          Order * const ordeptr,
                          pastix_int_t *treetab );

int pastix_task_order(pastix_data_t *pastix_data,
                      pastix_int_t   n,
                      pastix_int_t  *colptr,
                      pastix_int_t  *row,
                      pastix_int_t  *loc2glob,
                      pastix_int_t  *perm,
                      pastix_int_t  *invp);

#endif /* _ORDER_H_ */

