/* Copyright INRIA 2004
**
** This file is part of the Scotch distribution.
**
** The Scotch distribution is libre/free software; you can
** redistribute it and/or modify it under the terms of the
** GNU Lesser General Public License as published by the
** Free Software Foundation; either version 2.1 of the
** License, or (at your option) any later version.
**
** The Scotch distribution is distributed in the hope that
** it will be useful, but WITHOUT ANY WARRANTY; without even
** the implied warranty of MERCHANTABILITY or FITNESS FOR A
** PARTICULAR PURPOSE. See the GNU Lesser General Public
** License for more details.
**
** You should have received a copy of the GNU Lesser General
** Public License along with the Scotch distribution; if not,
** write to the Free Software Foundation, Inc.,
** 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
**
** $Id: order.h 2 2004-06-02 14:05:03Z ramet $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : order.h                                 **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                Mathieu FAVERGE                         **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the graph ordering routine.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 aug 1998     **/
/**                                 to     01 may 1999     **/
/**                # Version 2.0  : from : 25 oct 2003     **/
/**                                 to     22 apr 2004     **/
/**                                                        **/
/************************************************************/
#ifndef ORDER_H
#define ORDER_H

/*
**  The type and structure definitions.
*/

/*+ Ordering structure. vnodbas holds the base
    value for node indexings. vnodbas is equal
    to baseval for graphs, and to vnodbas for
    meshes. The same holds for rangtab, with
    rangtab[0] = vnodbas.                      +*/

typedef struct Order_ {
  pastix_int_t  cblknbr;              /*+ Number of column blocks             +*/
  pastix_int_t *rangtab;              /*+ Column block range array [based,+1] +*/
  pastix_int_t *permtab;              /*+ Permutation array [based]           +*/
  pastix_int_t *peritab;              /*+ Inverse permutation array [based]   +*/
} Order;

/*
**  The function prototypes.
*/

int  orderInit (      Order *          const ordeptr, pastix_int_t cblknbr, pastix_int_t vertnbr);
void orderExit (      Order *          const ordeptr);
int  orderLoad (      Order *          const ordeptr, FILE * const stream);
int  orderSave (const Order *          const ordeptr, FILE * const stream);
void orderBase (      Order * restrict const ordeptr, const pastix_int_t baseval);
int  orderCheck(const Order *          const ordeptr);

int orderComputeScotch(pastix_data_t *pastix_data);
int orderComputeMetis( pastix_data_t *pastix_data);
int orderLoadFiles(    pastix_data_t *pastix_data);
int orderSaveFiles(    pastix_data_t *pastix_data);

int pastix_task_order(pastix_data_t *pastix_data,
                      pastix_int_t  *perm,
                      pastix_int_t  *invp);



int pastix_task_scotch(pastix_data_t **pastix_data,
                       MPI_Comm        pastix_comm,
                       pastix_int_t             n,
                       pastix_int_t            *colptr,
                       pastix_int_t            *row,
                       pastix_int_t            *perm,
                       pastix_int_t            *invp);
int dpastix_task_scotch(pastix_data_t ** pastix_data,
                        MPI_Comm         pastix_comm,
                        pastix_int_t              n,
                        pastix_int_t            * colptr,
                        pastix_int_t            * row,
                        pastix_int_t            * perm,
                        pastix_int_t            * invp,
                        pastix_int_t            * loc2glob);

#endif /* ORDER_H */

