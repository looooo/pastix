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
** $Id: symbol_faxi_graph.c 285 2005-03-10 10:25:31Z pelegrin $
*/
/************************************************************/
/**                                                        **/
/**   NAME       : symbol_faxi_graph.c                     **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This is the block symbolic incomplete   **/
/**                factorization routine for generic       **/
/**                graphs.                                 **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 05 jun 2002     **/
/**                                 to     05 jun 2002     **/
/**                # Version 1.1  : from : 26 jun 2002     **/
/**                                 to     26 jun 2002     **/
/**                # Version 3.0  : from : 02 mar 2004     **/
/**                                 to     15 dec 2004     **/
/**                                                        **/
/**   NOTES      : # symbolFaxiGraph() could have called   **/
/**                  symbolFaxi() in the regular way.      **/
/**                  However, for efficiency reasons, we   **/
/**                  have decided to inline symbolFaxi(),  **/
/**                  to avoid a function call for every    **/
/**                  arc.                                  **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#define SYMBOL_FAXI
#define SYMBOL_FAXI_INCLUDED

#include "common.h"
#ifdef WITH_SCOTCH
#ifdef DISTRIBUTED
#include "ptscotch.h"
#else
#include "scotch.h"
#endif
#include "symbol.h"
#include "order.h"
#include "fax.h"
#include "symbol_faxi.h"

/***********************************/
/*                                 */
/* Symbolic factorization routine. */
/*                                 */
/***********************************/

/*+ This routine computes the block symbolic
*** factorization of the given matrix graph
*** according to the given vertex ordering.
*** It returns:
*** - 0   : on success.
*** - !0  : on error.
+*/

int
symbolFaxiGraph (
SymbolMatrix       * const  symbptr, /*+ Symbolic block matrix [based]              +*/
const SCOTCH_Graph * const  grafptr, /*+ Matrix adjacency structure [based]         +*/
const Order        * const  ordeptr, /*+ Matrix ordering                            +*/
const pastix_int_t                   levfmax) /*+ Inclusive maximum level of fill for blocks +*/
{
  pastix_int_t                   baseval;
  pastix_int_t                   vertnbr;
  pastix_int_t *                 verttab;
  const pastix_int_t * restrict  verttax;
  pastix_int_t                   edgenbr;
  pastix_int_t                   edgenum;
  pastix_int_t *                 edgetab;
  const pastix_int_t * restrict  edgetax;

  SCOTCH_graphData (grafptr, 
		    (SCOTCH_Num *)&baseval, 
		    (SCOTCH_Num *)&vertnbr, 
		    (SCOTCH_Num **)&verttab, 
		    NULL, NULL, NULL, 
		    (SCOTCH_Num *)&edgenbr, 
		    (SCOTCH_Num **)&edgetab, 
		    NULL);
  verttax = verttab - baseval;
  edgetax = edgetab - baseval;

#define SYMBOL_FAXI_ITERATOR(ngbdptr, vertnum, vertend) \
                                    for (edgenum = verttax[vertnum];     \
                                         edgenum < verttax[vertnum + 1]; \
                                         edgenum ++) {                   \
                                      vertend = edgetax[edgenum];

#define SYMBOL_FAXI_VERTEX_DEGREE(ngbdptr, vertnum) \
                                    (verttax[(vertnum) + 1] - verttax[(vertnum)])

  {
#include "symbol_faxi.c"
  }
}
#endif /* WITH_SCOTCH */
