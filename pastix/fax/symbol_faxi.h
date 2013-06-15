/************************************************************/
/**                                                        **/
/**   NAME       : symbol_faxi.h                           **/
/**                                                        **/
/**   AUTHORS    : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the incomplete symbolic             **/
/**                factorization routine.                  **/
/**                                                        **/
/**   DATES      : # Version 1.0  : from : 03 jun 2002     **/
/**                                 to     03 jun 2002     **/
/**                                                        **/
/************************************************************/

#define SYMBOL_FAXI_H

#ifdef CXREF_DOC
#include "common_pastix.h"
#ifndef GRAPH_H
#include "graph.h"
#endif /* GRAPH_H */
#ifndef SYMBOL_H
#include "symbol.h"
#endif /* SYMBOL_H */
#ifndef _ORDER_H_
#include "order.h"
#endif /* ORDER_H */
#ifndef FAX_H
#include "fax.h"
#endif /* FAX_H */
#endif /* CXREF_DOC */

/*
**  The defines.
*/

/* Prime number for hashing vertex numbers. */

#define SYMBOL_FAXI_HASHPRIME       17            /*+ Prime number for hashing +*/

/*
**  The type and structure definitions.
*/

/*+ The chained column block structure. These
    blocks are chained in a single linked list
    for block merge with blocks of left columns. +*/

typedef struct SymbolFaxiTlok_ {
  pastix_int_t                       frownum;              /*+ First row index            +*/
  pastix_int_t                       lrownum;              /*+ Last row index (inclusive) +*/
  pastix_int_t                       cblknum;              /*+ Facing column block        +*/
  pastix_int_t                       levfval;              /*+ Level of fill of block     +*/
  pastix_int_t                       nextnum;              /*+ Index of next block        +*/
} SymbolFaxiTlok;
