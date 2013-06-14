/************************************************************/
/**                                                        **/
/**   NAME       : symbol_cost.h                           **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the symbolic matrix cost computing  **/
/**                routine.                                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 14 oct 1998     **/
/**                                 to     16 oct 1998     **/
/**                                                        **/
/************************************************************/

#define SYMBOL_COST_H

#ifdef CXREF_DOC
#include "common_pastix.h"
#ifndef DOF_H
#include "dof.h"
#endif /* DOF_H */
#ifndef SYMBOL_H
#include "symbol.h"
#endif /* SYMBOL_H */
#endif /* CXREF_DOC */

/*
**  The function prototypes.
*/

#ifndef SYMBOL_COST
#define static
#endif

static void                 symbolCost2         (const SymbolCblk * const cblktax, const SymbolBlok * const bloktax, const Dof * const deofptr, double * const nnzptr, double * const opcptr, const pastix_int_t cblkmin, const pastix_int_t cblknbr);

#undef static
