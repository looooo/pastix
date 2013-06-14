/*+********************************************************+*/
/*+                                                        +*/
/*+   NAME       : symbol_cost.h                           +*/
/*+                                                        +*/
/*+   AUTHORS    : Pascal HENON                            +*/
/*+                                                        +*/
/*+   FUNCTION   : Part of a parallel direct block solver. +*/
/*+                Cost compute functions                  +*/
/*+                                                        +*/
/*+   DATES      : # Version 0.0  : from : 01 jan 1999     +*/
/*+                                 to     20 jan 1999     +*/
/*+                                                        +*/
/*+********************************************************+*/

#ifndef SYMBOL_COST_H
#define static
#endif

void           symbCost      (PASTIX_INT *, double *, const SymbolMatrix *, const Dof *);
double         recursive_sum (PASTIX_INT, PASTIX_INT, double (*fval)(PASTIX_INT, const SymbolMatrix *, const Dof *), 
			      const SymbolMatrix *, const Dof *);
double         cholesky      (PASTIX_INT, const SymbolMatrix *, const Dof *);
double         crout_hyb     (PASTIX_INT, const SymbolMatrix *, const Dof *);
double         crout_2t      (PASTIX_INT, const SymbolMatrix *, const Dof *);
double         crout_3t      (PASTIX_INT, const SymbolMatrix *, const Dof *);
double         nnz           (PASTIX_INT, const SymbolMatrix *, const Dof *);
double         crout_blok    (PASTIX_INT, const SymbolMatrix *, const Dof *);

#undef static
