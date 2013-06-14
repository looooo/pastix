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

void           symbCost      (pastix_int_t *, double *, const SymbolMatrix *, const Dof *);
double         recursive_sum (pastix_int_t, pastix_int_t, double (*fval)(pastix_int_t, const SymbolMatrix *, const Dof *), 
			      const SymbolMatrix *, const Dof *);
double         cholesky      (pastix_int_t, const SymbolMatrix *, const Dof *);
double         crout_hyb     (pastix_int_t, const SymbolMatrix *, const Dof *);
double         crout_2t      (pastix_int_t, const SymbolMatrix *, const Dof *);
double         crout_3t      (pastix_int_t, const SymbolMatrix *, const Dof *);
double         nnz           (pastix_int_t, const SymbolMatrix *, const Dof *);
double         crout_blok    (pastix_int_t, const SymbolMatrix *, const Dof *);

#undef static
