/*+********************************************************+*/
/*+                                                        +*/
/*+   NAME       : splifunc.h                              +*/
/*+                                                        +*/
/*+   AUTHORS    : Pascal HENON                            +*/
/*+                                                        +*/
/*+   FUNCTION   : Part of a parallel direct block solver. +*/
/*+                Compute optimal split                   +*/
/*+                                                        +*/
/*+   DATES      : # Version 0.0  : from : 01 jan 1999     +*/
/*+                                 to     20 jan 1999     +*/
/*+                                                        +*/
/*+********************************************************+*/

#ifndef SPLITFUNC_H
#define static
#endif
pastix_int_t           splitSeqCblk2D    (pastix_int_t, pastix_int_t, const SymbolMatrix *, const ExtraSymbolMatrix *, const Dof *, const BlendCtrl *,
			         pastix_int_t (*P)(pastix_int_t , pastix_int_t, const SymbolMatrix *, const ExtraSymbolMatrix *, const Dof *, pastix_int_t, const BlendCtrl *));
static void   virtualSplit      (pastix_int_t, pastix_int_t, const SymbolBlok *, pastix_int_t *, SymbolBlok *);
static double  cblkCost          (pastix_int_t, const SymbolBlok *, const Dof *);
static void   build_cblk        (pastix_int_t, const SymbolMatrix *, const ExtraSymbolMatrix *, SymbolBlok *);
static pastix_int_t    cblkNbr           (pastix_int_t,  const SymbolMatrix *, const ExtraSymbolMatrix *);
pastix_int_t           P1D               (pastix_int_t, pastix_int_t, const SymbolMatrix *, const ExtraSymbolMatrix *const extrasymbptr, const Dof *, pastix_int_t, const BlendCtrl *);
pastix_int_t           P2D               (pastix_int_t, pastix_int_t, const SymbolMatrix *, const ExtraSymbolMatrix *const extrasymbptr, const Dof *, pastix_int_t, const BlendCtrl *);
#undef static


