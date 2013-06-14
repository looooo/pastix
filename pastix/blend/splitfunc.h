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
PASTIX_INT           splitSeqCblk2D    (PASTIX_INT, PASTIX_INT, const SymbolMatrix *, const ExtraSymbolMatrix *, const Dof *, const BlendCtrl *,
			         PASTIX_INT (*P)(PASTIX_INT , PASTIX_INT, const SymbolMatrix *, const ExtraSymbolMatrix *, const Dof *, PASTIX_INT, const BlendCtrl *));
static void   virtualSplit      (PASTIX_INT, PASTIX_INT, const SymbolBlok *, PASTIX_INT *, SymbolBlok *);
static double  cblkCost          (PASTIX_INT, const SymbolBlok *, const Dof *);
static void   build_cblk        (PASTIX_INT, const SymbolMatrix *, const ExtraSymbolMatrix *, SymbolBlok *);
static PASTIX_INT    cblkNbr           (PASTIX_INT,  const SymbolMatrix *, const ExtraSymbolMatrix *);
PASTIX_INT           P1D               (PASTIX_INT, PASTIX_INT, const SymbolMatrix *, const ExtraSymbolMatrix *const extrasymbptr, const Dof *, PASTIX_INT, const BlendCtrl *);
PASTIX_INT           P2D               (PASTIX_INT, PASTIX_INT, const SymbolMatrix *, const ExtraSymbolMatrix *const extrasymbptr, const Dof *, PASTIX_INT, const BlendCtrl *);
#undef static


