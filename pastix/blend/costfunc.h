/*+********************************************************+*/
/*+                                                        +*/
/*+   NAME       : costfunc.h                              +*/
/*+                                                        +*/
/*+   AUTHORS    : Pascal HENON                            +*/
/*+                                                        +*/
/*+   FUNCTION   : Part of a parallel direct block solver. +*/
/*+                Cost compute functions                  +*/
/*+                                                        +*/
/*+   DATES      : # Version 0.0  : from : 27 sep 1998     +*/
/*+                                 to     03 oct 1998     +*/
/*+                                                        +*/
/*+********************************************************+*/

#ifndef COSTFUNC_H
#define static
#endif


void            costMatrixBuild       (CostMatrix *, const SymbolMatrix *, const Dof *);
void            costMatrixCorrect     (CostMatrix *, const SymbolMatrix *, Cand * candtab,  const Dof *);
double          subtreeUpdateCost     (PASTIX_INT, CostMatrix *, const EliminTree *);
double          subtreeUpdateCostLocal(PASTIX_INT, const BlendCtrl *, const SymbolMatrix *, const SimuCtrl *, const Dof *, PASTIX_INT);
double          cblkComputeCost       (PASTIX_INT, CostMatrix *, const SymbolMatrix *, const Dof *);
double          cblkComputeCost2D     (PASTIX_INT, CostMatrix *, const SymbolMatrix *, const Dof *);
				 		   
/** 2D **/			 	   
double          DIAGCost              (PASTIX_INT);
double          E1Cost                (PASTIX_INT, PASTIX_INT);
double          E2Cost                (PASTIX_INT, PASTIX_INT, PASTIX_INT);
				 		   
static double   computeCost           (PASTIX_INT, PASTIX_INT);
static double   contribCompCost       (PASTIX_INT, PASTIX_INT, PASTIX_INT);
static double   contribAddCost        (PASTIX_INT, PASTIX_INT);
double          costFtgtSend          (PASTIX_INT, PASTIX_INT, FanInTarget *, BlendCtrl *,  const Dof *);
				 		   
double          costFtgtAdd           (FanInTarget *, const Dof *);
double          cblkMaxCost           (PASTIX_INT, const CostMatrix *);
double          totalCost             (PASTIX_INT, const CostMatrix *);
void            printSolverInfo       (FILE *, const SolverMatrix *, const SymbolMatrix *, const Dof * const dofptr);
double          memorySpaceCost       (const SolverMatrix *);
static double   solverSpaceCost       (const SolverMatrix *);
static double   symbolSpaceCost       (const SymbolMatrix *);

#undef static




