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


void          costMatrixCorrect     (CostMatrix *, const SymbolMatrix *, Cand * candtab,  const Dof *);
double          subtreeUpdateCost     (pastix_int_t, CostMatrix *, const EliminTree *);
double          subtreeUpdateCostLocal(pastix_int_t, const BlendCtrl *, const SymbolMatrix *, const SimuCtrl *, const Dof *, pastix_int_t);
double          cblkComputeCost2D     (pastix_int_t, CostMatrix *, const SymbolMatrix *, const Dof *);

/** 2D **/
double          DIAGCost              (pastix_int_t);
double          E1Cost                (pastix_int_t, pastix_int_t);
double          E2Cost                (pastix_int_t, pastix_int_t, pastix_int_t);

static double   computeCost           (pastix_int_t, pastix_int_t);
static double   contribCompCost       (pastix_int_t, pastix_int_t, pastix_int_t);
static double   contribAddCost        (pastix_int_t, pastix_int_t);
double costFtgtSend( const BlendCtrl   *ctrl,
                     const Dof         *dofptr,
                     const d_FanInTarget *ftgt,
                     pastix_int_t clustsrc,
                     pastix_int_t sync_comm_nbr );


double          costFtgtAdd           (d_FanInTarget *, const Dof *);
double          cblkMaxCost           (pastix_int_t, const CostMatrix *);
double          totalCost             (pastix_int_t, const CostMatrix *);
void            printSolverInfo       (FILE *, const d_SolverMatrix *, const SymbolMatrix *, const Dof * const dofptr);
double          memorySpaceCost       (const d_SolverMatrix *);
static double   solverSpaceCost       (const d_SolverMatrix *);
static double   symbolSpaceCost       (const SymbolMatrix *);

#undef static
