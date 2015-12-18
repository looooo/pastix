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


void            costMatrixCorrect     (CostMatrix *, const SymbolMatrix *, Cand *);
double          subtreeUpdateCost     (pastix_int_t, CostMatrix *, const EliminTree *);
double          subtreeUpdateCostLocal(pastix_int_t, const BlendCtrl *, const SymbolMatrix *, const SimuCtrl *, pastix_int_t);
double          cblkComputeCost2D     (pastix_int_t, CostMatrix *, const SymbolMatrix *);

/** 2D **/
double          DIAGCost              (pastix_int_t);
double          E1Cost                (pastix_int_t, pastix_int_t);
double          E2Cost                (pastix_int_t, pastix_int_t, pastix_int_t);

static double   computeCost           (pastix_int_t, pastix_int_t);
static double   contribCompCost       (pastix_int_t, pastix_int_t, pastix_int_t);
static double   contribAddCost        (pastix_int_t, pastix_int_t);
double costFtgtSend( const BlendCtrl   *ctrl,
                     const FanInTarget *ftgt,
                     pastix_int_t clustsrc,
                     pastix_int_t sync_comm_nbr );


double          costFtgtAdd           (FanInTarget *);
double          cblkMaxCost           (pastix_int_t, const CostMatrix *);
double          totalCost             (pastix_int_t, const CostMatrix *);
void            printSolverInfo       (FILE *, const SolverMatrix *);
double          memorySpaceCost       (const SolverMatrix *);
static double   solverSpaceCost       (const SolverMatrix *);
static double   symbolSpaceCost       (const SymbolMatrix *);

#undef static
