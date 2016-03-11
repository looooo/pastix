/************************************************************/
/**                                                        **/
/**   NAME       : blend.h                                 **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Partition and distribute data           **/
/**                for an optimal parallel resolution      **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to                     **/
/**                                                        **/
/************************************************************/
#ifndef BLEND_H
#define BLEND_H

#include "solver.h"

/*
 * Function: solverBlend
 *
 * Main blend function
 *
 * Build the elimination graph from the symbolic partition.
 *
 * Build the cost matrix from the symbolic partition.
 *
 * Build the elimination tree from the symbolic partition.
 *
 * Distribute each column bloc on candidate processors.
 *
 * Build a new symbol matrix...
 *
 * Parameters:
 *   solvmtx    - Solver matrix structure.
 *   symbmtx    - Symbol matrix
 *   assemb1D   -
 *   assemb2D   -
 *   clustnbr   - Number of MPI processes.
 *   local_nbthrds - Number of threads.
 *   clustnum   - Processor ID number.
 *   option     - Blend parameters.
 */
void solverBlend(BlendCtrl    *ctrl,
                 SolverMatrix *solvmtx,
                 SymbolMatrix *symbmtx);

void splitSymbol( BlendCtrl    *ctrl,
                  SymbolMatrix *symbmtx );

void propMappTree( Cand               *candtab,
                   const EliminTree   *etree,
                   pastix_int_t        candnbr,
                   int nocrossproc, int allcand );

#endif /* BLEND_H */
