/**
 *
 * @file blend.h
 *
 * PaStiX analyse functions to manipulate candidates on the elimination tree
 * structure.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @addtogroup blend_dev_cand
 * @{
 *    This module contains all subroutines to initialize the candidates array
 *    for each supernode, as well as supernode properties that are defined by
 *    level such as 2D layouts and 2D tasks.
 *
 **/
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
void splitSymbol( BlendCtrl    *ctrl,
                  SymbolMatrix *symbmtx );

void propMappTree( Cand               *candtab,
                   const EliminTree   *etree,
                   pastix_int_t        candnbr,
                   int nocrossproc, int allcand );

int  solverMatrixGen (const pastix_int_t,
                      SolverMatrix *,
                      const SymbolMatrix *,
                      const SimuCtrl *,
                      const BlendCtrl *);

#endif /* BLEND_H */
