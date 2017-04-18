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
#ifndef BLEND_H
#define BLEND_H

#include "solver.h"

void splitSymbol    ( BlendCtrl    *ctrl,
                      SymbolMatrix *symbmtx );

void propMappTree   ( Cand             *candtab,
                      const EliminTree *etree,
                      pastix_int_t      candnbr,
                      int               nocrossproc,
                      int               allcand );

int  solverMatrixGen( const pastix_int_t,
                            SolverMatrix *,
                      const SymbolMatrix *,
                      const SimuCtrl *,
                      const BlendCtrl *);

#endif /* BLEND_H */

/**
 * @}
 */
