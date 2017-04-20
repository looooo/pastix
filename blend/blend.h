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
 * @addtogroup pastix_blend
 * @{
 *    This module contains all the subroutines and structures to perform the
 *    analyze step and prepare the numerical factorization and solve. It is
 *    composed of four main steps. The first one is the computation of the
 *    proportional mapping based on the elimination tree to attribute worker
 *    candidates to all nodes in the tree. The second one is the cutting step of
 *    the symbol matrix to generate more parallelism. Then, the simulation
 *    predict the best mapping out of the candidates and returns the associated
 *    static scheduling for the tasks. Finally, the local solver structure is
 *    created to control the numerical factorization and solve and store the
 *    problem data.
 *
 **/
#ifndef BLEND_H
#define BLEND_H

void propMappTree   ( Cand             *candtab,
                      const EliminTree *etree,
                      pastix_int_t      candnbr,
                      int               nocrossproc,
                      int               allcand );
void splitSymbol    ( BlendCtrl    *ctrl,
                      SymbolMatrix *symbmtx );
void simuRun        ( SimuCtrl *,
                      const BlendCtrl *,
                      const SymbolMatrix * );
int  solverMatrixGen( const pastix_int_t,
                            SolverMatrix *,
                      const SymbolMatrix *,
                      const SimuCtrl *,
                      const BlendCtrl *);

#endif /* BLEND_H */

/**
 * @}
 */
