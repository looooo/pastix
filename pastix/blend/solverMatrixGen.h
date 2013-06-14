/************************************************************/
/**                                                        **/
/**   NAME       : solverMatrixGen.h                       **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Genere the local solver matrix          **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 Oct 1998     **/
/**                                 to     15 Oct 1998     **/
/**                                                        **/
/************************************************************/
#ifndef SOLVERMATRIXGEN_H
#define SOLVERMATRIXGEN_H
PASTIX_INT *               solverMatrixGen (const PASTIX_INT, SolverMatrix *, const SymbolMatrix *, const SimuCtrl *, const BlendCtrl *, const Dof *);
void                allSolverMatrixSave(const char *, const SymbolMatrix *, const SimuCtrl *, const BlendCtrl *, const Dof *);
#endif /* SOLVERMATRIXGEN_H */
