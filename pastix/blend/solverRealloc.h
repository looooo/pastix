/************************************************************/
/**                                                        **/
/**   NAME       : solveRealloc.h                          **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Realloc in a contiguous way             **/
/**                the final solver  matrix                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 27 Nov 1998     **/
/**                                 to     12 Dec 1998     **/
/**                                                        **/
/************************************************************/

#ifndef SOLVER_REALLOC_H
#define SOLVER_REALLOC_H
#include "solver.h"
/* hack to call the function in a precision generated file */
void                     solverRealloc        (SolverMatrix *);

#endif /* SOLVER_REALLOC_H */
