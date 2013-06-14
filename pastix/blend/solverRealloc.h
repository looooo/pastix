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

#ifndef SOLVER_REALLOC
#define static
#endif

void                     solverRealloc        (SolverMatrix *, pastix_int_t*);
void                     solverExit           (SolverMatrix *);
void                     solverInit           (SolverMatrix *);
void                     setBcofPtr           (SolverMatrix *, const pastix_int_t *);
void                     setLocalBtagPtr      (SolverMatrix *);
#undef static
