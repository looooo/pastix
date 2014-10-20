/************************************************************/
/**                                                        **/
/**   NAME       : solver_io.c                             **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the input/output        **/
/**                routines for symbolic matrices.         **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 23 aug 1998     **/
/**                                 to     08 dec 1998     **/
/**                                                        **/
/************************************************************/

#define static

pastix_int_t solverLoad(SolverMatrix *, FILE *);
pastix_int_t solverSave(const SolverMatrix *, FILE *);


#undef static







