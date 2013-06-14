/*+********************************************************+*/
/*+                                                        +*/
/*+   NAME       : task.h                                  +*/
/*+                                                        +*/
/*+   AUTHORS    : Pascal HENON                            +*/
/*+                                                        +*/
/*+   FUNCTION   : Part of a parallel direct block solver. +*/
/*+                Structure used in the distribution      +*/
/*+                lead by simulation                      +*/
/*+   DATES      : # Version 0.0  : from : 28 sep 1998     +*/
/*+                                 to     05 oct 1998     +*/
/*+                                                        +*/
/*+********************************************************+*/


#include "common.h"

/*
**  The type and structure definitions.
*/
#ifndef TASK
#define static
#endif
void        taskBuild          (SimuCtrl *, SymbolMatrix *, Cand *, const Dof *, EliminGraph *, BlendCtrl *);
pastix_int_t         getFaceBlockE2     (pastix_int_t, pastix_int_t, pastix_int_t, const SymbolMatrix *, int);
double       taskSendCost       (SimuTask *, const pastix_int_t, const pastix_int_t, BlendCtrl *);  
#undef static
