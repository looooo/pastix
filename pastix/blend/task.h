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


#include "common_pastix.h"

/*
**  The type and structure definitions.
*/
#ifndef TASK
#define static
#endif
void        taskBuild          (SimuCtrl *, SymbolMatrix *, Cand *, const Dof *, EliminGraph *, BlendCtrl *);
PASTIX_INT         getFaceBlockE2     (PASTIX_INT, PASTIX_INT, PASTIX_INT, const SymbolMatrix *, int);
double       taskSendCost       (SimuTask *, const PASTIX_INT, const PASTIX_INT, BlendCtrl *);  
#undef static
