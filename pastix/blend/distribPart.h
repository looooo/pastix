/************************************************************/
/**                                                        **/
/**   NAME       : distribPart.h                           **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Distribution of the symbolic matrix     **/
/**                lead by simulation.                     **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 27 sep 1998     **/
/**                                 to     04 oct 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The function prototypes.
*/

#ifndef DISTRIB
#define static
#endif

void              distribPart              (SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       taskExec_DIAG            (PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       taskExec_E1              (PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       taskExec_E2              (PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       taskExec_COMP1D          (PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       computeTaskReceiveTime   (const PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static double      computeTaskOnProcReadyTime(PASTIX_INT, PASTIX_INT, SimuCtrl *, BlendCtrl *, Queue *);
static PASTIX_INT        getNextTaskNextProc      (SimuCtrl *, BlendCtrl *, PASTIX_INT *);
static PASTIX_INT        comp_int                 (const PASTIX_INT *, const PASTIX_INT *);
static PASTIX_INT        getNextProc              (SimuProc *, PASTIX_INT);
static PASTIX_INT        getTaskUnmapped          (Queue *, Queue *, SimuCtrl *);
static PASTIX_INT        getCostLowerTask         (SimuCtrl *, BlendCtrl *, PASTIX_INT *);
static PASTIX_INT        chooseCand               (PASTIX_INT, SimuCtrl *, BlendCtrl *);
static void       computeBlockCtrbNbr      (SimuCtrl *, SymbolMatrix *, BlendCtrl *);
static void       updateFtgtStruct         (PASTIX_INT, PASTIX_INT, PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *);
static void       queueReorder             (PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *ctrl);
static void       queueCostReorder         (PASTIX_INT, SimuCtrl *, BlendCtrl *ctrl);
static void       setClusterTime           (PASTIX_INT, SimuCtrl *);
static void       putInReadyQueue          (PASTIX_INT, PASTIX_INT, SymbolMatrix *, SimuCtrl *, BlendCtrl *);
#undef static
