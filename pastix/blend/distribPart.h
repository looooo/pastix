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
static void       taskExec_DIAG            (pastix_int_t, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       taskExec_E1              (pastix_int_t, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       taskExec_E2              (pastix_int_t, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       taskExec_COMP1D          (pastix_int_t, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static void       computeTaskReceiveTime   (const pastix_int_t, SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);
static double      computeTaskOnProcReadyTime(pastix_int_t, pastix_int_t, SimuCtrl *, BlendCtrl *, Queue *);
static pastix_int_t        getNextTaskNextProc      (SimuCtrl *, BlendCtrl *, pastix_int_t *);
static pastix_int_t        comp_int                 (const pastix_int_t *, const pastix_int_t *);
static pastix_int_t        getNextProc              (SimuProc *, pastix_int_t);
static pastix_int_t        getTaskUnmapped          (Queue *, Queue *, SimuCtrl *);
static pastix_int_t        getCostLowerTask         (SimuCtrl *, BlendCtrl *, pastix_int_t *);
static pastix_int_t        chooseCand               (pastix_int_t, SimuCtrl *, BlendCtrl *);
static void       computeBlockCtrbNbr      (SimuCtrl *, SymbolMatrix *, BlendCtrl *);
static void       updateFtgtStruct         (pastix_int_t, pastix_int_t, pastix_int_t, SymbolMatrix *, SimuCtrl *, BlendCtrl *);
static void       queueReorder             (pastix_int_t, SymbolMatrix *, SimuCtrl *, BlendCtrl *ctrl);
static void       queueCostReorder         (pastix_int_t, SimuCtrl *, BlendCtrl *ctrl);
static void       setClusterTime           (pastix_int_t, SimuCtrl *);
static void       putInReadyQueue          (pastix_int_t, pastix_int_t, SymbolMatrix *, SimuCtrl *, BlendCtrl *);
#undef static
