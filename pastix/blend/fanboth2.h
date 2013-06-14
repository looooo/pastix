/************************************************************/
/**                                                        **/
/**   NAME       : fanboth.h                               **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                Compute a moderate amalgamation in the  **/
/**                fan-in target to save memory            **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 01 April 2000   **/
/**                                 to     20 April 2000   **/
/**                                                        **/
/************************************************************/
#ifndef FANBOTH
#define static
typedef struct {
  PASTIX_INT ctrbnbr; /*+ number of contribution of the partial ftgt +*/ 
  PASTIX_INT ctrbcnt; 
  PASTIX_INT prionum; /*+ Priority of the partial ftgt +*/
  PASTIX_INT indnum;  /*+ index where the ftgt must be insert in the initial ftgttab +*/
  PASTIX_INT ftgtnum; /*+ Index of the initial ftgt from which is issue the partial ftgt +*/
  PASTIX_INT ftgtnewnum; /*+ index of the ftgt in the final ftgttab +*/
  PASTIX_INT next;    /*+ Chain to the next partial ftgt of the initial ftgt +*/
} ExtraFtgt;
PASTIX_INT Malt2(SolverMatrix *, double);
PASTIX_INT getFtgtInd2(SolverMatrix *, PASTIX_INT *, Queue *, Queue *);
PASTIX_INT getFtgtNextAccess(PASTIX_INT ind, PASTIX_INT ftgtaccessnbr, PASTIX_INT *ftgtaccesstab);
#undef static
#endif
