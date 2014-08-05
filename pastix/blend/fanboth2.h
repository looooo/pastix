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
  pastix_int_t ctrbnbr; /*+ number of contribution of the partial ftgt +*/ 
  pastix_int_t ctrbcnt; 
  pastix_int_t prionum; /*+ Priority of the partial ftgt +*/
  pastix_int_t indnum;  /*+ index where the ftgt must be insert in the initial ftgttab +*/
  pastix_int_t ftgtnum; /*+ Index of the initial ftgt from which is issue the partial ftgt +*/
  pastix_int_t ftgtnewnum; /*+ index of the ftgt in the final ftgttab +*/
  pastix_int_t next;    /*+ Chain to the next partial ftgt of the initial ftgt +*/
} ExtraFtgt;
pastix_int_t Malt2(d_SolverMatrix *, double);
pastix_int_t getFtgtInd2(d_SolverMatrix *, pastix_int_t *, Queue *, Queue *);
pastix_int_t getFtgtNextAccess(pastix_int_t ind, pastix_int_t ftgtaccessnbr, pastix_int_t *ftgtaccesstab);
#undef static
#endif
