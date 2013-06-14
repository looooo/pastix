/************************************************************/
/**                                                        **/
/**   NAME       : cand.h                                  **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the candidate group of a cblk       **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 sep 1998     **/
/**                                 to     08 oct 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/
#define D1 0
#define D2 1
#define DENSE 3

/*+ Processor candidate group to own a column blok      +*/
typedef struct Cand_{
  PASTIX_INT treelevel;    /*+ Level of the cblk in the elimination tree (deepness from the root) +*/
  double costlevel; /*+ Cost from root to node +*/
  PASTIX_INT fcandnum;     /*+ first processor number of this candidate group  +*/
  PASTIX_INT lcandnum;     /*+ last processor number of this candidate group   +*/
  PASTIX_INT fccandnum;    /*+ first cluster number of the cluster candidate group +*/
  PASTIX_INT lccandnum;    /*+ last cluster number of the cluster candidate group +*/
  PASTIX_INT distrib;      /*+ type of the distribution +*/
  PASTIX_INT cluster;      /*+ TRUE if cand are clusters (bubble number) +*/
#if defined(TRACE_SOPALIN) || defined(PASTIX_DYNSCHED)
  PASTIX_INT cand;         /*+ TRUE if cand are clusters (bubble number) +*/
#endif
} Cand;

