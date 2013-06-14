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
  pastix_int_t treelevel;    /*+ Level of the cblk in the elimination tree (deepness from the root) +*/
  double costlevel; /*+ Cost from root to node +*/
  pastix_int_t fcandnum;     /*+ first processor number of this candidate group  +*/
  pastix_int_t lcandnum;     /*+ last processor number of this candidate group   +*/
  pastix_int_t fccandnum;    /*+ first cluster number of the cluster candidate group +*/
  pastix_int_t lccandnum;    /*+ last cluster number of the cluster candidate group +*/
  pastix_int_t distrib;      /*+ type of the distribution +*/
  pastix_int_t cluster;      /*+ TRUE if cand are clusters (bubble number) +*/
#if defined(TRACE_SOPALIN) || defined(PASTIX_DYNSCHED)
  pastix_int_t cand;         /*+ TRUE if cand are clusters (bubble number) +*/
#endif
} Cand;

