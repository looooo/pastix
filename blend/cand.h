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

#ifndef _CAND_H_
#define _CAND_H_

/*
**  The type and structure definitions.
*/
#define D1 0
#define D2 1
#define DENSE 3

typedef enum cblktype_e {
    CBLK_1D,
    CBLK_SPLIT,
    CBLK_H,
    CBLK_DENSE,
    CBLK_SHUR
} cblktype_t;

#define CLUSTER   1
#define NOCLUSTER 0

/*+ Processor candidate group to own a column blok      +*/
typedef struct Cand_{
    double       costlevel;    /*+ Cost from root to node +*/
    pastix_int_t treelevel;    /*+ Level of the cblk in the elimination tree (deepness from the root) +*/
    pastix_int_t fcandnum;     /*+ first processor number of this candidate group  +*/
    pastix_int_t lcandnum;     /*+ last processor number of this candidate group   +*/
    pastix_int_t fccandnum;    /*+ first cluster number of the cluster candidate group +*/
    pastix_int_t lccandnum;    /*+ last cluster number of the cluster candidate group +*/
    pastix_int_t cluster;      /*+ Cluster id on which the task will be executed +*/
    cblktype_t   cblktype;     /*+ type of the distribution +*/
} Cand;

void candInit           ( Cand *candtab,
                          pastix_int_t cblknbr );
void candSetSubCandidate( Cand *candtab,
                          const EliminTree *etree,
                          pastix_int_t rootnum,
                          pastix_int_t procnum );
int  candCheck          ( Cand *candtab,
                          SymbolMatrix *symbmtx );
void candSetClusterCand ( Cand *candtab,
                          pastix_int_t  cblknbr,
                          pastix_int_t *core2clust,
                          pastix_int_t  coresnbr );
void candSave( const Cand *candtab,
               pastix_int_t cblknbr );

void candBuild( pastix_int_t autolevel, pastix_int_t level2D, double ratiolimit,
                Cand               *candtab,
                EliminTree         *etree,
                const SymbolMatrix *symbmtx,
                const CostMatrix   *costmtx );

#endif
