/************************************************************/
/**                                                        **/
/**   NAME       : param_blend.h                           **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the blend control structure.        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     08 sep 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ The parameters structure definition +*/

typedef struct netperf_ {
  double startup;
  double bandwidth;
} netperf;


/*+ Structure containing the structure passed through the blend primitives +*/
typedef struct BlendCtrl_ {
  BlendParam        *option;
  netperf           *perfptr;

  PASTIX_INT                clustnum;      /*+ Local cluster ID                       +*/
  PASTIX_INT                clustnbr;      /*+ Number of MPI process                  +*/
  PASTIX_INT                procnbr;       /*+ Number total of processors             +*/
  PASTIX_INT                proclocnbr;    /*+ Number of processors for one clustnum  +*/
  PASTIX_INT                thrdnbr;       /*+ Number total of threads                +*/
  PASTIX_INT                thrdlocnbr;    /*+ Number of threads for one clustnum     +*/
  PASTIX_INT                cudanbr;       /*+ Number of cuda device for one clustnum +*/
  PASTIX_INT                bublnbr;       /*+ Number of threads for one clustnum     +*/
  PASTIX_INT               *proc2clust;    /*+ proc2clust[i] = cluster of proc i      +*/
  BubbleTree        *btree;         /*+ arbre de bulles +*/
  EliminGraph       *egraph;        /*+ the elimination graph (only in vertex) +*/
  EliminTree        *etree;         /*+ the elimination tree                   +*/
  CostMatrix        *costmtx;       /*+ the cost bounded to each cblk and blok +*/
  Cand              *candtab;       /*+ processor candidate tab                +*/
  Queue             *lheap;         /*+ Use to order leaves                    +*/
  ExtendVectorINT   *intvec;        /*+ vector of PASTIX_INT used by several routines.
                                     The aim of this variable is to avoid
                                     repetedly memAlloc and memFree call      +*/
  ExtendVectorINT   *intvec2;       /*+ Another one                            +*/
  FILE              *tracefile;
} BlendCtrl;



/* Calcul le numero du process MPI */
#define CLUSTNUM(proc) (proc/(ctrl->procnbr/ctrl->clustnbr))
/* Calcul le noeud SMP sur lequel se trouve le process MPI si plusieurs par noeud SMP */
#define SMPNUM(clust) ( (clust*(ctrl->procnbr/ctrl->clustnbr)) / ctrl->option->procnbr)


PASTIX_INT      blendCtrlInit (BlendCtrl *, PASTIX_INT, PASTIX_INT, PASTIX_INT, PASTIX_INT, BlendParam *);
void     blendCtrlExit (BlendCtrl *);
void     perfcluster2  (PASTIX_INT procsrc, PASTIX_INT procdst, PASTIX_INT sync_comm_nbr, 
			netperf *np, BlendCtrl *ctrl);

