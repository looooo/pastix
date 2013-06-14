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

  pastix_int_t                clustnum;      /*+ Local cluster ID                       +*/
  pastix_int_t                clustnbr;      /*+ Number of MPI process                  +*/
  pastix_int_t                procnbr;       /*+ Number total of processors             +*/
  pastix_int_t                proclocnbr;    /*+ Number of processors for one clustnum  +*/
  pastix_int_t                thrdnbr;       /*+ Number total of threads                +*/
  pastix_int_t                thrdlocnbr;    /*+ Number of threads for one clustnum     +*/
  pastix_int_t                cudanbr;       /*+ Number of cuda device for one clustnum +*/
  pastix_int_t                bublnbr;       /*+ Number of threads for one clustnum     +*/
  pastix_int_t               *proc2clust;    /*+ proc2clust[i] = cluster of proc i      +*/
  BubbleTree        *btree;         /*+ arbre de bulles +*/
  EliminGraph       *egraph;        /*+ the elimination graph (only in vertex) +*/
  EliminTree        *etree;         /*+ the elimination tree                   +*/
  CostMatrix        *costmtx;       /*+ the cost bounded to each cblk and blok +*/
  Cand              *candtab;       /*+ processor candidate tab                +*/
  Queue             *lheap;         /*+ Use to order leaves                    +*/
  ExtendVectorINT   *intvec;        /*+ vector of pastix_int_t used by several routines.
                                     The aim of this variable is to avoid
                                     repetedly memAlloc and memFree call      +*/
  ExtendVectorINT   *intvec2;       /*+ Another one                            +*/
  FILE              *tracefile;
} BlendCtrl;



/* Calcul le numero du process MPI */
#define CLUSTNUM(proc) (proc/(ctrl->procnbr/ctrl->clustnbr))
/* Calcul le noeud SMP sur lequel se trouve le process MPI si plusieurs par noeud SMP */
#define SMPNUM(clust) ( (clust*(ctrl->procnbr/ctrl->clustnbr)) / ctrl->option->procnbr)


pastix_int_t      blendCtrlInit (BlendCtrl *, pastix_int_t, pastix_int_t, pastix_int_t, pastix_int_t, BlendParam *);
void     blendCtrlExit (BlendCtrl *);
void     perfcluster2  (pastix_int_t procsrc, pastix_int_t procdst, pastix_int_t sync_comm_nbr, 
			netperf *np, BlendCtrl *ctrl);

