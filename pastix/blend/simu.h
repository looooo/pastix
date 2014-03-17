/*+********************************************************+*/
/*+                                                        +*/
/*+   NAME       : simu.h                                  +*/
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


/*
**  The type and structure definitions.
*/
typedef struct SimuTimer_ {
  double s;
  double ms;
} SimuTimer;

typedef struct SimuCluster_ {
  pastix_int_t       fprocnum;
  pastix_int_t       lprocnum;
  ExtendVectorINT  * ftgtsend;         /*+ ftgt send by this proc (one vector per processor       +*/
  pastix_int_t       prionum;
} SimuCluster;

typedef struct SimuProc_ {
    SimuTimer          timer;            /*+ Simulated clock of the processor                       +*/

    pastix_queue_t   *readytask;   /* Heap of tasks ready to be executed */
    pastix_queue_t   *futuretask;  /* Heap of tasks ready to be executed in a near future */

    pastix_int_t       prionum;    /*+ Current priority to assign to a cblk mapp on this proc +*/
    ExtendVectorINT  * tasktab;    /* Vector to store tasks affected to the candidate */
} SimuProc;

typedef struct SimuFtgt_ {
  FanInTarget           ftgt;             /*+ the ftgt structures info                              +*/
  pastix_int_t          clustnum;         /*+ the cluster that send the ftgt                        +*/
  SimuTimer             timerecv;         /*+ time the ftgt will be receive                         +*/
  double                costsend;         /*+ time for the ftgt go from procsrc to procdest         +*/
  double                costadd;          /*+ cost to add the ftgt to its cblk destination          +*/
} SimuFtgt;

typedef struct SimuBlockTarget_ {
  pastix_int_t           prionum;
  pastix_int_t           taskid;
  pastix_int_t           bloksrc;
  pastix_int_t           cblkdest;
  pastix_int_t           frownum;
  pastix_int_t           lrownum;
} SimuBlockTarget;

typedef struct SimuBlok_ {
  pastix_int_t           tasknum;          /*+ Number of the task                                     +*/
  pastix_int_t           ftgtnum;          /*+ index of the first ftgt destinated to this cblk
                                                    in the ftgttab. This index is also used to find
                                                    the first cblk timer (one per cand proc) in the
                                                    timetab                                                +*/
  pastix_int_t           ctrbcnt;          /*+ counter for contributions remaining                    +*/
  pastix_int_t           ctrbnbr;          /*+ OIMBE temporaire sert juste pour DRUNK                 +*/

} SimuBlok;


typedef struct SimuTask_ {
  pastix_int_t  taskid;           /*+ identification of the task type                        +*/
  pastix_int_t  prionum;          /*+ priority of the task                                   +*/
  pastix_int_t  cblknum;          /*+ Number of the cblknum the task deals with              +*/
  pastix_int_t  bloknum;          /*+ number of the block that the task deals with           +*/
  pastix_int_t  bloknum2;
  pastix_int_t  facebloknum;      /*+ Number of the facing block for E2                      +*/
  SimuTimer     time;             /*+ Time the task is ready if it doesn't need message      +*/
  pastix_int_t  mesglen;          /*+ Time to send the block target                          +*/
  double        cost;             /*+ Cost of the task                                       +*/
  pastix_int_t  ctrbcnt;          /* nbr ftgt + le btgt (+ E1 pret si E2) */
  pastix_int_t  ftgtcnt;          /* nbr of contrib from fan-in target                       +*/
  pastix_int_t  tasknext;         /* chainage des E1 ou E2, si fin = -1 => liberer les btagptr */
} SimuTask;

/* task type allowed*/
#ifndef COMP_1D
#define COMP_1D    0
#define DIAG       1
#define E1         2
#define E2         3
#endif

typedef struct SimuCblk_ {
  pastix_int_t                   ctrbcnt;          /*+ counter for contributions remaining                    +*/
} SimuCblk;


typedef struct SimuCtrl_ {
  pastix_int_t          cblknbr;          /*+ Number of cblk                                          +*/
  pastix_int_t          ftgtprio;         /*+ Priority to assign to current ftgts                     +*/
  pastix_int_t          tasknbr;          /*+ Number of tasks                                         +*/
  pastix_int_t          ftgtcnt;          /*+ Number of received communication                        +*/
  SimuTask         *    tasktab;          /*+ SimuTask vector                                         +*/
  SimuProc         *    proctab;          /*+ Virtual processor tab                                   +*/
  SimuCluster      *    clustab;          /*+ Virtual cluster tab                                     +*/
  pastix_int_t     *    ownetab;          /*+ Vector containing the distribution of the diagonal blok +*/
  pastix_int_t     *    blprtab;          /*+ Vector containing the distribution of the blok          +*/
  SimuCblk         *    cblktab;          /*+ SimuCblk vector                                         +*/
  SimuBlok         *    bloktab;          /*+ SimuBlok vector                                         +*/
  SimuFtgt         *    ftgttab;          /*+ Vector containing the fan in target                     +*/
  pastix_int_t          ftgtnbr;
  SimuTimer        *    tasktimetab;      /*+ Vector containing a timer for each cand proc on each task +*/
  SimuTimer        *    ftgttimetab;      /*+ Vector containing a timer for each cluster on each ftgt  +*/
} SimuCtrl;

/*
**  The function prototypes.
*/

#ifndef SIMU
#define static
#endif

pastix_int_t        simuInit        (SimuCtrl *, SymbolMatrix *, pastix_int_t, pastix_int_t, pastix_int_t, pastix_int_t, Cand *);
pastix_int_t        simuRealloc     (SimuCtrl *, pastix_int_t, pastix_int_t);
void                simuExit        (SimuCtrl *, pastix_int_t, pastix_int_t, pastix_int_t);
pastix_int_t        compTimer       (SimuTimer *, SimuTimer *);
void                timerAdd        (SimuTimer *, double);
double              timerVal        (const SimuTimer *);
void                timerSet        (SimuTimer *, double);
void                timerSetMax(SimuTimer *timer, double t);
#undef static

void              simuRun              (SymbolMatrix *, SimuCtrl *, BlendCtrl *, const Dof *);


#define CLUST2INDEX(n,c) ((c) + simuctrl->bloktab[n].ftgtnum - ctrl->candtab[ctrl->egraph->ownetab[n]].fccandnum)
#define INDEX2CLUST(r,s) ((r) - simuctrl->bloktab[s].ftgtnum + ctrl->candtab[ctrl->egraph->ownetab[s]].fccandnum)
#define TIMER(pr)        (&(simuctrl->proctab[pr].timer))
