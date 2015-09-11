/**
 *
 * @file simu.h
 *
 *  PaStiX simulation routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#ifndef _SIMU_H_
#define _SIMU_H_

#include "simu_timer.h"

/*
 **  The type and structure definitions.
 */
typedef struct SimuCluster_ {
    pastix_int_t     fprocnum;   /*> Global index of the first processor belonging to the cluster (Check is it is not redundant) */
    pastix_int_t     lprocnum;   /*> Global index of the last processor belonging to the cluster (Check is it is not redundant)  */
    ExtendVectorINT *ftgtsend;   /*> Arrays of ftgt sent by this proc (one vector per processor)                                 */
    pastix_int_t     prionum;    /*> Counter to order tasks on one cluster                                                       */
} SimuCluster;

typedef struct SimuProc_ {
    SimuTimer        timer;      /*> Simulated clock of the processor                                  */
    pastix_queue_t  *readytask;  /*> Heap of tasks ready to be executed                                */
    pastix_queue_t  *futuretask; /*> Heap of tasks ready to be executed in a near future (after timer) */
    ExtendVectorINT *tasktab;    /*> Vector to store tasks affected to the candidate                   */
} SimuProc;

typedef struct SimuFtgt_ {
    FanInTarget  ftgt;       /*> Fan-in informations                            */
    pastix_int_t clustnum;   /*> Cluster sending the contribution               */
    SimuTimer    timerecv;   /*> Simulated clock of the reception time          */
    double       costsend;   /*> Cost to send the contribution                  */
    double       costadd;    /*> Cost to add the contribution to its final cblk */
} SimuFtgt;

typedef struct SimuCblk_ {
    pastix_int_t ctrbcnt;    /*> Counter of remaining contributions for the cblk */
} SimuCblk;

typedef struct SimuBlok_ {
    pastix_int_t tasknum;    /*> Task index opeating on this block (stored per block for 2D computations)   */
    pastix_int_t ftgtnum;    /*> Index of the first fanin destinated to this
                              *  block in the ftgttab. This index is also used to find the first cblk timer
                              *  (one per cand proc) in the timetab array                                   */
    pastix_int_t ctrbcnt;    /*> Counter of remaining contributions                                         */
    int          fccandnum;  /*> First candidate that is attributed to the cblk of the block                */
    int          ownerclust; /*> Processor on which the block is distributed                                */
} SimuBlok;

typedef struct SimuTask_ {
    pastix_int_t taskid;      /*> identification of the task type                        +*/
    pastix_int_t prionum;     /*> priority of the task                                   +*/
    pastix_int_t cblknum;     /*> Number of the cblknum the task deals with              +*/
    pastix_int_t bloknum;     /*> number of the block that the task deals with           +*/
    pastix_int_t bloknum2;
    pastix_int_t facebloknum; /*> Number of the facing block for E2                      +*/
    SimuTimer    time;        /*> Time the task is ready if it doesn't need message      +*/
    pastix_int_t mesglen;     /*> Time to send the block target                          +*/
    double       cost;        /*> Cost of the task                                       +*/
    pastix_int_t ctrbcnt;     /*> nbr ftgt + le btgt (+ E1 pret si E2) */
    pastix_int_t ftgtcnt;     /*> nbr of contrib from fan-in target                       +*/
    pastix_int_t tasknext;    /*> chainage des E1 ou E2, si fin = -1 => liberer les btagptr */
} SimuTask;

/* task type allowed*/
#ifndef COMP_1D
#define COMP_1D    0
#define DIAG       1
#define E1         2
#define E2         3
#endif

typedef struct SimuCtrl_ {
    pastix_int_t  cblknbr;     /*+ Number of cblk                                          +*/
    pastix_int_t  ftgtprio;    /*+ Priority to assign to current ftgts                     +*/
    pastix_int_t  tasknbr;     /*+ Number of tasks                                         +*/
    pastix_int_t  ftgtcnt;     /*+ Number of received communication                        +*/
    SimuTask     *tasktab;     /*+ SimuTask vector                                         +*/
    SimuProc     *proctab;     /*+ Virtual processor tab                                   +*/
    SimuCluster  *clustab;     /*+ Virtual cluster tab                                     +*/
    pastix_int_t *ownetab;     /*+ Vector containing the distribution of the diagonal blok +*/
    SimuCblk     *cblktab;     /*+ SimuCblk vector                                         +*/
    SimuBlok     *bloktab;     /*+ SimuBlok vector                                         +*/
    SimuFtgt     *ftgttab;     /*+ Vector containing the fan in target                     +*/
    pastix_int_t  ftgtnbr;
    SimuTimer    *ftgttimetab; /*+ Vector containing a timer for each cluster on each ftgt  +*/
} SimuCtrl;

/*
 **  The function prototypes.
 */
pastix_int_t simuInit   ( SimuCtrl *, const SymbolMatrix *, const Cand *, pastix_int_t, pastix_int_t );
pastix_int_t simuRealloc( SimuCtrl *, pastix_int_t, pastix_int_t );
void         simuExit   ( SimuCtrl *, pastix_int_t, pastix_int_t, pastix_int_t );
void         simuRun    ( SimuCtrl *, const BlendCtrl *, const SymbolMatrix * );

#define CLUST2INDEX(n,c) ((c) + simuctrl->bloktab[n].ftgtnum - simuctrl->bloktab[n].fccandnum)
#define INDEX2CLUST(r,s) ((r) - simuctrl->bloktab[s].ftgtnum + simuctrl->bloktab[s].fccandnum)
#define TIMER(pr)        (&(simuctrl->proctab[pr].timer))

#endif /* _SIMU_H_ */
