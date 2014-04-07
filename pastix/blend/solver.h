/************************************************************/
/**                                                        **/
/**   NAME       : solver.h                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the solver matrix.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     28 oct 1998     **/
/**                # Version 1.0  : from : 06 jun 2002     **/
/**                                 to     06 jun 2002     **/
/**                                                        **/
/************************************************************/

#ifndef SOLVER_H
#define SOLVER_H

/*
**  The type and structure definitions.
*/

#define COMP_1D                     0
#define DIAG                        1
#define E1                          2
#define E2                          3
#define DRUNK                       4

typedef struct Task_ {
  pastix_int_t                       taskid;               /*+ COMP_1D DIAG E1 E2                                     +*/
  pastix_int_t                       prionum;              /*+ Priority value for the factorization                   +*/
  pastix_int_t                       cblknum;              /*+ Attached column block                                  +*/
  pastix_int_t                       bloknum;              /*+ Attached block                                         +*/
  pastix_int_t volatile              ftgtcnt;              /*+ Number of fan-in targets                               +*/
  pastix_int_t volatile              ctrbcnt;              /*+ Total number of contributions                          +*/
  pastix_int_t                       indnum;               /*+ For E2 (COMP_1D), index of ftgt (>0) else if local = -taskdest
                                                            For DIAG and E1 , index of btag (>0) if there is a
                                                            local one it must be the first of the chain of local E1   +*/
#if (defined PASTIX_DYNSCHED) || (defined TRACE_SOPALIN)
  pastix_int_t                       threadid;             /*+ Index of the bubble which contains the task +*/
#endif
#ifdef TRACE_SOPALIN
  pastix_int_t                       fcandnum;             /*+ First thread candidate                      +*/
  pastix_int_t                       lcandnum;		  /*+ Last thread candidate                       +*/
  pastix_int_t                       id;                   /*+ Global cblknum of the attached column block +*/
#endif
} Task;

/*+ Solver column block structure. +*/

typedef struct SolverCblk_  {
  pastix_int_t                       fcolnum;              /*+ First column index                     +*/
  pastix_int_t                       lcolnum;              /*+ Last column index (inclusive)          +*/
  pastix_int_t                       bloknum;              /*+ First block in column (diagonal)       +*/
  pastix_int_t                       stride;               /*+ Column block stride                    +*/
  pastix_int_t                       procdiag;             /*+ Processor owner of diagonal block      +*/
  pastix_float_t * restrict          coeftab;              /*+ Coefficients access vector             +*/
  pastix_float_t * restrict          ucoeftab;             /*+ Coefficients access vector             +*/
} SolverCblk;

/*+ Solver block structure. +*/

typedef struct SolverBlok_ {
  pastix_int_t                       frownum;              /*+ First row index            +*/
  pastix_int_t                       lrownum;              /*+ Last row index (inclusive) +*/
  pastix_int_t                       cblknum;              /*+ Facing column block        +*/
  pastix_int_t                       levfval;              /*+ Level-of-fill value        +*/
  pastix_int_t                       coefind;              /*+ Index in coeftab           +*/
} SolverBlok;

/*+ Solver matrix structure. +*/

/* All data are local to one cluster */
typedef struct SolverMatrix_ {
  pastix_int_t              baseval;              /*+ Base value for numberings                         +*/

  pastix_int_t              nodenbr;              /*+ Number of nodes before dof extension              +*/
  pastix_int_t              coefnbr;              /*+ Number of coefficients (node after dof extension) +*/
  pastix_int_t              gcblknbr;             /*+ Global number of column blocks                    +*/
  pastix_int_t              cblknbr;              /*+ Number of column blocks                   +*/
  pastix_int_t              bloknbr;              /*+ Number of blocks                          +*/
  SolverCblk * restrict     cblktab;              /*+ Array of solver column blocks             +*/
  SolverBlok * restrict     bloktab;              /*+ Array of solver blocks                    +*/

#ifdef PASTIX_WITH_STARPU
  pastix_int_t              hcblknbr;
  pastix_int_t *            gcblk2halo;
  SolverHaloCblk * restrict hcblktab;
  SolverHaloBlok * restrict hbloktab;
#endif

  pastix_int_t              ftgtnbr;              /*+ Number of fanintargets                    +*/
  pastix_int_t              ftgtcnt;              /*+ Number of fanintargets to receive         +*/
  FanInTarget * restrict    ftgttab;              /*+ Fanintarget access vector                 +*/

  pastix_int_t              coefmax;              /*+ Working block max size (cblk coeff 1D)    +*/
  pastix_int_t              nbftmax;              /*+ Maximum block number in ftgt              +*/
  pastix_int_t              arftmax;              /*+ Maximum block area in ftgt                +*/

  pastix_int_t              clustnum;             /*+ current processor number                  +*/
  pastix_int_t              clustnbr;             /*+ number of processors                      +*/
  pastix_int_t              procnbr;              /*+ Number of physical processor used         +*/
  pastix_int_t              thrdnbr;              /*+ Number of local computation threads       +*/
  pastix_int_t              bublnbr;              /*+ Number of local computation threads       +*/
  BubbleTree  * restrict    btree;                /*+ Bubbles tree                              +*/

  pastix_int_t              indnbr;
  pastix_int_t * restrict   indtab;
  Task * restrict           tasktab;              /*+ Task access vector                        +*/
  pastix_int_t              tasknbr;              /*+ Number of Tasks                           +*/
  pastix_int_t **           ttsktab;              /*+ Task access vector by thread              +*/
  pastix_int_t *            ttsknbr;              /*+ Number of tasks by thread                 +*/

  pastix_int_t *            proc2clust;           /*+ proc -> cluster                           +*/
  pastix_int_t              gridldim;             /*+ Dimensions of the virtual processors      +*/
  pastix_int_t              gridcdim;             /*+ grid if dense end block                   +*/
  UpDownVector              updovct;              /*+ UpDown vector                             +*/
} SolverMatrix;

pastix_int_t
sizeofsolver(const SolverMatrix *solvptr,
             pastix_int_t *iparm );

#endif /* SOLVER_H */
