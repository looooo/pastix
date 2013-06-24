#ifndef _PASTIX_DATA_H_
#define _PASTIX_DATA_H_

#if defined(HAVE_SCOTCH)
#include <scotch.h>
#endif

#include "sopalin_define.h"
#include "dof.h"
#include "ftgt.h"
#include "symbol.h"
#include "csc.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"
/* #include "assembly.h" */
/* #include "param_blend.h" */
/* #include "order.h" */
/* #include "fax.h" */
/* #include "kass.h" */
/* #include "blend.h" */
/* #include "solverRealloc.h" */
/* #include "sopalin_thread.h" */
/* #include "stack.h" */
/* #include "sopalin3d.h" */
/* #include "sopalin_init.h" */
/* #include "sopalin_option.h" */
/* #include "csc_intern_updown.h" */
/* #include "csc_intern_build.h" */
/* #include "coefinit.h" */

/*
  struct: SopalinParam_

  Parameters for factorisation, updown and reffinement.
 */
typedef struct SopalinParam_ {
  CscMatrix      *cscmtx;          /*+ Compress Sparse Column matrix                    +*/
  double          epsilonraff;     /*+ epsilon to stop reffinement                      +*/
  double          rberror;         /*+ ||r||/||b||                                      +*/
  double          espilondiag;     /*+ epsilon critere for diag control                 +*/
  pastix_float_t *b;               /*+ b vector (RHS and solution)                      +*/
  pastix_float_t *transcsc;        /*+ transpose csc                                    +*/
  pastix_int_t    itermax;         /*+ max number of iteration                          +*/
  pastix_int_t    diagchange;      /*+ number of change of diag                         +*/
  pastix_int_t    gmresim;         /*+ Krylov subspace size for GMRES                   +*/
  pastix_int_t    fakefact;        /*+ Flag indicating if we want fake factorisation    +*/
  pastix_int_t    usenocsc;        /*+ Flag indicating if we want to use the intern CSC +*/
  int             factotype;       /*+ Type of factorization                            +*/
  int             symmetric;       /*+ Symmetric                                        +*/
  MPI_Comm        pastix_comm;     /*+ MPI communicator                                 +*/
  int             type_comm;       /*+ Communication mode                               +*/
  int             nbthrdcomm;      /*+ Communication's thread number                    +*/
  pastix_int_t   *iparm;           /*+ In/Out integer parameters                        +*/
  double         *dparm;           /*+ In/Out float parameters                          +*/
  int            *bindtab;         /*+ Define where to bin threads                      +*/
  int             stopthrd;        /*+ Boolean for communication thread controlling     +*/
  int             schur;           /*+ If API_YES won't compute last diag               +*/
  pastix_int_t    n;               /*+ size of the matrix                               +*/
  pastix_int_t    gN;
} SopalinParam;

/*
   struct: pastix_data_t

   Structure used to store datas for a step by step execution.
*/

struct pastix_data_s {
  SolverMatrix     solvmatr;           /*+ Matrix informations                                                 +*/
  CscMatrix        cscmtx;             /*+ Compress Sparse Column matrix                                       +*/
  SymbolMatrix    *symbmtx;            /*+ Symbol Matrix                                                       +*/
  SopalinParam     sopar;              /*+ Sopalin parameters                                                  +*/
  Order            ordemesh;           /*+ Order                                                               +*/
#ifdef WITH_SCOTCH
  SCOTCH_Graph     grafmesh;           /*+ Graph                                                               +*/
  int              malgrf;             /*+ boolean indicating if grafmesh has been allocated                   +*/
#endif /* WITH_SCOTCH */
#ifdef DISTRIBUTED
#ifdef WITH_SCOTCH
  SCOTCH_Dordering ordedat;            /*+ distributed scotch order                                            +*/
  SCOTCH_Dgraph    dgraph;
  pastix_int_t    *PTS_permtab;
  pastix_int_t    *PTS_peritab;
#endif /* WITH_SCOTCH */
  pastix_int_t    *glob2loc;           /*+ local column number of global column, or -(owner+1) is not local    +*/
  pastix_int_t     ncol_int;           /*+ Number of local columns in internal CSCD                            +*/
  pastix_int_t    *l2g_int;            /*+ Local to global column numbers in internal CSCD                     +*/
  int              malrhsd_int;        /*+ Indicates if internal distributed rhs has been allocated            +*/
  int              mal_l2g_int;
  pastix_float_t  *b_int;              /*+ Local part of the right-hand-side                                   +*/
#endif /* DISTRIBUTED */
  pastix_int_t    *loc2glob2;          /*+ local2global column number                                          +*/
  pastix_int_t     gN;                 /*+ global column number                                                +*/
  pastix_int_t     n;                  /*+ local column number                                                 +*/
  pastix_int_t    *iparm;              /*+ Vecteur de parametres entiers                                       +*/
  double          *dparm;              /*+ Vecteur de parametres floattant                                     +*/
  pastix_int_t     n2;                 /*+ Number of local columns                                             +*/
  pastix_int_t    *col2;               /*+ column tabular for the CSC matrix                                   +*/
                                       /*+ (index of first element of each col in row and values tabulars)     +*/
  pastix_int_t    *row2;               /*+ tabular containing row number of each element of                    +*/
                                       /*+  the CSC matrix, ordered by column.                                 +*/
  int              bmalcolrow;         /*+ boolean indicating if col2 ans row2 have been allocated             +*/
  int              malord;             /*+ boolean indicating if ordemesh has been allocated                   +*/
  int              malcsc;             /*+ boolean indicating if solvmatr->cscmtx has beek allocated           +*/
  int              malsmx;             /*+ boolean indicating if solvmatr->updovct.sm2xtab has been allocated  +*/
  int              malslv;             /*+ boolean indicating if solvmatr has been allocated                   +*/
  int              malcof;             /*+ boolean indicating if coeficients tabular(s) has(ve) been allocated +*/
  MPI_Comm         pastix_comm;        /*+ PaStiX MPI communicator                                             +*/
  MPI_Comm         intra_node_comm;    /*+ PaStiX intra node MPI communicator                                  +*/
  MPI_Comm         inter_node_comm;    /*+ PaStiX inter node MPI communicator                                  +*/
  int              procnbr;            /*+ Number of MPI tasks                                                 +*/
  int              procnum;            /*+ Local MPI rank                                                      +*/
  int              intra_node_procnbr; /*+ Number of MPI tasks in node_comm                                    +*/
  int              intra_node_procnum; /*+ Local MPI rank in node_comm                                         +*/
  int              inter_node_procnbr; /*+ Number of MPI tasks in node_comm                                    +*/
  int              inter_node_procnum; /*+ Local MPI rank in node_comm                                         +*/
  int             *bindtab;            /*+ Tabular giving for each thread a CPU to bind it too                 +*/
  pastix_int_t     nschur;             /*+ Number of entries for the Schur complement.                         +*/
  pastix_int_t    *listschur;          /*+ List of entries for the schur complement.                           +*/
  pastix_float_t  *schur_tab;
  pastix_int_t     schur_tab_set;
  int              cscInternFilled;
  int              scaling;            /*+ Indicates if the matrix has been scaled                             +*/
  pastix_float_t  *scalerowtab;        /*+ Describes how the matrix has been scaled                            +*/
  pastix_float_t  *iscalerowtab;
  pastix_float_t  *scalecoltab;
  pastix_float_t  *iscalecoltab;
#ifdef WITH_SEM_BARRIER
  sem_t           *sem_barrier;        /*+ Semaphore used for AUTOSPLIT_COMM barrier                           +*/
#endif
  pastix_int_t     pastix_id;          /*+ Id of the pastix instance (PID of first MPI task)                   +*/
};

#endif /* _PASTIX_DATA_H_ */
