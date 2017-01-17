/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.2.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#ifndef PASTIX_DATA_H_
#define PASTIX_DATA_H_

#include "isched.h"
#if defined(PASTIX_WITH_PARSEC)
#include <parsec.h>
#endif

#include "ftgt.h"
#include "symbol.h"
#include "queue.h"
#include "bulles.h"

/*
 * Steps of the pastix solver
 */
#define STEP_INIT      (1 << 0)
#define STEP_ORDERING  (1 << 1)
#define STEP_SYMBFACT  (1 << 2)
#define STEP_ANALYSE   (1 << 3)
#define STEP_CSC2BCSC  (1 << 4)
#define STEP_BCSC2CTAB (1 << 5)
#define STEP_NUMFACT   (1 << 6)
#define STEP_SOLVE     (1 << 7)
#define STEP_REFINE    (1 << 8)

struct pastix_bcsc_s;
typedef struct pastix_bcsc_s pastix_bcsc_t;

/* /\* */
/*   struct: SopalinParam_ */

/*   Parameters for factorisation, updown and refinement. */
/*  *\/ */
/* typedef struct SopalinParam_ { */
/*     pastix_bcsc_t  *bcsc;          /\*+ Compress Sparse Column matrix                    *\/ */
/*     double          epsilonraff;     /\*+ epsilon to stop refinement                      *\/ */
/*     double          rberror;         /\*+ ||r||/||b||                                      *\/ */
/*     double          espilondiag;     /\*+ epsilon critere for diag control                 *\/ */
/*     void *b;               /\*+ b vector (RHS and solution)                      *\/ */
/*     void *transcsc;        /\*+ transpose csc                                    *\/ */
/*     pastix_int_t    itermax;         /\*+ max number of iteration                          *\/ */
/*     pastix_int_t    diagchange;      /\*+ number of change of diag                         *\/ */
/*     pastix_int_t    gmresim;         /\*+ Krylov subspace size for GMRES                   *\/ */
/*     pastix_int_t    fakefact;        /\*+ Flag indicating if we want fake factorisation    *\/ */
/*     pastix_int_t    usenocsc;        /\*+ Flag indicating if we want to use the intern CSC *\/ */
/*     int             factotype;       /\*+ Type of factorization                            *\/ */
/*     int             symmetric;       /\*+ Symmetric                                        *\/ */
/*     MPI_Comm        pastix_comm;     /\*+ MPI communicator                                 *\/ */
/*     int             type_comm;       /\*+ Communication mode                               *\/ */
/*     int             nbthrdcomm;      /\*+ Communication's thread number                    *\/ */
/*     pastix_int_t   *iparm;           /\*+ In/Out integer parameters                        *\/ */
/*     double         *dparm;           /\*+ In/Out float parameters                          *\/ */
/*     int            *bindtab;         /\*+ Define where to bin threads                      *\/ */
/*     int             stopthrd;        /\*+ Boolean for communication thread controlling     *\/ */
/*     int             schur;           /\*+ If API_YES won't compute last diag               *\/ */
/*     pastix_int_t    n;               /\*+ size of the matrix                               *\/ */
/*     pastix_int_t    gN; */
/* } SopalinParam; */

/*
 struct: pastix_data_t

 Structure used to store datas for a step by step execution.
 */

struct pastix_data_s {
    pastix_int_t    *iparm;              /*< Store integer parameters (input/output)                             */
    double          *dparm;              /*< Store floating parameters (input/output)                            */

    pastix_int_t     steps;              /*< Bitmask of the steps performed or not                               */

    MPI_Comm         pastix_comm;        /*< PaStiX MPI communicator used for the ordering step                  */
    MPI_Comm         intra_node_comm;    /*< PaStiX intra node MPI communicator used for synchronizations        */
    MPI_Comm         inter_node_comm;    /*< PaStiX inter node MPI communicator used for the factorization       */
    int              initmpi;            /*< MPI Initialized by PaStiX                                           */
    int              procnbr;            /*< Total number of MPI processes                                       */
    int              procnum;            /*< Local MPI rank                                                      */
    int              intra_node_procnbr; /*< Number of MPI tasks in intra node communicator                      */
    int              intra_node_procnum; /*< Local MPI rank in intra node communicator                           */
    int              inter_node_procnbr; /*< Number of MPI tasks in inter node communicator                      */
    int              inter_node_procnum; /*< Local MPI rank in inter node communicator                           */

    isched_t        *isched;             /*< Internal scheduler structure that is always available               */
#if defined(PASTIX_WITH_PARSEC)
    parsec_context_t *parsec;             /*< PaRSEC Context if available                                         */
#endif
#if defined(PASTIX_WITH_STARPU)
    starpu_ctxt_t   *starpu;             /*< StarPU Context if available                                         */
#endif

    const pastix_csc_t *csc;             /*< Pointer to the user csc structure used as input                     */

    pastix_graph_t  *graph;              /*< Symmetrized graph of the problem used within ordering
                                             and symbolic factorization steps.                                   */
    pastix_int_t     schur_n;            /*< Number of entries for the Schur complement                          */
    pastix_int_t    *schur_list;         /*< List of entries for the schur complement                            */
    pastix_int_t     zeros_n;            /*< Number of diagonal entries considered as zeros                      */
    pastix_int_t    *zeros_list;         /*< List of diagonal entries considered as zeros                        */
    Order           *ordemesh;           /*< Ordering structure                                                  */

    SymbolMatrix    *symbmtx;            /*< Symbol Matrix                                                       */

    pastix_bcsc_t   *bcsc;               /*< Csc after reordering grouped by cblk                                */
    SolverMatrix    *solvmatr;           /*< Solver informations associted to the matrix problem                 */

    /* Backup for old pqstix interface */
    void            *b;
    void            *x0;

    /**
     * Former fields that are no longer used for now
     */
    /* SopalinParam     sopar;              /\* Sopalin parameters                                                  *\/ */
#ifdef PASTIX_DISTRIBUTED
#if defined(PASTIX_ORDERING_SCOTCH)
    pastix_int_t    *PTS_permtab;
    pastix_int_t    *PTS_peritab;
#endif /* PASTIX_ORDERING_SCOTCH */
    pastix_int_t    *glob2loc;           /*+ local column number of global column, or -(owner+1) is not local    */
    pastix_int_t     ncol_int;           /*+ Number of local columns in internal CSCD                            */
    pastix_int_t    *l2g_int;            /*+ Local to global column numbers in internal CSCD                     */
    void  *b_int;              /*+ Local part of the right-hand-side                                   */
#endif /* PASTIX_DISTRIBUTED */

    int             *bindtab;            /*+ Tabular giving for each thread a CPU to bind it too                 */
    void            *schur_tab;
    pastix_int_t     schur_tab_set;
    int              cscInternFilled;
    int              scaling;            /*+ Indicates if the matrix has been scaled                             */
    void  *scalerowtab;        /*+ Describes how the matrix has been scaled                            */
    void  *iscalerowtab;
    void  *scalecoltab;
    void  *iscalecoltab;
#ifdef WITH_SEM_BARRIER
    sem_t           *sem_barrier;        /*+ Semaphore used for AUTOSPLIT_COMM barrier                           */
#endif
    pastix_int_t     pastix_id;          /*+ Id of the pastix instance (PID of first MPI task)                   */
};

#endif /* PASTIX_DATA_H_ */
