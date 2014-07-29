/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
/*
 * File: z_pastix.c
 *
 * PaStiX external functions implementations.
 *
 * Authors:
 *   Mathieu FAVERGE  - faverge@labri.fr
 *   Xavier  LACOSTE  - lacoste@labri.fr
 *   Pierre  RAMET    - ramet@labri.fr
 */

#include "common.h"
#include <string.h>
#include <pthread.h>
#ifdef WITH_SEM_BARRIER
#  include <semaphore.h>
#  include <fcntl.h>
#  include <sys/ipc.h>
#  include <sys/shm.h>
#endif
#include <sys/stat.h>

#include "z_tools.h"
#include "sopalin_define.h"

#ifdef HAVE_SCOTCH
#  ifdef    PASTIX_DISTRIBUTED
#    include <ptscotch.h>
#  else
#    include <scotch.h>
#  endif /* PASTIX_DISTRIBUTED */
#endif /* WITH_SCOTCH */

#if defined(HAVE_METIS)
#include <metis.h>
#endif

#include "dof.h"
#include "z_ftgt.h"
#include "symbol.h"
#include "z_csc.h"
#include "z_updown.h"
#include "queue.h"
#include "bulles.h"
#include "z_solver.h"
#include "order.h"
#include "fax.h"
#include "kass.h"
#include "solverRealloc.h"
#include "sopalin_thread.h"
#include "stack.h"
#include "z_sopalin3d.h"
#include "z_sopalin_init.h"
#include "z_sopalin_option.h"
#include "z_csc_intern_updown.h"
#include "z_csc_intern_build.h"
#include "z_coefinit.h"
#include "out.h"
#include "z_pastix_internal.h"

#include "z_csc_utils.h"
#include "z_cscd_utils.h"
#include "z_cscd_utils_intern.h"
#include "sopalin_acces.h"
#include "z_scaling.c"
#include "perf.h"
#include "pastix_task_order.h"
#include "pastix_task_symbfact.h"

/*******************************************************************************
 * Section: Defines and Macros
 */

/*
 * Macro: FORTRAN_CALL
 *
 * Call a fortran function.
 *
 * Parameters:
 * name - Fortran function name.
 *
 */
#ifndef USE_NOFORTRAN
#  if (defined X_ARCHpower_ibm_aix)
#    define FORTRAN_CALL(name) name
#  else
#    define FORTRAN_CALL(name) name ## _
#  endif
#else
#  define FORTRAN_CALL(name)
#endif

#undef THREAD_FUNNELED_ON
#undef THREAD_FUNNELED_OFF
#define THREAD_FUNNELED_ON (                    \
        sopar->iparm[IPARM_THREAD_COMM_MODE] &  \
        API_THREAD_FUNNELED)
#define THREAD_FUNNELED_OFF (!THREAD_FUNNELED_ON)

#undef THREAD_COMM_ON
#undef THREAD_COMM_OFF
#define THREAD_COMM_ON  (                                       \
        sopar->iparm[IPARM_THREAD_COMM_MODE] &                  \
        ( API_THREAD_FUNNELED|API_THREAD_COMM_ONE|              \
          API_THREAD_COMM_DEFINED|API_THREAD_COMM_NBPROC ) )
#define THREAD_COMM_OFF (!THREAD_COMM_ON)

#define OPEN_SEM(sem, name, value) do {                                 \
        sem = sem_open(name, O_CREAT, S_IRUSR | S_IWUSR, value);        \
        if (sem == SEM_FAILED)                                          \
        {                                                               \
            errorPrint("sem_open failed\n");                            \
            perror("sem_open");                                         \
        }                                                               \
    } while(0)



/*
 * Defines: z_pastix.c defines
 *
 *   PASTIX_LOG            - If defined, start and end of this file
 *                           functions will be printed on stdout.
 *   COMPUTE               - If not defined, PaStiX will not use user's
 *                           coefficient.
 *   FORGET_PARTITION      - If defined, PaStiX won't use Scotch partition.
 *   DUMP_SYMBOLMATRIX     - Write the symbol matrix in a postscript format.
 *   STR_SIZE              - The default size of a string.
 *   TAG_RHS               - MPI tag used to communicate right-hand-side member.
 *   SCOTCH_STRAT_DIRECT   - Default Scotch strategy for the direct z_solver.
 *   SCOTCH_STRAT_INCOMP   - Default Scotch strategy for the incomplete z_solver.
 *   SCOTCH_STRAT_PERSO    - Parametrisable Scotch strategy for the direct
 *                           z_solver, can be set using several <IPARM_ACCESS>.
 *   PTSCOTCH_STRAT_DIRECT - Default PT-Scotch strategy for the direct z_solver.
 *   PTSCOTCH_STRAT_INCOMP - Default PT-Scotch strategy for the incomplete
 *                           z_solver.
 *   PTSCOTCH_STRAT_PERSO  - Parametrisable PT-Scotch strategy for the direct
 *                           z_solver, can be set using several <IPARM_ACCESS>.
 */
/*#define PASTIX_LOG*/
#define COMPUTE
/* #define FORGET_PARTITION  */
/* #define DUMP_SYMBOLMATRIX */

#define STR_SIZE               60
#define TAG_RHS                7
#define TAG_RHS2               8

/*******************************************************************************
 *  Section: Macros
 */

/*
 macro: print_onempi

 Print a string using processor 0.
 Uses printf syntax.

 Parameters:
 fmt - Format string (see printf manual).
 ... - Arguments depending on the format string.
 */
#define print_onempi(fmt, ...) if(procnum == 0) fprintf(stdout, fmt, __VA_ARGS__)

/*******************************************************************************
 * Section: Functions
 */

/*
 * Function: z_pastix_initParam
 *
 * sets default parameters for iparm and dparm
 *
 * Parameters:
 * iparm - tabular of IPARM_SIZE integer parameters.
 * dparm - tabular of DPARM_SIZE double parameters.
 */
void z_pastix_initParam(pastix_int_t    *iparm,
                      double *dparm)
{
#ifdef PASTIX_LOG
    fprintf(stdout, "-> api_init\n");
#endif

    pastix_int_t i;
    for (i=0; i<IPARM_SIZE; i++)
        iparm[i] = 0;

    for (i=0; i<DPARM_SIZE; i++)
        dparm[i] = 0;

    iparm[IPARM_MODIFY_PARAMETER]      = API_YES;             /* Indicate if parameters have been set by user         */
    iparm[IPARM_START_TASK]            = API_TASK_ORDERING;   /* Indicate the first step to execute (see PaStiX steps)*/
    iparm[IPARM_END_TASK]              = API_TASK_CLEAN;      /* Indicate the last step to execute (see PaStiX steps) */
    iparm[IPARM_VERBOSE]               = API_VERBOSE_NO;      /* Verbose mode (see Verbose modes)                     */
    iparm[IPARM_DOF_NBR]               = 1;                   /* Degree of freedom per node                           */
    iparm[IPARM_DOF_COST]              = 0;                   /* Degree of freedom for cost computation
                                                               (If different from IPARM_DOF_NBR) */
    iparm[IPARM_ITERMAX]               = 250;                 /* Maximum iteration number for refinement              */
    iparm[IPARM_MATRIX_VERIFICATION]   = API_YES;             /* Check the input matrix                               */
    iparm[IPARM_MC64]                  = 0;                   /* MC64 operation <z_pastix.h> IGNORE                     */
    iparm[IPARM_ONLY_RAFF]             = API_NO;              /* Refinement only                                      */
    iparm[IPARM_TRACEFMT]              = API_TRACE_PAJE;      /* Trace format (see Trace modes)                       */
    iparm[IPARM_GRAPHDIST]             = API_YES;             /* Specify if the given graph is distributed or not     */
    iparm[IPARM_AMALGAMATION_LEVEL]    = 5;                   /* Amalgamation level                                   */
    iparm[IPARM_ORDERING]              = API_ORDER_SCOTCH;    /* Choose ordering                                      */
    iparm[IPARM_ORDERING_DEFAULT]      = API_YES;             /* Use default ordering parameters with scotch or metis */

    /* Scotch default */
    iparm[IPARM_ORDERING_SWITCH_LEVEL] = 120;                 /* Ordering switch level    (see Scotch User's Guide)   */
    iparm[IPARM_ORDERING_CMIN]         = 0;                   /* Ordering cmin parameter  (see Scotch User's Guide)   */
    iparm[IPARM_ORDERING_CMAX]         = 100000;              /* Ordering cmax parameter  (see Scotch User's Guide)   */
    iparm[IPARM_ORDERING_FRAT]         = 8;                   /* Ordering frat parameter  (see Scotch User's Guide)   */

    /* Metis default */
#if defined(HAVE_METIS)
    iparm[IPARM_METIS_CTYPE   ] = METIS_CTYPE_SHEM;
    iparm[IPARM_METIS_RTYPE   ] = METIS_RTYPE_SEP1SIDED;
#endif
    iparm[IPARM_METIS_NO2HOP  ] = 0;
    iparm[IPARM_METIS_NSEPS   ] = 1;
    iparm[IPARM_METIS_NITER   ] = 10;
    iparm[IPARM_METIS_UFACTOR ] = 200;
    iparm[IPARM_METIS_COMPRESS] = 1;
    iparm[IPARM_METIS_CCORDER ] = 0;
    iparm[IPARM_METIS_PFACTOR ] = 0;
    iparm[IPARM_METIS_SEED    ] = 3452;
    iparm[IPARM_METIS_DBGLVL  ] = 0;

    iparm[IPARM_STATIC_PIVOTING]       = 0;                   /* number of control of diagonal magnitude              */
    iparm[IPARM_NNZEROS]               = 0;                   /* memory space for coefficients                        */
    iparm[IPARM_ALLOCATED_TERMS]       = 0;                   /* number of non zero in factorized sparse matrix       */
    iparm[IPARM_MIN_BLOCKSIZE]         = 60;                  /* min blocksize                                        */
    iparm[IPARM_MAX_BLOCKSIZE]         = 120;                 /* max blocksize                                        */
    iparm[IPARM_SCHUR]                 = API_NO;              /* Schur mode */
    iparm[IPARM_ISOLATE_ZEROS]         = API_NO;              /* Isolate null diagonal terms at the end of the matrix */
    iparm[IPARM_FACTORIZATION]         = API_FACT_LDLT;       /* LdLt     */
    iparm[IPARM_CPU_BY_NODE]           = 0;                   /* cpu/node */
    iparm[IPARM_BINDTHRD]              = API_BIND_AUTO;       /* Default binding method */
    iparm[IPARM_THREAD_NBR]            = 1;                   /* thread/mpi */
    iparm[IPARM_CUDA_NBR]              = 0;                   /* CUDA devices */
    iparm[IPARM_DISTRIBUTION_LEVEL]    = 0;                   /* 1d / 2d */
    iparm[IPARM_LEVEL_OF_FILL]         = 1;                   /* level of fill */
    iparm[IPARM_IO_STRATEGY]           = API_IO_NO;           /* I/O */
    iparm[IPARM_RHS_MAKING]            = API_RHS_B;           /* generate rhs */
    iparm[IPARM_REFINEMENT]            = API_RAF_GMRES;       /* gmres */
    iparm[IPARM_SYM]                   = API_SYM_YES;         /* Symmetric */
    iparm[IPARM_INCOMPLETE]            = API_NO;              /* direct */
    iparm[IPARM_ABS]                   = 1;                   /* ABS level to 1 */
    iparm[IPARM_ESP]                   = API_NO;              /* no esp */
#ifdef OOC
    iparm[IPARM_GMRES_IM]              = 1;                   /* gmres_im */
    iparm[IPARM_ITERMAX]               = 1;
#else
    iparm[IPARM_GMRES_IM]              = 25;                  /* gmres_im */
#endif
    iparm[IPARM_FREE_CSCUSER]          = API_CSC_PRESERVE;    /* Free user csc after coefinit */
    iparm[IPARM_FREE_CSCPASTIX]        = API_CSC_PRESERVE;    /* Free internal csc after coefinit */
    iparm[IPARM_OOC_LIMIT]             = 2000;                /* memory limit */
    iparm[IPARM_OOC_THREAD]            = 1;                   /* ooc thrdnbr */
    iparm[IPARM_OOC_ID]                = -1;                  /* Out of core run ID */
    iparm[IPARM_NB_SMP_NODE_USED]      = 0;                   /* Nb SMP node used (0 for 1 per MPI process) */
    iparm[IPARM_MURGE_REFINEMENT]      = API_YES;
    iparm[IPARM_TRANSPOSE_SOLVE]       = API_NO;
    /* Mode pour thread_comm :
     0 -> inutilisé
     1 -> 1 seul
     2 -> iparm[IPARM_NB_THREAD_COMM]
     3 -> Nbproc(iparm[IPARM_THREAD_NBR]))
     */
#ifdef PASTIX_FUNNELED
    iparm[IPARM_THREAD_COMM_MODE]      = API_THREAD_FUNNELED;
#else
    iparm[IPARM_THREAD_COMM_MODE]      = API_THREAD_MULTIPLE;
#endif
    iparm[IPARM_NB_THREAD_COMM]        = 1;                   /* Nb thread quand iparm[IPARM_THREAD_COMM_MODE] == API_THCOMM_DEFINED */
    iparm[IPARM_FILL_MATRIX]           = API_NO;              /* fill matrix */
    iparm[IPARM_INERTIA]               = -1;
    iparm[IPARM_ESP_NBTASKS]           = -1;
    iparm[IPARM_ESP_THRESHOLD]         = 16384;               /* Taille de bloc minimale pour passer en esp (2**14) = 128 * 128 */
    iparm[IPARM_ERROR_NUMBER]          = PASTIX_SUCCESS;
    iparm[IPARM_RHSD_CHECK]            = API_YES;
    iparm[IPARM_STARPU]                = API_NO;
    iparm[IPARM_AUTOSPLIT_COMM]        = API_NO;
#ifdef TYPE_COMPLEX
#  ifdef PREC_DOUBLE
    iparm[IPARM_FLOAT]                 = API_COMPLEXDOUBLE;
#  else
    iparm[IPARM_FLOAT]                 = API_COMPLEXSINGLE;
#  endif
#else
#  ifdef PREC_DOUBLE
    iparm[IPARM_FLOAT]                 = API_REALDOUBLE;
#  else
    iparm[IPARM_FLOAT]                 = API_REALSINGLE;
#  endif
#endif
    iparm[IPARM_STARPU_CTX_DEPTH]      = 3;
    iparm[IPARM_STARPU_CTX_NBR]        = -1;
    iparm[IPARM_PRODUCE_STATS]         = API_NO;
#ifdef PREC_DOUBLE
    dparm[DPARM_EPSILON_REFINEMENT] = 1e-12;
#else
    dparm[DPARM_EPSILON_REFINEMENT] = 1e-6;
#endif
    dparm[DPARM_RELATIVE_ERROR]     = -1;
    dparm[DPARM_SCALED_RESIDUAL]    = -1;
    dparm[DPARM_EPSILON_MAGN_CTRL]  = 1e-31;
    dparm[DPARM_FACT_FLOPS]         = 0;
    dparm[DPARM_SOLV_FLOPS]         = 0;
    dparm[63]                       = 0;

#ifdef PASTIX_LOG
    fprintf(stdout, "<- api_init\n");
#endif
}


/*
 Function: z_redispatch_rhs

 Redistribute right-hand-side member from *l2g*
 to *newl2g*

 Parameters:
 n        - Size of first right-hand-side member
 rhs      - Right-hand-side member
 l2g      - local to global column numbers
 newn     - New right-hand-side member size
 newrhs   - New right-hand-side member
 newl2g   - New local to global column numbers
 commSize - MPI communicator size
 commRank - MPI rank
 comm     - MPI communicator
 */
int z_redispatch_rhs(pastix_int_t      n,
                   pastix_complex64_t   *rhs,
                   pastix_int_t      nrhs,
                   pastix_int_t     *l2g,
                   pastix_int_t      newn,
                   pastix_complex64_t   *newrhs,
                   pastix_int_t     *newl2g,
                   int      commSize,
                   int      commRank,
                   MPI_Comm comm,
                   pastix_int_t dof)
{
    pastix_int_t           i,j, irhs;
#ifndef FORCE_NOMPI
    MPI_Status    status;
#endif
    pastix_int_t           gN = -1;
    pastix_int_t          *g2l;
    pastix_int_t         **count;
    MPI_Request * requests = NULL;
    pastix_int_t   ** toSendIdx;
    pastix_complex64_t ** toSendValues;

    z_cscd_build_g2l(newn,
                     newl2g,
                     comm,
                     &gN,
                     &g2l);

    /* allocate counter */
    MALLOC_INTERN(count, commSize, pastix_int_t*);
    for (i = 0; i < commSize; i++)
    {
        MALLOC_INTERN(count[i], commSize, pastix_int_t);
        for (j = 0; j < commSize; j++)
            count[i][j] = 0;
    }

    /* count how many entries we have to send
     to each processor
     complexity : n
     */
    for (i = 0; i < n; i++)
    {
        pastix_int_t globalidx = l2g[i];
        pastix_int_t localidx  = g2l[globalidx-1];
        if ( localidx > 0) {
            count[commRank][commRank] ++;
        }
        else {
            count[commRank][-localidx]++;
        }
    }

    /* Broadcast counters */
    for (i = 0; i < commSize; i++)
        MPI_Bcast(count[i], commSize, PASTIX_MPI_INT, i, comm);

    MALLOC_INTERN(toSendIdx, commSize, pastix_int_t *);
    MALLOC_INTERN(toSendValues, commSize, pastix_complex64_t *);
    for (irhs = 0; irhs < nrhs; irhs++)
    {
        /* Send data */
        for (i = 0; i < commSize; i++)
        {
            if (i != commRank) {
                MALLOC_INTERN(toSendIdx[i], count[commRank][i], pastix_int_t);
                MALLOC_INTERN(toSendValues[i], count[commRank][i]*dof, pastix_complex64_t);
                memset(toSendIdx[i], 0, count[commRank][i]*sizeof(pastix_int_t));
                memset(toSendValues[i], 0, count[commRank][i]*sizeof(pastix_complex64_t)*dof);
            }
        }

        for (i = 0; i < commSize; i++)
            count[commRank][i] = 0;

        for (i = 0; i < n; i++)
        {
            pastix_int_t globalidx = l2g[i];
            pastix_int_t localidx  = g2l[globalidx-1];
            if ( localidx > 0) {
                pastix_int_t d;
                for (d = 0; d < dof; d++)
                    newrhs[dof*(irhs*newn + localidx-1)+d] = rhs[dof*(irhs*n+i)+d];
            }
            else {
                pastix_int_t d;
                toSendIdx[-localidx][count[commRank][-localidx]] = globalidx;
                for (d = 0; d < dof; d++)
                    toSendValues[-localidx][dof*count[commRank][-localidx]+d] =
                        rhs[dof*(irhs*n+i)+d];
                count[commRank][-localidx]++;
            }
        }

        MALLOC_INTERN(requests, commSize, MPI_Request);

        for (i = 0; i < commSize; i++)
        {
            if (commRank != i)
            {
                MPI_Isend(toSendIdx[i], count[commRank][i],
                          PASTIX_MPI_INT, i, TAG_RHS, comm, &requests[i]);
                MPI_Isend(toSendValues[i], dof*count[commRank][i],
                          COMM_FLOAT, i, TAG_RHS2, comm, &requests[i]);
            }
        }

        for (i = 0; i < commSize; i++)
        {
            if (commRank != i)
            {
                pastix_int_t   * tmpIdx;
                pastix_complex64_t * tmpValues;
                MALLOC_INTERN(tmpIdx, count[i][commRank], pastix_int_t);
                MALLOC_INTERN(tmpValues, count[i][commRank]*dof, pastix_complex64_t);

                MPI_Recv(tmpIdx, count[i][commRank], PASTIX_MPI_INT,
                         i, TAG_RHS, comm, &status);
                MPI_Recv(tmpValues, dof*count[i][commRank], COMM_FLOAT,
                         i, TAG_RHS2, comm, &status);
                for (j= 0; j < count[i][commRank]; j++)
                {
                    pastix_int_t d;
                    for (d = 0; d < dof; d++)
                        newrhs[dof*(irhs*newn+g2l[tmpIdx[j]-1]-1)+d] = tmpValues[dof*j+d];
                }
                memFree_null(tmpIdx);
                memFree_null(tmpValues);
            }
        }

        for (i = 0; i < commSize; i++)
        {
            if (i != commRank)
            {
                MPI_Wait(&requests[i],&status);
                memFree_null(toSendIdx[i]);
                memFree_null(toSendValues[i]);
            }
        }
        memFree_null(requests);

    }
    memFree_null(toSendIdx);
    memFree_null(toSendValues);

    for (i = 0; i < commSize; i++)
        memFree_null(count[i]);
    memFree_null(count);
    memFree_null(g2l);
    return PASTIX_SUCCESS;
}

/*
 Function: z_buildUpDoVect

 Build z_UpDownVector from user B vector or
 computes its to obtain $$ \forall i X[i] = 1 $$ or $$ \forall i X[i] = i $$
 as the solution. (depending on iparm)

 Parameters:
 pastix_data - PaStiX global data structure.
 loc2glob2   - Global  column number of local columns.
 b           - User right-hand-side member.
 pastix_comm - MPI communicator.
 */
int z_buildUpdoVect(z_pastix_data_t    *pastix_data,
                    pastix_int_t       *loc2glob,
                    pastix_complex64_t *b,
                    MPI_Comm            pastix_comm)
{
    pastix_int_t           *iparm    = pastix_data->iparm;
    z_SolverMatrix         *solvmatr = &(pastix_data->solvmatr);
    Order                  *ordemesh = pastix_data->ordemesh;
    pastix_int_t            procnum  = pastix_data->procnum;
    pastix_int_t           *invp     = ordemesh->peritab;
    (void)loc2glob;

    /* Rhs taking from b */
    if (iparm[IPARM_RHS_MAKING]==API_RHS_B)
    {
        /* Using b */
        if (solvmatr->updovct.sm2xsze > 0 &&  b==NULL)
        {
            errorPrint("b must be allocated.");
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }
        else
        {
            /* Permuter b avec la permutation inverse */
            if (iparm[IPARM_GRAPHDIST] == API_NO )
            {
                z_CscUpdownRhs(&(solvmatr->updovct),
                             solvmatr,
                             b,
                             invp,
                             (int)iparm[IPARM_DOF_NBR]);
            }
#ifdef PASTIX_DISTRIBUTED
            else
            {
                z_CscdUpdownRhs(&(solvmatr->updovct),
                              solvmatr,
                              b,
                              invp,
                              pastix_data->glob2loc,
                              pastix_data->n2,
                              (int)iparm[IPARM_DOF_NBR]);
            }
#endif /* PASTIX_DISTRIBUTED */
        }
    }
    /* Generate rhs */
    else
    {
        if (procnum == 0)
        {
            if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
            {
                if (iparm[IPARM_RHS_MAKING]==API_RHS_1)
                    fprintf(stdout,GEN_RHS_1);
                else
                    fprintf(stdout,GEN_RHS_I);
            }
        }

        /* In updo step, if we only want to use reffinement.
         Set first solution.
         */
        if ((iparm[IPARM_ONLY_RAFF] == API_YES) && (iparm[IPARM_START_TASK] < API_TASK_REFINE))
        {
            fprintf(stdout,GEN_SOL_0);
            z_Csc2updown_X0(&(solvmatr->updovct),
                          solvmatr,
                          API_RHS_0,
                          pastix_comm);
        }
        else
        {
            z_Csc2updown(&(pastix_data->cscmtx),
                         &(solvmatr->updovct),
                         solvmatr,
                         iparm[IPARM_RHS_MAKING],
                         pastix_comm);
        }
    }
    return PASTIX_SUCCESS;
}

/*
 Function: z_global2localrhs

 Converts global right hand side to local right hand side.

 Parameters:
 lN       - local number of columns.
 lrhs     - local right hand side.
 grhs     - global right hand side.
 loc2glob - global index of local columns.
 */

void z_global2localrhs(pastix_int_t     lN,
                     pastix_complex64_t *lrhs,
                     pastix_complex64_t *grhs,
                     pastix_int_t   *loc2glob)
{
    pastix_int_t   i;

    for (i = 0; i < lN; i++)
        lrhs[i] = grhs [loc2glob[i]-1];
}

/*
 Function: z_global2localperm

 Converts global permutation (resp. reverse permutation) tabular to local
 permutation (resp. reverse permutation) tabular.

 Parameters:
 lN       - local number of columns.
 lperm    - local permutation tabular.
 gperm    - global permutation tabular.
 loc2glob - global index of local columns.
 */
void z_global2localperm(pastix_int_t  lN,
                      pastix_int_t *lperm,
                      pastix_int_t *gperm,
                      pastix_int_t *loc2glob)
{
    pastix_int_t   i;
    for (i = 0; i < lN; i++)
    {
        lperm[i] = gperm[loc2glob[i]-1];
    }
}

/*
 * Function: z_pastix_task_init
 *
 * Allocate and fill-in pastix_data
 *
 * Parameters:
 *   pastix_data - structure to build
 *   pastix_comm - PaStiX MPI communicator
 *   iparm       - integer parameters, to fill-in pastix_data
 *   dparm       - floating parameters, to fill-in pastix_data
 */
void z_pastix_task_init(z_pastix_data_t **pastix_data,
                      MPI_Comm        pastix_comm,
                      pastix_int_t   *iparm,
                      double         *dparm)
{

    /* Allocate pastix_data structure when we enter PaStiX for the first time.
     */
    MALLOC_INTERN(*pastix_data, 1, z_pastix_data_t);
    memset( *pastix_data, 0, sizeof(z_pastix_data_t) );

    /* Initialisation des champs de la structure */
    (*pastix_data)->n2               = -1;
    (*pastix_data)->malcsc           = 0;
    (*pastix_data)->malsmx           = 0;
    (*pastix_data)->malslv           = 0;
    (*pastix_data)->malcof           = 0;
    (*pastix_data)->iparm            = iparm;
    (*pastix_data)->dparm            = dparm;
    (*pastix_data)->pastix_comm      = pastix_comm;
    if (iparm != NULL && iparm[IPARM_AUTOSPLIT_COMM] == API_YES)
    {
        int i,len, mpisize;
        char procname[MPI_MAX_PROCESSOR_NAME];
        int color, key;

        MPI_Get_processor_name(procname,&len);
        MPI_Comm_rank(pastix_comm, &(key));
        color = 0;
        for (i = 0; i < len; i++)
            color = color*256*sizeof(char) + procname[i];
        MPI_Comm_split(pastix_comm, color, key, &(*pastix_data)->intra_node_comm);
        MPI_Comm_rank((*pastix_data)->intra_node_comm, &color);
        MPI_Comm_rank((*pastix_data)->intra_node_comm, &(key));
        MPI_Comm_split(pastix_comm, color, key, &(*pastix_data)->inter_node_comm);
        MPI_Comm_size((*pastix_data)->intra_node_comm, &mpisize);
        iparm[IPARM_THREAD_NBR] = mpisize;
    }
    else
    {
        (*pastix_data)->inter_node_comm      = pastix_comm;
        (*pastix_data)->intra_node_comm      = MPI_COMM_SELF;
    }

    (*pastix_data)->sopar.bindtab    = NULL;
    (*pastix_data)->sopar.b          = NULL;
    (*pastix_data)->sopar.transcsc   = NULL;
    (*pastix_data)->sopar.stopthrd   = API_NO;
    (*pastix_data)->bindtab          = NULL;
    (*pastix_data)->cscInternFilled  = API_NO;

#ifdef PASTIX_DISTRIBUTED
    (*pastix_data)->malrhsd_int      = API_NO;
    (*pastix_data)->l2g_int          = NULL;
    (*pastix_data)->mal_l2g_int      = API_NO;
    (*pastix_data)->glob2loc         = NULL;
    (*pastix_data)->PTS_permtab      = NULL;
    (*pastix_data)->PTS_peritab      = NULL;
#endif
    (*pastix_data)->schur_n          = 0;
    (*pastix_data)->schur_list       = NULL;
    (*pastix_data)->schur_tab        = NULL;
    (*pastix_data)->schur_tab_set    = API_NO;

    (*pastix_data)->solvmatr.updovct.cblktab    = NULL;
    (*pastix_data)->solvmatr.updovct.lblk2gcblk = NULL;
    (*pastix_data)->solvmatr.updovct.listblok   = NULL;
    (*pastix_data)->solvmatr.updovct.listcblk   = NULL;
    (*pastix_data)->solvmatr.updovct.gcblk2list = NULL;
    (*pastix_data)->solvmatr.updovct.loc2glob   = NULL;
    (*pastix_data)->solvmatr.updovct.cblktab    = NULL;
    (*pastix_data)->solvmatr.updovct.listptr    = NULL;

    (*pastix_data)->scaling  = API_NO;
    (*pastix_data)->scalerowtab = NULL;
    (*pastix_data)->iscalerowtab = NULL;
    (*pastix_data)->scalecoltab = NULL;
    (*pastix_data)->iscalecoltab = NULL;

    /* Récupération du nombre de proc */
    MPI_Comm_size(pastix_comm, &((*pastix_data)->procnbr));
    MPI_Comm_rank(pastix_comm, &((*pastix_data)->procnum));
    MPI_Comm_size((*pastix_data)->inter_node_comm, &((*pastix_data)->inter_node_procnbr));
    MPI_Comm_rank((*pastix_data)->inter_node_comm, &((*pastix_data)->inter_node_procnum));
    MPI_Comm_size((*pastix_data)->intra_node_comm, &((*pastix_data)->intra_node_procnbr));
    MPI_Comm_rank((*pastix_data)->intra_node_comm, &((*pastix_data)->intra_node_procnum));
    if ((*pastix_data)->procnum == 0)
    {
        (*pastix_data)->pastix_id = getpid();
    }
    MPI_Bcast(&((*pastix_data)->pastix_id), 1, PASTIX_MPI_INT, 0, pastix_comm);

#ifdef WITH_SEM_BARRIER
    if ((*pastix_data)->intra_node_procnbr > 1)
    {
        char sem_name[256];
        sprintf(sem_name, "/pastix_%d", (*pastix_data)->pastix_id);
        OPEN_SEM((*pastix_data)->sem_barrier, sem_name, 0);
    }
#endif

    if (iparm != NULL)
    {
        if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        {
            fprintf(stdout, "AUTOSPLIT_COMM : global rank : %d,"
                    " inter node rank %d,"
                    " intra node rank %d, threads %d\n",
                    (int)((*pastix_data)->procnum),
                    (int)((*pastix_data)->inter_node_procnum),
                    (int)((*pastix_data)->intra_node_procnum),
                    (int)iparm[IPARM_THREAD_NBR]);
        }

        iparm[IPARM_PID] = (*pastix_data)->pastix_id;
    }
    /* DIRTY Initialization for Scotch */
    srand(1);

    /* Environement variables */
    /* On Mac set VECLIB_MAXIMUM_THREADS if not setted */
    setenv("VECLIB_MAXIMUM_THREADS", "1", 0);

    (*pastix_data)->steps = STEP_INIT;
}

#ifdef MEMORY_USAGE
/*
 Function: z_pastix_print_memory_usage

 print memory usage during z_pastix.

 Parameters:
 iparm       - integer paréameters.
 pastix_comm - PaStiX MPI communicator
 */
void z_pastix_print_memory_usage(pastix_int_t      *iparm,
                                 MPI_Comm  pastix_comm)
{
    unsigned long smem[2], rmem[2];
    int           procnum;

    MPI_Comm_rank(pastix_comm,&procnum);
    smem[0] = memAllocGetMax();
    smem[1] = memAllocGetCurrent();
    MPI_Reduce(smem,rmem,2,MPI_LONG,MPI_MAX,0,pastix_comm);
    iparm[DPARM_MEM_MAX] = rmem[0];
    if (procnum == 0)
        if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        {
            fprintf(stdout,MAX_MEM_AF_CL,
                    MEMORY_WRITE(rmem[0]),
                    MEMORY_UNIT_WRITE(rmem[0]));

            fprintf(stdout,MEM_USED_AF_CL,
                    MEMORY_WRITE(rmem[1]),
                    MEMORY_UNIT_WRITE(rmem[1]));
        }
}
#else
#  ifdef z_pastix_print_memory_usage
#    undef z_pastix_print_memory_usage
#  endif
#  define z_pastix_print_memory_usage(pastix_data,pastix_comm)
#endif /* MEMORY_USAGE */

/*
 * Function: z_pastix_welcome_print
 *
 * Will print welcome message, options and parameters.
 *
 * Parameters:
 * pastix_data - PaStiX data structure
 * colptr      - starting index of each column in the CSC.
 * n           - number of columns.
 *
 */
void z_pastix_welcome_print(z_pastix_data_t *pastix_data,
                          pastix_int_t           *colptr,
                          pastix_int_t            ln)
{
    pastix_int_t    * iparm = pastix_data->iparm;
    double * dparm = pastix_data->dparm;
    pastix_int_t      gN    = 0;
    pastix_int_t      lNnz  = 0;
    pastix_int_t      gNnz  = 0;

    if (colptr != NULL)
        lNnz = colptr[ln]-colptr[0];

    if (iparm[IPARM_GRAPHDIST] == API_YES)
    {
        gN = 0;
        MPI_Allreduce(&ln, &gN, 1, PASTIX_MPI_INT, MPI_SUM,
                      pastix_data->pastix_comm);
        MPI_Reduce(&lNnz, &gNnz, 1, PASTIX_MPI_INT, MPI_SUM, 0,
                   pastix_data->pastix_comm);
    }
    else
    {
        gN   = ln;
        gNnz = lNnz;
    }
    dparm[DPARM_FILL_IN]       = 1.0/(double)gNnz;

    if ((iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) &&
        (pastix_data->procnum == 0 ))
    {
        fprintf(stdout, OUT_ENTETE_LINE1);
        fprintf(stdout, OUT_ENTETE_LINE2);
        fprintf(stdout, OUT_ENTETE_LINE3);

        /* TODO : en distribué afficher pour chaque proc... */
        if ( gN != 0)
        {
            fprintf(stdout, OUT_MATRIX_SIZE, (long)gN, (long)gN);
            if (gNnz != 0) fprintf(stdout, OUT_NNZ, (long)gNnz);
        }
        z_sopalin_option();
    }
}

/* Function: z_sopalin_check_param

 Check parameters consistency.

 Parameters:
 pastix_data - PaStiX data structure.

 Return:
 PASTIX_SUCCESS           - if no error occured
 BADPARAMETER_ERR - if Parameters are not correct on one proc.
 */
#define z_sopalin_check_param PASTIX_PREFIX_F(z_sopalin_check_param)
static inline
int z_sopalin_check_param(z_pastix_data_t *pastix_data)
{
    pastix_int_t           * iparm    = pastix_data->iparm;
    int             ret      = PASTIX_SUCCESS;
    int             ret_recv = PASTIX_SUCCESS;

    if ((iparm[IPARM_SYM]           == API_SYM_NO) &&
        (iparm[IPARM_FACTORIZATION] != API_FACT_LU))
    {
        errorPrint("With unsymmetric patterns LU decomposition should be used");
        ret = BADPARAMETER_ERR;
    }
#ifndef TYPE_COMPLEX
    if (iparm[IPARM_FACTORIZATION] == API_FACT_LDLH)
    {
        errorPrint("LDLH only available with complex");
        ret = BADPARAMETER_ERR;
    }
#endif
    MPI_Allreduce(&ret, &ret_recv, 1, MPI_INT, MPI_MAX, pastix_data->inter_node_comm);
    return ret_recv;
}


/* Function: z_pastix_check_param

 Check parameters consistency.

 Parameters:
 pastix_data - PaStiX data structure.

 Return:
 PASTIX_SUCCESS           - if no error occured
 BADPARAMETER_ERR - if Parameters are not correct on one proc.
 */
#define z_pastix_check_param PASTIX_PREFIX_F(z_pastix_check_param)
static inline
int z_pastix_check_param(z_pastix_data_t * pastix_data, int rhsnbr)
{
    int ret = PASTIX_SUCCESS, ret_recv;
    pastix_int_t * iparm = pastix_data->iparm;

    if (iparm[IPARM_THREAD_NBR] < 1) {
        errorPrint("iparm[IPARM_THREAD_NBR] must be strictly positive.");
        ret = BADPARAMETER_ERR;
    }

    if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD_CSC) ||
        PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD_GRAPH))
        iparm[IPARM_IO_STRATEGY] |= API_IO_LOAD;

    if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE_CSC) ||
        PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE_GRAPH))
        iparm[IPARM_IO_STRATEGY] |= API_IO_SAVE;

    if ((PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD)) &&
        (iparm[IPARM_ORDERING]    != API_ORDER_LOAD))
        iparm[IPARM_ORDERING]  = API_ORDER_LOAD;

    if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD) &&
        PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE))
    {
        errorPrint("Save and load strategy are not compatible\n");
        ret = BADPARAMETER_ERR;
    }


#ifndef MULT_SMX
    if (iparm[IPARM_START_TASK] < API_TASK_CLEAN &&
        iparm[IPARM_END_TASK] > API_TASK_NUMFACT &&
        rhsnbr > 1)
    {
        errorPrint("Need to rebuild with -DMULT_SMX to use multiple RHS\n");
        ret = BADPARAMETER_ERR;
    }
#endif
    if (iparm[IPARM_START_TASK] < API_TASK_CLEAN &&
        iparm[IPARM_END_TASK] > API_TASK_NUMFACT &&
        rhsnbr < 1)
    {
        errorPrint("Number of RHS must be a positive number\n");
        ret = BADPARAMETER_ERR;
    }

    if (iparm[IPARM_STARPU] == API_YES)
    {
#ifndef PASTIX_WITH_STARPU
        errorPrintW("To use StarPU scheduler please build"
                    " PaStiX with -DPASTIX_WITH_STARPU=ON");
        ret = BADPARAMETER_ERR;
#endif
    }

    if (iparm[IPARM_START_TASK] <= API_TASK_SOLVE  &&
        iparm[IPARM_END_TASK]   >= API_TASK_REFINE &&
        iparm[IPARM_ONLY_RAFF]  == API_YES) {
        errorPrintW("IPARM_ONLY_RAFF ignored, only possible if UPDO and RAFF are called in 2 steps");
        iparm[IPARM_ONLY_RAFF] = API_NO;
    }


    if (iparm[IPARM_FACTORIZATION] == API_FACT_LU &&
        iparm[IPARM_ONLY_RAFF]     == API_YES     &&
        iparm[IPARM_REFINEMENT]    == API_RAF_PIVOT)
    {
        errorPrint("Reffinement only is not compatible with the"
                   " simple iterative refinement.");
        ret = BADPARAMETER_ERR;
    }
    MPI_Allreduce(&ret, &ret_recv, 1, MPI_INT, MPI_MAX, pastix_data->pastix_comm);
    return ret_recv;
}

#ifdef COMPUTE
/*
 * Function: z_pastix_fake_fillin_csc
 *
 * Fill in the internal csc based on the user csc and fill in the coeftab structure
 *
 * Parameters:
 * pastix_data - PaStiX data structure.
 * pastix_comm - PaStiX MPI communicator.
 * n           - Size of the matrix.
 * colptr      - starting index of each column in row and avals.
 * row         - row of each element of the matrix.
 * avals       - value of each element of the matrix.
 * b           - Right hand side.
 * nrhs        - Number of right-hand-sides.
 * loc2glob    - global number of local columns, NULL if not distributed.
 */
int z_pastix_fake_fillin_csc( z_pastix_data_t *pastix_data,
                            MPI_Comm       pastix_comm,
                            pastix_int_t     n,
                            pastix_int_t    *colptr,
                            pastix_int_t    *row,
                            pastix_complex64_t  *avals,
                            pastix_complex64_t  *b,
                            pastix_int_t     nrhs,
                            pastix_int_t    *loc2glob)
{
#  ifdef PASTIX_DISTRIBUTED
    pastix_int_t     *iparm    = pastix_data->iparm;
    pastix_int_t     *l_colptr = NULL;
    pastix_int_t     *l_row    = NULL;
    pastix_complex64_t   *l_val    = NULL;
    pastix_int_t     *l_l2g    = NULL;
    pastix_complex64_t   *l_b      = NULL;
    pastix_int_t      l_n      = 0;
    int             mal_l_l2g = API_NO;
    int             mal_l_b   = API_NO;
    int             OK       = 0;
    int             OK_RECV  = 0;
    int             retval = PASTIX_SUCCESS;
    int             retval_recv;
    pastix_int_t      iter;
    pastix_int_t      gN = -1;
    (void)pastix_comm; (void)nrhs;

    if (iparm[IPARM_GRAPHDIST] == API_YES)
    {
        /* If user has not specified that he is
         absolutly certain that is CSCd is
         correctly distributed */
        if (iparm[IPARM_CSCD_CORRECT] == API_NO)
        {
            /* Test que la cscd utilisateur correspond a la cscd z_pastix */
            l_n = 0;
            l_l2g = NULL;

            OK = 0;
            if (l_n != n)
            {
                OK = 1;
            }
            else
            {
                for (iter = 0; iter < l_n; iter++)
                {
                    if (l_l2g[iter] != loc2glob[iter])
                    {
                        OK = 1;
                        break;
                    }
                }
            }
            MPI_Allreduce(&OK, &OK_RECV, 1, MPI_INT, MPI_SUM, pastix_comm);
        }

        if (NULL != pastix_data->glob2loc)
            memFree_null(pastix_data->glob2loc);

        /* Building global to local column number correspondance */
        z_cscd_build_g2l(l_n,
                         l_l2g,
                         pastix_comm,
                         &gN,
                         &(pastix_data->glob2loc));

    }

    /* Si la cscd correspond pas, on la corrige */
    if (OK_RECV > 0)
    {
        if (pastix_data->procnum == 0 &&
            iparm[IPARM_VERBOSE] >= API_VERBOSE_YES)
            fprintf(stdout,OUT_REDIS_CSC);

        /* redistributing cscd with correct local to global correspondance */
        retval = z_cscd_redispatch_int(  n,    colptr,    row,  avals,    b, nrhs, loc2glob,
                                         l_n, &l_colptr, &l_row, &l_val, &l_b, l_l2g,
                                         API_YES, pastix_comm, iparm[IPARM_DOF_NBR]);
        memFree_null(l_colptr); /* in fake fillin we do nothing with that */
        MPI_Allreduce(&retval, &retval_recv, 1, MPI_INT, MPI_MAX, pastix_comm);
        if (retval_recv != PASTIX_SUCCESS)
            return retval_recv;

        mal_l_b = API_YES;
    }
    else
    {
        /* the user CSCd is correctly distributed,
         only the pointers to the user CSCd arrays are used.
         */
        if ((iparm[IPARM_CSCD_CORRECT] == API_NO)
            && (l_l2g != NULL))
            memFree_null(l_l2g);

        l_n      = n;
        l_colptr = colptr;
        l_row    = row;
        l_val    = avals;
        l_l2g    = loc2glob;
        mal_l_l2g = API_NO;
        mal_l_b   = API_NO;
        l_b      = b;
    }

    if (pastix_data->mal_l2g_int == API_YES)
        memFree_null(pastix_data->l2g_int);
    if (pastix_data->malrhsd_int == API_YES)
        memFree_null(pastix_data->b_int);

    pastix_data->ncol_int = l_n;
    pastix_data->l2g_int  = l_l2g;
    pastix_data->b_int    = l_b;
    pastix_data->mal_l2g_int = mal_l_l2g;
    pastix_data->malrhsd_int = mal_l_b;
#else
    (void)pastix_data; (void)pastix_comm; (void)n; (void)colptr; (void)row; (void)avals; (void)b;
    (void)nrhs; (void)loc2glob;
#endif

    return PASTIX_SUCCESS;
}

/*
 *  Function: z_pastix_fillin_csc
 *
 *  Fill in the internal csc based on the user csc and fill in the coeftab structure
 *
 *  Parameters:
 *  pastix_data - PaStiX data structure.
 *  pastix_comm - PaStiX MPI communicator.
 *  n           - Size of the matrix.
 *  colptr      - starting index of each column in row and avals.
 *  row         - row of each element of the matrix.
 *  avals       - value of each element of the matrix.
 *  b           - Right hand side.
 *  nrhs        - Number of right-hand-sides.
 *  loc2glob    - global number of local columns, NULL if not distributed.
 */
int z_pastix_fillin_csc( z_pastix_data_t *pastix_data,
                       MPI_Comm       pastix_comm,
                       pastix_int_t            n,
                       pastix_int_t           *colptr,
                       pastix_int_t           *row,
                       pastix_complex64_t         *avals,
                       pastix_complex64_t         *b,
                       pastix_int_t            nrhs,
                       pastix_int_t           *loc2glob)
{
    pastix_int_t            *iparm    = pastix_data->iparm;
    z_SolverMatrix   *solvmatr = &(pastix_data->solvmatr);
    Order          *ordemesh = pastix_data->ordemesh;
    pastix_int_t            *l_colptr = NULL;
    pastix_int_t            *l_row    = NULL;
    pastix_complex64_t          *l_val    = NULL;
    pastix_int_t            *l_l2g    = loc2glob;
    pastix_complex64_t          *l_b      = NULL;
    pastix_int_t             l_n      = n;
    pastix_complex64_t         **transcsc = NULL;
    pastix_int_t             procnum  = pastix_data->procnum;
    int             malcsc   = 0;
    int             forcetr  = 0;
    char            Type[4];
    double          clk;
    int             mal_l_l2g = API_NO;
    int             mal_l_b   = API_NO;
    (void)pastix_comm;
    (void)nrhs;
#  ifdef PASTIX_DISTRIBUTED
    int             OK       = 0;
    int             OK_RECV  = 0;
    int             retval = PASTIX_SUCCESS;
    int             retval_recv;
    pastix_int_t             iter;
    pastix_int_t             gN = -1;

    if (iparm[IPARM_GRAPHDIST] == API_YES)
    {
        /* If user has not specified that he is
         absolutly certain that is CSCd is
         correctly distributed */
        if (iparm[IPARM_CSCD_CORRECT] == API_NO)
        {
            /* Test que la cscd utilisateur correspond a la cscd z_pastix */
            l_n = z_pastix_getLocalNodeNbr(&pastix_data);
            MALLOC_INTERN(l_l2g, l_n, pastix_int_t);
            if (l_l2g != NULL) mal_l_l2g = API_YES;
            z_pastix_getLocalNodeLst(&pastix_data, l_l2g);

            OK = 0;
            if (l_n != n)
            {
                OK = 1;
            }
            else
            {
                for (iter = 0; iter < l_n; iter++)
                {
                    if (l_l2g[iter] != loc2glob[iter])
                    {
                        OK = 1;
                        break;
                    }
                }
            }
            MPI_Allreduce(&OK, &OK_RECV, 1, MPI_INT, MPI_SUM, pastix_comm);
        }

        if (NULL != pastix_data->glob2loc)
            memFree_null(pastix_data->glob2loc);

        /* Building global to local column number correspondance */
        z_cscd_build_g2l(l_n,
                         l_l2g,
                         pastix_comm,
                         &gN,
                         &(pastix_data->glob2loc));

    }

    /* Si la cscd correspond pas, on la corrige */
    if (OK_RECV > 0)
    {
        if (procnum == 0 &&
            iparm[IPARM_VERBOSE] >= API_VERBOSE_YES)
            fprintf(stdout,OUT_REDIS_CSC);

        /* redistributing cscd with correct local to global correspondance */
        clockInit(clk);
        clockStart(clk);
        retval = z_cscd_redispatch_int(  n,    colptr,    row,  avals,    b, nrhs, loc2glob,
                                       l_n, &l_colptr, &l_row, &l_val, &l_b, l_l2g,
                                       API_YES, pastix_comm, iparm[IPARM_DOF_NBR]);
        clockStop((clk));

        if (iparm[IPARM_VERBOSE] >= API_VERBOSE_YES)
            print_onempi(OUT_REDISCSCDTIME, (double)clockVal(clk));

        MPI_Allreduce(&retval, &retval_recv, 1, MPI_INT, MPI_MAX, pastix_comm);
        if (retval_recv != PASTIX_SUCCESS)
            return retval_recv;

        malcsc = API_YES;
        mal_l_b = API_YES;
    }
    else
#  endif /* PASTIX_DISTRIBUTED */
    {
        /* the user CSCd is correctly distributed,
         only the pointers to the user CSCd arrays are used.
         */
        malcsc = API_NO;
        if ((iparm[IPARM_CSCD_CORRECT] == API_NO)
            && (l_l2g != NULL))
            memFree_null(l_l2g);

        l_n      = n;
        l_colptr = colptr;
        l_row    = row;
        l_val    = avals;
        l_l2g    = loc2glob;
        mal_l_l2g = API_NO;
        mal_l_b   = API_NO;
        l_b      = b;
    }

#  ifdef PASTIX_DISTRIBUTED
    if (pastix_data->mal_l2g_int == API_YES)
        memFree_null(pastix_data->l2g_int);
    if (pastix_data->malrhsd_int == API_YES)
        memFree_null(pastix_data->b_int);

    pastix_data->ncol_int = l_n;
    pastix_data->l2g_int  = l_l2g;
    pastix_data->b_int    = l_b;
    pastix_data->mal_l2g_int = mal_l_l2g;
    pastix_data->malrhsd_int = mal_l_b;
#  endif

    Type[0] = 'R';
    Type[1] = 'S';
    Type[2] = 'A';
    Type[3] = '\0';

    /* Fill in of the internal CSC */
    if (iparm[IPARM_FILL_MATRIX] == API_NO) /* not false factorisation */
    {
        if (pastix_data->malcsc)
        {
            z_CscExit(&(pastix_data->cscmtx));
            pastix_data->malcsc=0;
        }

        clockInit(clk);
        clockStart(clk);

        /* Choix des parametres pour z_CscOrdistrib */
        if (iparm[IPARM_SYM] == API_SYM_YES || iparm[IPARM_SYM] == API_SYM_HER) /* symmetric mtx */
        {
            if (iparm[IPARM_FACTORIZATION] == API_FACT_LU) /* LU */
            {
                /* Lu on RSA */
                forcetr = 1;
                Type[1] = 'S';
            }
            else
            {
                forcetr = 0;
                Type[1] = 'S';
            }
            if (iparm[IPARM_SYM] == API_SYM_HER)
                Type[1] = 'H';
        }
        else
        {
            forcetr = 0;
            Type[1] = 'U';
        }
        transcsc = &(pastix_data->sopar.transcsc);

        if (iparm[IPARM_ONLY_RAFF] == API_YES)
            transcsc = NULL;

        /* Build internal CSCD from user CSC */
        if (iparm[IPARM_GRAPHDIST] == API_NO)
        {
            z_CscOrdistrib(&(pastix_data->cscmtx), Type,
                         transcsc, ordemesh,
                         l_n, l_n, l_colptr[l_n]-1, l_colptr,
                         l_row, l_val, forcetr,
                         solvmatr, procnum, iparm[IPARM_DOF_NBR]);
        }
#  ifdef PASTIX_DISTRIBUTED
        else
        {
            /* Build internal CSCD from user CSCD */
            z_CscdOrdistrib(&(pastix_data->cscmtx), Type,
                          transcsc, ordemesh,
                          l_n, l_colptr,
                          l_row, l_val,
                          l_l2g,
                          gN, pastix_data->glob2loc,
                          forcetr, solvmatr, procnum,
                          iparm[IPARM_DOF_NBR], pastix_data->inter_node_comm);
        }
#  endif /* PASTIX_DISTRIBUTED */
        pastix_data->malcsc = 1;

        /* Libération de la csc interne temporaire */
        if (malcsc == API_YES)
        {
            /* l2g and rhs are not freed here beacause
             they will be used for solution filling */
            memFree_null(l_colptr);
            memFree_null(l_row);
            memFree_null(l_val);
        }

        clockStop(clk);

        if (iparm[IPARM_VERBOSE] >= API_VERBOSE_YES)
            print_onempi(OUT_FILLCSCTIME, (double)clockVal(clk));

        /* User Csc is useless after cscordistrib */
        if (iparm[IPARM_FREE_CSCUSER] == API_CSC_FREE)
        {
            free(colptr);
            free(row);
            free(avals);
            colptr = NULL;
            row    = NULL;
            avals  = NULL;
        }
    }

    if (pastix_data->malcof)
    {
        if (iparm[IPARM_SCHUR] == API_YES && pastix_data->schur_tab_set == API_YES)
        {
            z_SolverMatrix * datacode = &(pastix_data->solvmatr);
            pastix_int_t            cblk;

            if (SOLV_TASKNBR > 0)
            {
                cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
                if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
                {
                    SOLV_COEFTAB(cblk) = NULL;
                }
            }
        }
        z_CoefMatrix_Free( &(pastix_data->sopar), solvmatr, iparm[IPARM_FACTORIZATION]);
        pastix_data->malcof=0;
    }

    /* On alloue par bloc colonne et pas tte la matrice */
    /* L'initialisation de la matrice est faite dans z_sopalin_init_smp              */
    pastix_data->sopar.iparm = iparm;
    if (iparm[IPARM_ONLY_RAFF] == API_NO)
    {
        if (iparm[IPARM_SCHUR] == API_YES && pastix_data->schur_tab_set == API_YES)
        {
            z_SolverMatrix * datacode = &(pastix_data->solvmatr);
            pastix_int_t            cblk;

            if (SOLV_TASKNBR > 0)
            {
                cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
                if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
                {
                    SOLV_COEFTAB(cblk) = pastix_data->schur_tab;
                }
            }
        }

        pastix_data->malcof = 1;
    }

    return PASTIX_SUCCESS;
}
#else
#  define z_pastix_fillin_csc(ptx_data, comm, n, col, row, a, b, nrhs, l2g) PASTIX_SUCCESS
#endif /* COMPUTE */

/*
 Function: z_pastix_task_sopalin

 Factorisation, updown and raffinement tasks.

 Parameters:
 pastix_data - PaStiX data structure.
 pastix_comm - PaStiX MPI communicator.
 n           - Size of the matrix.
 colptr      - starting index of each column in row and avals.
 row         - row of each element of the matrix.
 avals       - value of each element of the matrix.
 b           - Right hand side.
 loc2glob    - global number of local columns, NULL if not distributed.
 */
int z_pastix_task_sopalin( z_pastix_data_t *pastix_data,
                         MPI_Comm       pastix_comm,
                         pastix_int_t            n,
                         pastix_int_t           *colptr,
                         pastix_int_t           *row,
                         pastix_complex64_t         *avals,
                         pastix_complex64_t         *b,
                         pastix_int_t            rhsnbr,
                         pastix_int_t           *loc2glob)
{
    long            spivot, rpivot;
    double          sfacttime, rfacttime;
    double          ssolvtime,rsolvtime;
    double          srafftime,rrafftime;
    pastix_int_t           * iparm    = pastix_data->iparm;
    double        * dparm    = pastix_data->dparm;
    z_SolverMatrix  * solvmatr = &(pastix_data->solvmatr);
    Order         * ordemesh = pastix_data->ordemesh;
    z_SopalinParam  * sopar    = &(pastix_data->sopar);
    pastix_int_t             procnum  = pastix_data->inter_node_procnum;
    int             ret = 0;

    print_debug(DBG_STEP, "->API_TASK_NUMFACT\n");
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    {
        switch(iparm[IPARM_FACTORIZATION])
        {
        case API_FACT_LU:
            print_onempi("%s", OUT_STEP_NUMFACT_LU);
            break;
        case API_FACT_LLT:
            print_onempi("%s", OUT_STEP_NUMFACT_LLT);
            break;
        case API_FACT_LDLH:
            print_onempi("%s", OUT_STEP_NUMFACT_LDLH);
            break;
        case API_FACT_LDLT:
        default:
            print_onempi("%s", OUT_STEP_NUMFACT_LDLT);
        }
    }

    /* the error value has been reduced */
    if (PASTIX_SUCCESS != (ret = z_sopalin_check_param(pastix_data)))
        return ret;

    /* Remplissage de la csc interne */
    if (pastix_data->cscInternFilled == API_NO)
    {
        z_pastix_fillin_csc(pastix_data, (pastix_data)->pastix_comm, n,
                          colptr, row, avals, b, rhsnbr, loc2glob);
    }
#ifdef PASTIX_DISTRIBUTED
    else
    {
        pastix_data->l2g_int = loc2glob;
        pastix_data->b_int   = b;
    }
#endif
    /* sopalin */
    sopar->cscmtx      = &(pastix_data->cscmtx);
    sopar->itermax     = iparm[IPARM_ITERMAX];
    sopar->diagchange  = 0;
    sopar->epsilonraff = dparm[DPARM_EPSILON_REFINEMENT];
    sopar->rberror     = 0;
    sopar->espilondiag = dparm[DPARM_EPSILON_MAGN_CTRL];
    sopar->fakefact    = (iparm[IPARM_FILL_MATRIX] == API_YES) ? API_YES : API_NO;
    sopar->usenocsc    = 0;
    sopar->factotype   = iparm[IPARM_FACTORIZATION];
    sopar->symmetric   = iparm[IPARM_SYM];
    sopar->pastix_comm = pastix_comm;
    sopar->iparm       = iparm;
    sopar->dparm       = dparm;
    sopar->schur       = iparm[IPARM_SCHUR];
#ifdef PASTIX_DISTRIBUTED
    sopar->n           = pastix_data->ncol_int;
    MPI_Allreduce(&(pastix_data->ncol_int), &sopar->gN, 1, PASTIX_MPI_INT, MPI_SUM, pastix_comm);
#else
    sopar->n           = n;
    sopar->gN          = n;
#endif


    if (sopar->b != NULL)
        memFree_null(sopar->b);
    sopar->bindtab     = pastix_data->bindtab;

#ifdef FORCE_NOMPI
    iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
#elif (defined PASTIX_DYNSCHED)
    if ((iparm[IPARM_THREAD_COMM_MODE] != API_THREAD_COMM_ONE))
    {
        iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_COMM_ONE;
        if (procnum == 0 && iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
            errorPrintW("Dynsched require API_THREAD_COMM_ONE, forced.");
    }
#endif


    switch (iparm[IPARM_THREAD_COMM_MODE])
    {
    case API_THREAD_COMM_ONE:
    case API_THREAD_FUNNELED:
        iparm[IPARM_NB_THREAD_COMM] = 1;

    case API_THREAD_COMM_DEFINED:
        iparm[IPARM_NB_THREAD_COMM] = MAX(iparm[IPARM_NB_THREAD_COMM],1);
        break;

    case API_THREAD_COMM_NBPROC:
        iparm[IPARM_NB_THREAD_COMM] = iparm[IPARM_THREAD_NBR];
        break;

    default:
        iparm[IPARM_NB_THREAD_COMM] = 0;
    }
    sopar->type_comm  = iparm[IPARM_THREAD_COMM_MODE];
    sopar->nbthrdcomm = iparm[IPARM_NB_THREAD_COMM];

    switch (iparm[IPARM_END_TASK])
    {
    case API_TASK_NUMFACT: /* Only sopalin */

        print_debug(DBG_STEP,"FACTO SEULE\n");

        /* no facto if only raff */
        if (iparm[IPARM_ONLY_RAFF] == API_NO)
        {
            switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
                ge_z_sopalin_thread(solvmatr, sopar);
                break;
            case API_FACT_LLT:
                po_z_sopalin_thread(solvmatr, sopar);
                break;
            case API_FACT_LDLH:
                he_z_sopalin_thread(solvmatr, sopar);
                break;
            case API_FACT_LDLT:
            default:
                sy_z_sopalin_thread(solvmatr, sopar);
            }
        }
        break;

    case API_TASK_SOLVE: /* Sopalin and updown */

        print_debug(DBG_STEP,"FACTO + UPDO\n");
#ifndef PASTIX_FUNNELED
        if (THREAD_FUNNELED_ON)
        {
            if (procnum == 0)
                errorPrintW("API_THREAD_FUNNELED require -DPASTIX_FUNNELED,"
                            " force API_THREAD_MULTIPLE");
            sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
        }
#endif /* PASTIX_FUNNELED */
#ifndef PASTIX_THREAD_COMM
        if (THREAD_COMM_ON)
        {
            if (procnum == 0)
                errorPrintW("API_THREAD_COMM_* require -DPASTIX_THREAD_COMM,"
                            " force API_THREAD_MULTIPLE");
            sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
        }
#endif /* PASTIX_THREAD_COMM */

        /* Pour l'instant uniquement si on est en 1d */
        if (iparm[IPARM_DISTRIBUTION_LEVEL] == 0)
        {
            /* attention MRHS_ALLOC */
            if (pastix_data->malsmx)
            {
                memFree_null(solvmatr->updovct.sm2xtab);
                pastix_data->malsmx=0;
            }
            solvmatr->updovct.sm2xnbr = rhsnbr;
            MALLOC_INTERN(solvmatr->updovct.sm2xtab,
                          solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                          pastix_complex64_t);

            pastix_data->malsmx=1;

            z_buildUpdoVect(pastix_data,
#ifdef PASTIX_DISTRIBUTED
                          pastix_data->l2g_int,
                          pastix_data->b_int,
#else
                          NULL,
                          b,
#endif
                          pastix_comm);
        }

        if (iparm[IPARM_ONLY_RAFF] == API_NO)
        {
            /* Pour l'instant uniquement si on est en 1d */
            if (iparm[IPARM_DISTRIBUTION_LEVEL] == 0)
            {
                /* setting sopar->b for reffinement */
                /* Only 1 rhs is saved in sopar->b */
                if (sopar->b == NULL)
                {
                    MALLOC_INTERN(sopar->b, solvmatr->updovct.sm2xsze, pastix_complex64_t);
                }
                memcpy(sopar->b, solvmatr->updovct.sm2xtab,
                       solvmatr->updovct.sm2xsze*sizeof(pastix_complex64_t));
            }

            switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
                ge_z_sopalin_updo_thread(solvmatr, sopar);
                break;
            case API_FACT_LLT:
                po_z_sopalin_updo_thread(solvmatr, sopar);
                break;
            case API_FACT_LDLH:
                he_z_sopalin_updo_thread(solvmatr, sopar);
                break;
            case API_FACT_LDLT:
            default:
                sy_z_sopalin_updo_thread(solvmatr, sopar);
            }
        }

        iparm[IPARM_START_TASK]++;
        break;
    case API_TASK_REFINE: /* Sopalin, updown and raff */
    case API_TASK_CLEAN:

        print_debug(DBG_STEP,"FACTO + UPDO + RAFF (+ CLEAN)\n");
#ifndef PASTIX_UPDO_ISEND
        if (THREAD_COMM_ON)
        {
            if (procnum == 0)
                errorPrintW("THREAD_COMM require -DPASTIX_UPDO_ISEND,"
                            " force API_THREAD_MULTIPLE");
            sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
        }
#endif /* PASTIX_UPDO_ISEND */
#ifndef STORAGE
        if (THREAD_COMM_ON)
        {
            if (procnum == 0)
                errorPrintW("THREAD_COMM require -DSTORAGE,"
                            " force API_THREAD_MULTIPLE");
            sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
        }
#endif /* STORAGE */
        /* Pour l'instant uniquement si on est en 1d */
        if (iparm[IPARM_DISTRIBUTION_LEVEL] == 0)
        {
            /* attention MRHS_ALLOC */
            if (pastix_data->malsmx)
            {
                memFree_null(solvmatr->updovct.sm2xtab);
                pastix_data->malsmx=0;
            }
            if (rhsnbr > 1)
                errorPrintW("Reffinement works only with 1 rhs, please call them one after the other.");
            solvmatr->updovct.sm2xnbr = 1;
            MALLOC_INTERN(solvmatr->updovct.sm2xtab,
                          solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                          pastix_complex64_t);
            pastix_data->malsmx=1;
            z_buildUpdoVect(pastix_data,
#ifdef PASTIX_DISTRIBUTED
                          pastix_data->l2g_int,
                          pastix_data->b_int,
#else
                          NULL,
                          b,
#endif
                          pastix_comm);
            /* setting sopar->b for reffinement */
            if (sopar->b == NULL)
            {
                MALLOC_INTERN(sopar->b,
                              solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                              pastix_complex64_t);
            }
            memcpy(sopar->b, solvmatr->updovct.sm2xtab,
                   solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze*sizeof(pastix_complex64_t));
        }
        sopar->itermax     = iparm[IPARM_ITERMAX];
        sopar->epsilonraff = dparm[DPARM_EPSILON_REFINEMENT];
#ifdef OOC
        if (iparm[IPARM_GMRES_IM] != 1)
        {
            iparm[IPARM_GMRES_IM] = 1;
            if (procnum == 0)
                errorPrintW("IPARM_GMRES_IM force to 1 when using OOC");
        }
#endif
        sopar->gmresim = iparm[IPARM_GMRES_IM];
        switch (iparm[IPARM_REFINEMENT])
        {
        case API_RAF_GMRES:
            switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
                ge_z_sopalin_updo_gmres_thread(solvmatr, sopar);
                break;
            case API_FACT_LLT:
                po_z_sopalin_updo_gmres_thread(solvmatr, sopar);
                break;
            case API_FACT_LDLH:
                he_z_sopalin_updo_gmres_thread(solvmatr, sopar);
                break;
            case API_FACT_LDLT:
                sy_z_sopalin_updo_gmres_thread(solvmatr, sopar);
                break;
            default:
                errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
                return BADPARAMETER_ERR;
            }
            break;
        case API_RAF_PIVOT:
            switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
                ge_z_sopalin_updo_pivot_thread(solvmatr, sopar);
                break;
            case API_FACT_LLT:
            case API_FACT_LDLH:
            case API_FACT_LDLT:
                errorPrint("Refinement method and factorization type are incompatibles");
                return BADPARAMETER_ERR;
            default:
                errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
                return BADPARAMETER_ERR;
            }
            break;
        case API_RAF_GRAD:
            switch(iparm[IPARM_FACTORIZATION])
            {
            case API_FACT_LU:
                errorPrint("Refinement method and factorization type are incompatibles");
                return BADPARAMETER_ERR;
            case API_FACT_LLT:
                po_z_sopalin_updo_grad_thread(solvmatr, sopar);
                break;
            case API_FACT_LDLH:
                he_z_sopalin_updo_grad_thread(solvmatr, sopar);
                break;
            case API_FACT_LDLT:
                sy_z_sopalin_updo_grad_thread(solvmatr, sopar);
                break;
            default:
                errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
                return BADPARAMETER_ERR;
            }
            break;
            /* case API_RAF_BICGSTAB: */
            /*   switch(iparm[IPARM_FACTORIZATION]) */
            /*     { */
            /*     case API_FACT_LU: */
            /*       ge_z_sopalin_updo_bicgstab_thread(solvmatr, sopar); */
            /*       break; */
            /*     case API_FACT_LLT: */
            /*     case API_FACT_LDLH: */
            /*     case API_FACT_LDLT: */
            /*       errorPrint("Refinement method and factorization type are incompatibles"); */
            /*       return BADPARAMETER_ERR; */
            /*     default: */
            /*       errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]); */
            /*       return BADPARAMETER_ERR; */
            /*     } */
            /*   break; */
        default:
            errorPrint("Undefined refinement method : %ld", (long)iparm[IPARM_REFINEMENT]);
            return BADPARAMETER_ERR;
        }
        /* sopar->b was only needed for raff */
        memFree_null(sopar->b);
        iparm[IPARM_START_TASK]++;
        iparm[IPARM_START_TASK]++;
        iparm[IPARM_NBITER]         = sopar->itermax;
        dparm[DPARM_RELATIVE_ERROR] = sopar->rberror;
    }

    if ((iparm[IPARM_END_TASK] > API_TASK_NUMFACT) /* Not only sopalin */
        && (iparm[IPARM_DISTRIBUTION_LEVEL] == 0))
    {
        /* b <- solution */
        if (iparm[IPARM_GRAPHDIST] == API_NO)
        {
            if (iparm[IPARM_ONLY_RAFF] == API_NO)
            {
                z_CscRhsUpdown(&(solvmatr->updovct),
                             solvmatr,
                             b, n, ordemesh->peritab,
                             iparm[IPARM_DOF_NBR],
                             iparm[IPARM_RHS_MAKING],
                             pastix_comm);
            }
        }
#ifdef PASTIX_DISTRIBUTED
        else
        {
            z_CscdRhsUpdown(&(solvmatr->updovct),
                          solvmatr,
                          pastix_data->b_int,
                          pastix_data->ncol_int,
                          pastix_data->glob2loc,
                          ordemesh->peritab,
                          (int)iparm[IPARM_DOF_NBR],
                          pastix_comm);
        }
#endif
    }

    iparm[IPARM_STATIC_PIVOTING] = sopar->diagchange;

    /*
     * Memory statistics
     */
#ifdef MEMORY_USAGE
    {
        unsigned long smem[2], rmem[2];

        smem[0] = memAllocGetMax();
        smem[1] = memAllocGetCurrent();

        dparm[DPARM_MEM_MAX] = (double)smem[0];
        MPI_Reduce(smem, rmem, 2, MPI_LONG, MPI_MAX, 0, pastix_comm);

        if (procnum == 0)
        {
            dparm[DPARM_MEM_MAX] = (double)rmem[0];

            if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
            {
                fprintf(stdout, OUT_MAX_MEM_AF_SOP,  MEMORY_WRITE(rmem[0]), MEMORY_UNIT_WRITE(rmem[0]));
                fprintf(stdout, OUT_MEM_USED_AF_SOP, MEMORY_WRITE(rmem[1]), MEMORY_UNIT_WRITE(rmem[1]));
            }
        }
    }
#endif /* MEMORY_USAGE */


    spivot    = (long)  iparm[IPARM_STATIC_PIVOTING];
    sfacttime = (double)dparm[DPARM_FACT_TIME];
    MPI_Reduce(&spivot,   &rpivot,   1,MPI_LONG,  MPI_SUM,0,pastix_comm);
    MPI_Reduce(&sfacttime,&rfacttime,1,MPI_DOUBLE,MPI_MAX,0,pastix_comm);

    if (iparm[IPARM_ONLY_RAFF] == API_NO)
    {
        /*
         * Factorization Time
         */
        if ((procnum == 0) && (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT))
        {
            fprintf(stdout, OUT_STATIC_PIVOTING, rpivot);
            if (sopar->iparm[IPARM_INERTIA] != -1)
            {
                if (rpivot == 0)
                {
                    fprintf(stdout, OUT_INERTIA, (long)sopar->iparm[IPARM_INERTIA]);
                }
                else
                {
                    fprintf(stdout, OUT_INERTIA_PIVOT, (long)sopar->iparm[IPARM_INERTIA]);
                }
            }
            if (sopar->iparm[IPARM_ESP_NBTASKS] != -1)
                fprintf(stdout, OUT_ESP_NBTASKS,     (long)sopar->iparm[IPARM_ESP_NBTASKS]);
            fprintf(stdout, OUT_TIME_FACT,       rfacttime);
            fprintf(stdout, OUT_FLOPS_FACT,
                    PRINT_FLOPS((dparm[DPARM_FACT_FLOPS]/rfacttime)),
                    PRINT_FLOPS_UNIT((dparm[DPARM_FACT_FLOPS]/rfacttime)));
        }

        /*fprintf(stdout," Terms allocated during factorization %ld\n",(long)iparm[IPARM_ALLOCATED_TERMS]);*/

        /*
         * Solve Time
         */
        if (iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
        {
            ssolvtime = dparm[DPARM_SOLV_TIME];
            MPI_Reduce(&ssolvtime,&rsolvtime,1,MPI_DOUBLE,MPI_MAX,0,pastix_comm);

            if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
                print_onempi(OUT_TIME_SOLV, rsolvtime);
        }

        /*
         * Refinement Time
         */
        if (iparm[IPARM_END_TASK] > API_TASK_SOLVE)
        {
            srafftime = dparm[DPARM_RAFF_TIME];
            MPI_Reduce(&srafftime, &rrafftime, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);

            if ((procnum == 0) && (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT))
            {
                fprintf(stdout, OUT_RAFF_ITER_NORM,
                        (long)  iparm[IPARM_NBITER],
                        (double)dparm[DPARM_RELATIVE_ERROR]);
                fprintf(stdout, OUT_TIME_RAFF, rrafftime);
                if (iparm[IPARM_PRODUCE_STATS] == API_YES) {
                    if (dparm[DPARM_RELATIVE_ERROR] > 0)
                        print_onempi(OUT_PREC1, dparm[DPARM_RELATIVE_ERROR]);
                    if (dparm[DPARM_SCALED_RESIDUAL] > 0)
                        print_onempi(OUT_PREC2, dparm[DPARM_SCALED_RESIDUAL]);
                }

            }
        }
    }

    iparm[IPARM_START_TASK]++;
    return PASTIX_SUCCESS;
}



/*
 Function: z_pastix_task_updown

 Updown task.

 Parameters:
 pastix_data - PaStiX data structure.
 pastix_comm - PaStiX MPI communicator.
 n           - Matrix size.
 b           - Right hand side.
 loc2glob    - local to global column number.

 */
void z_pastix_task_updown(z_pastix_data_t *pastix_data,
                        MPI_Comm       pastix_comm,
                        pastix_int_t            n,
                        pastix_complex64_t         *b,
                        pastix_int_t            rhsnbr,
                        pastix_int_t           *loc2glob)
{
    pastix_int_t           * iparm    = pastix_data->iparm;
    double        * dparm    = pastix_data->dparm;
    z_SolverMatrix  * solvmatr = &(pastix_data->solvmatr);
    z_SopalinParam  * sopar    = &(pastix_data->sopar);
    Order         * ordemesh = pastix_data->ordemesh;
    pastix_int_t             procnum  = pastix_data->procnum;
    double          ssolvtime,rsolvtime;

    print_debug(DBG_STEP,"-> z_pastix_task_updown\n");
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        print_onempi("%s", OUT_STEP_SOLVE);

    if (sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
    {
        if (procnum == 0)
            errorPrintW("Updown step incompatible with 2D distribution");
        return;
    }

#ifndef PASTIX_UPDO_ISEND
    if (THREAD_COMM_ON)
    {
        if (procnum == 0)
            errorPrintW("THREAD_COMM require -DPASTIX_UPDO_ISEND,"
                        " force API_THREAD_MULTIPLE");
        sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
    }
#endif /* PASTIX_UPDO_ISEND */
#ifndef STORAGE
    if (THREAD_COMM_ON)
    {
        if (procnum == 0)
            errorPrintW("THREAD_COMM require -DSTORAGE,"
                        " force API_THREAD_MULTIPLE");
        sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
    }
#endif /* STORAGE */


    if ((iparm[IPARM_ONLY_RAFF] == API_YES) && (iparm[IPARM_END_TASK] > API_TASK_SOLVE))
    {
        errorPrintW("IPARM_ONLY_RAFF ignored, only possible if UPDO and RAFF are called in 2 steps");
        iparm[IPARM_ONLY_RAFF] = API_NO;
    }

    /* attention MRHS_ALLOC */
    if (pastix_data->malsmx)
    {
        memFree_null(solvmatr->updovct.sm2xtab);
        pastix_data->malsmx=0;
    }

    solvmatr->updovct.sm2xnbr = rhsnbr;
    MALLOC_INTERN(solvmatr->updovct.sm2xtab,
                  solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                  pastix_complex64_t);

    pastix_data->malsmx=1;

    z_buildUpdoVect(pastix_data,
                  loc2glob,
                  b,
                  pastix_comm);

    if (iparm[IPARM_ONLY_RAFF] == API_NO)
    {
        /* setting sopar->b for reffinement */
        /* Only 1 rhs is saved in sopar->b */
        if (sopar->b == NULL)
        {
            MALLOC_INTERN(sopar->b, solvmatr->updovct.sm2xsze, pastix_complex64_t);
        }
        memcpy(sopar->b, solvmatr->updovct.sm2xtab,
               solvmatr->updovct.sm2xsze*sizeof(pastix_complex64_t));

        sopar->iparm = iparm;
        sopar->dparm = dparm;

        switch(iparm[IPARM_FACTORIZATION])
        {
        case API_FACT_LU:
            ge_z_updo_thread(solvmatr, sopar);
            break;
        case API_FACT_LLT:
            po_z_updo_thread(solvmatr, sopar);
            break;
        case API_FACT_LDLH:
            he_z_updo_thread(solvmatr, sopar);
            break;
        case API_FACT_LDLT:
        default:
            sy_z_updo_thread(solvmatr, sopar);
        }

        /*
         if ((procnum == 0) && (iparm[IPARM_END_TASK] < API_TASK_REFINE))
         errorPrintW("Need a call to step 6 (refinement) to put the solution in the user vector.");
         */
        if ((iparm[IPARM_END_TASK] < API_TASK_REFINE ||
             iparm[IPARM_TRANSPOSE_SOLVE] == API_YES)
            && (iparm[IPARM_DISTRIBUTION_LEVEL] == 0))
        {
            /* b <- solution */
            if (iparm[IPARM_GRAPHDIST] == API_NO)
            {
                z_CscRhsUpdown(&(solvmatr->updovct),
                             solvmatr,
                             b, n, ordemesh->peritab,
                             iparm[IPARM_DOF_NBR],
                             iparm[IPARM_RHS_MAKING],
                             pastix_comm);
            }
#ifdef PASTIX_DISTRIBUTED
            else
            {
                z_CscdRhsUpdown(&(solvmatr->updovct),
                              solvmatr,
                              b, n,
                              pastix_data->glob2loc,
                              ordemesh->peritab,
                              iparm[IPARM_DOF_NBR], pastix_comm);
            }
#endif
        }

        if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) {
            ssolvtime = (double)dparm[DPARM_SOLV_TIME];
            MPI_Reduce(&ssolvtime,&rsolvtime,1,MPI_DOUBLE,MPI_MAX,0,pastix_comm);
            print_onempi(OUT_TIME_SOLV,rsolvtime);
            if (iparm[IPARM_PRODUCE_STATS] == API_YES) {
                if (dparm[DPARM_RELATIVE_ERROR] > 0)
                    print_onempi(OUT_PREC1, dparm[DPARM_RELATIVE_ERROR]);
                if (dparm[DPARM_SCALED_RESIDUAL] > 0)
                    print_onempi(OUT_PREC2, dparm[DPARM_SCALED_RESIDUAL]);
            }
        }
    }
    iparm[IPARM_START_TASK]++;
}

/*
 Function: z_pastix_task_raff

 Reffinement task

 Parameters:
 pastix_data - PaStiX data structure.
 pastix_comm - PaStiX MPI communicator.
 n           - Matrix size.
 b           - Right hand side.
 loc2glob    - local to global column number.
 */
void z_pastix_task_raff(z_pastix_data_t *pastix_data,
                      MPI_Comm       pastix_comm,
                      pastix_int_t            n,
                      pastix_complex64_t         *b,
                      pastix_int_t            rhsnbr,
                      pastix_int_t           *loc2glob)
{
    pastix_int_t           * iparm    = pastix_data->iparm;
    double        * dparm    = pastix_data->dparm;
    z_SopalinParam  * sopar    = &(pastix_data->sopar);
    z_SolverMatrix  * solvmatr = &(pastix_data->solvmatr);
    Order         * ordemesh = pastix_data->ordemesh;
    double          srafftime,rrafftime;
    pastix_int_t             procnum  = pastix_data->procnum;;
    pastix_complex64_t         * tmp;

    print_debug(DBG_STEP, "->z_pastix_task_raff\n");
#ifndef PASTIX_UPDO_ISEND
    if (THREAD_COMM_ON)
    {
        if (procnum == 0)
            errorPrintW("THREAD_COMM require -DPASTIX_UPDO_ISEND,"
                        " force API_THREAD_MULTIPLE");
        sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
    }
#endif /* PASTIX_UPDO_ISEND */
#ifndef STORAGE
    if (THREAD_COMM_ON)
    {
        if (procnum == 0)
            errorPrintW("THREAD_COMM require -DSTORAGE,"
                        " force API_THREAD_MULTIPLE");
        sopar->iparm[IPARM_THREAD_COMM_MODE] = API_THREAD_MULTIPLE;
    }
#endif /* STORAGE */

    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        print_onempi("%s", OUT_STEP_REFF);

    if (sopar->iparm[IPARM_DISTRIBUTION_LEVEL] != 0)
    {
        if (procnum == 0)
            errorPrintW("Refinment step incompatible with 2D distribution");
        return;
    }

    if (rhsnbr > 1)
    {
        errorPrintW("Reffinement works only with 1 rhs, please call them one after the other.");
        solvmatr->updovct.sm2xnbr = 1;
    }

    if (iparm[IPARM_ONLY_RAFF] == API_YES )
    {

        /* setting sopar->b for reffinement */
        if (sopar->b == NULL)
        {
            MALLOC_INTERN(sopar->b,
                          solvmatr->updovct.sm2xnbr*solvmatr->updovct.sm2xsze,
                          pastix_complex64_t);
        }

        tmp = solvmatr->updovct.sm2xtab;
        solvmatr->updovct.sm2xtab = sopar->b;

        z_buildUpdoVect(pastix_data,
                      loc2glob,
                      b,
                      pastix_comm);

        sopar->b = solvmatr->updovct.sm2xtab;
        solvmatr->updovct.sm2xtab = tmp;

    }

    sopar->itermax     = iparm[IPARM_ITERMAX];
    sopar->epsilonraff = dparm[DPARM_EPSILON_REFINEMENT];
#ifdef OOC
    if (iparm[IPARM_GMRES_IM] != 1)
    {
        iparm[IPARM_GMRES_IM] = 1;
        if (procnum == 0)
            errorPrintW("IPARM_GMRES_IM force to 1 when using OOC");
    }
#endif
    sopar->gmresim = iparm[IPARM_GMRES_IM];

    switch (iparm[IPARM_REFINEMENT])
    {
    case API_RAF_GMRES:
        switch(iparm[IPARM_FACTORIZATION])
        {
        case API_FACT_LU:
            ge_z_gmres_thread(solvmatr, sopar);
            break;
        case API_FACT_LLT:
            po_z_gmres_thread(solvmatr, sopar);
            break;
        case API_FACT_LDLH:
            he_z_gmres_thread(solvmatr, sopar);
            break;
        case API_FACT_LDLT:
            sy_z_gmres_thread(solvmatr, sopar);
            break;
        default:
            errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
            iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
            return;

        }
        break;
    case API_RAF_PIVOT:
        switch(iparm[IPARM_FACTORIZATION])
        {
        case API_FACT_LU:
            ge_z_pivot_thread(solvmatr, sopar);
            break;
        case API_FACT_LLT:
        case API_FACT_LDLH:
        case API_FACT_LDLT:
            errorPrint("Refinement method and factorization type are incompatibles");
            iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
            return;
        default:
            errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
            iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
            return;
        }
        break;
    case API_RAF_GRAD:
        switch(iparm[IPARM_FACTORIZATION])
        {
        case API_FACT_LU:
            errorPrint("Refinement method and factorization type are incompatibles");
            iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
            return;
        case API_FACT_LLT:
            po_z_grad_thread(solvmatr, sopar);
            break;
        case API_FACT_LDLH:
            he_z_grad_thread(solvmatr, sopar);
            break;
        case API_FACT_LDLT:
            sy_z_grad_thread(solvmatr, sopar);
            break;
        default:
            errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]);
            iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
            return;
        }
        break;
        /* case API_RAF_BICGSTAB: */
        /*   switch(iparm[IPARM_FACTORIZATION]) */
        /*     { */
        /*     case API_FACT_LU: */
        /*       ge_z_bicgstab_thread(solvmatr, sopar); */
        /*       break; */
        /*     case API_FACT_LLT: */
        /*     case API_FACT_LDLH: */
        /*     case API_FACT_LDLT: */
        /*       errorPrint("Refinement method and factorization type are incompatibles"); */
        /*       iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR; */
        /*       return; */
        /*     default: */
        /*       errorPrint("Undefined factorization type : %ld", (long)iparm[IPARM_FACTORIZATION]); */
        /*       iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR; */
        /*       return; */
        /*     } */
        /*   break; */
    default:
        errorPrint("Undefined refinement method : %ld", (long)iparm[IPARM_REFINEMENT]);
        iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
        return;
    }

    dparm[DPARM_RELATIVE_ERROR] = sopar->rberror;
    iparm[IPARM_NBITER]         = sopar->itermax;

    /* sopar->b was only needed for raff */
    memFree_null(sopar->b);

    /* b <- solution */
    if (iparm[IPARM_GRAPHDIST] == API_NO)
    {
        z_CscRhsUpdown(&(solvmatr->updovct),
                     solvmatr,
                     b, n, ordemesh->peritab,
                     iparm[IPARM_DOF_NBR],
                     iparm[IPARM_RHS_MAKING],
                     pastix_comm);
    }
#ifdef PASTIX_DISTRIBUTED
    else
    {
        z_CscdRhsUpdown(&(solvmatr->updovct),
                      solvmatr,
                      b, n,
                      pastix_data->glob2loc,
                      ordemesh->peritab,
                      iparm[IPARM_DOF_NBR], pastix_comm);
    }
#endif

    /* Fin du roerdering */

    srafftime = (double)dparm[DPARM_RAFF_TIME];
    MPI_Reduce(&srafftime,&rrafftime,1,MPI_DOUBLE,MPI_MAX,0,pastix_comm);

    if ((procnum == 0) && (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT))
    {
        fprintf(stdout, OUT_RAFF_ITER_NORM, (long)iparm[IPARM_NBITER], (double)dparm[DPARM_RELATIVE_ERROR]);
        if (iparm[IPARM_PRODUCE_STATS] == API_YES) {
            if (dparm[DPARM_RELATIVE_ERROR] > 0)
                print_onempi(OUT_PREC1, dparm[DPARM_RELATIVE_ERROR]);
            if (dparm[DPARM_SCALED_RESIDUAL] > 0)
                print_onempi(OUT_PREC2, dparm[DPARM_SCALED_RESIDUAL]);
        }

        fprintf(stdout, OUT_TIME_RAFF, rrafftime);
    }
    iparm[IPARM_START_TASK]++;

    return;
}

/*
 * Function: z_pastix_task_clean
 *
 * Cleaning task
 *
 * Parameters:
 *
 */
void z_pastix_task_clean(z_pastix_data_t **pastix_data,
                       MPI_Comm        pastix_comm)
{
    pastix_int_t             i;
    pastix_int_t           * iparm    = (*pastix_data)->iparm;
    z_SopalinParam  * sopar    = &((*pastix_data)->sopar);
    z_SolverMatrix  * solvmatr = &((*pastix_data)->solvmatr);
#ifdef PASTIX_DEBUG
    int             procnum  = (*pastix_data)->procnum;
    double        * dparm    = (*pastix_data)->dparm;
    FILE          * stream;
#endif
    (void)pastix_comm;

    print_debug(DBG_STEP, "->z_pastix_task_clean\n");
#ifdef PASTIX_DISTRIBUTED
    if ((*pastix_data)->mal_l2g_int == API_YES)
        memFree_null((*pastix_data)->l2g_int   );

    if ((*pastix_data)->malrhsd_int == API_YES)
    {
        if (iparm[IPARM_RHSD_CHECK] == API_YES)
        {
            memFree_null((*pastix_data)->b_int);
        }

        (*pastix_data)->malrhsd_int = API_NO;
    }
#endif

    if ((*pastix_data)->malcof)
    {
        if (iparm[IPARM_SCHUR] == API_YES && (*pastix_data)->schur_tab_set == API_YES)
        {
            z_SolverMatrix * datacode = &((*pastix_data)->solvmatr);
            pastix_int_t            cblk;

            if (SOLV_TASKNBR > 0)
            {
                cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
                if (SYMB_LCOLNUM(cblk) == (*pastix_data)->n2*(*pastix_data)->iparm[IPARM_DOF_NBR]-1)
                {
                    SOLV_COEFTAB(cblk) = NULL;
                }
            }
        }
    }

#ifdef OOC
    {

        char str[STR_SIZE];
        struct stat stFileInfo;

#  ifndef OOC_DIR
#    define OOC_DIR                "/tmp/z_pastix"
#  endif

        for (i = 0; i <solvmatr->cblknbr; i++)
        {
            sprintf(str,"%s/pastix_coef_%d/%d",OOC_DIR,(int)iparm[IPARM_OOC_ID],(int)i);
            if (-1 != stat(str,&stFileInfo) && (-1 == remove(str)))
            {
                perror("remove");
                EXIT(MOD_SOPALIN,UNKNOWN_ERR);
            }
        }
        sprintf(str,"%s/pastix_coef_%d",OOC_DIR,(int)iparm[IPARM_OOC_ID]);
        if (-1 == remove(str))
        {
            perror("remove");
            EXIT(MOD_SOPALIN,UNKNOWN_ERR);
        }

        if (iparm[IPARM_FACTORIZATION] == API_FACT_LU)
        {
            for (i = 0; i <solvmatr->cblknbr; i++)
            {
                sprintf(str,"%s/pastix_ucoef_%d/%d",OOC_DIR, (int)iparm[IPARM_OOC_ID], (int)i);
                if (-1 != stat(str,&stFileInfo) && (-1 == remove(str)))
                {
                    perror("remove");
                    EXIT(MOD_SOPALIN,UNKNOWN_ERR);
                }
            }
            sprintf(str,"%s/pastix_ucoef_%d",OOC_DIR, (int)iparm[IPARM_OOC_ID]);
            if (-1 == remove(str))
            {
                perror("remove");
                EXIT(MOD_SOPALIN,UNKNOWN_ERR);
            }
        }

#  ifdef OOC_FTGT
        for (i = 0; i <solvmatr->ftgtnbr; i++)
        {
            sprintf(str,"%s/pastix_ftgt_%d/%d",OOC_DIR,(int)iparm[IPARM_OOC_ID],(int)i);
            if (-1 != stat(str,&stFileInfo) && -1 == remove(str))
            {
                perror("remove");
                EXIT(MOD_SOPALIN,UNKNOWN_ERR);
            }
        }

        sprintf(str,"%s/pastix_ftgt_%d",OOC_DIR,(int)iparm[IPARM_OOC_ID]);
        if (-1 == remove(str))
        {
            perror("remove");
            EXIT(MOD_SOPALIN,UNKNOWN_ERR);
        }
#  endif /* OOC_FTGT */
    }
#endif /* OOC */

#ifdef PASTIX_DEBUG
    {
        char filename[256];
        sprintf(filename, "parm%ld.dump", (long) procnum);
        PASTIX_FOPEN(stream, filename, "w");
        api_dumparm(stream,iparm,dparm);
        fclose(stream);
    }
#endif

    if ((*pastix_data)->ordemesh != NULL)
    {
        orderExit((*pastix_data)->ordemesh);
        memFree_null( (*pastix_data)->ordemesh );
    }

    if ((*pastix_data)->malcsc)
    {
        z_CscExit(&((*pastix_data)->cscmtx));
        (*pastix_data)->malcsc=0;
    }

    if ((*pastix_data)->malsmx)
    {
        memFree_null(solvmatr->updovct.sm2xtab);
        (*pastix_data)->malsmx=0;
    }

    /* Pour l'instant uniquement si on est en 1d */
    if (iparm[IPARM_DISTRIBUTION_LEVEL] == 0)
    {
        if (solvmatr->updovct.cblktab)
            for (i=0; i<solvmatr->cblknbr; i++)
            {
                if (solvmatr->updovct.cblktab[i].browcblktab)
                    memFree_null(solvmatr->updovct.cblktab[i].browcblktab);

                if (solvmatr->updovct.cblktab[i].browproctab)
                    memFree_null(solvmatr->updovct.cblktab[i].browproctab);
            }

        memFree_null(solvmatr->updovct.lblk2gcblk);
        memFree_null(solvmatr->updovct.listblok);
        memFree_null(solvmatr->updovct.listcblk);
        memFree_null(solvmatr->updovct.gcblk2list);
        memFree_null(solvmatr->updovct.loc2glob);
        memFree_null(solvmatr->updovct.cblktab);
        memFree_null(solvmatr->updovct.listptr);
    }

    if (NULL != sopar->b)
        memFree_null(sopar->b);

    if ((*pastix_data)->malslv)
    {
        solverExit((d_SolverMatrix*)solvmatr);
#ifdef PASTIX_DYNSCHED
        Bubble_Free(solvmatr->btree);
        memFree_null(solvmatr->btree);
#endif
        (*pastix_data)->malslv=0;
    }

    if ((*pastix_data)->sopar.bindtab != NULL)
        memFree_null((*pastix_data)->sopar.bindtab);

    if ((*pastix_data)->schur_list != NULL)
        memFree_null((*pastix_data)->schur_list);
#ifdef PASTIX_DISTRIBUTED
    if ((*pastix_data)->glob2loc != NULL)
        memFree_null((*pastix_data)->glob2loc);
#endif
#ifdef WITH_SEM_BARRIER
    if ((*pastix_data)->intra_node_procnbr > 1)
    {
        if (sem_close((*pastix_data)->sem_barrier) < 0)
        {
            perror("sem_close");
        }
        if ((*pastix_data)->intra_node_procnum == 0)
        {
            char sem_name[256];
            sprintf(sem_name, "/pastix_%d", (*pastix_data)->pastix_id);
            if (sem_unlink(sem_name) < 0)
            {
                perror("sem_unlink");
            }
        }
    }
#endif

    if (*pastix_data != NULL)
        memFree_null(*pastix_data);

    z_FreeMpiType();
    z_FreeMpiSum();
    z_pastix_print_memory_usage(iparm,pastix_comm);

}

void z_pastix_unscale(z_pastix_data_t *pastix_data, pastix_int_t sym) {
#ifndef FORCE_MPI
    int size;
    MPI_Comm_size(pastix_data->pastix_comm, &size);
    if(size > 1) {
        errorPrint("pastix_task_unscale is not implemented in the distributed case (scaletabs must be distributed).");
        exit(1);
    }
#endif
    if(pastix_data->scaling == API_YES) {
        if(sym == API_YES)
            z_Matrix_Unscale_Sym(pastix_data, &pastix_data->solvmatr, pastix_data->scalerowtab, pastix_data->iscalerowtab);
        else
            z_Matrix_Unscale_Unsym(pastix_data, &pastix_data->solvmatr, pastix_data->scalerowtab, pastix_data->iscalerowtab, pastix_data->scalecoltab, pastix_data->iscalecoltab);
    }
}

#ifdef WITH_SEM_BARRIER
#  define SEM_BARRIER do {                                              \
        if ((*pastix_data)->intra_node_procnum == 0)                    \
        {                                                               \
            int si_iter;                                                \
            for (si_iter = 0;                                           \
                 si_iter < (*pastix_data)->intra_node_procnbr-1;        \
                 si_iter++)                                             \
            {                                                           \
                sem_post((*pastix_data)->sem_barrier);                  \
            }                                                           \
        }                                                               \
        else                                                            \
        {                                                               \
            sem_wait((*pastix_data)->sem_barrier);                      \
        }                                                               \
    } while(0)
#else
#  define SEM_BARRIER do {} while (0)
#endif
#define SYNC_IPARM do {                                         \
        if ((*pastix_data)->intra_node_procnbr > 1)             \
        {                                                       \
            SEM_BARRIER;                                        \
            MPI_Bcast(iparm,                                    \
                      IPARM_SIZE, PASTIX_MPI_INT,               \
                      0, (*pastix_data)->intra_node_comm);      \
        }                                                       \
    } while (0)

#define WAIT_AND_RETURN do {                                    \
        if ( *pastix_data != NULL )                             \
        {                                                       \
            SYNC_IPARM;                                         \
            if (iparm[IPARM_START_TASK] > API_TASK_SOLVE)       \
            {                                                   \
                MPI_Bcast(b, n, COMM_FLOAT, 0,                  \
                          (*pastix_data)->intra_node_comm );    \
            }                                                   \
        }                                                       \
        else {                                                  \
            MPI_Barrier((*pastix_data)->intra_node_comm);       \
        }                                                       \
        return;                                                 \
    } while (0)

/** ****************************************************************************************
 *
 *  Function: z_pastix
 *
 *  Computes one to all steps of the resolution of Ax=b linear system, using direct methods.
 *
 *  The matrix is given in CSC format.
 *
 *  Parameters:
 *  pastix_data - Data used for a step by step execution.
 *  pastix_comm - MPI communicator which compute the resolution.
 *  n           - Size of the system.
 *  colptr      - Tabular containing the start of each column in row and avals tabulars.
 *  row         - Tabular containing the row number for each element sorted by column.
 *  avals       - Tabular containing the values of each elements sorted by column.
 *  perm        - Permutation tabular for the renumerotation of the unknowns.
 *  invp        - Reverse permutation tabular for the renumerotation of the unknowns.
 *  b           - Right hand side vector(s).
 *  rhs         - Number of right hand side vector(s).
 *  iparm       - Integer parameters given to z_pastix.
 *  dparm       - Double parameters given to pâstix.
 *
 *  About: Example
 *
 *  from file <simple.c> :
 *
 *  > /\*******************************************\/
 *  > /\*    Check Matrix format                  *\/
 *  > /\*******************************************\/
 *  > /\*
 *  >  * Matrix needs :
 *  >  *    - to be in fortran numbering
 *  >  *    - to have only the lower triangular part in symmetric case
 *  >  *    - to have a graph with a symmetric structure in unsymmetric case
 *  >  *\/
 *  > mat_type = API_SYM_NO;
 *  > if (MTX_ISSYM(type)) mat_type = API_SYM_YES;
 *  > if (MTX_ISHER(type)) mat_type = API_SYM_HER;
 *  > z_pastix_checkMatrix(MPI_COMM_WORLD, verbosemode,
 *  >                    mat_type,  API_YES,
 *  >                    ncol, &colptr, &rows, &values, NULL);
 *  >
 *  > /\*******************************************\/
 *  > /\* Initialize parameters to default values *\/
 *  > /\*******************************************\/
 *  > iparm[IPARM_MODIFY_PARAMETER] = API_NO;
 *  > z_pastix(&pastix_data, MPI_COMM_WORLD,
 *  >        ncol, colptr, rows, values,
 *  >        perm, invp, rhs, 1, iparm, dparm);
 *  >
 *  > /\*******************************************\/
 *  > /\*       Customize some parameters         *\/
 *  > /\*******************************************\/
 *  > iparm[IPARM_THREAD_NBR] = nbthread;
 *  > iparm[IPARM_SYM] = mat_type;
 *  > switch (mat_type)
 *  >   {
 *  >     case API_SYM_YES:
 *  >       iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
 *  >       break;
 *  >     case API_SYM_HER:
 *  >       iparm[IPARM_FACTORIZATION] = API_FACT_LDLH;
 *  >       break;
 *  >     default:
 *  >       iparm[IPARM_FACTORIZATION] = API_FACT_LU;
 *  >    }
 *  > iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
 *  > iparm[IPARM_END_TASK]            = API_TASK_CLEAN;
 *  >
 *  > /\*******************************************\/
 *  > /\*           Save the rhs                  *\/
 *  > /\*    (it will be replaced by solution)    *\/
 *  > /\*******************************************\/
 *  > rhssaved = malloc(ncol*sizeof(pastix_complex64_t));
 *  > memcpy(rhssaved, rhs, ncol*sizeof(pastix_complex64_t));
 *  >
 *  > /\*******************************************\/
 *  > /\*           Call z_pastix                   *\/
 *  > /\*******************************************\/
 *  > perm = malloc(ncol*sizeof(pastix_int_t));
 *  > invp = malloc(ncol*sizeof(pastix_int_t));
 *  >
 *  > z_pastix(&pastix_data, MPI_COMM_WORLD,
 *  >  ncol, colptr, rows, values,
 *  >  perm, invp, rhs, 1, iparm, dparm);
 */
void z_pastix(z_pastix_data_t **pastix_data,
            MPI_Comm        pastix_comm,
            pastix_int_t             n,
            pastix_int_t            *colptr,
            pastix_int_t            *row,
            pastix_complex64_t          *avals,
            pastix_int_t            *perm,
            pastix_int_t            *invp,
            pastix_complex64_t          *b,
            pastix_int_t             rhs,
            pastix_int_t            *iparm,
            double         *dparm)
{
    int ret = PASTIX_SUCCESS;
#ifdef FIX_SCOTCH /* Pour le debug au cines */
    _SCOTCHintRandInit();
#endif

    iparm[IPARM_GRAPHDIST] = API_NO;
    if (iparm[IPARM_MODIFY_PARAMETER] == API_NO) /* init task */
    {
        /* init with default for iparm & dparm */
        z_pastix_initParam(iparm, dparm);
        iparm[IPARM_GRAPHDIST] = API_NO;
        return;
    }
    /*
     * Init : create pastix_data structure if it's first time
     */
    if (*pastix_data == NULL)
    {
        /* Need to be set to -1 in every cases */
        iparm[IPARM_OOC_ID] = -1;

        /* Allocation de la structure pastix_data qd on rentre dans
         z_pastix pour la première fois */
        z_pastix_task_init(pastix_data, pastix_comm, iparm, dparm);
        if ((*pastix_data)->intra_node_procnum == 0) {
            /* Affichage des options */
            z_pastix_welcome_print(*pastix_data, colptr, n);

            /* Matrix verification */
            if (iparm[IPARM_MATRIX_VERIFICATION] == API_YES)
                if ( PASTIX_SUCCESS != (ret = z_pastix_checkMatrix((*pastix_data)->inter_node_comm,
                                                                 iparm[IPARM_VERBOSE], iparm[IPARM_SYM],
                                                                 API_NO, n, &colptr, &row,
                                                                 (avals == NULL)?NULL:(&avals),
                                                                 NULL, iparm[IPARM_DOF_NBR])))
                {
                    errorPrint("The matrix is not in the correct format");
                    iparm[IPARM_ERROR_NUMBER] = ret;
                    return;
                }
        }
    }

    (*pastix_data)->n2 = n;

    if (PASTIX_SUCCESS != (ret = z_pastix_check_param(*pastix_data, rhs)))
    {
        iparm[IPARM_ERROR_NUMBER] = ret;
        return;
    }


    if ((*pastix_data)->intra_node_procnum == 0) {
        /* only master node do computations */
        if (iparm[IPARM_END_TASK]<API_TASK_ORDERING) {
            WAIT_AND_RETURN;
        }

        if (iparm[IPARM_ISOLATE_ZEROS] == API_YES &&
            iparm[IPARM_SCHUR] == API_YES)
        {
            errorPrint("Schur complement is incompatible with diagonal zeros isolation.");
            iparm[IPARM_ERROR_NUMBER] = BADPARAMETER_ERR;
            WAIT_AND_RETURN;
        }

        /*
         * Scotch : Ordering
         */
        if (iparm[IPARM_START_TASK] == API_TASK_ORDERING) /* scotch task */
        {
            if (iparm[IPARM_ISOLATE_ZEROS] == API_YES)
            {
                pastix_int_t itercol;
                pastix_int_t iterrow;
                pastix_int_t iterzeros = 0;
                int found;

                (*pastix_data)->zeros_n = 0;
                for (itercol = 0; itercol < n; itercol++)
                {
                    found = API_NO;
                    for (iterrow = colptr[itercol]-1; iterrow <  colptr[itercol+1]-1; iterrow++)
                    {
                        if (row[iterrow]-1 == itercol)
                        {
                            if (ABS_FLOAT(avals[iterrow]) < dparm[DPARM_EPSILON_REFINEMENT])
                            {
                                (*pastix_data)->zeros_n++;
                            }
                            found = API_YES;
                            break;
                        }
                    }
                    if (found == API_NO)
                    {
                        (*pastix_data)->zeros_n++;
                    }
                }
                MALLOC_INTERN((*pastix_data)->zeros_list,
                              (*pastix_data)->zeros_n,
                              pastix_int_t);
                for (itercol = 0; itercol < n; itercol++)
                {
                    found = API_NO;
                    for (iterrow = colptr[itercol]-1; iterrow <  colptr[itercol+1]-1; iterrow++)
                    {
                        if (row[iterrow]-1 == itercol)
                        {
                            if (ABS_FLOAT(avals[iterrow]) < dparm[DPARM_EPSILON_REFINEMENT])
                            {
                                (*pastix_data)->zeros_list[iterzeros] = itercol+1;
                                iterzeros++;
                            }
                            found = API_YES;
                            break;
                        }
                    }
                    if (found == API_NO)
                    {
                        (*pastix_data)->zeros_list[iterzeros] = itercol+1;
                        iterzeros++;
                    }
                }
            }

            // TODO: (*pastix_data)->inter_node_comm,
            if (PASTIX_SUCCESS !=
                (ret = pastix_task_order( (d_pastix_data_t*)*pastix_data,
                                          n, colptr, row, NULL, perm, invp)))
            {
                iparm[IPARM_ERROR_NUMBER] = ret;
                WAIT_AND_RETURN;
            }
            if (iparm[IPARM_ISOLATE_ZEROS] == API_YES)
            {
                memFree_null((*pastix_data)->zeros_list);
            }
        }

        if (iparm[IPARM_END_TASK]<API_TASK_SYMBFACT) {
            WAIT_AND_RETURN;
        }

        /*
         * Fax : Facto symbolic
         */
        if (iparm[IPARM_START_TASK] == API_TASK_SYMBFACT) /* Fax task */
            pastix_task_symbfact( (d_pastix_data_t*)*pastix_data, perm, invp );

        if (iparm[IPARM_END_TASK] < API_TASK_ANALYSE) {
            WAIT_AND_RETURN;
        }

        /*
         * Blend : Scheduling
         */
        if (iparm[IPARM_START_TASK] == API_TASK_ANALYSE) /* Blend task */
            z_pastix_task_blend(*pastix_data);

        if (iparm[IPARM_END_TASK]<API_TASK_NUMFACT) {
            WAIT_AND_RETURN;
        }

#if defined(PROFILE) && defined(MARCEL)
        profile_activate(FUT_ENABLE, MARCEL_PROF_MASK, 0);
        marcel_printf("DEBUT profil marcel\n");
#endif

        /*
         * Sopalin : Factorisation
         */
        if (iparm[IPARM_START_TASK] == API_TASK_NUMFACT) /* Sopalin task */
        {
            ret = z_pastix_task_sopalin(*pastix_data,
                                      (*pastix_data)->inter_node_comm, n,
                                      colptr, row, avals, b, rhs, NULL);


            MPI_Bcast(&ret, 1, MPI_INT, 0, (*pastix_data)->inter_node_comm);
            if (PASTIX_SUCCESS != ret) {
                iparm[IPARM_ERROR_NUMBER] = ret;
                WAIT_AND_RETURN;
            }
        }
        if (iparm[IPARM_END_TASK]<iparm[IPARM_START_TASK]) {
            WAIT_AND_RETURN;
        }

        /*
         * Updo : solve
         */
        if (iparm[IPARM_START_TASK] == API_TASK_SOLVE) /* Updown task */
        {
            /* For thread comm */
            (*pastix_data)->sopar.stopthrd = API_YES;
            z_pastix_task_updown(*pastix_data, (*pastix_data)->inter_node_comm,
                               n, b, rhs, NULL);
            /* For thread comm */
            (*pastix_data)->sopar.stopthrd = API_NO;
        }
        if (iparm[IPARM_END_TASK]<API_TASK_REFINE) {
            WAIT_AND_RETURN;
        }

        /*
         * Raff
         */
        if (iparm[IPARM_START_TASK] == API_TASK_REFINE) /* Refinement task */
        {
            z_pastix_task_raff(*pastix_data, (*pastix_data)->inter_node_comm,
                             n, b, rhs, NULL);
        }
        if (iparm[IPARM_END_TASK]<API_TASK_CLEAN) {
            WAIT_AND_RETURN;
        }

        /*
         * Clean
         */
    } /* (*pastix_data)->intra_node_procnum == 0 */

    SYNC_IPARM;
    if (iparm[IPARM_END_TASK]<API_TASK_CLEAN)
        return;

    if (iparm[IPARM_START_TASK] == API_TASK_CLEAN)
        z_pastix_task_clean(pastix_data, pastix_comm);

#if defined(PROFILE) && defined(MARCEL)
    profile_stop();
    marcel_printf("FIN profil marcel\n");
#endif
}

#define REDISTRIBUTE_RHS                                                \
    {                                                                   \
                                                                        \
        if (b != NULL && rhsHasBeenRedistributed == API_NO)             \
        {                                                               \
            rhsHasBeenRedistributed = API_YES;                          \
            if (rhs_need_redispatch == API_YES)                         \
            {                                                           \
                if ((*pastix_data)->procnum == 0 &&                     \
                    iparm[IPARM_VERBOSE] >= API_VERBOSE_YES)            \
                    fprintf(stdout,OUT_REDIS_RHS);                      \
                                                                        \
                /* Distribute the user RHS into                         \
                 intern distribution */                                 \
                if (b_int == NULL)                                      \
                {                                                       \
                    MALLOC_INTERN(b_int, ncol_int*rhs, pastix_complex64_t); \
                    (*pastix_data)->b_int = b_int;                      \
                }                                                       \
                z_redispatch_rhs(n,                                       \
                               b,                                       \
                               rhs,                                     \
                               loc2glob,                                \
                               ncol_int,                                \
                               b_int,                                   \
                               l2g_int,                                 \
                               (*pastix_data)->procnbr,                 \
                               (*pastix_data)->procnum,                 \
                               pastix_comm,                             \
                               iparm[IPARM_DOF_NBR]);                   \
                                                                        \
            }                                                           \
            else                                                        \
            {                                                           \
                b_int = b;                                              \
            }                                                           \
        }                                                               \
    }
#define REDISTRIBUTE_SOL                                        \
    {                                                           \
        if (rhs_need_redispatch == API_YES)                     \
        {                                                       \
            if ((*pastix_data)->procnum == 0 &&                 \
                iparm[IPARM_VERBOSE] >= API_VERBOSE_YES)        \
                fprintf(stdout,OUT_REDIS_SOL);                  \
                                                                \
            z_redispatch_rhs(ncol_int,                            \
                           b_int,                               \
                           rhs,                                 \
                           l2g_int,                             \
                           n,                                   \
                           b,                                   \
                           loc2glob,                            \
                           (*pastix_data)->procnbr,             \
                           (*pastix_data)->procnum,             \
                           pastix_comm,                         \
                           iparm[IPARM_DOF_NBR]);               \
        }                                                       \
    }

/** Function: z_dpastix

 Computes one to all steps of the resolution of
 Ax=b linear system, using direct methods.
 Here the matrix is given distributed.

 The matrix is given in CSCD format.

 Parameters:
 pastix_data - Data used for a step by step execution.
 pastix_comm - MPI communicator which compute the resolution.
 n           - Size of the system.
 colptr      - Tabular containing the start of each column in row and avals tabulars.
 row         - Tabular containing the row number for each element sorted by column.
 avals       - Tabular containing the values of each elements sorted by column.
 loc2glob    - Global column number of the local columns.
 perm        - Permutation tabular for the renumerotation of the unknowns.
 invp        - Reverse permutation tabular for the renumerotation of the unknowns.
 b           - Right hand side vector(s).
 rhs         - Number of right hand side vector(s).
 iparm       - Integer parameters given to z_pastix.
 dparm       - Double parameters given to pâstix.
 */
void z_dpastix(z_pastix_data_t    **pastix_data,
               MPI_Comm             pastix_comm,
               pastix_int_t         n,
               pastix_int_t        *colptr,
               pastix_int_t        *row,
               pastix_complex64_t  *avals,
               pastix_int_t        *loc2glob,
               pastix_int_t        *perm,
               pastix_int_t        *invp,
               pastix_complex64_t  *b,
               pastix_int_t         rhs,
               pastix_int_t        *iparm,
               double              *dparm)
{
#ifdef PASTIX_DISTRIBUTED
    pastix_int_t    ncol_int                = 0;
    pastix_int_t   *l2g_int                 = NULL;
    pastix_complex64_t *b_int                   = NULL;
    int    ret                     = PASTIX_SUCCESS;
    int    ret_rcv                 = PASTIX_SUCCESS;
    pastix_int_t    gN                      = -1;
    int    rhs_need_redispatch     = API_NO;
    int    rhsHasBeenRedistributed = API_NO;
    int    mayNeedReturnSol        = API_NO;
#  ifdef FIX_SCOTCH /* Pour le debug au cines */
    _SCOTCHintRandInit();
#  endif


    if (iparm[IPARM_MODIFY_PARAMETER] == API_NO) /* init task */
    {
        /* init with default for iparm & dparm */
        z_pastix_initParam(iparm, dparm);
        return;
    }

    /* Si pastix_data est nul, c'est qu'on rentre dans
     la fonction pour la première fois */
    if (*pastix_data == NULL)
    {
        iparm[IPARM_OOC_ID]         = -1;
        /* initialisation et allocation de pastix_data */
        z_pastix_task_init(pastix_data,pastix_comm,iparm,dparm);

        /* Affichage des options */
        z_pastix_welcome_print(*pastix_data, colptr, n);

        /* multiple RHS see MRHS_ALLOC */
        if ( ((*pastix_data)->procnum == 0)  && (rhs!=1) )
            errorPrintW("multiple right-hand-side not tested...");

        /* Matrix verification */
        if (iparm[IPARM_MATRIX_VERIFICATION] == API_YES)
            if (PASTIX_SUCCESS != (ret = z_pastix_checkMatrix(pastix_comm, iparm[IPARM_VERBOSE],
                                                            iparm[IPARM_SYM], API_NO,
                                                            n, &colptr, &row , (avals == NULL)?NULL:(&avals),
                                                            ( (iparm[IPARM_GRAPHDIST] == API_NO)?
                                                              NULL:(&loc2glob) ),
                                                            iparm[IPARM_DOF_NBR])))
            {
                errorPrint("The matrix is not in the correct format");
                iparm[IPARM_ERROR_NUMBER] = ret;
                return;
            }

    }
    (*pastix_data)->n2 = n;


    if (PASTIX_SUCCESS != (ret = z_pastix_check_param(*pastix_data, rhs)))
    {
        iparm[IPARM_ERROR_NUMBER] = ret;
        return;
    }

    /*
     * WARNING: tant que tt n'est pas parallele IPARM_FREE_CSCUSER
     * est incompatible avec z_dpastix car il est appliqué sur la
     * csc interne globable et non sur la csc utilisateur
     */
    if (iparm[IPARM_FREE_CSCUSER] == API_CSC_FREE)
    {
        iparm[IPARM_FREE_CSCUSER] = API_CSC_PRESERVE;
        if ((*pastix_data)->procnum == 0)
            errorPrintW("Free CSC user is forbiden with z_dpastix for now");
    }

    /*
     * Scotch : Ordering
     */
    if (iparm[IPARM_START_TASK] == API_TASK_ORDERING) /* scotch task */
    {
        if ((*pastix_data)->intra_node_procnum == 0) {
            if (PASTIX_SUCCESS != (ret = pastix_task_order(*pastix_data,
                                                           n, colptr, row, loc2glob,
                                                           perm, invp)))
            {
                errorPrint("Error in ordering task\n");
                iparm[IPARM_ERROR_NUMBER] = ret;
            }
        }
        SYNC_IPARM;
        if (iparm[IPARM_ERROR_NUMBER] != PASTIX_SUCCESS)
            return;
    }
    if (iparm[IPARM_END_TASK]<API_TASK_SYMBFACT)
        return;

    /*
     * Fax : Facto symbolic
     */
    if (iparm[IPARM_START_TASK] == API_TASK_SYMBFACT) /* Fax task */
    {
        pastix_task_symbfact( *pastix_data,
                              perm, invp );
        SYNC_IPARM;
    }

    if (iparm[IPARM_END_TASK]<API_TASK_ANALYSE)
        return;

    if ((*pastix_data)->intra_node_procnum == 0)
    {
        if (iparm[IPARM_START_TASK] == API_TASK_ANALYSE) /* Blend task */
        {
            z_pastix_task_blend(*pastix_data);
        }
    }
    SYNC_IPARM;
    if (iparm[IPARM_END_TASK]<API_TASK_NUMFACT)
        return;

#  if defined(PROFILE) && defined(MARCEL)
    profile_activate(FUT_ENABLE, MARCEL_PROF_MASK, 0);
    marcel_printf("DEBUT profil marcel\n");
#  endif

    if (iparm[IPARM_START_TASK] == API_TASK_NUMFACT) /* Sopalin task */
    {
        mayNeedReturnSol = API_YES;
        if ((*pastix_data)->intra_node_procnum == 0)
        {
            ret = z_pastix_task_sopalin(*pastix_data,
                                      (*pastix_data)->inter_node_comm,
                                      n,
                                      colptr,
                                      row,
                                      avals,
                                      b,
                                      rhs,
                                      loc2glob);

            MPI_Allreduce(&ret, &ret_rcv, 1, MPI_INT, MPI_MAX,
                          (*pastix_data)->inter_node_comm);
            if (PASTIX_SUCCESS != ret_rcv)
            {
                errorPrint("Error in numeric factorisation task\n");
                iparm[IPARM_ERROR_NUMBER] = ret_rcv;
            }
        }
        else
        {
            /* Remplissage de la csc interne */
            if ((*pastix_data)->cscInternFilled == API_NO)
            {
                z_pastix_fake_fillin_csc(*pastix_data, pastix_comm, n,
                                       colptr, row, avals, b, rhs, loc2glob);
            }
#  ifdef PASTIX_DISTRIBUTED
            else
            {
                (*pastix_data)->l2g_int = loc2glob;
                (*pastix_data)->b_int   = b;
            }
#  endif

        }
        SYNC_IPARM;

        if (iparm[IPARM_ERROR_NUMBER] != PASTIX_SUCCESS)
            return;

        if (iparm[IPARM_START_TASK] > API_TASK_SOLVE)
        {
            ncol_int = (*pastix_data)->ncol_int;
            b_int    = (*pastix_data)->b_int;
            l2g_int  = (*pastix_data)->l2g_int;
            if ((*pastix_data)->malrhsd_int)
                rhs_need_redispatch = API_YES;
            REDISTRIBUTE_SOL;
        }

        if (iparm[IPARM_END_TASK] < API_TASK_CLEAN)
            return;

    }



    if ((*pastix_data)->intra_node_procnum == 0)
    {
        ncol_int = (*pastix_data)->ncol_int;
        b_int    = (*pastix_data)->b_int;
        l2g_int  = (*pastix_data)->l2g_int;

        /* User can change CSCD after blend */
        if ((iparm[IPARM_GRAPHDIST] == API_YES) &&
            ((*pastix_data)->glob2loc == NULL))
        {
            z_cscd_build_g2l(ncol_int,
                           l2g_int,
                           (*pastix_data)->inter_node_comm,
                           &gN,
                           &((*pastix_data)->glob2loc));
        }


        /* Updown task */

        /* If user has not specified that he is
         absolutly certain that is CSCd is
         correctly distributed */
        if ( ( iparm[IPARM_START_TASK] == API_TASK_SOLVE ||
               iparm[IPARM_START_TASK] == API_TASK_REFINE ) &&
             iparm[IPARM_CSCD_CORRECT] == API_NO )
        {
            pastix_int_t my_n;
            pastix_int_t * my_l2g = NULL;
            int OK, OK_RECV;
            pastix_int_t iter;
            /* Test que la cscd utilisateur correspond a la cscd z_pastix */
            my_n = z_pastix_getLocalNodeNbr(pastix_data);

            OK = 0;
            if (my_n != n)
            {
                OK = 1;
            }
            else
            {
                if ((*pastix_data)->l2g_int) {
                    my_l2g = (*pastix_data)->l2g_int;
                } else {
                    MALLOC_INTERN(my_l2g, my_n, pastix_int_t);
                    z_pastix_getLocalNodeLst(pastix_data, my_l2g);
                }
                for (iter = 0; iter < my_n; iter++)
                {
                    if (my_l2g[iter] != loc2glob[iter])
                    {
                        OK = 1;
                        break;
                    }
                }
                if (!(*pastix_data)->l2g_int) {
                    memFree_null(my_l2g);
                }
            }
            MPI_Allreduce(&OK, &OK_RECV, 1, MPI_INT, MPI_SUM, pastix_comm);
            if (OK_RECV != 0)
                rhs_need_redispatch = API_YES;
        }

        if (iparm[IPARM_START_TASK] == API_TASK_SOLVE)
        {
            mayNeedReturnSol = API_YES;
            REDISTRIBUTE_RHS;
            /* For thread comm */
            (*pastix_data)->sopar.stopthrd = API_YES;
            z_pastix_task_updown(*pastix_data,
                               (*pastix_data)->inter_node_comm,
                               ncol_int, b_int, rhs, l2g_int);
            /* For thread comm */
            (*pastix_data)->sopar.stopthrd = API_NO;
        }
    }
    else
    {
        /* If user has not specified that he is
         absolutly certain that is CSCd is
         correctly distributed */
        if ( ( iparm[IPARM_START_TASK] == API_TASK_SOLVE ||
               iparm[IPARM_START_TASK] == API_TASK_REFINE ) &&
             iparm[IPARM_CSCD_CORRECT] == API_NO )
        {
            int OK, OK_RECV;
            OK = 0;

            MPI_Allreduce(&OK, &OK_RECV, 1, MPI_INT, MPI_SUM, pastix_comm);
            if (OK_RECV != 0)
                rhs_need_redispatch = API_YES;
        }

        if (iparm[IPARM_START_TASK] == API_TASK_SOLVE)
        {
            mayNeedReturnSol = API_YES;
            REDISTRIBUTE_RHS;
        }
    }

    if ( mayNeedReturnSol &&
         iparm[IPARM_END_TASK]<API_TASK_REFINE)
    {
        REDISTRIBUTE_SOL;
        return;
    }

    if (iparm[IPARM_START_TASK] == API_TASK_REFINE) /* Refinement task */
    {
        /* If it wasn't done just after solve */
        REDISTRIBUTE_RHS;
        if ((*pastix_data)->intra_node_procnum == 0)
        {
            z_pastix_task_raff(*pastix_data, (*pastix_data)->inter_node_comm,
                             ncol_int, b_int, rhs, l2g_int);
        }
        REDISTRIBUTE_SOL;
    }

    SYNC_IPARM;
    if (iparm[IPARM_END_TASK]<API_TASK_CLEAN)
        return;

    if (iparm[IPARM_START_TASK] == API_TASK_CLEAN)
        z_pastix_task_clean(pastix_data, pastix_comm);

#  if defined(PROFILE) && defined(MARCEL)
    profile_stop();
    marcel_printf("FIN profil marcel\n");
#  endif
#else
    (void)pastix_data; (void)pastix_comm; (void)n; (void)colptr; (void)row;
    (void)avals; (void)loc2glob; (void)perm; (void)invp; (void)b; (void)rhs;
    (void)iparm; (void)dparm;
    errorPrint("To use z_dpastix please compile with -DPASTIX_DISTRIBUTED");
    iparm[IPARM_ERROR_NUMBER] = BAD_DEFINE_ERR;
#endif /* PASTIX_DISTRIBUTED */
}



/*
 Function: z_pastix_bindThreads

 Set bindtab in pastix_data, it gives for each thread the CPU to bind in to.
 bindtab follows this organisation :

 bindtab[threadnum] = cpu to set thread threadnum.

 Parameters:
 pastix_data - Structure de donnée pour l'utilisation step by step
 thrdnbr     - Nombre de threads / Taille du tableau
 bindtab     - Tableau de correspondance entre chaque thread et coeur de la machine
 */

void z_pastix_bindThreads ( z_pastix_data_t *pastix_data, pastix_int_t thrdnbr, pastix_int_t *bindtab)
{
    int i;

    if ( pastix_data == NULL )
    {
        errorPrint("Pastix_data need to be initialized before to try to set bindtab.");
        EXIT(MOD_SOPALIN, BADPARAMETER_ERR);
    }

    /* Copy association tab between threads and cores */
    MALLOC_INTERN(pastix_data->sopar.bindtab, thrdnbr, int);
    for (i = 0; i < thrdnbr; i++)
    {
        pastix_data->sopar.bindtab[i] = bindtab[i];
    }
    /* Check values in bindtab */
    {
        int nbproc;
#ifdef MARCEL
        nbproc = marcel_nbvps();
#else
        nbproc = sysconf(_SC_NPROCESSORS_ONLN);
#endif

        for (i=0; i< thrdnbr; i++)
            if (!(pastix_data->sopar.bindtab[i] < nbproc))
            {
                errorPrint("Try to bind thread on an unavailable core.");
                EXIT(MOD_SOPALIN, BADPARAMETER_ERR);
            }

    }

    pastix_data->bindtab = pastix_data->sopar.bindtab;
    return;
}
/*
 * Function: z_pastix_checkMatrix_int
 *
 * Check the matrix :
 * - Renumbers in Fortran numerotation (base 1) if needed (base 0)
 * - Check that the matrix contains no doubles,  with flagcor == API_YES,
 *   correct it.
 * - Can scale the matrix if compiled with -DMC64 -DSCALING (untested)
 * - Checks the symetry of the graph in non symmetric mode.
 *   With non distributed matrices, with flagcor == API_YES,
 *   correct the matrix.
 * - sort the CSC.
 *
 * Parameters:
 *   pastix_comm - PaStiX MPI communicator
 *   verb        - Level of prints (API_VERBOSE_[NOT|NO|YES])
 *   flagsym     - Indicate if the given matrix is symetric
 *                 (API_SYM_YES or API_SYM_NO)
 *   flagcor     - Indicate if we permit the function to reallocate the matrix.
 *   n           - Number of local columns.
 *   colptr      - First element of each row in *row* and *avals*.
 *   row         - Row of each element of the matrix.
 *   avals       - Value of each element of the matrix.
 *   loc2glob    - Global column number of local columns
 *                 (NULL if not distributed).
 *   dof         - Number of degrees of freedom.
 *   flagalloc   - indicate if allocation on CSC uses internal malloc.
 */
pastix_int_t z_pastix_checkMatrix_int(MPI_Comm pastix_comm,
                                    pastix_int_t      verb,
                                    pastix_int_t      flagsym,
                                    pastix_int_t      flagcor,
                                    pastix_int_t      n,
                                    pastix_int_t    **colptr,
                                    pastix_int_t    **row,
                                    pastix_complex64_t  **avals,
                                    pastix_int_t    **loc2glob,
                                    pastix_int_t      dof,
                                    pastix_int_t      flagalloc)
{
    int  procnum;
    int  ret;
    int  OK;
    int  OK_RECV;
    pastix_int_t  old;
    pastix_int_t  i;
    pastix_int_t  l2g_sum_n[2];
    pastix_int_t  l2g_sum_n_reduced[2];

    MPI_Comm_rank(pastix_comm, &procnum);

    if (verb > API_VERBOSE_NOT)
        print_onempi("%s","Check : Numbering");

    if (!(*colptr)[0])
    {
        /* fortran-style numbering */
        if (verb > API_VERBOSE_NOT)
            print_onempi("%s", "\n\tC numbering to Fortran Numbering\tOK\n");
        z_csc_Cnum2Fnum(*row,*colptr,n);
        if (loc2glob != NULL)
            for (i = 0; i <  n; i++)
                (*loc2glob)[i]++;
    }
    else
    {
        if (verb > API_VERBOSE_NOT)
            print_onempi("%s", "\t\tOK\n");
    }

    if (loc2glob != NULL)
    {
        pastix_int_t l2g_OK;
        pastix_int_t l2g_OK_rcv;
        l2g_sum_n[0] = 0;
        l2g_sum_n[1] = n;
        l2g_OK       = 0;

        for (i = 0; i < n ; i++)
        {
            l2g_sum_n[0] += (*loc2glob)[i];
            if ( i > 0 && (*loc2glob)[i] <= (*loc2glob)[i-1] )
            {
                l2g_OK = 1;
            }
        }
        MPI_Allreduce(&l2g_OK,  &l2g_OK_rcv,        1,
                      PASTIX_MPI_INT, MPI_SUM, pastix_comm);
        if (l2g_OK_rcv > 0)
        {
            print_onempi("%s", "Local column must be ordered increasingly\n");
            return BADPARAMETER_ERR;
        }
        MPI_Allreduce(l2g_sum_n, l2g_sum_n_reduced, 2, PASTIX_MPI_INT, MPI_SUM, pastix_comm);
        /* All column have been distributed */
        if (2*l2g_sum_n_reduced[0] != (l2g_sum_n_reduced[1]*(l2g_sum_n_reduced[1]+1)))
        {
            print_onempi("%s", "All column must be destributed once and only once\n");
            return BADPARAMETER_ERR;
        }
    }


    /* sorting */
    if (verb > API_VERBOSE_NOT)
        print_onempi("%s", "Check : Sort CSC");

    if (avals != NULL)
        z_csc_sort(n,*colptr,*row,*avals, dof);
    else
        z_csc_sort(n,*colptr,*row,NULL, 0);
    if (verb > API_VERBOSE_NOT)
        print_onempi("%s","\t\tOK\n");


    if (verb > API_VERBOSE_NOT)
        print_onempi("%s", "Check : Duplicates");

    old = (*colptr)[n]-1;
    /* Preserve sorting */
    ret = z_csc_check_doubles(n,
                              *colptr,
                              row,
                              avals,
                              dof,
                              flagcor,
                              flagalloc);

    if (loc2glob != NULL)
    {
        if (ret == API_YES) {OK = 0;}
        MPI_Allreduce(&OK, &OK_RECV, 1, MPI_INT, MPI_SUM, pastix_comm);
        if (OK_RECV>0) {ret = API_NO;}
    }

    if (ret == API_YES)
    {
        if (verb > API_VERBOSE_NOT)
            print_onempi("%s", "\t\tOK\n");
    }
    else
    {
        if (verb > API_VERBOSE_NOT)
            print_onempi("%s", "\t\t\tKO\n");
        return MATRIX_ERR;
    }

    if (verb > API_VERBOSE_NOT)
    {
        if (old != ((*colptr)[n] - 1))
        {
            if (loc2glob != NULL)
            {
                fprintf(stdout,
                        "\n\t%ld double terms merged on proc %ld \n",
                        (long)(old - (*colptr)[n]-1),(long)procnum);
            }
            else
            {
                print_onempi("\t%ld double terms merged\n",
                             (long)(old - ((*colptr)[n]-1)));
            }
        }
    }
    {
        pastix_int_t cnt_lower     = 0;
        pastix_int_t cnt_upper     = 0;
        pastix_int_t cnt_diag      = 0;
        pastix_int_t cnt_num_zeros = 0;
        pastix_int_t globn         = n;
        pastix_int_t itercol;
        pastix_int_t iterrow;

        for (itercol = 0; itercol < n; itercol ++)
        {
            for ( iterrow = (*colptr)[itercol] - 1;
                  iterrow < (*colptr)[itercol+1] - 1;
                  iterrow++)
            {
                if ((*row)[iterrow]-1 >
                    ((loc2glob == NULL)?itercol:((*loc2glob)[itercol]-1)))
                {
                    cnt_lower++;
                }
                else
                {
                    if ((*row)[iterrow]-1 <
                        ((loc2glob == NULL)?itercol:((*loc2glob)[itercol]-1)))
                    {
                        cnt_upper++;
                    }
                    else
                    {
                        cnt_diag++;
                        if (avals != NULL && (*avals)[iterrow*dof*dof] == 0.)
                            cnt_num_zeros++;
                    }
                }
            }
        }


        if (loc2glob != NULL)
        {
            pastix_int_t send_data[5];
            pastix_int_t recv_data[5];

            send_data[0] = cnt_lower;
            send_data[1] = cnt_upper;
            send_data[2] = cnt_diag;
            send_data[3] = cnt_num_zeros;
            send_data[4] = globn;
            MPI_Allreduce(send_data, recv_data, 5, PASTIX_MPI_INT, MPI_SUM, pastix_comm);
            cnt_lower      = recv_data[0];
            cnt_upper      = recv_data[1];
            cnt_diag       = recv_data[2];
            cnt_num_zeros  = recv_data[3];
            globn          = recv_data[4];
        }

        if (cnt_diag != globn)
            if (verb > API_VERBOSE_NOT)
                errorPrintW("%d/%d structural zeros found on the diagonal.",
                            globn-cnt_diag, globn);

        if (cnt_num_zeros != 0)
            if (verb > API_VERBOSE_NOT)
                errorPrintW("%d numerical zeros found on the diagonal.", cnt_num_zeros);

        if (cnt_upper == cnt_lower &&
            cnt_lower != 0 &&
            ( flagsym == API_SYM_YES  ||
              flagsym == API_SYM_HER ) &&
            flagcor == API_YES)
        {
            pastix_int_t index = 0;
            pastix_int_t lastindex = 0;
            pastix_int_t   * tmprows;
            pastix_complex64_t * tmpvals;
            errorPrintW("Upper and lower part given on a symmetric matrix, dropping upper");
            for (itercol = 0; itercol < n; itercol++)
            {

                for (iterrow = (*colptr)[itercol]-1;
                     iterrow < (*colptr)[itercol+1]-1;
                     iterrow++)
                {
                    if ((*row)[iterrow]-1 >=
                        (loc2glob == NULL)?itercol:((*loc2glob)[itercol]-1))
                    {
                        (*row)[index] = (*row)[iterrow];
                        if (avals != NULL)
                            (*avals)[index] = (*avals)[iterrow];
                        index++;
                    }
                }
                (*colptr)[itercol] = lastindex+1;
                lastindex = index;
            }
            (*colptr)[n] = lastindex+1;
            MALLOC_EXTERN(tmprows, lastindex, pastix_int_t);
            memcpy(tmprows, (*row),   lastindex*sizeof(pastix_int_t));
            free((*row));
            (*row) = tmprows;
            if (avals != NULL)
            {
                MALLOC_EXTERN(tmpvals, lastindex, pastix_complex64_t);
                memcpy(tmpvals, (*avals), lastindex*sizeof(pastix_complex64_t));
                free((*avals));
                (*avals) = tmpvals;
            }
        }
        else
        {
            if ( ( flagsym == API_SYM_YES ||
                   flagsym == API_SYM_HER ) &&
                 cnt_lower != 0 && cnt_upper != 0 )
            {
                errorPrint("Only lower or upper part should be given (lower %d upper %d diag %d)",
                           cnt_lower, cnt_upper, cnt_diag);
                return MATRIX_ERR;
            }
        }
    }


    /* Pre-conditionnement mc64 */
#ifdef SCALING
#  ifdef MC64
    if (sizeof(int) != sizeof(pastix_int_t))
    {
        errorPrint("MC64 only works with classical integers\n");
        return INTEGER_TYPE_ERR;
    }

    errorPrintW("NOT TESTED");
    for (i = 0; i > (*colptr)[n]-1; i++)
        if ((*row)[i] == 0)
            errorPrint("Et Merde\n");

    if ((flagcor == API_YES) && (iparm[IPARM_MC64] == 1))
    {
        pastix_int_t    job;
        pastix_int_t    m      = n;
        pastix_int_t    ne     = (*colptr)[n]-1;
        pastix_int_t    num;
        pastix_int_t   *p, *ip;
        pastix_int_t    liw    = 3*m+2*n+ne;
        pastix_int_t   *iw;
        pastix_int_t    ldw    = n+3*m+ne;
        pastix_int_t    nicntl = 10;
        pastix_int_t    ncntl  = 10;
        pastix_int_t   *icntl;
        pastix_int_t    info;
        pastix_complex64_t *dw;
        pastix_complex64_t *cntl;

        print_onempi("%s", "Preconditioning...\n");

        MALLOC_INTERN(p,     m,      pastix_int_t);
        MALLOC_INTERN(ip,    m,      pastix_int_t);
        MALLOC_INTERN(iw,    liw,    pastix_int_t);
        MALLOC_INTERN(dw,    ldw,    pastix_complex64_t);
        MALLOC_INTERN(icntl, nicntl, pastix_int_t);
        MALLOC_INTERN(cntl,  ncntl,  pastix_complex64_t);

        for (i=0;i<m;i++)
        {
            p[i]  = i+1;
            ip[i] = i+1;
        }
        for (i=0;i<ldw;i++)
            dw[i] = 0.0;
        for (i=0;i<liw;i++)
            iw[i] = 0;
        for (i=0;i<nicntl;i++)
            icntl[i] = 0;
        for (i=0;i<ncntl;i++)
            cntl[i] = 0.0;

        /* compute scaling and unsymmetric column permutation */
        print_onempi("%s", "compute scaling and unsymmetric column permutation...\n");

        FORTRAN_CALL(mc64id)(icntl,cntl);

        cntl[1] = DBL_MAX;
        job     = 6;
        num     = n;

        printf("adresse1 job=%p m=%p n=%p ne=%p ia=%p\n",&job,&m,&n,&ne,*colptr);
        FORTRAN_CALL(mc64ad)(&job,&m,&n,&ne,*colptr,*row,*avals,&num,p,&liw,iw,
                             &ldw,dw,icntl,cntl,&info);
        fprintf(stdout,"return info=%ld (num=%ld n=%ld)\n",info,num,n);
        printf("adresse1 job=%p m=%p n=%p ne=%p ia=%p\n",&job,&m,&n,&ne,*colptr);

        if (num<0 || info<0)
        {
            errorPrint("Error in MC64AD !!!");
            memFree_null(p);
            memFree_null(ip);
            memFree_null(iw);
            memFree_null(dw);
            memFree_null(icntl);
            memFree_null(cntl);
        }
        else
        {
            /* scaling */
            for (i=0;i<m+n;i++)
                dw[i]=exp(dw[i]); /* a_ij := aij * exp(u_i + u_j) */

            print_onempi("%s", "scaling rows...\n");
            z_csc_rowScale(n,*colptr,*row,*avals,dw);

            print_onempi("%s", "z_scaling.columns...\n");
            z_csc_colScale(n,*colptr,*row,*avals,dw+m);

            /* apply unsymmetric column permutation */
            for (i=0;i<m;i++)
                ip[p[i]-1]=i+1; /* inverse permutation */

            print_onempi("s%", "apply unsymmetric column permutation...\n");
            z_csc_colPerm(n,*colptr,*row,*avals,ip);

            memFree_null(p);
            memFree_null(ip);
            memFree_null(iw);
            memFree_null(dw);
            memFree_null(icntl);
            memFree_null(cntl);
        }
    }
#  endif /* MC64 */
#endif /* SCALING */

    /* Symmetrisation du graphe des matrices non-symmétriques */
    if (flagsym == API_SYM_NO)
    {
        if (verb > API_VERBOSE_NOT)
            print_onempi("%s", "Check : Graph symmetry");

        old = (*colptr)[n]-1;

        /* Version distribuée */
        if ((loc2glob != NULL))
        {
            /* Preserve sorting */
            if (EXIT_SUCCESS == z_cscd_checksym(n, *colptr, row, avals, *loc2glob,
                                                flagcor,
                                                flagalloc, dof,  pastix_comm))
            {
                if (verb > API_VERBOSE_NOT)
                {
                    print_onempi("%s", "\t\tOK\n");
                }
                if (verb > API_VERBOSE_NOT)
                {
                    if (old != ((*colptr)[n] - 1))
                    {
                        fprintf(stdout, "\tAdd %ld null terms on proc %ld\n",
                                (long)((*colptr)[n]-1-old), (long)procnum);
                    }
                }

            }
            else
            {
                if (verb > API_VERBOSE_NOT)
                {
                    print_onempi("%s", "\t\tKO\n");
                }
                return MATRIX_ERR;
            }
        }
        /* Version non distribuée */
        else
        {
            /* Preserve sorting */
            if (EXIT_SUCCESS == z_csc_checksym(n, *colptr, row, avals, flagcor, flagalloc, dof))
            {
                if (verb > API_VERBOSE_NOT)
                {
                    print_onempi("%s", "\t\tOK\n");
                }
            }
            else
            {
                if (verb > API_VERBOSE_NOT)
                {
                    print_onempi("%s", "\t\tKO\n");
                }
                return MATRIX_ERR;
            }
            if (verb > API_VERBOSE_NOT)
            {
                if (old != ((*colptr)[n] - 1))
                {
                    print_onempi("\tAdd %ld null terms\t\n",(long)((*colptr)[n]-1-old));
                }
            }


        }
    }

    return PASTIX_SUCCESS;
}

/*
 Function: z_pastix_getLocalUnknownNbr

 Return the node number in the new distribution computed by blend.
 Needs blend to be runned with pastix_data before.

 Parameters:
 pastix_data - Data used for a step by step execution.

 Returns:
 Number of local nodes/columns in new distribution.
 */
pastix_int_t z_pastix_getLocalUnknownNbr(z_pastix_data_t ** pastix_data)
{
    z_SolverMatrix  * solvmatr = &((*pastix_data)->solvmatr);
    pastix_int_t index;
    pastix_int_t nodenbr;

    nodenbr = 0;
    if ((*pastix_data)->intra_node_procnum == 0)
        for (index=0; index<solvmatr->cblknbr; index++)
        {
            nodenbr += solvmatr->cblktab[index].lcolnum-solvmatr->cblktab[index].fcolnum+1;
        }
    return nodenbr;
}

/*
 Function: z_pastix_getLocalNodeNbr

 Return the node number in the new distribution computed by blend.
 Needs blend to be runned with pastix_data before.

 Parameters:
 pastix_data - Data used for a step by step execution.

 Returns:
 Number of local nodes/columns in new distribution.
 */
pastix_int_t z_pastix_getLocalNodeNbr(z_pastix_data_t ** pastix_data)
{
    z_SolverMatrix  * solvmatr = &((*pastix_data)->solvmatr);
    pastix_int_t index;
    pastix_int_t nodenbr;

    nodenbr = 0;
    if ((*pastix_data)->intra_node_procnum == 0)
        for (index=0; index<solvmatr->cblknbr; index++)
        {
            nodenbr += solvmatr->cblktab[index].lcolnum-solvmatr->cblktab[index].fcolnum+1;
        }
    nodenbr = nodenbr/(*pastix_data)->iparm[IPARM_DOF_NBR];
    return nodenbr;
}
/* qsort int comparison function */
int
cmpint(const void *p1, const void *p2)
{
    const pastix_int_t *a = (const pastix_int_t *)p1;
    const pastix_int_t *b = (const pastix_int_t *)p2;

    return (int) *a - *b;
}
/*
 Function: z_pastix_getLocalUnknownLst

 Fill in unknowns with the list of local nodes/clumns.
 Needs nodelst to be allocated with nodenbr*sizeof(pastix_int_t),
 where nodenbr has been computed by <z_pastix_getLocalUnknownNbr>.

 Parameters:
 pastix_data - Data used for a step by step execution.
 nodelst     - An array where to write the list of local nodes/columns.
 */
pastix_int_t z_pastix_getLocalUnknownLst(z_pastix_data_t **pastix_data,
                                       pastix_int_t            *nodelst)
{

    z_SolverMatrix  * solvmatr = &((*pastix_data)->solvmatr);
    Order         * ordemesh = (*pastix_data)->ordemesh;
    pastix_int_t index;
    pastix_int_t index2;
    pastix_int_t index3;
    int dof = (*pastix_data)->iparm[IPARM_DOF_NBR];

    index3 = 0;
    if ((*pastix_data)->intra_node_procnum == 0)
        for (index=0; index<solvmatr->cblknbr; index++)
        {
            for (index2 = solvmatr->cblktab[index].fcolnum;
                 index2 < solvmatr->cblktab[index].lcolnum + 1;
                 index2++)
                nodelst[index3++] = dof*ordemesh->peritab[(index2-index2%dof)/dof]+1+index2%dof;
        }
    qsort(nodelst, index3, sizeof(pastix_int_t), cmpint);

    return PASTIX_SUCCESS;
}


/*
 Function: z_pastix_getLocalNodeLst

 Fill in nodelst with the list of local nodes/clumns.
 Needs nodelst to be allocated with nodenbr*sizeof(pastix_int_t),
 where nodenbr has been computed by <z_pastix_getLocalNodeNbr>.

 Parameters:
 pastix_data - Data used for a step by step execution.
 nodelst     - An array where to write the list of local nodes/columns.
 */
pastix_int_t z_pastix_getLocalNodeLst(z_pastix_data_t **pastix_data,
                                    pastix_int_t            *nodelst)
{

    z_SolverMatrix  * solvmatr = &((*pastix_data)->solvmatr);
    Order         * ordemesh = (*pastix_data)->ordemesh;
    pastix_int_t index;
    pastix_int_t index2;
    pastix_int_t index3;
    int dof = (*pastix_data)->iparm[IPARM_DOF_NBR];

    index3 = 0;
    if ((*pastix_data)->intra_node_procnum == 0)
        for (index=0; index<solvmatr->cblknbr; index++)
        {
            for (index2 = solvmatr->cblktab[index].fcolnum;
                 index2 < solvmatr->cblktab[index].lcolnum + 1;
                 index2+=dof)
                nodelst[index3++] = ordemesh->peritab[index2/dof]+1;
        }
    qsort(nodelst, index3, sizeof(pastix_int_t), cmpint);

    return PASTIX_SUCCESS;
}

/*
 Function: z_pastix_setSchurUnknownList

 Set the list of unknowns to isolate at the end
 of the matrix via permutations.

 Parameters:
 pastix_data - Data used for a step by step execution.
 n           - Number of unknowns.
 list        - List of unknowns.
 */
pastix_int_t z_pastix_setSchurUnknownList(z_pastix_data_t * pastix_data,
                                        pastix_int_t  n,
                                        pastix_int_t *list)
{
    if (pastix_data == NULL)
        return STEP_ORDER_ERR;

    if (n == 0 || list == NULL)
        return BADPARAMETER_ERR;

    pastix_data->schur_n = n;
    MALLOC_INTERN(pastix_data->schur_list, n, pastix_int_t);
    memcpy(pastix_data->schur_list, list, n*sizeof(pastix_int_t));
    return PASTIX_SUCCESS;
}
/*
 Function: z_pastix_getSchurLocalNodeNbr

 Compute the number of nodes in the local part of the Schur.

 Parameters:
 pastix_data - Common data structure for PaStiX calls.
 nodeNbr     - (out) Number of nodes in schur (local).

 Returns:
 PASTIX_SUCCESS      - For the moment

 TODO: Error management.
 */
pastix_int_t z_pastix_getSchurLocalNodeNbr(z_pastix_data_t * pastix_data, pastix_int_t * nodeNbr)
{
    z_SolverMatrix * datacode = &(pastix_data->solvmatr);
    int            owner = API_NO;
    pastix_int_t            cblk;

    if (SOLV_TASKNBR > 0)
    {
        cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
        if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
        {
            owner = API_YES;
        }

    }

    if (owner == API_YES)
    {
        *nodeNbr = SYMB_LCOLNUM(cblk) - SYMB_FCOLNUM(cblk) + 1;
    }
    else
    {
        *nodeNbr = 0;
    }
    return PASTIX_SUCCESS;
}

/*
 Function: z_pastix_getSchurLocalUnkownNbr

 Compute the number of unknowns in the local part of the Schur.

 Parameters:
 pastix_data - Common data structure for PaStiX calls.
 unknownNbr  - (out) Number of unknowns in schur (local).

 Returns:
 PASTIX_SUCCESS      - For the moment

 TODO: Error management.
 */
pastix_int_t z_pastix_getSchurLocalUnkownNbr(z_pastix_data_t * pastix_data, pastix_int_t * unknownNbr)
{
    z_SolverMatrix * datacode = &(pastix_data->solvmatr);
    int            owner = API_NO;
    pastix_int_t            cblk;

    if (SOLV_TASKNBR > 0)
    {
        cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
        if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
        {
            owner = API_YES;
        }

    }

    if (owner == API_YES)
    {
        fprintf(stdout, "SYMB_LCOLNUM(cblk) %ld\n", (long)SYMB_LCOLNUM(cblk));
        fprintf(stdout, "SYMB_FCOLNUM(cblk) %ld\n", (long)SYMB_FCOLNUM(cblk));
        *unknownNbr = (SYMB_LCOLNUM(cblk) - SYMB_FCOLNUM(cblk) + 1)*pastix_data->iparm[IPARM_DOF_NBR];
    }
    else
    {
        *unknownNbr = 0;
    }
    return PASTIX_SUCCESS;
}

/*
 Function: z_pastix_getSchurLocalNodeList

 Compute the list of nodes in the local part of the Schur.

 Parameters:
 pastix_data - Common data structure for PaStiX calls.
 nodes     - (out) Nodes in schur (local).

 Returns:
 PASTIX_SUCCESS      - For the moment

 TODO: Error management.
 */
pastix_int_t z_pastix_getSchurLocalNodeList(z_pastix_data_t * pastix_data, pastix_int_t * nodes)
{
    z_SolverMatrix * datacode = NULL;
    Order        * ordemesh = NULL;
    int            owner = API_NO;
    pastix_int_t            cblk;
    pastix_int_t            intern_index;
    pastix_int_t            dof;
    pastix_int_t            intern_index_dof;
    datacode = &(pastix_data->solvmatr);
    ordemesh = pastix_data->ordemesh;

    if (SOLV_TASKNBR > 0)
    {
        cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
        if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
        {
            owner = API_YES;
        }

    }

    if (owner == API_YES)
    {
        pastix_int_t iter;
        for (iter = 0; iter < SYMB_LCOLNUM(cblk) - SYMB_FCOLNUM(cblk) + 1; iter+=pastix_data->iparm[IPARM_DOF_NBR])
        {
            intern_index = iter + SYMB_FCOLNUM(cblk);
            dof = intern_index % pastix_data->iparm[IPARM_DOF_NBR];
            intern_index_dof = (intern_index - dof) / pastix_data->iparm[IPARM_DOF_NBR];
            nodes[iter/pastix_data->iparm[IPARM_DOF_NBR]] = ordemesh->peritab[intern_index_dof];
        }
    }

    return PASTIX_SUCCESS;
}


/*
 Function: z_pastix_getSchurLocalUnkownList

 Compute the list of unknowns in the local part of the Schur.

 Parameters:
 pastix_data - Common data structure for PaStiX calls.
 unknowns    - (out) Unknowns in schur (local).

 Returns:
 PASTIX_SUCCESS      - For the moment

 TODO: Error management.
 */
pastix_int_t pastix_getSchurLocalUnknownList(z_pastix_data_t * pastix_data, pastix_int_t * unknowns)
{
    z_SolverMatrix * datacode = NULL;
    Order        * ordemesh = NULL;
    int            owner = API_NO;
    pastix_int_t            cblk;
    pastix_int_t            intern_index;
    pastix_int_t            dof;
    pastix_int_t            intern_index_dof;
    datacode = &(pastix_data->solvmatr);
    ordemesh = pastix_data->ordemesh;

    if (SOLV_TASKNBR > 0)
    {
        cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
        if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
        {
            owner = API_YES;
        }

    }

    if (owner == API_YES)
    {
        pastix_int_t iter;
        for (iter = 0; iter < SYMB_LCOLNUM(cblk) - SYMB_FCOLNUM(cblk) + 1; iter++)
        {
            intern_index = iter + SYMB_FCOLNUM(cblk);
            dof = intern_index % pastix_data->iparm[IPARM_DOF_NBR];
            intern_index_dof = (intern_index - dof) / pastix_data->iparm[IPARM_DOF_NBR];
            unknowns[iter] = (ordemesh->peritab[intern_index_dof]*pastix_data->iparm[IPARM_DOF_NBR])+dof;
        }
    }

    return PASTIX_SUCCESS;
}


/*
 Function: z_pastix_getSchurLocalUnkownList

 Give user memory area to store schur in PaStiX.

 Parameters:
 pastix_data - Common data structure for PaStiX calls.
 array       - Memory area to store the schur.

 Returns:
 PASTIX_SUCCESS      - For the moment

 TODO: Error management.
 */
pastix_int_t z_pastix_setSchurArray(z_pastix_data_t * pastix_data, pastix_complex64_t * array)
{
    pastix_data->schur_tab = array;
    pastix_data->schur_tab_set = API_YES;
    return PASTIX_SUCCESS;
}
/*
 Function: z_pastix_getSchur

 Get the Schur complement from PaStiX.

 Schur complement is a dense block in a
 column scheme.

 Parameters:
 pastix_data - Data used for a step by step execution.
 schur - Array to fill-in with Schur complement.

 */
pastix_int_t z_pastix_getSchur(z_pastix_data_t * pastix_data,
                             pastix_complex64_t * schur)
{
    z_SolverMatrix * datacode = &(pastix_data->solvmatr);
    int            owner = API_NO;
    pastix_int_t            send[2];
    pastix_int_t            recv[2];
    pastix_int_t            cblk;

    if (SOLV_TASKNBR > 0)
    {
        cblk = TASK_CBLKNUM(SOLV_TASKNBR-1);
        if (SYMB_LCOLNUM(cblk) == pastix_data->n2*pastix_data->iparm[IPARM_DOF_NBR]-1)
        {
            owner = API_YES;
        }

    }

    if (owner == API_YES)
    {
        pastix_int_t coefnbr  = SOLV_STRIDE(cblk) * (SYMB_LCOLNUM(cblk) - SYMB_FCOLNUM(cblk) + 1);
        memcpy(schur, SOLV_COEFTAB(cblk), coefnbr*sizeof(pastix_complex64_t));
        send[0] = coefnbr;
        send[1] = SOLV_PROCNUM;

        MPI_Allreduce(&send, &recv, 2, PASTIX_MPI_INT, MPI_SUM, pastix_data->pastix_comm);
    }
    else
    {
        send[0] = 0;
        send[1] = 0;

        MPI_Allreduce(&send, &recv, 2, PASTIX_MPI_INT, MPI_SUM, pastix_data->pastix_comm);
    }
    MPI_Bcast(schur, recv[0], COMM_FLOAT, recv[1], pastix_data->pastix_comm);
    return PASTIX_SUCCESS;
}
/*
 * Function: z_pastix_checkMatrix
 *
 * Check the matrix :
 * - Renumbers in Fortran numerotation (base 1) if needed (base 0)
 * - Check that the matrix contains no doubles,  with flagcor == API_YES,
 *   correct it.
 * - Can scale the matrix if compiled with -DMC64 -DSCALING (untested)
 * - Checks the symetry of the graph in non symmetric mode.
 *   With non distributed matrices, with flagcor == API_YES,
 *   correct the matrix.
 * - sort the CSC.
 *
 * Parameters:
 *   pastix_comm - PaStiX MPI communicator
 *   verb        - Level of prints (API_VERBOSE_[NOT|NO|YES])
 *   flagsym     - Indicate if the given matrix is symetric
 *                 (API_SYM_YES or API_SYM_NO)
 *   flagcor     - Indicate if we permit the function to reallocate the matrix.
 *   n           - Number of local columns.
 *   colptr      - First element of each row in *row* and *avals*.
 *   row         - Row of each element of the matrix.
 *   avals       - Value of each element of the matrix.
 *   loc2glob    - Global column number of local columns
 *                 (NULL if not distributed).
 *   dof         - Number of degrees of freedom.
 */
pastix_int_t z_pastix_checkMatrix(MPI_Comm pastix_comm,
                                pastix_int_t      verb,
                                pastix_int_t      flagsym,
                                pastix_int_t      flagcor,
                                pastix_int_t      n,
                                pastix_int_t    **colptr,
                                pastix_int_t    **row,
                                pastix_complex64_t  **avals,
                                pastix_int_t    **loc2glob,
                                pastix_int_t      dof)
{
    return z_pastix_checkMatrix_int(pastix_comm,
                                  verb,
                                  flagsym,
                                  flagcor,
                                  n,
                                  colptr,
                                  row,
                                  avals,
                                  loc2glob,
                                  dof,
                                  API_NO);
}

#define z_pastix_getMemoryUsage PASTIX_EXTERN_F(z_pastix_getMemoryUsage)
unsigned long z_pastix_getMemoryUsage() {
#ifdef MEMORY_USAGE
    return memAllocGetCurrent();
#else
    return -1;
#endif
}

#define pastix_getMaxMemoryUsage PASTIX_EXTERN_F(pastix_getMaxMemoryUsage)
unsigned long pastix_getMaxMemoryUsage() {
#ifdef MEMORY_USAGE
    return memAllocGetMax();
#else
    return -1;
#endif
}
