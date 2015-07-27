/*
  header: api.h

  Header file containing constants used in PaStiX and provided to users.

  Authors:
    Mathieu Faverge - faverge@labri.fr
    Xavier   Lacoste - lacoste@labri.fr
    Pierre Ramet     - ramet@labri.fr

 */

#ifndef _PASTIX_API_H_
#define _PASTIX_API_H_

#define PASTIX_MASK_ISTRUE(var, mask) (var == (var | mask))

/* Acces au tableau iparm*/
/*
   enum: IPARM_ACCESS

   Integer parameters tabular accessors

   IPARM_MODIFY_PARAMETER      - Indicate if parameters have been set by user             Default: API_YES             IN
   IPARM_START_TASK            - Indicate the first step to execute (see PaStiX steps)    Default: API_TASK_ORDERING   IN
   IPARM_END_TASK              - Indicate the last step to execute (see PaStiX steps)     Default: API_TASK_CLEAN      IN
   IPARM_VERBOSE               - Verbose mode (see Verbose modes)                         Default: API_VERBOSE_NO      IN
   IPARM_DOF_NBR               - Degree of freedom per node                               Default: 1                   IN
   IPARM_ITERMAX               - Maximum iteration number for refinement                  Default: 250                 IN
   IPARM_MATRIX_VERIFICATION   - Check the input matrix                                   Default: API_NO              IN
   IPARM_MC64                  - MC64 operation <pastix.h> IGNORE                         Default: 0                   IN
   IPARM_ONLY_RAFF             - Refinement only                                          Default: API_NO              IN
   IPARM_CSCD_CORRECT          - Indicate if the cscd has been redistributed after blend  Default: API_NO              IN
   IPARM_NBITER                - Number of iterations performed in refinement             Default: -                   OUT
   IPARM_TRACEFMT              - Trace format (see Trace modes)                           Default: API_TRACE_PAJE      IN

   IPARM_ORDERING              - Choose ordering                                          Default: API_ORDER_SCOTCH    IN
   IPARM_ORDERING_DEFAULT      - Use default ordering parameters with \scotch{} or \metis{} Default: API_YES           IN

   IPARM_SCOTCH_SWITCH_LEVEL   - Ordering switch level    (see \scotch{} User's Guide)    Default: 120                 IN
   IPARM_SCOTCH_CMIN           - Ordering cmin parameter  (see \scotch{} User's Guide)    Default: 0                   IN
   IPARM_SCOTCH_CMAX           - Ordering cmax parameter  (see \scotch{} User's Guide)    Default: 100000              IN
   IPARM_SCOTCH_FRAT           - Ordering frat parameter  (see \scotch{} User's Guide)    Default: 8                   IN

   IPARM_METIS_CTYPE           - \metis{} parameters (see \metis{} Manual), used only if IPARM_ORDERING_DEFAULT set to API_NO   Default: METIS_CTYPE_SHEM       IN
   IPARM_METIS_RTYPE           - \metis{} parameters (see \metis{} Manual), used only if IPARM_ORDERING_DEFAULT set to API_NO   Default: METIS_RTYPE_SEP1SIDED  IN
   IPARM_METIS_NO2HOP          - \metis{} parameters (see \metis{} Manual), used only if IPARM_ORDERING_DEFAULT set to API_NO   Default: 0                      IN
   IPARM_METIS_NSEPS           - \metis{} parameters (see \metis{} Manual), used only if IPARM_ORDERING_DEFAULT set to API_NO   Default: 1                      IN
   IPARM_METIS_NITER           - \metis{} parameters (see \metis{} Manual), used only if IPARM_ORDERING_DEFAULT set to API_NO   Default: 10                     IN
   IPARM_METIS_UFACTOR         - \metis{} parameters (see \metis{} Manual), used only if IPARM_ORDERING_DEFAULT set to API_NO   Default: 200                    IN
   IPARM_METIS_COMPRESS        - \metis{} parameters (see \metis{} Manual), used only if IPARM_ORDERING_DEFAULT set to API_NO   Default: 1                      IN
   IPARM_METIS_CCORDER         - \metis{} parameters (see \metis{} Manual), used only if IPARM_ORDERING_DEFAULT set to API_NO   Default: 0                      IN
   IPARM_METIS_PFACTOR         - \metis{} parameters (see \metis{} Manual), used only if IPARM_ORDERING_DEFAULT set to API_NO   Default: 0                      IN
   IPARM_METIS_SEED            - \metis{} parameters (see \metis{} Manual)                                                      Default: 3452                   IN
   IPARM_METIS_DBGLVL          - \metis{} parameters (see \metis{} Manual)                                                      Default: 0                      IN

   IPARM_SF_KASS               - Force the use of KASS instead of Fax to perform the symbolic factorization  Default: API_NO                   IN
   IPARM_AMALGAMATION_LVLCBLK  - Amalgamation level                                       Default: 5                   IN
   IPARM_AMALGAMATION_LVLBLAS  - Amalgamation level                                       Default: 5                   IN

   IPARM_STATIC_PIVOTING       - Static pivoting                                          Default: -                   OUT
   IPARM_NNZEROS               - Number of nonzero entries in the factorized matrix       Default: -                   OUT
   IPARM_ALLOCATED_TERMS       - Maximum memory allocated for matrix terms                Default: -                   OUT
   IPARM_BASEVAL               - Baseval used for the matrix                              Default: 0                   IN
   IPARM_MIN_BLOCKSIZE         - Minimum block size                                       Default: 60                  IN
   IPARM_MAX_BLOCKSIZE         - Maximum block size                                       Default: 120                 IN
   IPARM_SCHUR                 - Schur mode                                               Default: API_NO              IN
   IPARM_ISOLATE_ZEROS         - Isolate null diagonal terms at the end of the matrix     Default: API_NO              IN
   IPARM_RHSD_CHECK            - Set to API_NO to avoid RHS redistribution                Default: API_YES             IN
   IPARM_FACTORIZATION         - Factorization mode (see Factorization modes)             Default: API_FACT_LDLT       IN
   IPARM_NNZEROS_BLOCK_LOCAL   - Number of nonzero entries in the local block factorized matrix Default: -                   OUT
   IPARM_CPU_BY_NODE           - Number of CPUs per SMP node                              Default: 0                   IN
   IPARM_BINDTHRD              - Thread binding mode (see Thread binding modes)           Default: API_BIND_AUTO       IN
   IPARM_THREAD_NBR            - Number of threads per MPI process                        Default: 1                   IN
   IPARM_DISTRIBUTION_LEVEL    - Distribution level IGNORE                                Default:                     IN
   IPARM_LEVEL_OF_FILL         - Level of fill for incomplete factorization               Default: 1                   IN
   IPARM_IO_STRATEGY           - IO strategy (see Checkpoints modes)                      Default: API_IO_NO           IN
   IPARM_RHS_MAKING            - Right-hand-side making (see Right-hand-side modes)      Default: API_RHS_B           IN
   IPARM_REFINEMENT            - Refinement type (see Refinement modes)                   Default: API_RAF_GMRES       IN
   IPARM_SYM                   - Symmetric matrix mode (see Symmetric modes)              Default: API_SYM_YES         IN
   IPARM_INCOMPLETE            - Incomplete factorization                                 Default: API_NO              IN
   IPARM_ABS                   - ABS level (Automatic Blocksize Splitting)                Default: 1                   IN
   IPARM_ESP                   - ESP (Enhanced Sparse Parallelism)                        Default: API_NO              IN
   IPARM_GMRES_IM              - GMRES restart parameter                                  Default: 25                  IN
   IPARM_FREE_CSCUSER          - Free user CSC                                            Default: API_CSC_PRESERVE    IN
   IPARM_FREE_CSCPASTIX        - Free internal CSC (Use only without call to Refin. step) Default: API_CSC_PRESERVE    IN
   IPARM_OOC_LIMIT             - Out of core memory limit (Mo)                            Default: 2000                IN
   IPARM_OOC_THREAD            - Out of core thread number IGNORE                         Default: 1                   IN
   IPARM_OOC_ID                - Out of core run ID        IGNORE                         Default: -                   OUT
   IPARM_NB_SMP_NODE_USED      - Number of SMP node used   IGNORE                         Default:                     IN
   IPARM_THREAD_COMM_MODE      - Threaded communication mode (see Communication modes)    Default: API_THREAD_MULT     IN
   IPARM_NB_THREAD_COMM        - Number of thread(s) for communication                    Default: 1                   IN
   IPARM_FILL_MATRIX           - Initialize matrix coefficients (for test only)  IGNORE   Default:                     IN
   IPARM_INERTIA               - Return the inertia (symmetric matrix without pivoting)   Default: -                   OUT
   IPARM_ESP_NBTASKS           - Return the number of tasks generated by ESP              Default: -                   OUT
   IPARM_ESP_THRESHOLD         - Minimal block sizee to switch in ESP mode (128 * 128)    Default: 16384               IN
   IPARM_DOF_COST              - Degree of freedom for cost computation (If different from IPARM_DOF_NBR) Default: 0                    IN
   IPARM_MURGE_REFINEMENT      - Enable refinement in MURGE                               Default: API_YES             IN
   IPARM_STARPU                - Use StarPU runtime                                       Default: API_NO              IN
   IPARM_AUTOSPLIT_COMM        - Automaticaly split communicator to have one MPI task by node             Default: API_NO               IN
   IPARM_FLOAT                 - Indicate the floating point type  IGNORE                 Default: -                   INOUT
   IPARM_PID                   - Pid of the first process (used for naming the log directory) Default: -1                  OUT
   IPARM_ERROR_NUMBER          - Return value                                             Default: -                   OUT
   IPARM_CUDA_NBR              - Number of cuda devices                                   Default: 0                   IN
   IPARM_TRANSPOSE_SOLVE       - Use transposed matrix during solve                       Default: API_NO              IN
   IPARM_STARPU_CTX_DEPTH      - Tree depth of the contexts given to StarPU               Default:3                    IN
   IPARM_STARPU_CTX_NBR        - Number of contexts created                               Default:-1                   INOUT
   IPARM_PRODUCE_STATS         - Compute some statistiques (such as precision error)      Default:API_NO               IN
   IPARM_GPU_CRITERIUM         - Criterium for sorting GPU                                Default:0                    IN
   IPARM_MURGE_MAY_REFINE      - Enable refinement in MURGE                               Default: API_NO             IN
   IPARM_SIZE                  - Iparm Size                IGNORE                         Default:                     IN
*/
enum IPARM_ACCESS {
  IPARM_MODIFY_PARAMETER,
  IPARM_START_TASK,
  IPARM_END_TASK,
  IPARM_VERBOSE,
  IPARM_DOF_NBR,
  IPARM_ITERMAX,
  IPARM_MATRIX_VERIFICATION,
  IPARM_MC64,
  IPARM_ONLY_RAFF,
  IPARM_CSCD_CORRECT,
  IPARM_NBITER,
  IPARM_TRACEFMT,
  IPARM_GRAPHDIST,

  /* Ordering */
  IPARM_ORDERING,
  IPARM_ORDERING_DEFAULT,
  IPARM_SCOTCH_SWITCH_LEVEL,
  IPARM_SCOTCH_CMIN,
  IPARM_SCOTCH_CMAX,
  IPARM_SCOTCH_FRAT,

  IPARM_METIS_CTYPE,
  IPARM_METIS_RTYPE,
  IPARM_METIS_NO2HOP,
  IPARM_METIS_NSEPS,
  IPARM_METIS_NITER,
  IPARM_METIS_UFACTOR,
  IPARM_METIS_COMPRESS,
  IPARM_METIS_CCORDER,
  IPARM_METIS_PFACTOR,
  IPARM_METIS_SEED,
  IPARM_METIS_DBGLVL,

  /* Symbolic Factoization */
  IPARM_SF_KASS,
  IPARM_AMALGAMATION_LVLBLAS,
  IPARM_AMALGAMATION_LVLCBLK,

  IPARM_STATIC_PIVOTING,
  IPARM_NNZEROS,
  IPARM_ALLOCATED_TERMS,
  IPARM_BASEVAL,
  IPARM_MIN_BLOCKSIZE,
  IPARM_MAX_BLOCKSIZE,
  IPARM_SCHUR,
  IPARM_ISOLATE_ZEROS,
  IPARM_RHSD_CHECK,
  IPARM_FACTORIZATION,
  IPARM_NNZEROS_BLOCK_LOCAL,
  IPARM_CPU_BY_NODE,
  IPARM_BINDTHRD,
  IPARM_THREAD_NBR,
  IPARM_DISTRIBUTION_LEVEL,
  IPARM_LEVEL_OF_FILL,
  IPARM_IO_STRATEGY,
  IPARM_RHS_MAKING,
  IPARM_REFINEMENT,
  IPARM_SYM,
  IPARM_INCOMPLETE,
  IPARM_ABS,
  IPARM_ESP,
  IPARM_GMRES_IM,
  IPARM_FREE_CSCUSER,
  IPARM_FREE_CSCPASTIX,
  IPARM_OOC_LIMIT,
  IPARM_OOC_THREAD,
  IPARM_OOC_ID,
  IPARM_NB_SMP_NODE_USED,
  IPARM_THREAD_COMM_MODE,
  IPARM_NB_THREAD_COMM,
  IPARM_FILL_MATRIX,
  IPARM_INERTIA,
  IPARM_ESP_NBTASKS,
  IPARM_ESP_THRESHOLD,
  IPARM_DOF_COST,
  IPARM_MURGE_REFINEMENT,
  IPARM_STARPU,
  IPARM_AUTOSPLIT_COMM,
  IPARM_FLOAT,
  IPARM_PID,
  IPARM_ERROR_NUMBER,
  IPARM_CUDA_NBR,
  IPARM_TRANSPOSE_SOLVE,
  IPARM_STARPU_CTX_DEPTH,
  IPARM_STARPU_CTX_NBR,
  IPARM_PRODUCE_STATS,
  IPARM_GPUS_NBR,
  IPARM_GPU_CRITERIUM,

  IPARM_MURGE_MAY_REFINE,

  IPARM_SIZE
};


/*
 * Backward compatibility
 */
enum IPARM_ACCESS_DEPRECATED {
    IPARM_DEFAULT_ORDERING      = IPARM_ORDERING_DEFAULT,
    IPARM_ORDERING_SWITCH_LEVEL = IPARM_SCOTCH_SWITCH_LEVEL,
    IPARM_ORDERING_CMIN         = IPARM_SCOTCH_CMIN,
    IPARM_ORDERING_CMAX         = IPARM_SCOTCH_CMAX,
    IPARM_ORDERING_FRAT         = IPARM_SCOTCH_FRAT,
    IPARM_AMALGAMATION_LEVEL    = IPARM_AMALGAMATION_LVLCBLK
};

/* Acces au tableau dparm */
/*
   Enum: DPARM_ACCESS

   Floating point parameters tabular accossors

   DPARM_FILL_IN            - Fill-in                                           Default: -                OUT
   DPARM_MEM_MAX            - Maximum memory (-DMEMORY_USAGE)                   Default: -                OUT
   DPARM_EPSILON_REFINEMENT - Epsilon for refinement                            Default: 1e^{-12}         IN
   DPARM_RELATIVE_ERROR     - Relative backward error                           Default: -                OUT
   DPARM_EPSILON_MAGN_CTRL  - Epsilon for magnitude control                     Default: 1e^{-31}         IN
   DPARM_ANALYZE_TIME       - Time for Analyse step (wallclock)                 Default: -                OUT
   DPARM_PRED_FACT_TIME     - Predicted factorization time                      Default: -                OUT
   DPARM_FACT_TIME          - Time for Numerical Factorization step (wallclock) Default: -                OUT
   DPARM_SOLV_TIME          - Time for Solve step (wallclock)                   Default: -                OUT
   DPARM_FACT_FLOPS         - Numerical Factorization flops (rate!)             Default: -                OUT
   DPARM_SOLV_FLOPS         - Solve flops (rate!)                               Default: -                OUT
   DPARM_RAFF_TIME          - Time for Refinement step (wallclock)              Default: -                OUT
   DPARM_SIZE               - Dparm Size         IGNORE                         Default: -                IN
 */
enum DPARM_ACCESS {
  DPARM_FILL_IN                 = 1,
  DPARM_MEM_MAX                 = 2,
  DPARM_EPSILON_REFINEMENT      = 5,
  DPARM_RELATIVE_ERROR          = 6,
  DPARM_SCALED_RESIDUAL         = 7,
  DPARM_EPSILON_MAGN_CTRL       = 10,
  DPARM_ANALYZE_TIME            = 18,
  DPARM_PRED_FACT_TIME          = 19,
  DPARM_FACT_TIME               = 20,
  DPARM_SOLV_TIME               = 21,
  DPARM_FACT_THFLOPS            = 22,
  DPARM_FACT_RLFLOPS            = 25,
  DPARM_SOLV_FLOPS              = 23,
  DPARM_RAFF_TIME               = 24,
  DPARM_A_NORM                  = 25,
  DPARM_SIZE                    = 64 /* Need to be greater or equal to 64 for backward compatibility */
};

#define DPARM_FACT_FLOPS DPARM_FACT_THFLOPS

/** Etapes de résolution de PaStiX */
/*
  Enum: API_TASK

  PaStiX step modes (index IPARM_START_TASK and IPARM_END_TASK)

  API_TASK_INIT       - Set default parameters
  API_TASK_ORDERING   - Ordering
  API_TASK_SYMBFACT   - Symbolic factorization
  API_TASK_ANALYSE    - Tasks mapping and scheduling
  API_TASK_NUMFACT    - Numerical factorization
  API_TASK_SOLVE      - Numerical solve
  API_TASK_REFINE     - Numerical refinement
  API_TASK_CLEAN      - Clean
 */
/* _POS_ 1 */
enum API_TASK {
  API_TASK_INIT       = 0,
  API_TASK_ORDERING   = 1,
  API_TASK_SYMBFACT   = 2,
  API_TASK_ANALYSE    = 3,
  API_TASK_NUMFACT    = 4,
  API_TASK_SOLVE      = 5,
  API_TASK_REFINE     = 6,
  API_TASK_CLEAN      = 7
};

/** Etapes de résolution de PaStiX pour compatibilte avec les anciennes version */
/*
  Enum: API_TASK_OLD

  API_TASK_SCOTCH     - Ordering
  API_TASK_FAX        - Symbolic factorization
  API_TASK_BLEND      - Tasks mapping and scheduling
  API_TASK_SOPALIN    - Numerical factorization
  API_TASK_UPDOWN     - Numerical solve
  API_TASK_REFINEMENT - Numerical refinement
 */
/* _POS_ -1 */
enum API_TASK_OLD {
  API_TASK_SCOTCH     = 1,
  API_TASK_FAX        = 2,
  API_TASK_BLEND      = 3,
  API_TASK_SOPALIN    = 4,
  API_TASK_UPDOWN     = 5,
  API_TASK_REFINEMENT = 6
};

/** Affichage de PaStiX */
/*
   Enum: API_VERBOSE

   Verbose modes (index IPARM_VERBOSE)

   API_VERBOSE_NOT        - Silent mode, no messages
   API_VERBOSE_NO         - Some messages
   API_VERBOSE_YES        - Many messages
   API_VERBOSE_CHATTERBOX - Like a gossip
   API_VERBOSE_UNBEARABLE - Really talking too much...
*/
/* _POS_ 5 */
enum API_VERBOSE {
  API_VERBOSE_NOT        = 0, /* Nothing  */
  API_VERBOSE_NO         = 1, /* Default  */
  API_VERBOSE_YES        = 2, /* Extended */
  API_VERBOSE_CHATTERBOX = 3,
  API_VERBOSE_UNBEARABLE = 4
};

/** Load strategy for graph and ordering */
/*
  Enum: API_IO

  Check-points modes (index IPARM_IO)

  API_IO_NO         - No output or input
  API_IO_LOAD       - Load ordering during ordering step and symbol matrix instead of symbolic factorisation.
  API_IO_SAVE       - Save ordering during ordering step and symbol matrix instead of symbolic factorisation.
  API_IO_LOAD_GRAPH - Load graph during ordering step.
  API_IO_SAVE_GRAPH - Save graph during ordering step.
  API_IO_LOAD_CSC   - Load CSC(d) during ordering step.
  API_IO_SAVE_CSC   - Save CSC(d) during ordering step.
 */
/* _POS_ 6 */
enum API_IO {
  API_IO_NO         = 0,
  API_IO_LOAD       = 1,
  API_IO_SAVE       = 2,
  API_IO_LOAD_GRAPH = 4,
  API_IO_SAVE_GRAPH = 8,
  API_IO_LOAD_CSC   = 16,
  API_IO_SAVE_CSC   = 32
};

/** Genération du second membre */
/*
  Enum: API_RHS

  Right-hand-side modes (index IPARM_RHS)

  API_RHS_B - User's right hand side
  API_RHS_1 - $ \forall i, X_i = 1 $
  API_RHS_I - $ \forall i, X_i = i $

 */
/* _POS_ 7 */
enum API_RHS {
  API_RHS_B = 0, /* Utilisation du second membre fournit */
  API_RHS_1 = 1, /* Utilisation d'un second membre dont tous les coefficients valent 1 */
  API_RHS_I = 2, /* Utilisation d'un second membre tel que RHS(i) = i */
  API_RHS_0 = 3  /* Initialisation en mode ONLY_RAFF d'une solution X0(i) = 0 */
};

/**
 * Type of RHS generated to test the solver
 */
typedef enum pastix_rhstype_e {
  PastixRhsOne,
  PastixRhsI,
  PastixRhsRndX,
  PastixRhsRndB
} pastix_rhstype_t;


/** Type de raffinement utilisé */
/*
  Enum: API_RAF

  Refinement modes (index IPARM_REFINEMENT)

  API_RAF_GMRES   - GMRES
  API_RAF_GRAD    - Conjugate Gradient ($LL^t$ or $LDL^t$ factorization)
  API_RAF_PIVOT   - Iterative Refinement (only for $LU$ factorization)
  API_RAF_BICGSTAB - BICGSTAB
 */
/* _POS_ 8 */
enum API_RAF {
  API_RAF_GMRES   = 0, /* Utilisation de GMRES */
  API_RAF_GRAD    = 1, /* Utilisation du gradient conjugue */
  API_RAF_PIVOT   = 2, /* Utilisation de la methode du pivot */
  API_RAF_BICGSTAB = 3
};

/**
 * Factorization algorithms available for IPARM_FACTORIZATION parameter
 */
typedef enum pastix_factotype_e {
  PastixFactLLT  = 0,
  PastixFactLDLT = 1,
  PastixFactLU   = 2,
  PastixFactLDLH = 3
} pastix_factotype_t;

#define PastixGeneral       111
#define PastixSymmetric     112
#define PastixHermitian     113

/** ****************************************************************************
 *
 *  PaStiX constants - Compatible with CBLAS & LAPACK
 *  The naming and numbering is consistent with:
 *
 *    1) CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz)
 *    2) C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/)
 *    3) Plasma (http://icl.cs.utk.edu/plasma/index.html)
 *
 **/
typedef enum pastix_trans_e {
    PastixNoTrans   = 111,
    PastixTrans     = 112,
    PastixConjTrans = 113
} pastix_trans_t;

typedef enum pastix_uplo_e {
    PastixUpper      = 121,
    PastixLower      = 122,
    PastixUpperLower = 123
} pastix_uplo_t;

typedef enum pastix_diag_e {
    PastixNonUnit = 131,
    PastixUnit    = 132
} pastix_diag_t;

typedef enum pastix_side_e {
    PastixLeft  = 141,
    PastixRight = 142
} pastix_side_t;

typedef enum pastix_normtype_e {
    PastixOneNorm       = 171,
    PastixFrobeniusNorm = 174,
    PastixInfNorm       = 175,
    PastixMaxNorm       = 177
} pastix_normtype_t;

/** Supressing user CSC(D) when not usefull anymore */
/*
  Enum: API_ERASE_CSC

  CSC Management modes (index IPARM_FREE_CSCUSER and IPARM_FREE_CSCPASTIX)

  API_CSC_PRESERVE - Do not free the CSC
  API_CSC_FREE     - Free the CSC when no longer needed
 */
/* _POS_ 11 */
enum API_ERASE_CSC{
  API_CSC_PRESERVE = 0,
  API_CSC_FREE     = 1
};

/** DMP communication mode */
/*
  Enum: API_THREAD_MODE

  Comunication modes (index IPARM_THREAD_COMM_MODE)

  API_THREAD_MULTIPLE      - All threads communicate.
  API_THREAD_FUNNELED      - One thread perform all the MPI Calls.
  API_THREAD_COMM_ONE      - One dedicated communication thread will receive messages.
  API_THREAD_COMM_DEFINED  - Then number of threads receiving the messages is given by IPARM_NB_THREAD_COMM.
  API_THREAD_COMM_NBPROC   - One communication thread per computation thread will receive messages.
 */
/* _POS_ 9 */
enum API_THREAD_MODE {
  API_THREAD_MULTIPLE      = 1,
  API_THREAD_FUNNELED      = 2,
  API_THREAD_COMM_ONE      = 4,
  API_THREAD_COMM_DEFINED  = 8,
  API_THREAD_COMM_NBPROC   = 16
};

/** Thread binding */
/*
  Enum: API_BIND_MODE

  Thread-binding modes (index IPARM_BINTHRD)

  API_BIND_NO   - Do not bind thread
  API_BIND_AUTO - Default binding
  API_BIND_TAB  - Use vector given by pastix_setBind
*/
/* _POS_ 12 */
enum API_BIND_MODE {
  API_BIND_NO      = 0, /* Do not bind threads                            */
  API_BIND_AUTO    = 1, /* Default thread binding                         */
  API_BIND_TAB     = 2  /* Use tabular given by pastix_setBind to bind */
};

/** Boolean */
/*
  Enum: API_BOOLEAN

  Boolean modes (All boolean except IPARM_SYM)

  API_NO  - No
  API_YES - Yes
 */
/* _POS_ 2 */
enum API_BOOLEAN {
  API_NO  = 0,
  API_YES = 1
};

/** Trace format */
/*
  Enum: API_TRACEFMT

  Trace modes (index IPARM_TRACEFMT)

  API_TRACE_PAJE       - Use Paje trace format
  API_TRACE_HUMREAD    - Use human-readable text trace format
  API_TRACE_UNFORMATED - Unformated trace format
 */
/* _POS_ 10 */
enum API_TRACEFMT {
  API_TRACE_PAJE       = 1, /* Use Paje trace format       */
  API_TRACE_HUMREAD    = 2, /* Use text trace format       */
  API_TRACE_UNFORMATED = 3  /* Use unformated trace format */
};

/*
  Enum: API_ORDER

  Ordering modes (index IPARM_ORDERING)

  API_ORDER_SCOTCH   - Use \scotch{} ordering
  API_ORDER_METIS    - Use \metis{} ordering
  API_ORDER_PERSONAL - Apply user's permutation
  API_ORDER_LOAD     - Load ordering from disk
 */
/* _POS_ 11 */
enum API_ORDER {
  API_ORDER_SCOTCH    = 0,
  API_ORDER_METIS     = 1,
  API_ORDER_PERSONAL  = 2,
  API_ORDER_LOAD      = 3,
  API_ORDER_PTSCOTCH  = 4
};

/**
 * Type of elements used in the value array of the matrix and defined th IPARM_FLOAT parameter.
 * (Start at 2 for compatibility with Plasma in case of kernels use)
 */
typedef enum pastix_coeftype_e {
    PastixPattern   = 0,
    PastixFloat     = 2,
    PastixDouble    = 3,
    PastixComplex32 = 4,
    PastixComplex64 = 5
} pastix_coeftype_t;

/*
 * Enum: API_GPU_CRITERIUM
 *
 * Criterium used to decide to put tasks on GPUs.
 *
 * API_GPU_CRITERION_UPDATES  - Number of updates on the panel.
 * API_GPU_CRITERION_CBLKSIZE - Size of the target panel.
 * API_GPU_CRITERION_FLOPS    - Number of FLOP involved in updates.
 * API_GPU_CRITERION_PRIORITY - Priority computed in static scheduler.
 */
enum API_GPU_CRITERIUM {
  API_GPU_CRITERION_UPDATES  = 0,
  API_GPU_CRITERION_CBLKSIZE = 1,
  API_GPU_CRITERION_FLOPS    = 2,
  API_GPU_CRITERION_PRIORITY = 3
};
/*
  Enum: MODULES

  Module Identification number.

  If an error occurs, error value is set to
  MODULE + EER_NUMBER.

  User can catch error by computing iparm[IPARM_ERROR_NUMBER]%100.

  MODULE can be catch by computing iparm[IPARM_ERROR_NUMBER] - iparm[IPARM_ERROR_NUMBER]%100.

  MOD_UNKNOWN - Unknown module
  MOD_SOPALIN - Numerical factorisation module
  MOD_BLEND   - Analysing module
  MOD_SCOTCH  - Scotch module
  MOD_FAX     - Symbolic factorisation module
  MOD_ORDER   - Order module
  MOD_COMMON  - Common module
  MOD_SI      -
  MOD_GRAPH   - Graph module
  MOD_SYMBOL  - Symbol structure module
  MOD_KASS    - Kass module
  MOD_BUBBLE  - Bubble
  MOD_MURGE   - Murge

*/
enum MODULES {
  MOD_UNKNOWN =    0,
  MOD_SOPALIN =  100,
  MOD_BLEND   =  200,
  MOD_SCOTCH  =  300,
  MOD_FAX     =  400,
  MOD_ORDER   =  500,
  MOD_COMMON  =  600,
  MOD_SI      =  700,
  MOD_GRAPH   =  800,
  MOD_SYMBOL  =  900,
  MOD_KASS    = 1000,
  MOD_BUBBLE  = 1100,
  MOD_MURGE   = 1200
};

/* Enum: ERR_NUMBERS

   Error Numbers

   PASTIX_SUCCESS            - No error
   PASTIX_ERR_UNKNOWN        - Unknown error
   PASTIX_ERR_ALLOC          - Allocation error
   PASTIX_ERR_ASSERT         - Error in one assertion
   PASTIX_ERR_NOTIMPLEMENTED - Not implemented feature
   PASTIX_ERR_OUTOFMEMORY    - Not enough memory (OOC)
   PASTIX_ERR_THREAD         - Error with threads
   PASTIX_ERR_INTERNAL       - Internal error
   PASTIX_ERR_BADPARAMETER   - Bad parameters given
   PASTIX_ERR_FILE           - Error in In/Out operations
   PASTIX_ERR_BAD_DEFINEOR   - Error with defines during compilation
   PASTIX_ERR_INTEGER_TYPE   - Error with integer types
   PASTIX_ERR_IO             - Error with input/output
   PASTIX_ERR_MATRIX         - Wrongly defined matrix
   PASTIX_ERR_FLOAT_TYPE     - Wrong type of floating point values
   PASTIX_ERR_STEP_ORDER     - Error in ordering
   PASTIX_ERR_MPI            - Error with MPI calls
*/
/* Need to conserve it MURGE compliant */
enum pastix_error_e {
  PASTIX_SUCCESS            = 0,
  PASTIX_ERR_UNKNOWN        = 1,
  PASTIX_ERR_ALLOC          = 2,
  PASTIX_ERR_ASSERT         = 3,
  PASTIX_ERR_NOTIMPLEMENTED = 4,
  PASTIX_ERR_OUTOFMEMORY    = 5,
  PASTIX_ERR_THREAD         = 6,
  PASTIX_ERR_INTERNAL       = 7,
  PASTIX_ERR_BADPARAMETER   = 8,
  PASTIX_ERR_FILE           = 9,
  PASTIX_ERR_BAD_DEFINE     = 10,
  PASTIX_ERR_INTEGER_TYPE   = 11,
  PASTIX_ERR_IO             = 12,
  PASTIX_ERR_MATRIX         = 13,
  PASTIX_ERR_FLOAT_TYPE     = 14,
  PASTIX_ERR_STEP_ORDER     = 15,
  PASTIX_ERR_MPI            = 16
};

#endif /* _PASTIX_API_H_ */
