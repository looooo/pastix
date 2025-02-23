/**
 *
 * @file api.c
 *
 * PaStiX API routines
 *
 * @copyright 2004-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Esragul Korkmaz
 * @author Gregoire Pichon
 * @author Matias Hastaran
 * @author Tony Delarue
 * @author Alycia Lisito
 * @author Brieuc Nicolas
 * @author Tom Moenne-Loccoz
 * @date 2024-07-09
 *
 **/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define _GNU_SOURCE 1
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#include "common.h"
#if defined(PASTIX_ORDERING_METIS)
#include <metis.h>
#endif
#include "graph/graph.h"
#include "pastix/order.h"
#include "blend/solver.h"
#include "bcsc/bcsc.h"
#include "isched.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "models.h"
#if defined(PASTIX_WITH_PARSEC)
#include "sopalin/parsec/pastix_parsec.h"
#endif
#if defined(PASTIX_WITH_STARPU)
#include "sopalin/starpu/pastix_starpu.h"
#endif
#if defined(HAVE_BLAS_SET_NUM_THREADS)
#if defined(HAVE_BLI_THREAD_SET_NUM_THREADS)
#include <blis.h>
#elif defined(HAVE_MKL_SET_NUM_THREADS)
#include <mkl_service.h>
#elif defined(HAVE_OPENBLAS_SET_NUM_THREADS)
#include "cblas.h"
#endif
#endif

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#if defined(PASTIX_OS_WINDOWS)
#define pastix_mkdir( __str ) mkdir( (__str) )
#else
#define pastix_mkdir( __str ) mkdir( (__str), 0700 )
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

#if defined(PASTIX_WITH_MPI)
static int pastix_mpi_in_use = 0;  /**< Counter of the number of Pastix instances using MPI          */
static int pastix_mpi_init   = 0;  /**< Boolean to know if MPI has been initialized by pastix or not */
static volatile pastix_atomic_lock_t pastix_mpi_lock = PASTIX_ATOMIC_UNLOCKED; /**< Lock to protect the MPI initialization */
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Generate a unique temporary directory to store output files
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          On entry, the pastix_data structure associated to the pastix instance.
 *          On exit, if not already initialized, pastix_data->dir_global and
 *          pastix_data->dir_local are initialized.
 *
 *******************************************************************************/
void
pastix_gendirectories( pastix_data_t *pastix_data )
{
    char **dir_global = &(pastix_data->dir_global);
    mode_t old_mask;
    int rc, len;

    if ( *dir_global != NULL ) {
        return;
    }

    if ( pastix_data->procnum == 0 )
    {
        char *name = pastix_getenv_get_value_str( "PASTIX_OUTPUT_DIR", "pastix" );
        rc = asprintf( dir_global, "%s-XXXXXX", name );
        assert( rc != -1 );
        free( name );

#if !defined(HAVE_MKDTEMP)
        {
            *dir_global = mktemp( *dir_global );
            if ( (*dir_global)[0] == '\0' ) {
                perror( "pastix_gendirectories/global/mktemp" );
                free( *dir_global );
                *dir_global = NULL;
                return;
            }

            old_mask = umask(S_IWGRP | S_IWOTH);
            rc = pastix_mkdir( *dir_global );
            (void)umask(old_mask);

            if ( rc == -1 ) {
                perror( "pastix_gendirectories/global/mkdir" );
                free( *dir_global );
                *dir_global = NULL;
                return;
            }
        }
#else
        {
            old_mask = umask(S_IWGRP | S_IWOTH);
            *dir_global = mkdtemp( *dir_global );
            (void)umask(old_mask);
            if ( *dir_global == NULL ) {
                perror( "pastix_gendirectories/global/mkdtemp" );
                return;
            }
        }
#endif
        /* Broadcast the main directory to everyone */
        len = strlen( pastix_data->dir_global );

        MPI_Bcast( &len, 1, MPI_INT,
                   0, pastix_data->inter_node_comm );
        MPI_Bcast( pastix_data->dir_global, len+1, MPI_CHAR,
                   0, pastix_data->inter_node_comm );

        fprintf( stdout, "OUTPUTDIR: %s\n", pastix_data->dir_global );
    }
    else {
        len = 0;
        MPI_Bcast( &len, 1, MPI_INT,
                   0, pastix_data->inter_node_comm );
        pastix_data->dir_global = malloc( (len+1) * sizeof(char) );
        MPI_Bcast( pastix_data->dir_global, len+1, MPI_CHAR,
                   0, pastix_data->inter_node_comm );
    }

    assert( pastix_data->dir_global != NULL );

    /*
     * Create the local directory name
     */
#if defined(PASTIX_WITH_MPI)
    if (pastix_data->procnbr > 1)
    {
        char *localdir;
        rc = asprintf( &localdir, "%s/%0*d",
                       pastix_data->dir_global,
                       (int)pastix_iceil( pastix_data->procnbr, 10 ),
                       pastix_data->procnum );

        old_mask = umask(S_IWGRP | S_IWOTH);
        rc = pastix_mkdir( localdir );
        (void)umask(old_mask);

        if ( rc == -1 ) {
            perror( "pastix_gendirectories/local/mkdir" );
            free( localdir );
            return;
        }
        pastix_data->dir_local = localdir;
    }
    else
#endif
    {
        pastix_data->dir_local = strdup( pastix_data->dir_global );
    }
    (void)rc;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Generate the full filename within local or global directory.
 *
 *******************************************************************************
 *
 * @param[in] dirname
 *          The pointer to the directory string associated to the instance.
 *          It must have been initialized before calling.
 *
 * @param[in] filename
 *          The filename to create in the unique directory.
 *
 *******************************************************************************
 *
 * @return The filename to use. NULL if failed.
 *
 *******************************************************************************/
char *
pastix_fname( const char *dirname,
              const char *filename )
{
    char *fullname = NULL;
    int rc;

    /* Make sure dirname has been initialized before calling fopen */
    assert( dirname );

    rc = asprintf( &fullname, "%s/%s", dirname, filename );
    if ( rc <= 0 ) {
        pastix_print_error( "pastix_fname: Couldn't not generate the tempory filename for the output file" );
        return NULL;
    }

    return fullname;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Open a file in the unique directory of the pastix instance
 *
 *******************************************************************************
 *
 * @param[in] dirname
 *          The pointer to the directory string associated to the instance.
 *          It must have been initialized before calling.
 *
 * @param[in] filename
 *          The filename to create in the unique directory.
 *
 * @param[in] mode
 *          Opening mode of the file as described by the fopen() function.
 *
 *******************************************************************************
 *
 * @return The FILE structure of the opened file. NULL if failed.
 *
 *******************************************************************************/
FILE *
pastix_fopenw( const char *dirname,
               const char *filename,
               const char *mode )
{
    char *fullname;
    FILE *f = NULL;

    /* Make sure dirname has been initialized before calling fopen */
    assert( dirname );

    fullname = pastix_fname( dirname, filename );
    if ( !fullname ) {
        return NULL;
    }

    f = fopen( fullname, mode );
    free( fullname );

    if ( NULL == f )
    {
        perror( "pastix_fopenw" );
        pastix_print_error( "pastix_fopenw: Couldn't open file: %s with mode %s\n",
                    filename, mode );
        return NULL;
    }

    return f;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Open a file in the current directory in read only mode.
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The filename of the file to open.
 *
 *******************************************************************************
 *
 * @return The FILE structure of the opened file. NULL if failed.
 *
 *******************************************************************************/
FILE *
pastix_fopen( const char *filename )
{
    FILE *f = fopen(filename, "r");

    if (NULL == f)
    {
        perror("pastix_fopen");
        pastix_print_error( "pastix_fopen: Couldn't open file: %s with mode r\n",
                    filename );
        return NULL;
    }

    return f;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Print information about the solver configuration
 *
 *******************************************************************************
 *
 * @param[in] pastix
 *          The main data structure.
 *
 *******************************************************************************/
void
pastixWelcome( const pastix_data_t *pastix )
{
    pastix_print( pastix->procnum, 0, OUT_HEADER,
                  /* Version    */ PASTIX_VERSION_MAJOR, PASTIX_VERSION_MINOR, PASTIX_VERSION_MICRO,
                  /* Sched. seq */ "Enabled",
                  /* Sched. sta */ (pastix->isched ? "Started" : "Disabled"),
                  /* Sched. dyn */ ( (pastix->iparm[IPARM_SCHEDULER] == PastixSchedDynamic) ? "Started" : "Disabled" ),
                  /* Sched. PaR */
#if defined(PASTIX_WITH_PARSEC)
                  (pastix->parsec == NULL ? "Enabled" : "Started" ),
#else
                  "Disabled",
#endif
                  /* Sched. SPU */
#if defined(PASTIX_WITH_STARPU)
                  (pastix->starpu == NULL ? "Enabled" : "Started" ),
#else
                  "Disabled",
#endif
                  /* MPI nbr    */ pastix->procnbr,
                  /* Thrd nbr   */ (int)(pastix->iparm[IPARM_THREAD_NBR]),
                  /* GPU nbr    */ (int)(pastix->iparm[IPARM_GPU_NBR]),
                  /* MPI mode   */ pastix_mpithreadmode_getstr( pastix->iparm[IPARM_MPI_THREAD_LEVEL] ),
                  /* Distrib    */ ((pastix->iparm[IPARM_TASKS2D_LEVEL] == 0) ? "1D" : "2D"),
                                   ((pastix->iparm[IPARM_TASKS2D_LEVEL] <  0) ? ((long)pastix->iparm[IPARM_TASKS2D_WIDTH]) :
                                                                               -((long)pastix->iparm[IPARM_TASKS2D_LEVEL])),
                  /* Block size */ (long)pastix->iparm[IPARM_MIN_BLOCKSIZE],
                                   (long)pastix->iparm[IPARM_MAX_BLOCKSIZE],
                  /* Model CPU  */ pastix->cpu_models->name,
                  /* Model GPU  */ pastix->gpu_models->name,
                  /* Strategy   */ ((pastix->iparm[IPARM_COMPRESS_WHEN] == PastixCompressNever) ? "No compression" :
                                    (pastix->iparm[IPARM_COMPRESS_WHEN] == PastixCompressWhenBegin) ? "Memory Optimal" : "Just-In-Time") );


    if ( pastix->iparm[IPARM_COMPRESS_WHEN] != PastixCompressNever ) {
        pastix_print( pastix->procnum, 0, OUT_HEADER_LR,
                      /* Tolerance       */ (double)pastix->dparm[DPARM_COMPRESS_TOLERANCE],
                      /* Compress method */ compmeth_shnames[pastix->iparm[IPARM_COMPRESS_METHOD]],
                      /* Compress width  */ (long)pastix->iparm[IPARM_COMPRESS_MIN_WIDTH],
                      /* Compress height */ (long)pastix->iparm[IPARM_COMPRESS_MIN_HEIGHT],
                      /* Min ratio       */ (double)pastix->dparm[DPARM_COMPRESS_MIN_RATIO],
                      /* Tolerance used  */ pastix->iparm[IPARM_COMPRESS_RELTOL] ? "Relative" : "Absolute",
                      /* Ortho method    */ ((pastix->iparm[IPARM_COMPRESS_ORTHO] == PastixCompressOrthoCGS) ? "CGS" :
                                             (pastix->iparm[IPARM_COMPRESS_ORTHO] == PastixCompressOrthoQR)  ? "QR" : "partialQR"),
                      /* Splitting strategy    */ ((pastix->iparm[IPARM_SPLITTING_STRATEGY] == PastixSplitNot)  ? "Not used" :
                                                   (pastix->iparm[IPARM_SPLITTING_STRATEGY] == PastixSplitKway) ? "KWAY" : "KWAY and projections"),
                      /* Levels of projections */ (long)pastix->iparm[IPARM_SPLITTING_LEVELS_PROJECTIONS],
                      /* Levels of kway        */ (long)pastix->iparm[IPARM_SPLITTING_LEVELS_KWAY],
                      /* Projections distance  */ (long)pastix->iparm[IPARM_SPLITTING_PROJECTIONS_DISTANCE],
                      /* Projections depth     */ (long)pastix->iparm[IPARM_SPLITTING_PROJECTIONS_DEPTH],
                      /* Projections width     */ (long)pastix->iparm[IPARM_SPLITTING_PROJECTIONS_WIDTH] );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Print summary information
 *
 *******************************************************************************
 *
 * @param[in] pastix
 *          The main data structure.
 *
 *******************************************************************************/
void
pastixSummary( const pastix_data_t *pastix )
{
    (void)pastix;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Initialize the iparm and dparm arrays to their default
 * values.
 *
 * This is performed only if iparm[IPARM_MODIFY_PARAMETER] is set to 0.
 *
 *******************************************************************************
 *
 * @param[inout] iparm
 *          The integer array of parameters to initialize.
 *
 * @param[inout] dparm
 *          The floating point array of parameters to initialize.
 *
 *******************************************************************************/
void
pastixInitParam( pastix_int_t *iparm,
                 double       *dparm )
{
    memset( iparm, 0, IPARM_SIZE * sizeof(pastix_int_t) );
    memset( dparm, 0, DPARM_SIZE * sizeof(double) );

    iparm[IPARM_VERBOSE]               = PastixVerboseNo;
    iparm[IPARM_IO_STRATEGY]           = PastixIONo;

    /* Stats */
    iparm[IPARM_NNZEROS]               = 0;
    iparm[IPARM_NNZEROS_BLOCK_LOCAL]   = 0;
    iparm[IPARM_ALLOCATED_TERMS]       = 0;
    iparm[IPARM_PRODUCE_STATS]         = 0;
    iparm[IPARM_TRACE]                 = PastixTraceNot;

    /* Scaling */
    iparm[IPARM_MC64]                  = 0;

    /*
     * Ordering parameters
     */
#if defined(PASTIX_ORDERING_SCOTCH)
    iparm[IPARM_ORDERING]              = PastixOrderScotch;
#elif defined(PASTIX_ORDERING_METIS)
    iparm[IPARM_ORDERING]              = PastixOrderMetis;
#else
    iparm[IPARM_ORDERING]              = -1;
#endif
    iparm[IPARM_ORDERING_DEFAULT]      = 1;

    /* Scotch */
    {
        iparm[IPARM_SCOTCH_MT]           = 1;
        iparm[IPARM_SCOTCH_SWITCH_LEVEL] = 120;
        iparm[IPARM_SCOTCH_CMIN]         = 0;
        iparm[IPARM_SCOTCH_CMAX]         = 100000;
        iparm[IPARM_SCOTCH_FRAT]         = 8;
    }

    /* Metis */
    {
#if defined(PASTIX_ORDERING_METIS)
        iparm[IPARM_METIS_CTYPE   ] = METIS_CTYPE_SHEM;
        iparm[IPARM_METIS_RTYPE   ] = METIS_RTYPE_SEP1SIDED;
#else
        iparm[IPARM_METIS_CTYPE   ] = -1;
        iparm[IPARM_METIS_RTYPE   ] = -1;
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
    }

    /* Symbolic factorization */
    iparm[IPARM_AMALGAMATION_LVLCBLK]  = 5;
    iparm[IPARM_AMALGAMATION_LVLBLAS]  = 5;

    /* Reordering */
    iparm[IPARM_REORDERING_SPLIT]      = 0;
    iparm[IPARM_REORDERING_STOP]       = PASTIX_INT_MAX;

    /* Splitting */
    iparm[IPARM_SPLITTING_STRATEGY]             = PastixSplitKway;
    iparm[IPARM_SPLITTING_LEVELS_PROJECTIONS]   = 0;
    iparm[IPARM_SPLITTING_LEVELS_KWAY]          = PASTIX_INT_MAX;
    iparm[IPARM_SPLITTING_PROJECTIONS_DEPTH]    = 3;
    iparm[IPARM_SPLITTING_PROJECTIONS_DISTANCE] = 3;
    iparm[IPARM_SPLITTING_PROJECTIONS_WIDTH]    = 1;

    /* Analyze */
    iparm[IPARM_MIN_BLOCKSIZE]         = 160;
    iparm[IPARM_MAX_BLOCKSIZE]         = 320;
    iparm[IPARM_TASKS2D_LEVEL]         = -1;
    iparm[IPARM_TASKS2D_WIDTH]         = iparm[IPARM_MIN_BLOCKSIZE];
    iparm[IPARM_ALLCAND]               = 0;

    /* Incomplete */
    iparm[IPARM_INCOMPLETE]            = 0;
    iparm[IPARM_LEVEL_OF_FILL]         = 0;

    /* Factorization */
    iparm[IPARM_FACTORIZATION]         = PastixFactLU;
    iparm[IPARM_FACTO_LOOK_SIDE]       = PastixFactRightLooking;
    iparm[IPARM_STATIC_PIVOTING]       = 0;
    iparm[IPARM_FREE_CSCUSER]          = 0;
    iparm[IPARM_SCHUR_FACT_MODE]       = PastixFactModeLocal;

    /* Solve */
    iparm[IPARM_TRANSPOSE_SOLVE]       = PastixNoTrans;
    iparm[IPARM_SCHUR_SOLV_MODE]       = PastixSolvModeLocal;
    iparm[IPARM_APPLYPERM_WS]          = 1;

    /* Refinement */
    iparm[IPARM_REFINEMENT]            = PastixRefineGMRES;
    iparm[IPARM_NBITER]                = 0;
    iparm[IPARM_ITERMAX]               = 250;
    iparm[IPARM_GMRES_IM]              = 25;

    /* Context */
    iparm[IPARM_SCHEDULER]             = PastixSchedDynamic;
    iparm[IPARM_THREAD_NBR]            = -1;
    iparm[IPARM_SOCKET_NBR]            = -1;
    iparm[IPARM_AUTOSPLIT_COMM]        = 0;

    /* GPU */
    iparm[IPARM_GPU_NBR]               = 0;
    iparm[IPARM_GPU_MEMORY_PERCENTAGE] = 95;
    iparm[IPARM_GPU_MEMORY_BLOCK_SIZE] = 32 * 1024;
    iparm[IPARM_GLOBAL_ALLOCATION]     = 0;

    /* Compression */
    iparm[IPARM_COMPRESS_MIN_WIDTH]    = 128;
    iparm[IPARM_COMPRESS_MIN_HEIGHT]   = 20;
    iparm[IPARM_COMPRESS_WHEN]         = PastixCompressNever;
    iparm[IPARM_COMPRESS_METHOD]       = PastixCompressMethodPQRCP;
    iparm[IPARM_COMPRESS_ORTHO]        = PastixCompressOrthoCGS;
    iparm[IPARM_COMPRESS_PRESELECT]    = 1;
    iparm[IPARM_COMPRESS_ILUK]         = -2;

    /* Mixed-Precision */
    iparm[IPARM_MIXED] = 0;

    /* MPI modes */
    iparm[IPARM_MPI_THREAD_LEVEL] = PastixMpiNone;

    /* Subset for old pastix interface  */
    iparm[IPARM_MODIFY_PARAMETER] = 1;
    iparm[IPARM_START_TASK] = PastixTaskOrdering;
    iparm[IPARM_END_TASK]   = PastixTaskClean;
    iparm[IPARM_FLOAT]      = -1;
    iparm[IPARM_MTX_TYPE]   = -1;
    iparm[IPARM_DOF_NBR]    = 1;

    dparm[DPARM_FILL_IN]            =  0.;
    dparm[DPARM_EPSILON_REFINEMENT] = -1.;
    dparm[DPARM_RELATIVE_ERROR]     = -1.;
    dparm[DPARM_EPSILON_MAGN_CTRL]  =  0.;
    dparm[DPARM_ANALYZE_TIME]       =  0.;
    dparm[DPARM_PRED_FACT_TIME]     =  0.;
    dparm[DPARM_FACT_TIME]          =  0.;
    dparm[DPARM_SOLV_TIME]          =  0.;
    dparm[DPARM_FACT_FLOPS]         =  0.;
    dparm[DPARM_FACT_THFLOPS]       =  0.;
    dparm[DPARM_FACT_RLFLOPS]       =  0.;
    dparm[DPARM_FACT_ENERGY]        =  0.;
    dparm[DPARM_MEM_FR]             =  0.;
    dparm[DPARM_MEM_LR]             =  0.;
    dparm[DPARM_SOLV_FLOPS]         =  0.;
    dparm[DPARM_SOLV_THFLOPS]       =  0.;
    dparm[DPARM_SOLV_RLFLOPS]       =  0.;
    dparm[DPARM_SOLV_ENERGY]        =  0.;
    dparm[DPARM_REFINE_TIME]        =  0.;
    dparm[DPARM_A_NORM]             = -1.;
    dparm[DPARM_COMPRESS_TOLERANCE] =  1e-8;
    dparm[DPARM_COMPRESS_MIN_RATIO] =  1.;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_internal
 *
 * @brief Internal function that setups the multiple communicators in order
 * to perform the ordering step in MPI only mode, and the factorization in
 * MPI+Thread mode with the same amount of ressources.
 *
 *******************************************************************************
 *
 * @param[inout] pastix
 *          The pastix data structure to initialize.
 *
 * @param[in] comm
 *          The MPI communicator associated to the pastix data structure
 *
 * @param[in] autosplit
 *          Enable automatic split when multiple processes run on the same node
 *
 *******************************************************************************/
static inline void
apiInitMPI( pastix_data_t *pastix,
            PASTIX_Comm    comm,
            int            autosplit )
{
    /*
     * Setup all communicators for autosplitmode and initialize number/rank of
     * processes.
     */
    if ( comm == 0 ) {
        comm = MPI_COMM_WORLD;
    }
    pastix->pastix_comm = comm;
    MPI_Comm_size(pastix->pastix_comm, &(pastix->procnbr));
    MPI_Comm_rank(pastix->pastix_comm, &(pastix->procnum));

#if defined(PASTIX_WITH_MPI)
    if ( autosplit )
    {
        int     i, len;
        char    procname[MPI_MAX_PROCESSOR_NAME];
        int     rc, key = pastix->procnum;
        int64_t color;
        (void)rc;

        /*
         * Get hostname to generate a hash that will be the color of each node
         * MPI_Get_processor_name is not used as it can returned different
         * strings for processes of a same physical node.
         */
        rc = gethostname(procname, MPI_MAX_PROCESSOR_NAME-1);
        assert(rc == 0);
        procname[MPI_MAX_PROCESSOR_NAME-1] = '\0';
        len = strlen( procname );

        /* Compute hash */
        color = 0;
        for (i = 0; i < len; i++) {
            color = color*256*sizeof(char) + procname[i];
        }

        /* Create intra-node communicator */
        MPI_Comm_split(pastix->pastix_comm, color, key, &(pastix->intra_node_comm));
        MPI_Comm_size(pastix->intra_node_comm, &(pastix->intra_node_procnbr));
        MPI_Comm_rank(pastix->intra_node_comm, &(pastix->intra_node_procnum));

        /* Create inter-node communicator */
        MPI_Comm_split(pastix->pastix_comm, pastix->intra_node_procnum, key, &(pastix->inter_node_comm));
        MPI_Comm_size(pastix->inter_node_comm, &(pastix->inter_node_procnbr));
        MPI_Comm_rank(pastix->inter_node_comm, &(pastix->inter_node_procnum));
    }
    else
#endif
    {
        pastix->intra_node_comm    = MPI_COMM_SELF;
        pastix->intra_node_procnbr = 1;
        pastix->intra_node_procnum = 0;
        pastix->inter_node_comm    = pastix->pastix_comm;
        pastix->inter_node_procnbr = pastix->procnbr;
        pastix->inter_node_procnum = pastix->procnum;
    }

    assert( pastix->inter_node_procnbr * pastix->intra_node_procnbr == pastix->procnbr );
    (void)autosplit;
    (void)comm;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Initialize the solver instance with a bintab array to specify the
 * thread binding.
 *
 * @remark You should always prefer the pastixInit() function when hwloc is
 * available, and use the pastixInitWithAffinity() function only if you know
 * what you want to do with your threads.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The main data structure.
 *
 * @param[in] pastix_comm
 *          The MPI communicator.
 *
 * @param[inout] iparm
 *          The integer array of parameters to initialize.
 *
 * @param[inout] dparm
 *          The floating point array of parameters to initialize.
 *
 * @param[in] bindtab
 *          Integer array of size iparm[IPARM_THREAD_NBR] that will specify the
 *          thread binding. NULL if let to the system.
 *          Each thread i will be bound to to the core bindtab[i] if
 *          bindtab[i] >= 0, or not bound if bindtab[i] < 0.
 *          If other libraries of the main application are spawning their own threads
 *          too (eg. OpenMP), we strongly recommend not to bind the main thread,
 *          and let bindtab[0] = -1 to avoid binding impact on other libraries.
 *
 *******************************************************************************/
void
pastixInitWithAffinity( pastix_data_t **pastix_data,
                        PASTIX_Comm     pastix_comm,
                        pastix_int_t   *iparm,
                        double         *dparm,
                        const int      *bindtab )
{
    pastix_data_t *pastix;

    /*
     * Allocate pastix_data structure when we enter PaStiX for the first time.
     */
    MALLOC_INTERN(pastix, 1, pastix_data_t);
    memset( pastix, 0, sizeof(pastix_data_t) );

    /*
     * Initialize iparm/dparm vectors and set them to default values if not set
     * by the user.
     */
    if ( iparm[IPARM_MODIFY_PARAMETER] == 0 ) {
        pastixInitParam( iparm, dparm );
    }

    /*
     * Check if MPI is initialized
     */
#if defined(PASTIX_WITH_MPI)
    {
        int provided = MPI_THREAD_SINGLE;
        int flag = 0;

        pastix_atomic_lock( &pastix_mpi_lock );
        pastix_mpi_in_use++;

        MPI_Initialized(&flag);
        if ( !flag ) {
            MPI_Init_thread( NULL, NULL, MPI_THREAD_MULTIPLE, &provided );
            pastix_mpi_init = 1;
        }
        else {
            MPI_Query_thread( &provided );
        }

        switch ( provided ) {
        case MPI_THREAD_MULTIPLE:
            iparm[IPARM_MPI_THREAD_LEVEL] = PastixMpiThreadMultiple;
            break;
        case MPI_THREAD_SERIALIZED:
            iparm[IPARM_MPI_THREAD_LEVEL] = PastixMpiThreadSerialized;
            break;
        case MPI_THREAD_FUNNELED:
            iparm[IPARM_MPI_THREAD_LEVEL] = PastixMpiThreadFunneled;
            break;
        case MPI_THREAD_SINGLE:
            iparm[IPARM_MPI_THREAD_LEVEL] = PastixMpiThreadSingle;
            break;
        default:
            iparm[IPARM_MPI_THREAD_LEVEL] = PastixMpiNone;
        }
        pastix_atomic_unlock( &pastix_mpi_lock );
    }
#endif

    /*
     * TODO : Replace this id by a list of ids to know
     *        which sub-problem this instance represents.
     */
    pastix->id    = 953833;
    pastix->iparm = iparm;
    pastix->dparm = dparm;

    /*
     * Set bits to flush subnormals in mixed-precision and avoid slow-downs
     */
    if ( pastix->iparm[IPARM_FTZ] ) {
        set_ftz();
    }

    pastix->steps = 0;
    pastix->sched = PastixSchedDynamic;

    pastix->isched = NULL;
#if defined(PASTIX_WITH_PARSEC)
    pastix->parsec = NULL;
#endif
#if defined(PASTIX_WITH_STARPU)
    pastix->starpu = NULL;
#endif

    apiInitMPI( pastix, pastix_comm, iparm[IPARM_AUTOSPLIT_COMM] );

    if ( (pastix->intra_node_procnbr > 1) &&
         (pastix->iparm[IPARM_THREAD_NBR] != -1 ) ) {
        pastix_print( pastix->procnum, 0,
                      "WARNING: Thread number forced by MPI autosplit feature\n" );
        iparm[IPARM_THREAD_NBR] = pastix->intra_node_procnbr;
    }

#if defined(PASTIX_GENERATE_MODEL)
    pastix_print( pastix->procnum, 0,
                  "WARNING: PaStiX compiled with -DPASTIX_GENERATE_MODEL forces single thread computations\n" );
    iparm[IPARM_THREAD_NBR] = 1;
#endif

    /*
     * Start the internal threads
     */
    pastix->isched = ischedInit( pastix->iparm[IPARM_THREAD_NBR], bindtab );
    pastix->iparm[IPARM_THREAD_NBR] = pastix->isched->world_size;

    if ( ( pastix->iparm[IPARM_SOCKET_NBR] == -1 ) ||
         ( pastix->iparm[IPARM_SOCKET_NBR] > pastix->isched->socketsnbr ) ) {
        pastix->iparm[IPARM_SOCKET_NBR] = pastix->isched->socketsnbr;
    }

    /*
     * Start PaRSEC if compiled with it and scheduler set to PaRSEC
     */
#if defined(PASTIX_WITH_PARSEC)
    if ( (pastix->parsec == NULL) &&
         (iparm[IPARM_SCHEDULER] == PastixSchedParsec) ) {
        int argc = 0;
        pastix_parsec_init( pastix, &argc, NULL, bindtab );
    }
#endif /* defined(PASTIX_WITH_PARSEC) */

    /*
     * Start StarPU if compiled with it and scheduler set to StarPU
     */
#if defined(PASTIX_WITH_STARPU)
    if ( (pastix->starpu == NULL) &&
         (iparm[IPARM_SCHEDULER] == PastixSchedStarPU) ) {
        int argc = 0;
        pastix_starpu_init( pastix, &argc, NULL, bindtab );
    }
#endif /* defined(PASTIX_WITH_STARPU) */

    pastix->graph      = NULL;
    pastix->schur_n    = 0;
    pastix->schur_list = NULL;
    pastix->zeros_n    = 0;
    pastix->zeros_list = NULL;
    pastix->ordemesh   = NULL;

    pastix->symbmtx    = NULL;

    pastix->bcsc       = NULL;
    pastix->solvmatr   = NULL;

    pastix->cpu_models = NULL;
    pastix->gpu_models = NULL;

    pastix->dir_global = NULL;
    pastix->dir_local  = NULL;

    pastixModelsLoad( pastix );

    if (iparm[IPARM_VERBOSE] > PastixVerboseNot) {
        pastixWelcome( pastix );
    }

    /* Initialization step done, overwrite anything done before */
    pastix->steps = STEP_INIT;

    papiEnergyInit( pastix->iparm[IPARM_SOCKET_NBR] );

    /* If PASTIX_WITH_EZTRACE starts (and pauses) trace. */
#if defined(PASTIX_WITH_EZTRACE)
    kernelsTraceInit( pastix, pastix->iparm[IPARM_TRACE] );
#endif

    *pastix_data = pastix;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Initialize the solver instance
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The main data structure.
 *
 * @param[in] pastix_comm
 *          The MPI communicator.
 *
 * @param[inout] iparm
 *          The integer array of parameters to initialize.
 *
 * @param[inout] dparm
 *          The floating point array of parameters to initialize.
 *
 *******************************************************************************/
void
pastixInit( pastix_data_t **pastix_data,
            PASTIX_Comm     pastix_comm,
            pastix_int_t   *iparm,
            double         *dparm )
{
    pastixInitWithAffinity( pastix_data, pastix_comm,
                            iparm, dparm, NULL );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Finalize the solver instance
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The main data structure.
 *
 *******************************************************************************/
void
pastixFinalize( pastix_data_t **pastix_data )
{
    pastix_data_t *pastix = *pastix_data;

    papiEnergyFinalize();

    /* If PASTIX_WITH_EZTRACE stops trace. */
#if defined(PASTIX_WITH_EZTRACE)
    kernelsTraceFinalize( *pastix_data );
#endif

    pastixSummary( *pastix_data );

    ischedFinalize( pastix->isched );

    if ( pastix->graph != NULL )
    {
        graphExit( pastix->graph );
        memFree_null( pastix->graph );
    }

    if ( pastix->ordemesh != NULL )
    {
        pastixOrderExit( pastix->ordemesh );
        memFree_null( pastix->ordemesh );
    }

    if ( pastix->symbmtx != NULL )
    {
        pastixSymbolExit( pastix->symbmtx );
        memFree_null( pastix->symbmtx );
    }

    if ( pastix->solvloc != NULL )
    {
        solverExit( pastix->solvloc );
        memFree_null( pastix->solvloc );
    }

    if ( pastix->solvglob != NULL )
    {
        solverExit( pastix->solvglob );
        memFree_null( pastix->solvglob );
    }

    if ( pastix->solvmatr != NULL )
    {
        pastix->solvmatr = NULL;
    }

    if ( pastix->bcsc != NULL )
    {
        bcscExit( pastix->bcsc );
        memFree_null( pastix->bcsc );
    }

    if (pastix->schur_list != NULL )
    {
        memFree_null( pastix->schur_list );
    }
#if defined(PASTIX_WITH_PARSEC)
    if (pastix->parsec != NULL) {
        pastix_parsec_finalize( pastix );
    }
#endif /* defined(PASTIX_WITH_PARSEC) */
#if defined(PASTIX_WITH_STARPU)
    if (pastix->starpu != NULL) {
        pastix_starpu_finalize( pastix );
    }
#endif /* defined(PASTIX_WITH_STARPU) */

#if defined(PASTIX_WITH_MPI)
    pastix_atomic_lock( &pastix_mpi_lock );
    pastix_mpi_in_use--;
    if ( (pastix_mpi_in_use == 0) && pastix_mpi_init ) {
        MPI_Finalize();
    }
    pastix_atomic_unlock( &pastix_mpi_lock );
#endif

    if ( pastix->cpu_models != NULL ) {
        pastixModelsFree( pastix->cpu_models );
        pastix->cpu_models = NULL;
    }
    if ( pastix->gpu_models != NULL ) {
        pastixModelsFree( pastix->gpu_models );
        pastix->gpu_models = NULL;
    }

    if ( pastix->dir_global != NULL ) {
        free( pastix->dir_global );
    }
    if ( pastix->dir_local != NULL ) {
        free( pastix->dir_local );
    }
    memFree_null(*pastix_data);
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Dump the iparm and dparm parameters in a CSV file.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The main data structure.
 *
 *******************************************************************************/
void
pastixDumpParam( const pastix_data_t *pastix_data )
{
    FILE   *csv      = NULL;
    char   *fullname = NULL;
    int     rc;
    int32_t lidx;
    static volatile int32_t id = 0;

    if( pastix_data->inter_node_procnum != 0 ) {
        return;
    }

    lidx = pastix_atomic_add_32b( &id, 1 );
    rc  = asprintf( &fullname, "idparam_%d.csv", lidx );

    if ( rc <= 0 ) {
        pastix_print_error( "pastixDumpParam: Couldn't not generate the filename for the output file" );
        return;
    }

    csv = pastix_fopenw( pastix_data->dir_global, fullname, "w" );
    pastix_param2csv( pastix_data, csv );

    fclose(csv);
    free(fullname);
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Check the values of iparm and dparm arrays.
 *
 *******************************************************************************
 *
 * @param[in] iparm
 *          The iparm options array. Can be NULL.
 *
 * @param[in] dparm
 *          The dparm options array. Can be NULL.
 *
 *******************************************************************************
 *
 * @return The amount of iparm/dparm with an incorrect value.
 *
 *******************************************************************************/
int
pastixCheckParam( const pastix_int_t *iparm,
                  const double       *dparm )
{
    int irc = 0;
    int drc = 0;

    if ( iparm != NULL ) {
        irc = iparm_check_values( iparm );
    }
    if ( dparm != NULL ) {
        drc = dparm_check_values( dparm );
    }

    return irc + drc;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Return the current number of threads used for blas calls.
 *
 *******************************************************************************
 *
 * @return The number of threads.
 *
 *******************************************************************************/
int
pastixBlasGetNumThreads( void )
{
#if defined(HAVE_BLI_THREAD_SET_NUM_THREADS)
    return bli_thread_get_num_threads();
#elif defined(HAVE_MKL_SET_NUM_THREADS)
    return mkl_get_max_threads();
#elif defined(HAVE_OPENBLAS_SET_NUM_THREADS)
    return openblas_get_num_threads();
#else
    return 1;
#endif
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Set the number of threads for blas calls (BLIS, MKL, OpenBLAS) and
 * return the previous number of blas threads.
 *
 * @param[in] nt
 *          Number of threads.
 *
 *******************************************************************************
 *
 * @return The previous number of blas threads.
 *
 *******************************************************************************/
int
pastixBlasSetNumThreads( int nt )
{
    int prevnt = pastixBlasGetNumThreads();
#if defined(HAVE_BLI_THREAD_SET_NUM_THREADS)
    bli_thread_set_num_threads( nt );
#elif defined(HAVE_MKL_SET_NUM_THREADS)
    mkl_set_num_threads( nt );
#elif defined(HAVE_OPENBLAS_SET_NUM_THREADS)
    openblas_set_num_threads( nt );
#endif
    (void)nt;
    return prevnt;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_api
 *
 * @brief Set the number of threads for blas calls (BLIS, MKL, OpenBLAS) to 1
 * and return the previous number of blas threads.
 *
 * @param[inout] pastix
 *          The pastix data structure where the original number of threads is
 *          saved, see blas_nt_origin, just before setting to 1.
 *
 *******************************************************************************
 *
 * @return The previous number of blas threads.
 *
 *******************************************************************************/
int
pastixBlasSetNumThreadsOne( void )
{
    return pastixBlasSetNumThreads( 1 );
}
