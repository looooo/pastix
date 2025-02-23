/**
 * @file reentrant.c
 *
 * @brief A reentrant example that runs two threads then run two instances of the solver in each thread.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Mathieu Faverge
 * @author Matias Hastaran
 * @author Pierre Ramet
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2024-07-05
 *
 * @ingroup pastix_examples
 * @code
 *
 */
#include <pthread.h>
#include <pastix.h>
#include <spm.h>

/**
 *  Struct: solv_param
 *
 *  Structure containing information to give
 *  to each thread.
 */
typedef struct solve_param {
    pastix_int_t  iparm[IPARM_SIZE];
    double        dparm[DPARM_SIZE];
    char         *filename;
    spm_driver_t  driver;
    int           check;
    int           scatter;
    int           id;
    int           rc;
} solve_param_t;

/**
 * Function: solve_smp
 *
 * Thread routine to launch the solver
 *
 * Parameters:
 *   arg - a pointer to a <solve_param> structure.
 */
static void *solve_smp(void *arg)
{
    pastix_data_t *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    spmatrix_t    *spm;
    spmatrix_t     spm2;
    void          *x, *b, *x0 = NULL;
    size_t         size;
    int            check;
    int            scatter;
    int            rc;
    int            nrhs = 1;
    solve_param_t  param = *(solve_param_t *)arg;

    if ( param.iparm[IPARM_THREAD_NBR] == -1) {
        param.iparm[IPARM_THREAD_NBR] = 2;
    }
    check   = param.check;
    scatter = param.scatter;

    /**
     * Startup PaStiX
     */
    {
        int *bindtab = malloc( param.iparm[IPARM_THREAD_NBR] * sizeof(int) );
        int i;

        for(i=0; i<param.iparm[IPARM_THREAD_NBR]; i++ ) {
            bindtab[i] = param.iparm[IPARM_THREAD_NBR] * param.id + i;
        }

        pastixInitWithAffinity( &pastix_data, MPI_COMM_WORLD,
                                param.iparm, param.dparm,
                                bindtab );
        free( bindtab );
    }

#if defined(PASTIX_WITH_MPI)
    {
        int size, rank;
        MPI_Comm_size( MPI_COMM_WORLD, &size );
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        if( size > 1 ) {
            if ( rank == 0 ) {
                fprintf( stderr, "\nWarning: Reentrant doesn't work with multiple MPI instances\n" );
                fprintf( stderr, "Quitting now\n" );
            }
            pastixFinalize( &pastix_data );
            exit(0);
        }
    }
#endif

    /**
     * Read the sparse matrix with the driver
     */
    spm = malloc( sizeof( spmatrix_t ) );
    if ( scatter ) {
        rc = spmReadDriverDist( param.driver, param.filename, spm, MPI_COMM_WORLD );
    }
    else {
        rc = spmReadDriver( param.driver, param.filename, spm );
    }
    if ( rc != SPM_SUCCESS ) {
        pastixFinalize( &pastix_data );
        return NULL;
    }

    spmPrintInfo( spm, stdout );

    rc = spmCheckAndCorrect( spm, &spm2 );
    if ( rc != 0 ) {
        spmExit( spm );
        *spm = spm2;
        rc = 0;
    }

    /**
     * Generate a Fake values array if needed for the numerical part
     */
    if ( spm->flttype == SpmPattern ) {
        spmGenFakeValues( spm );
    }

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_task_analyze( pastix_data, spm );

    /**
     * Normalize A matrix (optional, but recommended for low-rank functionality)
     */
    double normA = spmNorm( SpmFrobeniusNorm, spm );
    spmScal( 1./normA, spm );

    /**
     * Perform the numerical factorization
     */
    pastix_task_numfact( pastix_data, spm );

    /**
     * Generates the b and x vector such that A * x = b
     * Compute the norms of the initial vectors if checking purpose.
     */
    size = pastix_size_of( spm->flttype ) * spm->nexp * nrhs;
    x = malloc( size );
    b = malloc( size );

    if ( check )
    {
        if ( check > 1 ) {
            x0 = malloc( size );
        }
        spmGenRHS( SpmRhsRndX, nrhs, spm, x0, spm->nexp, b, spm->nexp );
        memcpy( x, b, size );
    }
    else {
        spmGenRHS( SpmRhsRndB, nrhs, spm, NULL, spm->nexp, x, spm->nexp );

        /* Apply also normalization to b vectors */
        spmScalMat( 1./normA, spm, nrhs, b, spm->nexp );

        /* Save b for refinement */
        memcpy( b, x, size );
    }

    /**
     * Solve the linear system (and perform the optional refinement)
     */
    pastix_task_solve( pastix_data, spm->nexp, nrhs, x, spm->nexp );
    pastix_task_refine( pastix_data, spm->nexp, nrhs, b, spm->nexp, x, spm->nexp );

    if ( check )
    {
        param.rc = spmCheckAxb( param.dparm[DPARM_EPSILON_REFINEMENT], nrhs,
                                spm, x0, spm->nexp, b, spm->nexp, x, spm->nexp );

        if ( x0 ) {
            free( x0 );
        }
    }

    spmExit( spm );
    free( x );
    free( b );
    free( spm );
    pastixFinalize( &pastix_data );

    return NULL;
}

int main (int argc, char **argv)
{
    pastix_int_t   iparm[IPARM_SIZE];      /*< Integer in/out parameters for pastix                */
    double         dparm[DPARM_SIZE];      /*< Floating in/out parameters for pastix               */
    spm_driver_t   driver;                 /*< Matrix driver(s) requested by user                  */
    char          *filename = NULL;               /*< Filename(s) given by user                           */
    int            nbcallingthreads = 2;
    solve_param_t *solve_param;
    pthread_t     *threads;
    int            scatter = 0;
    int            check   = 1;
    int            rc      = 0;
    int            i;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );

    /**
     * Get options from command line
     */
    pastixGetOptions( argc, argv,
                      iparm, dparm,
                      &check, &scatter, &driver, &filename );

    /**
     *    Set parameters for each thread
     */
    solve_param = (solve_param_t*) malloc(nbcallingthreads * sizeof(solve_param_t));
    threads     = (pthread_t*)     malloc(nbcallingthreads * sizeof(pthread_t));

    for (i = 0; i < nbcallingthreads; i++)
    {
        memcpy(solve_param[i].iparm, iparm, sizeof(solve_param[i].iparm));
        memcpy(solve_param[i].dparm, dparm, sizeof(solve_param[i].dparm));
        solve_param[i].check    = check;
        solve_param[i].scatter  = scatter;
        solve_param[i].id       = i;
        solve_param[i].driver   = driver;
        solve_param[i].filename = filename;
        solve_param[i].rc       = 0;

        /**
         *   Launch instance of solver
         */
        pthread_create(&threads[i], NULL, solve_smp, (void *)&solve_param[i]);
    }

    /**
     *     Wait for the end of thread
     */
    for (i = 0; i < nbcallingthreads; i++) {
        pthread_join(threads[i],(void**)NULL);
        rc += solve_param[i].rc;
    }

    free( filename );
    free( threads );
    free( solve_param );
    return rc;
}

/**
 * @endcode
 */
