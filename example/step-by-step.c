/**
 * @file step-by-step.c
 *
 * @brief A step-by-step example that runs one full analyze (ordering, symbolic
 *       factorization, analyze), then loops over 2 factorizations that are both
 *       used for 2 solves each.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Hastaran Matias
 * @date 2017-01-17
 *
 * @ingroup pastix_examples
 * @code
 *
 */
#include <pastix.h>
#include <spm.h>

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t    iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double          dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    pastix_driver_t driver;
    char           *filename;
    pastix_spm_t   *spm, *spm2;
    void           *x, *b, *x0 = NULL;
    size_t          size;
    int             check = 1;
    int             nrhs = 1;
    int             nfact = 2;
    int             nsolv = 2;
    long            i,j;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );

    /**
     * Get options from command line
     */
    pastix_getOptions( argc, argv,
                       iparm, dparm,
                       &check, &driver, &filename );

    /**
     * Startup PaStiX
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Read the sparse matrix with the driver
     */
    spm = malloc( sizeof( pastix_spm_t ) );
    spmReadDriver( driver, filename, spm, MPI_COMM_WORLD );
    free(filename);

    spmPrintInfo( spm, stdout );

    spm2 = spmCheckAndCorrect( spm );
    if ( spm2 != spm ) {
        spmExit( spm );
        free(spm);
        spm = spm2;
    }

    /**
     * Generate a Fake values array if needed for the numerical part
     */
    if ( spm->flttype == PastixPattern ) {
        spmGenFakeValues( spm );
    }

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_subtask_order( pastix_data, spm, NULL );
    pastix_subtask_symbfact( pastix_data );
    pastix_subtask_reordering( pastix_data );
    pastix_subtask_blend( pastix_data );

    size = pastix_size_of( spm->flttype ) * spm->n;
    x = malloc( size );
    b = malloc( size );
    if ( check > 1 ) {
        x0 = malloc( size );
    }

    /* Do nfact factorization */
    for (i = 0; i < nfact; i++)
    {
        /**
         * Perform the numerical factorization
         */
        pastix_subtask_spm2bcsc( pastix_data, spm );
        pastix_subtask_bcsc2ctab( pastix_data );
        pastix_subtask_sopalin( pastix_data );
        for (j = 0; j < nsolv; j++)
        {
            /**
             * Generates the b and x vector such that A * x = b
             * Compute the norms of the initial vectors if checking purpose.
             */
            if ( check )
            {
                spmGenRHS( PastixRhsRndX, nrhs, spm, x0, spm->n, b, spm->n );
                memcpy( x, b, size );
            }
            else {
                spmGenRHS( PastixRhsRndB, nrhs, spm, NULL, spm->n, x, spm->n );
                /* Save b for refinement: TODO: make 2 examples w/ or w/o refinement */
                memcpy( b, x, size );
            }

            /**
             * Solve the linear system
             */
            pastix_task_solve( pastix_data, nrhs, x, spm->n );
            pastix_task_refine( pastix_data, x, nrhs, b );

            if ( check ) {
                spmCheckAxb( nrhs, spm, x0, spm->n, b, spm->n, x, spm->n );
            }
        }
    }

    pastixFinalize( &pastix_data );

    spmExit( spm );
    free(spm);
    free(b);
    free(x);
    if (x0) free(x0);

    return EXIT_SUCCESS;
}

/**
 * @endcode
 */
