/**
 * @file ordering_grid.c
 *
 * @brief Partition by hand a regular Laplacian
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Pierre Ramet
 * @date 2017-05-17
 *
 * @ingroup pastix_examples
 * @code
 *
 */
#include <pastix.h>
#include <spm.h>
#include <pastix/order.h>
#include "drivers/laplacian.h"

int main (int argc, char **argv)
{
    pastix_data_t    *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t      iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double            dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    pastix_driver_t   driver;
    char             *filename;
    pastix_spm_t     *spm, *spm2;
    void             *x, *b, *x0 = NULL;
    size_t            size;
    int               check = 1;
    int               nrhs  = 1;
    int               rc    = 0;
    pastix_order_t   *ord;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );

    /**
     * Get options from command line
     */
    pastixGetOptions( argc, argv,
                      iparm, dparm,
                      &check, &driver, &filename );

    /**
     * Build optimal ordering (for laplacians only)
     */
    if ( driver != PastixDriverLaplacian ){
        fprintf(stderr, "Grid ordering can be used with PastixDriverLaplacian driver only\n");
        return EXIT_FAILURE;
    }

    /**
     * Read the sparse matrix with the driver
     */
    spm = malloc( sizeof( pastix_spm_t ) );
    spmReadDriver( driver, filename, spm, MPI_COMM_WORLD );

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
     * Parse Laplacian dimensions and generate associate manual ordering
     */
    {
        pastix_int_t      dim1, dim2, dim3;
        pastix_coeftype_t flttype;
        double            alpha, beta;

        laplacian_parse_info( filename, &flttype, &dim1, &dim2, &dim3, &alpha, &beta );
        ord = malloc(sizeof(pastix_order_t));
        pastixOrderGrid( &ord, dim1, dim2, dim3 );

        (void) flttype;
        (void) alpha;
        (void) beta;
    }
    free(filename);

    /**
     * Startup PaStiX
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    iparm[IPARM_ORDERING] = PastixOrderPersonal;
    pastix_subtask_order( pastix_data, spm, ord );
    pastix_subtask_symbfact( pastix_data );
    pastix_subtask_reordering( pastix_data );
    pastix_subtask_blend( pastix_data );

    /**
     * Normalize A matrix (optional, but recommended for low-rank functionality)
     */
    double normA = spmNorm( PastixFrobeniusNorm, spm );
    spmScalMatrix( 1./normA, spm );

    /**
     * Perform the numerical factorization
     */
    pastix_task_numfact( pastix_data, spm );

    /**
     * Generates the b and x vector such that A * x = b
     * Compute the norms of the initial vectors if checking purpose.
     */
    size = pastix_size_of( spm->flttype ) * spm->n * nrhs;
    x = malloc( size );
    b = malloc( size );

    if ( check )
    {
        if ( check > 1 ) {
            x0 = malloc( size );
        }
        spmGenRHS( PastixRhsRndX, nrhs, spm, x0, spm->n, b, spm->n );
        memcpy( x, b, size );
    }
    else {
        spmGenRHS( PastixRhsRndB, nrhs, spm, NULL, spm->n, x, spm->n );

        /* Apply also normalization to b vectors */
        spmScalRHS( spm->flttype, 1./normA, spm->n, nrhs, b, spm->n );

        /* Save b for refinement */
        memcpy( b, x, size );
    }

    /**
     * Solve the linear system (and perform the optional refinement)
     */
    pastix_task_solve( pastix_data, nrhs, x, spm->n );
    pastix_task_refine( pastix_data, spm->n, nrhs, b, spm->n, x, spm->n );

    if ( check )
    {
        rc = spmCheckAxb( nrhs, spm, x0, spm->n, b, spm->n, x, spm->n );

        if (x0) free(x0);
    }

    pastixOrderExit( ord );
    spmExit( spm );
    free( spm );
    free( x );
    free( b );
    pastixFinalize( &pastix_data );

    return rc;
}

/**
 * @endcode
 */
