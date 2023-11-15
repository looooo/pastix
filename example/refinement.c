/**
 * @file refinement.c
 *
 * @brief A refinement example that runs iterative methods without applying
 *        the preconditioner (factorization and solve steps are not called).
 *        Based on the step-by-step example.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.1
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2023-07-21
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
    spm_driver_t driver;
    char           *filename = NULL;
    pastix_spm_t   *spm, spm2;
    void           *x, *b, *x0 = NULL;
    size_t          size;
    int             scatter = 0;
    int             check   = 1;
    int             nrhs    = 1;
    int             rc      = 0;

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
     * Startup PaStiX
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Read the sparse matrix with the driver
     */
    spm = malloc( sizeof( pastix_spm_t ) );
    if ( scatter ) {
        rc = spmReadDriverDist( driver, filename, spm, MPI_COMM_WORLD );
    }
    else {
        rc = spmReadDriver( driver, filename, spm );
    }
    free( filename );
    if ( rc != SPM_SUCCESS ) {
        pastixFinalize( &pastix_data );
        return rc;
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

    /**
     * Normalize A matrix (optional, but recommended for low-rank functionality)
     */
    double normA = spmNorm( SpmFrobeniusNorm, spm );
    spmScal( 1./normA, spm );

    size = pastix_size_of( spm->flttype ) * spm->nexp * nrhs;
    x = malloc( size );
    b = malloc( size );
    if ( check > 1 ) {
        x0 = malloc( size );
    }

    pastix_subtask_spm2bcsc( pastix_data, spm );

    /**
     * Generates the b and x vector such that A * x = b
     * Compute the norms of the initial vectors if checking purpose.
     */
    if ( check )
      {
        spmGenRHS( SpmRhsRndX, nrhs, spm, x0, spm->nexp, b, spm->nexp );
        memcpy( x, b, size );
      }
    else
      {
        spmGenRHS( SpmRhsRndB, nrhs, spm, NULL, spm->nexp, x, spm->nexp );

        /* Apply also normalization to b vectors */
        spmScalMat( 1./normA, spm, nrhs, b, spm->nexp );

        /* Save b for refinement */
        memcpy( b, x, size );
      }

    /**
     * Refine the linear system
     */
    pastix_task_refine( pastix_data, spm->nexp, nrhs, b, spm->nexp, x, spm->nexp );

    if ( check ) {
      rc |= spmCheckAxb( dparm[DPARM_EPSILON_REFINEMENT], nrhs, spm, x0, spm->nexp, b, spm->nexp, x, spm->nexp );
    }
    if (iparm[IPARM_NBITER]>=iparm[IPARM_ITERMAX]) {
      rc |= 1;
    }

    spmExit( spm );
    free( spm );
    free( b );
    free( x );
    if ( x0 ) {
        free( x0 );
    }
    pastixFinalize( &pastix_data );

    return rc;
}

/**
 * @endcode
 */
