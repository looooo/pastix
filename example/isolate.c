/**
 * @file isolate.c
 *
 * @brief A simple example that shows how to isolate a set of problematic
 * unknowns (disconnected, 0 on the diagonal) to the end of the factorization to
 * fix numerical issues that can occrus when they are put too early in the
 * ordering.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2024-07-05
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
    spm_driver_t    driver;
    char           *filename = NULL;
    spmatrix_t     *spm, spm2;
    void           *x, *b, *x0 = NULL;
    size_t          size;
    int             check   = 1;
    int             scatter = 0;
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
    spm = malloc( sizeof( spmatrix_t ) );
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
    if ( spm->flttype == SpmPattern ) {
        spmGenFakeValues( spm );
    }

    /**
     * Initialize the list of elements with zeroes on the diagonal to force them
     * to the end of the ordering.  Here the example takes randomly a third of
     * the unknowns, but you should define this list based on zeroes diagonal
     * entries, and/or on disconnected unknowns.
     *
     * Note that it can be used jointly with Schur unknowns and they will be
     * placed just before the set of Schur unknowns in the ordering.
     */
    {
        pastix_int_t nzeros = spm->gN / 3;
        /* Set to a maximum to avoid memory problem with the test */
        nzeros = (nzeros > 5000) ? 5000 : nzeros;

        if ( nzeros > 0 ) {
            pastix_int_t i;
            pastix_int_t baseval = spmFindBase(spm);
            pastix_int_t *list = (pastix_int_t*)malloc(nzeros * sizeof(pastix_int_t));

            for (i=0; i<nzeros; i++) {
                list[i] = i+baseval;
            }
            pastixIsolateUnknowns( pastix_data, nzeros, list );
            free( list );
        }
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
        rc = spmCheckAxb( dparm[DPARM_EPSILON_REFINEMENT], nrhs, spm, x0, spm->nexp, b, spm->nexp, x, spm->nexp );

        if ( x0 ) {
            free( x0 );
        }
    }

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
