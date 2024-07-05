/**
 * @file simple_trans.c
 *
 * @brief This example extends the simple one to show how to avoid costly
 * CSR/CSC conversion.
 *
 * For users who have a CSR matrix and prefer to avoid the
 * potentially costly translation to CSC format, it is possible to
 * factorize A^t and then solve the problem with these results.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Tony Delarue
 * @author Mathieu Faverge
 * @author Pierre Ramet
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
    spm_driver_t    driver;
    char           *filename = NULL;
    spmatrix_t     *spm, spm2;
    void           *x, *b, *x0 = NULL;
    void           *tmp;
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

    /**
     * Check we are in the correct conditions
     */
    if ( (spm->mtxtype == SpmHermitian) &&
         ((spm->flttype == SpmComplex32) || (spm->flttype == SpmComplex64)) &&
         ((iparm[IPARM_FACTORIZATION] == PastixFactLLH) || (iparm[IPARM_FACTORIZATION] == PastixFactLDLH)) )
    {
        fprintf( stderr, "WARNING: This is not possible to use the CSR trick with Hermitian matrices and LL^h / LDL^h factorizations.\n" );

        spmExit( spm );
        free( spm );
        pastixFinalize( &pastix_data );
        return 0;
    }

    /**
     * Check the spm structure to symmetrize pattern and remove duplicated entries
     */
    rc = spmCheckAndCorrect( spm, &spm2 );
    if ( rc != 0 ) {
        spmExit( spm );
        *spm = spm2;
        rc = 0;
    }

    /**
     * In this example, we show how to work with a CSR matrix, so let's convert it
     */
    spmConvert( SpmCSR, spm );

    /**
     * Generate a fake values array if needed for the numerical part
     */
    if ( spm->flttype == SpmPattern ) {
        spmGenFakeValues( spm );
    }

    /**
     * Normalize A matrix (optional, but recommended for low-rank functionality)
     */
    double normA = spmNorm( SpmFrobeniusNorm, spm );
    spmScal( 1./normA, spm );

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
     * PaStiX works with CSC matrices only for now, so let's lure him by telling him
     * that our CSR matrix is a CSC matrix.
     * Note that this is not available yet for Hermitian matrices and LL^h,
     * LDL^h factorizations.
     */
    {
        spm->fmttype = SpmCSC;
        tmp = spm->colptr;
        spm->colptr = spm->rowptr;
        spm->rowptr = tmp;

        iparm[IPARM_TRANSPOSE_SOLVE] = PastixTrans;
    }

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_task_analyze( pastix_data, spm );

    /**
     * Perform the numerical factorization
     */
    pastix_task_numfact( pastix_data, spm );

    /**
     * Solve the linear system (and perform the optional refinement)
     */
    pastix_task_solve( pastix_data, spm->nexp, nrhs, x, spm->nexp );
    pastix_task_refine( pastix_data, spm->nexp, nrhs, b, spm->nexp, x, spm->nexp );

    /**
     * Restore the original CSR format of the matrix to perform the correct
     * operation in the check.
     */
    {
        spm->fmttype = SpmCSR;
        tmp = spm->colptr;
        spm->colptr = spm->rowptr;
        spm->rowptr = tmp;
    }

    /**
     * Check the results
     */
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
