/**
 * @file multidof.c
 *
 * @brief A simple example that checks the multiple dof implementation.
 *
 * This test takes any matrix as input, arbitrarily expand the matrix by adding
 * multiple degrees of freedom to each unknow. A constant number (4) if check <
 * 3, and a variadic number (\in [1,4]) if check >= 3. To do this test, the
 * values of the matrix are dropped, and random values are generated afterwards
 * because some spm functions are not implemented yet with the support for
 * multiple degrees of freedom.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2022-10-13
 *
 * @ingroup pastix_examples
 * @code
 *
 */
#include <spm.h>
#include <pastix.h>

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t    iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double          dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    spm_driver_t    driver;
    char           *filename = NULL;
    spmatrix_t     *spm, *spm2;
    void           *x, *b, *x0 = NULL;
    size_t          size;
    int             scatter  = 0;
    int             check    = 1;
    int             variadic = 0; /* 0 => Constant DoFs, 1 => variadic DoFs */
    int             dofmax   = 4; /* Maximal DoF per unknown */
    int             nrhs     = 10;
    int             rc       = 0;

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
     * Let's hack this test to use different dofs:
     *    if check < 3, then we check constant dof,
     *    otherwise, we check variadic dof
     */
    variadic = check / 3;
    check    = check % 3;

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

    spm2 = malloc( sizeof( spmatrix_t ) );
    {
        rc = spmDofExtend( spm, variadic, dofmax, spm2 );
        /* Switch spm and spm2 */
        if ( rc == SPM_SUCCESS ) {
            spmatrix_t *spmtmp = spm;
            spm  = spm2;
            spm2 = spmtmp;
        }
        else {
            fprintf( stderr, "multidof example: Failed to generate fake df for the matrix. Switch back to single dof mode.\n" );
            assert( 0 );
        }
    }

    /* Check the spm */
    rc = spmCheckAndCorrect( spm, spm2 );
    if ( rc != 0 ) {
        spmExit( spm );
        memcpy( spm, spm2, sizeof(pastix_spm_t) );
        rc = 0;
    }
    free( spm2 );

    spmPrintInfo( spm, stdout );

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_subtask_order( pastix_data, spm, NULL );
    pastix_subtask_symbfact( pastix_data );
    pastix_subtask_reordering( pastix_data );

    /**
     * Expand the ordering, symbol, and spm for the final analyse step, and for
     * numerical computations. This is required only if the dof is variadic as
     * this case is not yet fully supported.
     */
    if ( spm->dof < 1 ) {
        pastixExpand( pastix_data, spm );
    }

    pastix_subtask_blend( pastix_data );

    /**
     * Generate a Fake values array if needed for the numerical part
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
