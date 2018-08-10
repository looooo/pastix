/**
 * @file step-by-step.c
 *
 * @brief A simple example that reads the matrix and then runs pastix in one call.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Hastaran Matias
 * @date 2018-07-16
 *
 * @ingroup pastix_examples
 * @code
 *
 */
#include <pastix.h>
#include <spm.h>
#include "pastix/order.h"
#include "common/common.h"
#include "blend/solver.h"
#include "common/pastixdata.h"

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t    iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double          dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    spm_driver_t    driver;
    char           *filename;
    spmatrix_t     *spm, *spm2;
    void           *x, *b, *x0 = NULL;
    size_t          size;
    int             check = 1;
    int             variadic = 1; /* 0 => Constant DoFs, 1 => variadic DoFs */
    int             dofmax   = 4; /* Maximal DoF per unknown */
    int             nrhs  = 10;
    int             rc    = 0;

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
     * Read the sparse matrix with the driver
     */
    spm = malloc( sizeof( spmatrix_t ) );
    spmReadDriver( driver, filename, spm );
    free( filename );

    spmPrintInfo( spm, stdout );

    /*
     * For now many operations from the spm check do not work with values and
     * multiple dofs (variadic or constant), thus we remove the value for
     * testing
     */
    {
        if ( spm->values ) {
            free( spm->values );
            spm->values = NULL;
        }
        spm->flttype = SpmPattern;

        spm2 = spmDofExtend( spm, variadic, dofmax );
        spmExit( spm );
        free( spm );
        spm = spm2;
    }

    /* Check the spm */
    spm2 = malloc( sizeof( spmatrix_t ) );
    rc = spmCheckAndCorrect( spm, spm2 );
    if ( rc != 0 ) {
        spmExit( spm );
        memcpy( spm, spm2, sizeof(pastix_spm_t) );
    }
    free( spm2 );

    /**
     * Startup PaStiX
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_subtask_order( pastix_data, spm, NULL );
    pastix_subtask_symbfact( pastix_data );
    pastix_subtask_reordering( pastix_data );

    /**
     * Expand the ordering, symbol, and spm for the final analyse step, and for
     * numerical computations.
     */
    pastixSymbolExpand( pastix_data->symbmtx );
    pastixOrderExpand( pastix_data->ordemesh, spm );

    spm2 = spmExpand( spm );
    spmExit( spm );
    memcpy( spm, spm2, sizeof(pastix_spm_t) );
    free( spm2 );

#if !defined(NDEBUG)
    pastixOrderCheck( pastix_data->ordemesh );
#endif

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
        spmGenRHS( SpmRhsRndX, nrhs, spm, x0, spm->n, b, spm->n );
        memcpy( x, b, size );
    }
    else {
        spmGenRHS( SpmRhsRndB, nrhs, spm, NULL, spm->n, x, spm->n );

        /* Apply also normalization to b vectors */
        spmScalVector( spm->flttype, 1./normA, spm->n * nrhs, b, 1 );

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
        rc = spmCheckAxb( dparm[DPARM_EPSILON_REFINEMENT], nrhs, spm, x0, spm->n, b, spm->n, x, spm->n );

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
