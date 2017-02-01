/**
 *
 * @file schur.c
 *
 * PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 * LaBRI, University of Bordeaux 1 and IPB.
 * Construct the schur complement of the matrix.
 *
 * @version 5.1.0
 * @author  Hastaran Matias
 * @date    2017-01-17
 *
 **/
#include <pastix.h>
#include <spm.h>
#include "drivers.h"

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t    iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double          dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    pastix_driver_t driver;
    char           *filename;
    pastix_spm_t   *spm, *spm2;
    void           *x0, *x, *b;
    size_t          size;
    int             check = 2;
    int             nrhs = 1;
    double          normA;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );

    /**
     * Get options from command line
     */
    pastix_ex_getoptions( argc, argv,
                          iparm, dparm,
                          &driver, &filename );

    /**
     * Read the sparse matrix with the driver
     */
    spm = malloc( sizeof( pastix_spm_t ) );
    spmReadDriver( driver, filename, spm, MPI_COMM_WORLD );
    free(filename);

    spm2 = spmCheckAndCorrect( spm );
    if ( spm2 != spm ) {
        spmExit( spm );
        free(spm);
        spm = spm2;
    }

    /**
     * Scal the matrix to avoid unexpected rouding errors
     */
    normA = spmNorm( PastixFrobeniusNorm, spm );
    spmScal( 1./normA, spm );

    /**
     * Startup PaStiX
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Initialize the schur list with the first third of the unknowns
     */
    {
        pastix_int_t n = spm->gN / 3;

        if ( n > 0 ) {
            pastix_int_t i;
            pastix_int_t baseval = spmFindBase(spm);
            pastix_int_t *list = (pastix_int_t*)malloc(n * sizeof(pastix_int_t));

            for (i=0; i<n; i++) {
                list[i] = i+baseval;
            }
            pastix_setSchurUnknownList( pastix_data, n, list );
            free( list );

            iparm[IPARM_SCHUR] = API_YES;
        }
    }

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_task_order( pastix_data, spm, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );
    pastix_task_blend( pastix_data );

    /**
     * Perform the numerical factorization
     */
    pastix_task_sopalin( pastix_data, spm );

    /**
     * Generates the b and x vector such that A * x = b
     * Compute the norms of the initial vectors if checking purpose.
     */
    size = pastix_size_of( spm->flttype ) * spm->n;
    x = malloc( size );
    b = malloc( size );

    if ( check )
    {
        if ( check > 1 ) {
            x0 = malloc( size );
        } else {
            x0 = NULL;
        }
        spmGenRHS( PastixRhsRndX, nrhs, spm, x0, spm->n, b, spm->n );
        memcpy( x, b, size );
    }
    else {
        spmGenRHS( PastixRhsRndB, nrhs, spm, NULL, spm->n, x, spm->n );

        /* Save b for refinement: TODO: make 2 examples w/ or w/o refinement */
        b = malloc( size );
        memcpy( b, x, size );
    }

    /**
     * Solve the linear system
     */
    pastix_task_solve( pastix_data, spm, nrhs, x, spm->n );

    pastix_task_raff(pastix_data, x, nrhs, b);

    if ( check )
    {
        spmCheckAxb( nrhs, spm, x0, spm->n, b, spm->n, x, spm->n );

        if (x0) free(x0);

    }

    spmExit( spm );
    free( spm );
    free(x);
    free(b);
    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    return EXIT_SUCCESS;
}
