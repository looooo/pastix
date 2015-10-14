/**
 *  @file: analyze.c
 *
 *  A simple example that performs only the analyses steps onto the given graph.
 *  These tests doesn't require the values of the matrix.
 *
 */
#include <pastix.h>
#include <csc.h>
#include "../matrix_drivers/drivers.h"

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t    iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double          dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    pastix_driver_t driver;
    char           *filename;
    pastix_csc_t   *csc, *csc2;

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
     * Startup PaStiX
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Read the sparse matrix with the driver
     */
    csc = malloc( sizeof( pastix_csc_t ) );
    cscReadFromFile( driver, filename, csc, MPI_COMM_WORLD );
    free(filename);
    csc2 = spmCheckAndCorrect( csc );
    if ( csc2 != csc ) {
        spmExit( csc );
        free(csc);
        csc = csc2;
    }

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_task_order( pastix_data, csc, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );
    pastix_task_blend( pastix_data );

    spmExit( csc );
    free( csc );
    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    return EXIT_SUCCESS;
}
