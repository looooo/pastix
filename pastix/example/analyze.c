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
    pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix           */
    pastix_int_t    iparm[IPARM_SIZE];  /* integer parameters for pastix                             */
    double          dparm[DPARM_SIZE];  /* floating parameters for pastix                            */
    pastix_driver_t driver;        /* Matrix driver(s) requested by user                        */
    char           *filename;           /* Filename(s) given by user                                 */
    pastix_csc_t    csc;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );

    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    iparm[IPARM_IO_STRATEGY]   = API_IO_SAVE;
    iparm[IPARM_MIN_BLOCKSIZE] = 60;
    iparm[IPARM_MAX_BLOCKSIZE] = 120;

    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    split_level       = 0;
    stop_criteria     = INT_MAX;
    stop_when_fitting = 0;

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
    cscReadFromFile( driver, filename, &csc, MPI_COMM_WORLD );
    free(filename);

    printf("Reordering parameters are %d %d %d\n", split_level, stop_criteria, stop_when_fitting);
    if (stop_when_fitting != 0 && stop_when_fitting != 1){
        fprintf(stderr, "Fatal error in reordering parameters -R split_level:stop_criteria:stop_when_fitting\n");
        exit(1);
    }

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_task_order( pastix_data, &csc, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );
    pastix_task_reordering( pastix_data, split_level, stop_criteria, stop_when_fitting );
    pastix_task_blend( pastix_data );

    spmExit( &csc );
    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    return EXIT_SUCCESS;
}
