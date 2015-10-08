/**
 *  @file: clifexec.c
 *
 *  A simple example :
 *  read the matrix, check it is correct and correct it if needed,
 *  then run pastix in one call.
 *
 */
#include <pastix.h>
#include <csc.h>
#include <limits.h>
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

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

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

    /**
     * Perform ordering, symbolic factorization
     */
    pastix_task_order( pastix_data, &csc, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );

    spmExit( &csc );
    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    return EXIT_SUCCESS;
}
