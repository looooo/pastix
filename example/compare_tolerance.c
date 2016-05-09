/**
 *  @file: simple.c
 *
 *  A simple example :
 *  read the matrix, check it is correct and correct it if needed,
 *  then run pastix in one call.
 *
 */
#include <pastix.h>
#include <spm.h>
#include "../matrix_drivers/drivers.h"

// TODO: will need some header cleanup
#include "common.h"
#include "blend/solver.h"
#include "sopalin/coeftab.h"

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t    iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double          dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    pastix_driver_t driver;
    char           *filename;
    pastix_spm_t   *spm, *spm2;
    SolverMatrix   *dense, *lr;
    int             rc;
    char            tolerance[128];
    double          init_tolerance;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );

    /**
     * Check if parameters are well set
     */
    char *tol  = getenv("TOLERANCE");
    if (tol == NULL){
        printf("Set TOLERANCE variable\n");
        exit(1);
    }
    init_tolerance = atof(tol);

    iparm[IPARM_ITERMAX]          = 100;
    iparm[IPARM_REORDERING_SPLIT] = 0;

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
    spm = malloc( sizeof( pastix_spm_t ) );
    cscReadFromFile( driver, filename, spm, MPI_COMM_WORLD );
    free(filename);
    spm2 = spmCheckAndCorrect( spm );
    if ( spm2 != spm ) {
        spmExit( spm );
        free(spm);
        spm = spm2;
    }

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_task_order( pastix_data, spm, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );

    pastix_task_reordering( pastix_data );
    pastix_task_blend( pastix_data );


    /**
     * Copy the A solver structure
     */
    dense = solverCopy( pastix_data->solvmatr, spm->flttype );

    /**
     * Perform the low-rank numerical factorization
     */
    pastix_task_sopalin( pastix_data, spm );
    lr = pastix_data->solvmatr;

    /**
     * Perform the infinite precision factorization
     */
    sprintf(tolerance, "0.00000000000000000001\n");
    setenv("TOLERANCE", tolerance, 1);
    pastix_data->solvmatr = dense;
    pastix_task_sopalin( pastix_data, spm );

    /**
     * Compare LR factorization with infinite-precision factorization
     */
    rc = coeftabDiff[spm->flttype-2]( lr, dense );
    if (rc) {
        fprintf(stderr,
                " -- Test factors between tolerances %e and %s: FAILED !!! -- \n"
                "    %d cblk have not been correctly compressed\n",
                init_tolerance, tolerance, rc );
    }


    spmExit( spm );
    free( spm );
    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    return EXIT_SUCCESS;
}
