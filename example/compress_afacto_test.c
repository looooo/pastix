/**
 *  @file: analyze.c
 *
 *  A simple example that performs only the analyses steps onto the given graph.
 *  These tests doesn't require the values of the matrix.
 *
 */
#include <pastix.h>
#include <spm.h>
#include <limits.h>
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
    SolverMatrix   *copy;
    int rc;

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

    /* Compute the norm of A, to scale the epsilon parameter for pivoting */
    fprintf( stdout, "-- ||A||_2  =                                   %lg\n",
             spmNorm( PastixFrobeniusNorm, spm ) );

    rc = pastix_subtask_spm2bcsc( pastix_data, spm );
    if (rc != PASTIX_SUCCESS)
        return EXIT_FAILURE;

    rc = pastix_subtask_bcsc2ctab( pastix_data, spm );
    if (rc != PASTIX_SUCCESS)
        return EXIT_FAILURE;

    /* Let's copy the uncompressed coeftab */
    copy = solverCopy( pastix_data->solvmatr, spm->flttype );

    /* Compress/Uncompress the copy to see what is lost */
    coeftabCompress[spm->flttype-2]( copy );
    coeftabUncompress[spm->flttype-2]( copy );

    /* Let's check the difference between the two matrices */
    rc = coeftabDiff[spm->flttype-2]( pastix_data->solvmatr, copy );
    if (rc) {
        fprintf(stderr,
                " -- Test compression on matrix before factorization: FAILED !!! --\n"
                "    %d cblk have not been correctly compressed\n",
                rc );
    }

    spmExit( spm );
    free( spm );
    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    return EXIT_SUCCESS;
}
