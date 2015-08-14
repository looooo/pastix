/**
 *  @file: simple.c
 *
 *  A simple example :
 *  read the matrix, check it is correct and correct it if needed,
 *  then run pastix in one call.
 *
 *  @precisions normal z => s, d, c
 */
#include <pastix.h>
#include <csc.h>
#include <lapacke.h>
#include "../matrix_drivers/drivers.h"

void CORE_dplrnt( int m, int n, double *A, int lda,
                  int gM, int m0, int n0, unsigned long long int seed );

int core_dgeadd(int trans, int M, int N, double alpha,
                const double *A, int LDA,
                double *B, int LDB);

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t    iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double          dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    pastix_driver_t driver;
    char           *filename;
    pastix_csc_t    csc;
    void           *x0, *x, *b;
    size_t          size;
    int             check = 2;
    int             nrhs = 1;

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
    cscReadFromFile( driver, filename, &csc, MPI_COMM_WORLD );
    free(filename);

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    iparm[IPARM_FACTORIZATION] = API_FACT_LU;
    pastix_task_order( pastix_data, &csc, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );
    pastix_task_blend( pastix_data );

    /**
     * Perform the numerical factorization
     */
    pastix_task_sopalin( pastix_data, &csc );

    /**
     * Generates the b and x vector such that A * x = b
     * Compute the norms of the initial vectors if checking purpose.
     */
    size = pastix_size_of( csc.flttype ) * csc.n;
    x = malloc( size );

    if ( check )
    {
        b = malloc( size );

        if ( check > 1 ) {
            x0 = malloc( size );
        } else {
            x0 = NULL;
        }
        spmGenRHS( PastixRhsRndX, nrhs, &csc, x0, csc.n, b, csc.n );
        memcpy( x, b, size );
    }
    else {
        spmGenRHS( PastixRhsRndB, nrhs, &csc, NULL, csc.n, x, csc.n );
        memcpy( b, x, size );
    }

    /**
     * Solve the linear system
     */
    pastix_task_solve( pastix_data, &csc, nrhs, x, csc.n );

    pastix_task_raff(pastix_data, x, nrhs, b);

    if ( check )
    {
        spmCheckAxb( nrhs, &csc, x0, csc.n, b, csc.n, x, csc.n );

        if (x0) free(x0);

        free(x); free(b);
    }
    spmExit( &csc );

    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    return EXIT_SUCCESS;
}
