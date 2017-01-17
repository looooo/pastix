/* File: simple.c
 *
 *  A simple example :
 *  read the matrix, check it is correct and correct it if needed,
 *  then run pastix in one call.
 *
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#include <assert.h>

/* to access functions from the libpastix, respect this order */
#include <pastix.h>
#include <spm.h>
#include "../matrix_drivers/drivers.h"

int main (int argc, char **argv)
{
    pastix_data_t   *pastix_data = NULL; /* Pointer to a storage structure needed by pastix           */
    pastix_float_t  *b           = NULL; /* right hand side                                           */
    pastix_int_t     iparm[IPARM_SIZE]; /* integer parameters for pastix                             */
    double           dparm[DPARM_SIZE]; /* floating parameters for pastix                            */
    char            *filename;  /* Filename(s) given by user                                 */
    int              nrhs        = 1;
    pastix_spm_t    *spm, *spm2;
    double           normA;
    pastix_driver_t  driver;
    void            *x, *x0;
    size_t           size;
    int              check       = 2;
    int              ret         = PASTIX_SUCCESS;

    /**
     * Initialize parameters to default values
     */
    iparm[IPARM_MODIFY_PARAMETER] = API_NO;
    pastix( &pastix_data, MPI_COMM_WORLD,
            -1, NULL, NULL, NULL,
            NULL, NULL, NULL, 1, iparm, dparm );

    /**
     * Update options from command line, and get the matrix filename
     */
    pastix_ex_getoptions( argc, argv,
                          iparm, dparm,
                          &driver, &filename );

    /**
     * Read Matrice
     */
    spm = malloc( sizeof( pastix_spm_t ) );
    spmReadDriver( driver, filename, spm, MPI_COMM_WORLD );
    free(filename);

    /**
     * Check Matrix format
     */
    spm2 = spmCheckAndCorrect( spm );
    if ( spm2 != spm ) {
        spmExit( spm );
        free(spm);
        spm = spm2;
    }

    /**
     * Scale the matrix to avoid unexpected rouding errors
     */
    normA = spmNorm( PastixFrobeniusNorm, spm );
    spmScal( 1./normA, spm );

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
        memcpy( b, x, size );
    }

    /**
     * Call pastix
     */
    iparm[IPARM_START_TASK] = API_TASK_INIT;
    iparm[IPARM_END_TASK]   = API_TASK_CLEAN;
    pastix( &pastix_data, MPI_COMM_WORLD,
            spm->n, spm->colptr, spm->rowptr, spm->values,
            NULL, NULL, x, nrhs, iparm, dparm );
    if (iparm[IPARM_ERROR_NUMBER] != PASTIX_SUCCESS)
        return ret;

    /**
     * Check the solution
     */
    if ( check )
    {
        spmCheckAxb( nrhs, spm, x0, spm->n, b, spm->n, x, spm->n );
        if (x0) free(x0);
    }

    spmExit( spm );
    free(b);
    free(x);
    return EXIT_SUCCESS;
}
