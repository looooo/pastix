/* File: schur.c

 Construct the schur complement of the matrix.

 */

#include <pastix.h>
#include <spm.h>
#include "../matrix_drivers/drivers.h"
#include "../common/pastixdata.h"
#include "../order/order.h"

int main (int argc, char **argv)
{
    pastix_data_t   *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t     iparm[IPARM_SIZE]; /*< Integer in/out parameters for pastix                */
    double           dparm[DPARM_SIZE]; /*< Floating in/out parameters for pastix               */
    char            *filename;
    int              nrhs        = 1;
    pastix_spm_t    *spm, *spm2;
    double           normA;
    pastix_driver_t  driver;
    void            *x, *x0, *b;
    size_t           size;
    int              check       = 2;
    int              ret         = PASTIX_SUCCESS;

    /**
     * Initialize parameters to default values
     */
    iparm[IPARM_MODIFY_PARAMETER] = API_NO;
    pastix( &pastix_data, MPI_COMM_WORLD, -1, NULL, NULL, NULL,
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
        b = malloc( size );
        memcpy( b, x, size );
    }

    /**
     * Initialize pastix
     */
    iparm[IPARM_START_TASK] = API_TASK_INIT;
    iparm[IPARM_END_TASK]   = API_TASK_INIT;
    pastix( &pastix_data, MPI_COMM_WORLD,
            -1, NULL, NULL, NULL,
            NULL, NULL, NULL, 1, iparm, dparm );

    /**
     * Initialize the schur
     */
    int baseval=spmFindBase(spm);
    pastix_data->schur_n = spm->gN/3;
    pastix_data->schur_list = (pastix_int_t*)malloc(pastix_data->schur_n*sizeof(pastix_int_t));
    for (int i=0; i<pastix_data->schur_n; i++)
        pastix_data->schur_list[i]=i+baseval;
    iparm[IPARM_SCHUR]=API_YES;


    /**
     * Call pastix
     */
    iparm[IPARM_START_TASK] = API_TASK_ORDERING;
    iparm[IPARM_END_TASK]   = API_TASK_CLEAN;
    pastix(&pastix_data, MPI_COMM_WORLD,
           spm->n, spm->colptr, spm->rowptr, spm->values,
           NULL, NULL, x, nrhs, iparm, dparm );

    if ( check ) {
        spmCheckAxb( nrhs, spm, x0, spm->n, b, spm->n, x, spm->n );
        if (x0) free(x0);
    }

    free(x);
    free(b);
    free(spm);
    return EXIT_SUCCESS;
}
