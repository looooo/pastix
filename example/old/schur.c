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
#include <pastix/old_api.h>
#include <spm.h>

int main (int argc, char **argv)
{
    pastix_data_t   *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t     iparm[IPARM_SIZE]; /*< Integer in/out parameters for pastix                */
    double           dparm[DPARM_SIZE]; /*< Floating in/out parameters for pastix               */
    char            *filename;
    int              nrhs = 1;
    pastix_spm_t    *spm, *spm2;
    pastix_driver_t  driver;
    void            *x, *x0, *b;
    size_t           size;
    int              check = 1;

    /**
     * Initialize parameters to default values
     */
    iparm[IPARM_MODIFY_PARAMETER] = API_NO;
    pastix( &pastix_data, MPI_COMM_WORLD, -1, NULL, NULL, NULL,
            NULL, NULL, NULL, 1, iparm, dparm );

    /**
     * Update options from command line, and get the matrix filename
     */
    pastix_getOptions( argc, argv,
                       iparm, dparm,
                       &check, &driver, &filename );

    /**
     * Read Matrice
     */
    spm = malloc( sizeof( pastix_spm_t ) );
    spmReadDriver( driver, filename, spm, MPI_COMM_WORLD );
    free(filename);

    spmPrintInfo( spm, stdout );

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
     * Generates the b and x vector such that A * x = b
     * Compute the norms of the initial vectors if checking purpose.
     */
    size = pastix_size_of( spm->flttype ) * spm->n;
    x = malloc( size );
    b = malloc( size );

    x0 = NULL;
    if ( check )
    {
        if ( check > 1 ) {
            x0 = malloc( size );
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
        }
    }

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
