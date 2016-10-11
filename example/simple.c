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
#include "../order/order.h"
#include "../common/common.h"
#include "../blend/solver.h"

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

    printf("\tH-PaStiX parameters are\n");
    printf("\tSPLITSYMBOL %ld %ld\n", iparm[IPARM_MIN_BLOCKSIZE], iparm[IPARM_MAX_BLOCKSIZE]);
    printf("\tCOMPRESS_SIZE %ld\n", iparm[IPARM_COMPRESS_SIZE]);
    printf("\tTOLERANCE %.3g\n", dparm[DPARM_COMPRESS_TOLERANCE]);

    printf("\tCOMPRESS_WHEN %ld COMPRESS_METHOD %ld\n",
           iparm[IPARM_COMPRESS_WHEN],
           iparm[IPARM_COMPRESS_METHOD]);

    /* TO BE CLEAN !!! */
    compress_when      = iparm[IPARM_COMPRESS_WHEN];
    compress_method    = iparm[IPARM_COMPRESS_METHOD];
    compress_tolerance = dparm[DPARM_COMPRESS_TOLERANCE];

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


    /* Run another symbolic factorization to eliminate zero blocks created by splitting */
    if (0)
    {
        pastix_int_t new_cblknbr = pastix_data->symbmtx->cblknbr;
        pastix_int_t i;

        pastix_data->ordemesh->cblknbr = new_cblknbr;
        pastix_data->ordemesh->rangtab = malloc((new_cblknbr+1) * sizeof(pastix_int_t));
        pastix_data->ordemesh->treetab = malloc((new_cblknbr+1) * sizeof(pastix_int_t));

        SymbolCblk *cblk = pastix_data->symbmtx->cblktab;
        for (i=0; i<new_cblknbr; i++, cblk++){
            pastix_data->ordemesh->rangtab[i] = cblk->fcolnum;
            pastix_data->ordemesh->treetab[i] = i+1;
        }
        pastix_data->ordemesh->rangtab[new_cblknbr] = spm->n;
        pastix_data->ordemesh->treetab[new_cblknbr-1] = -1;

        symbolExit(pastix_data->symbmtx);
        memFree_null(pastix_data->symbmtx);
        pastix_data->symbmtx = NULL;

        /* The symbolic structure was already split, avoid another splitting */
        iparm[IPARM_MIN_BLOCKSIZE] = 100000000000;
        iparm[IPARM_MAX_BLOCKSIZE] = 100000000000;

        pastix_task_symbfact( pastix_data, NULL, NULL );
        pastix_task_reordering( pastix_data );
        pastix_task_blend( pastix_data );
    }

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

    if ( check )
    {
        b = malloc( size );

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

        free(x); free(b);
    }
    spmExit( spm );
    free( spm );
    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    return EXIT_SUCCESS;
}
