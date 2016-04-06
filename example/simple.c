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
     * Check if parameters are well set
     */
    char *tol  = getenv("TOLERANCE");
    if (tol == NULL){
        printf("Set TOLERANCE variable\n");
        exit(1);
    }
    double tolerance = atof(tol);

    char *surf = getenv("SURFACE");
    if (surf == NULL){
        printf("Set SURFACE\n");
        exit(1);
    }
    pastix_int_t surface = atoi(surf);
    if (surface < 0){
        printf("Set correctly SURFACE (integer >= 0)\n");
        exit(1);
    }

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

    printf("\tH-PaStiX parameters are\n");
    printf("\tSPLITSYMBOL %ld %ld\n", iparm[IPARM_MIN_BLOCKSIZE], iparm[IPARM_MAX_BLOCKSIZE]);
    printf("\tCOMPRESS_SIZE %ld\n", iparm[IPARM_COMPRESS_SIZE]);
    printf("\tTOLERANCE %.3g\n", tolerance);
    printf("\t SURFACE %ld\n", surface);

    current_cblk  = 0;
    total_memory  = 0.;
    total_memory2 = 0.;

    gain_L = 0 ;
    gain_D = 0 ;
    gain_U = 0 ;
    tot_surface = 0.;

    time_comp    = 0.0;
    time_uncomp  = 0.0;
    time_recomp  = 0.0;
    time_fact    = 0.0;
    time_trsm    = 0.0;
    time_update  = 0.0;

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
     * Perform reordering after splitting
     */
    /* pastix_task_symbfact( pastix_data, NULL, NULL ); */
    /* pastix_task_reordering( pastix_data ); */
    /* pastix_task_blend( pastix_data ); */

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     * Second run to build new split cblk
     */
    /* rangtab_new     = malloc(csc->n * sizeof(pastix_int_t)); */
    /* rangtab_current = 0; */
    /* memset(rangtab_new, 0, csc->n * sizeof(pastix_int_t)); */

    /* Perform ordering, symbolic factorization, and analyze steps */
    /* rangtab_new[rangtab_current] = csc->n; */

    /* permtab_saved = malloc(csc->n*sizeof(pastix_int_t)); */
    /* peritab_saved = malloc(csc->n*sizeof(pastix_int_t)); */
    /* memcpy(permtab_saved, pastix_data->ordemesh->permtab, csc->n*sizeof(pastix_int_t)); */
    /* memcpy(peritab_saved, pastix_data->ordemesh->peritab, csc->n*sizeof(pastix_int_t)); */

    /* pastix_task_order( pastix_data, csc, NULL, NULL ); */
    /* pastix_task_symbfact( pastix_data, NULL, NULL ); */

    /* if (reordering == 1) */
    /*     pastix_task_reordering( pastix_data ); */

    /* rangtab_current = 0; */
    /* pastix_task_blend( pastix_data ); */


    /**
     * Perform the numerical factorization
     */
    pastix_task_sopalin( pastix_data, spm );

    printf("Total memory of the solver %10f Mo (symmetric) %10f (unsymmetric)\n", total_memory, total_memory2);
    printf("Gain_L %10f Mo Gain_D %10f Mo, Gain_U %10f Mo\n", gain_L, gain_D, gain_U);
    printf("Extra memory %.3g Mo\n", tot_surface);

    printf("Time compression   %10.3g s\n", time_comp  );
    printf("Time uncompression %10.3g s\n", time_uncomp);
    printf("Time recompression %10.3g s\n", time_recomp);
    printf("Time factorization %10.3g s\n", time_fact  );
    printf("Time trsm panel    %10.3g s\n", time_trsm  );
    printf("Time update        %10.3g s\n", time_update);
    printf("\n");

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
