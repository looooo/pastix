/**
 *  @file: simple.c
 *
 *  A simple example :
 *  read the matrix, check it is correct and correct it if needed,
 *  then run pastix in one call.
 *
 */
#include <pastix.h>
#include <csc.h>
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
    pastix_csc_t   *csc, *csc2;
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
    char *lvls = getenv("LEVELS");
    if (lvls == NULL){
        printf("Set LEVELS variable\n");
        exit(1);
    }
    pastix_int_t levels = atoi(lvls);

    char *tol  = getenv("TOLERANCE");
    if (tol == NULL){
        printf("Set TOLERANCE variable\n");
        exit(1);
    }
    double tolerance = atof(tol);

    char *splt = getenv("SPLITSIZE");
    if (splt == NULL){
        printf("Set SPLITSIZE variable\n");
        exit(1);
    }
    pastix_int_t splitsize = atoi(splt);

    char *thr  = getenv("THRESHOLD");
    if (thr == NULL){
        printf("Set THRESHOLD variable\n");
        exit(1);
    }
    pastix_int_t threshold = atoi(thr);
    if (threshold > splitsize){
        printf("Threshold (%ld) should be lower than splitsize (%ld)\n",
               threshold, splitsize);
        exit(1);
    }

    char *tree = getenv("HODLRTREE");
    if (tree == NULL){
        printf("Set HODLRTREE (0- unactive, 1- active)\n");
        exit(1);
    }
    pastix_int_t hodlrtree = atoi(tree);
    if (hodlrtree != 0 && hodlrtree != 1){
        printf("Set correctly HODLRTREE (0- unactive, 1- active)\n");
        exit(1);
    }

    char *ordo = getenv("REORDERING");
    if (ordo == NULL){
        printf("Set REORDERING (0- unactive, 1- active)\n");
        exit(1);
    }
    pastix_int_t reordering = atoi(ordo);
    if (reordering != 0 && reordering != 1){
        printf("Set correctly REORDERING (0- unactive, 1- active)\n");
        exit(1);
    }

    char *comp = getenv("COMPRESS");
    if (comp == NULL){
        printf("Set COMPRESS (0- unactive, 1- active)\n");
        exit(1);
    }
    pastix_int_t compress = atoi(comp);
    if (compress != 0 && compress != 1){
        printf("Set correctly COMPRESS (0- unactive, 1- active)\n");
        exit(1);
    }

    char *face = getenv("FACING");
    if (face == NULL){
        printf("Set FACING (0- unactive, 1- active)\n");
        exit(1);
    }
    pastix_int_t facing = atoi(face);
    if (facing != 0 && facing != 1){
        printf("Set correctly FACING (0- unactive, 1- active)\n");
        exit(1);
    }

    printf("\tH-PaStiX parameters are\n");
    printf("\tSPLITSIZE %ld THRESHOLD %ld\n", splitsize, threshold);
    printf("\tTOLERANCE %.3g\n", tolerance);
    printf("\tHODLRTREE %ld (on last supernode right now)\n", hodlrtree);
    printf("\tREORDERING %ld\n", reordering);
    printf("\tCOMPRESS %ld\n", compress);
    printf("\tFACING %ld\n", facing);
    printf("\tLEVELS %ld (unused right now)\n", levels);

    iparm[IPARM_MIN_BLOCKSIZE] = 100;
    iparm[IPARM_MAX_BLOCKSIZE] = 150;

    printf("\tSPLITSYMBOL %ld %ld\n", iparm[IPARM_MIN_BLOCKSIZE], iparm[IPARM_MAX_BLOCKSIZE]);

    iparm[IPARM_ITERMAX]          = 10;
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


    current_cblk  = 0;
    total_memory  = 0.;
    total_memory2 = 0.;

    gain_L = 0 ;
    gain_D = 0 ;
    gain_U = 0 ;

    /**
     * Read the sparse matrix with the driver
     */
    csc = malloc( sizeof( pastix_csc_t ) );
    cscReadFromFile( driver, filename, csc, MPI_COMM_WORLD );
    free(filename);
    csc2 = spmCheckAndCorrect( csc );
    if ( csc2 != csc ) {
        spmExit( csc );
        free(csc);
        csc = csc2;
    }

    /**
     * Perform ordering, symbolic factorization, and analyze steps
     */
    pastix_task_order( pastix_data, csc, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );

    if (reordering == 1)
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
    pastix_task_sopalin( pastix_data, csc );

    printf("Total memory of the solver %10f Mo (symmetric) %10f (unsymmetric)\n", total_memory, total_memory2);
    printf("Gain_L %10f Mo Gain_D %10f Mo, Gain_U %10f Mo\n", gain_L, gain_D, gain_U);

    /**
     * Generates the b and x vector such that A * x = b
     * Compute the norms of the initial vectors if checking purpose.
     */
    size = pastix_size_of( csc->flttype ) * csc->n;
    x = malloc( size );

    if ( check )
    {
        b = malloc( size );

        if ( check > 1 ) {
            x0 = malloc( size );
        } else {
            x0 = NULL;
        }
        spmGenRHS( PastixRhsRndX, nrhs, csc, x0, csc->n, b, csc->n );
        memcpy( x, b, size );
    }
    else {
        spmGenRHS( PastixRhsRndB, nrhs, csc, NULL, csc->n, x, csc->n );

        /* Save b for refinement: TODO: make 2 examples w/ or w/o refinement */
        b = malloc( size );
        memcpy( b, x, size );
    }

    /**
     * Solve the linear system
     */
    pastix_task_solve( pastix_data, csc, nrhs, x, csc->n );

    pastix_task_raff(pastix_data, x, nrhs, b);

    if ( check )
    {
        spmCheckAxb( nrhs, csc, x0, csc->n, b, csc->n, x, csc->n );

        if (x0) free(x0);

        free(x); free(b);
    }
    spmExit( csc );
    free( csc );
    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    return EXIT_SUCCESS;
}
