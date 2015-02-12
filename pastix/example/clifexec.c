/**
 *  @file: clifexec.c
 *
 *  A simple example :
 *  read the matrix, check it is correct and correct it if needed,
 *  then run pastix in one call.
 *
 *  @precisions normal z => s, d, c
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "pastix.h"
#include "../matrix_drivers/drivers.h"
#include "../symbol/symbol.h"
#include <scotch.h>

#ifdef FORCE_NOMPI
#define MPI_COMM_WORLD 0
#endif

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix           */
    pastix_int_t    iparm[IPARM_SIZE];  /* integer parameters for pastix                             */
    double          dparm[DPARM_SIZE];  /* floating parameters for pastix                            */
#ifndef FORCE_NOMPI
    int             required;           /* MPI thread level required                                 */
    int             provided;           /* MPI thread level provided                                 */
#endif
    int             mpid;
    pastix_driver_t driver;        /* Matrix driver(s) requested by user                        */
    char           *filename;           /* Filename(s) given by user                                 */
    pastix_csc_t    csc;

    /*******************************************/
    /*          MPI initialisation             */
    /*******************************************/
#if defined(PASTIX_WITH_MPI)
    required = MPI_THREAD_MULTIPLE;
    provided = -1;
    MPI_Init_thread(&argc, &argv, required, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpid);
    if (mpid == 0)
    {
        switch (provided)
        {
        case MPI_THREAD_SINGLE:
            printf("MPI_Init_thread level = MPI_THREAD_SINGLE\n");
            break;
        case MPI_THREAD_FUNNELED:
            printf("MPI_Init_thread level = MPI_THREAD_FUNNELED\n");
            break;
        case MPI_THREAD_SERIALIZED:
            printf("MPI_Init_thread level = MPI_THREAD_SERIALIZED\n");
            break;
        case MPI_THREAD_MULTIPLE:
            printf("MPI_Init_thread level = MPI_THREAD_MULTIPLE\n");
            break;
        default:
            printf("MPI_Init_thread level = ???\n");
        }
    }
#else
    mpid = 0;
#endif
    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Get options from command line
     */
    pastix_ex_getoptions( argc, argv,
                          iparm, dparm,
                          &driver, &filename );

    cscReadFromFile( driver, filename, &csc, MPI_COMM_WORLD );
    free(filename);

    orderingGregoire = 0;
    pastix_task_order( pastix_data, csc.n, csc.colptr, csc.rows, NULL, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );

    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    orderingGregoire = 1;
    pastix_task_order( pastix_data, csc.n, csc.colptr, csc.rows, NULL, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );

    //pastix_task_blend( pastix_data );
    //pastix_task_sopalin( pastix_data, &csc );

    //cscExit( csc );
    free(csc.colptr);
    free(csc.rows);
    //free(csc.avals);

    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );
#ifndef FORCE_NOMPI
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
