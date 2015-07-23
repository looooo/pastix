/**
 *  @file: analyze.c
 *
 *  A simple example that performs only the analyses steps onto the given graph.
 *  These tests doesn't require the values of the matrix.
 *
 */
#include <pastix.h>
#include "../matrix_drivers/drivers.h"
#include <csc.h>

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix           */
    pastix_int_t    iparm[IPARM_SIZE];  /* integer parameters for pastix                             */
    double          dparm[DPARM_SIZE];  /* floating parameters for pastix                            */
#if defined(PASTIX_WITH_MPI)
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
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Get options from command line
     */
    pastix_ex_getoptions( argc, argv,
                          iparm, dparm,
                          &driver, &filename );

    cscReadFromFile( driver, filename, &csc, MPI_COMM_WORLD );
    free(filename);

    pastix_task_order( pastix_data, &csc, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );
    pastix_task_blend( pastix_data );

    spmExit( &csc );

    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

#if defined(PASTIX_WITH_MPI)
    MPI_Finalize();
#endif

    return EXIT_SUCCESS;
}
