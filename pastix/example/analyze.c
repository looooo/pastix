/**
 *  @file: analyze.c
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
#include <pastix.h>
#include "../matrix_drivers/drivers.h"
#include <csc.h>

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
    void           *rhs = NULL;

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

    cscReadFromFile( driver, filename, &csc, &rhs, MPI_COMM_WORLD );
    free(filename);

    pastix_task_order( pastix_data, csc.n, csc.colptr, csc.rows, NULL, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );
    pastix_task_blend( pastix_data );
    //pastix_task_sopalin( pastix_data, &csc );

    //cscExit( csc );
    free(csc.colptr);
    free(csc.rows);
    free(csc.avals);
    free(rhs);

    /* if (!PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD)) */
    /* { */
    /*     pastix_complex64_t *rhs     = NULL; */
    /*     char           *type    = NULL; */
    /*     char           *rhstype = NULL; */
    /*     pastix_int_t    mat_type; */

    /*     /\*******************************************\/ */
    /*     /\*      Read Matrice                       *\/ */
    /*     /\*******************************************\/ */
    /*     read_matrix(filename[0], &ncol, &colptr, &rows, &values, */
    /*                   &rhs, &type, &rhstype, driver_type[0], MPI_COMM_WORLD); */

    /*     free(filename[0]); */
    /*     free(filename); */
    /*     free(driver_type); */
    /*     free(rhs); */
    /*     free(rhstype); */

    /*     mat_type = API_SYM_NO; */
    /*     if (MTX_ISSYM(type)) mat_type = API_SYM_YES; */
    /*     if (MTX_ISHER(type)) mat_type = API_SYM_HER; */
    /*     iparm[IPARM_SYM] = mat_type; */
    /*     switch (mat_type) */
    /*     { */
    /*     case API_SYM_YES: */
    /*         iparm[IPARM_FACTORIZATION] = API_FACT_LDLT; */
    /*     break; */
    /*     case API_SYM_HER: */
    /*         iparm[IPARM_FACTORIZATION] = API_FACT_LDLH; */
    /*         break; */
    /*     default: */
    /*         iparm[IPARM_FACTORIZATION] = API_FACT_LU; */
    /*     } */

    /*     free(type); */
    /* } */
    /* else */
    /* { */
    /*     iparm[IPARM_START_TASK] = API_TASK_SYMBFACT; */
    /* } */

    /* iparm[IPARM_END_TASK]            = API_TASK_ANALYSE; */
    /* pastix(&pastix_data, MPI_COMM_WORLD, */
    /*          ncol, colptr, rows, values, */
    /*          NULL, NULL, NULL, 0, iparm, dparm); */

    /* /\* Clean *\/ */
    /* iparm[IPARM_START_TASK] = API_TASK_CLEAN; */
    /* iparm[IPARM_END_TASK]   = API_TASK_CLEAN; */
    /* pastix(&pastix_data, MPI_COMM_WORLD, */
    /*        0, NULL, NULL, NULL, */
    /*        NULL, NULL, NULL, 0, iparm, dparm); */

    /* if (colptr != NULL) free(colptr); */
    /* if (rows   != NULL) free(rows); */
    /* if (values != NULL) free(values); */

    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );
#if defined(PASTIX_WITH_MPI)
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
