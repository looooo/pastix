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
#include <limits.h>
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
    iparm[IPARM_IO_STRATEGY]   = API_IO_SAVE;
    iparm[IPARM_MIN_BLOCKSIZE] = 60;
    iparm[IPARM_MAX_BLOCKSIZE] = 120;

    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    split_level       = 0;
    stop_criteria     = INT_MAX;
    stop_when_fitting = 0;

    /**
     * Get options from command line
     */
    pastix_ex_getoptions( argc, argv,
                          iparm, dparm,
                          &driver, &filename );

    cscReadFromFile( driver, filename, &csc, MPI_COMM_WORLD );
    free(filename);

    printf("Reordering parameters are %d %d %d\n", split_level, stop_criteria, stop_when_fitting);
    if (stop_when_fitting != 0 && stop_when_fitting != 1){
        fprintf(stderr, "Fatal error in reordering parameters -R split_level:stop_criteria:stop_when_fitting\n");
        exit(1);
    }

    pastix_task_order( pastix_data, csc.n, csc.colptr, csc.rowptr, NULL, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );
    pastix_task_reordering( pastix_data, split_level, stop_criteria, stop_when_fitting );
    pastix_task_blend( pastix_data );

    /* pastix_task_sopalin( pastix_data, &csc ); */

    if (0)
    {
        double normA, normB, normX, normR, result, eps;
        double mdone = -1.;
        double done  = 1.;
        double dzero = 0.;
        int loud = 4;
        size_t size = pastix_size_of( csc.flttype ) * csc.gN;
        void *x0 = malloc( size );
        void *x  = malloc( size );
        void *b  = malloc( size );

        switch( csc.flttype ) {
        case PastixComplex64:
            CORE_zplrnt( csc.gN, 1, x0, 1, 1, 0, 0, 243759 );
            break;
        case PastixComplex32:
            CORE_cplrnt( csc.gN, 1, x0, 1, 1, 0, 0, 243759 );
            break;
        case PastixDouble:
            CORE_dplrnt( csc.gN, 1, x0, 1, 1, 0, 0, 243759 );
            break;
        case PastixFloat:
            CORE_splrnt( csc.gN, 1, x0, 1, 1, 0, 0, 243759 );
            break;
        default:
            ;
        }
        /* Copy b to x */
        //if (0)
        {
            pastix_int_t i;
            double *lb = (double*)x0;
            for(i=0; i<csc.gN; i++) {
                lb[i] = (double)i+1;
            }
        }
        spmMatVec( PastixNoTrans, &done, &csc, x0, &dzero, b );
        memcpy( x, b, size );

        {
            int PRHS_ii;
            fprintf(stdout,"%s (Proc %d) : ", "RHS", 0);
            for (PRHS_ii= 0; PRHS_ii<csc.gN; PRHS_ii++)
                fprintf(stdout,"%.3g ", ((double*)b)[PRHS_ii]);
            fprintf(stdout,"\n");
        }
        pastix_task_solve( pastix_data, &csc, 1, x, csc.gN );

        eps = LAPACKE_dlamch( 'e' );
        normB = LAPACKE_dlange( LAPACK_COL_MAJOR, 'I', csc.gN, 1, b, csc.gN );
        normX = LAPACKE_dlange( LAPACK_COL_MAJOR, 'I', csc.gN, 1, x, csc.gN );
        normA = dparm[DPARM_A_NORM];

        spmMatVec( PastixNoTrans, &mdone, &csc, x, &done, b );
        normR = LAPACKE_dlange( LAPACK_COL_MAJOR, 'I', csc.gN, 1, b, csc.gN );

        {
            int PRHS_ii;
            fprintf(stdout,"%s (Proc %d) : ", "SOL", 0);
            for (PRHS_ii= 0; PRHS_ii<csc.gN; PRHS_ii++)
                fprintf(stdout,"%.3g ",((double*)x)[PRHS_ii]);
            fprintf(stdout,"\n");
        }
        result = normR / ( ( normA * normX + normB ) * csc.gN * eps ) ;

        if ( loud > 2 ) {
            printf("============\n");
            printf("Checking the Residual of the solution \n");
            if ( loud > 3 )
                printf( "-- ||A||_oo = %e, ||x_0||_oo = %e, ||x||_oo = %e, ||b||_oo= %e, ||A x - b||_oo = %e\n",
                        normA,
                        LAPACKE_dlange( LAPACK_COL_MAJOR, 'I', csc.gN, 1, x0, csc.gN ),
                        normX, normB, normR );

            printf("-- ||Ax-b||_oo/((||A||_oo||x||_+||b||_oo).N.eps) = %e \n", result);
        }

        if (  isnan(normX) || isinf(normX) || isnan(result) || isinf(result) || (result > 60.0) ) {
            if( loud ) printf("-- Solution is suspicious ! \n");
        }
        else{
            if( loud ) printf("-- Solution is CORRECT ! \n");
        }

        core_dgeadd( PastixNoTrans, csc.gN, 1, -1., x0, 1, x, 1 );
        normR = LAPACKE_dlange( LAPACK_COL_MAJOR, 'M', csc.gN, 1, x, csc.gN );
        printf( "-- ||X0-X||_oo = %e\n", normR );

        free(x0); free(x); free(b);
    }
    spmExit( &csc );

    /* if (!PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_LOAD)) */
    /* { */
    /*     pastix_complex64_t *rhs     = NULL; */
    /*     char           *type    = NULL; */
    /*     char           *rhstype = NULL; */
    /*     pastix_int_t    mat_type; */

    /*     /\*******************************************\/ */
    /*     /\*      Read Matrice                       *\/ */
    /*     /\*******************************************\/ */
    /*     read_matrix(filename[0], &ncol, &colptr, &rowptr, &values, */
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
