/* File: schur.c

 Construct the schur complement of the matrix.

 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <complex.h>

/* to access functions from the libpastix, respect this order */
#include <pastix.h>
#include <spm.h>
#include "../matrix_drivers/drivers.h"
#include "../common/pastixdata.h"
#include "../order/order.h"
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
    /* pastix_data_t  *pastix_data = NULL; /\* Pointer to a storage structure needed by pastix           *\/ */
    /* pastix_int_t    ncol;               /\* Size of the matrix                                        *\/ */
    /* pastix_int_t   *colptr      = NULL; /\* Indexes of first element of each column in row and values *\/ */
    /* pastix_int_t   *rows        = NULL; /\* Row of each element of the matrix                         *\/ */
    /* pastix_float_t *values      = NULL; /\* Value of each element of the matrix                       *\/ */
    /* pastix_float_t *rhs         = NULL; /\* right hand side                                           *\/ */
    /* pastix_float_t *rhssaved    = NULL; /\* right hand side (save)                                    *\/ */
    /* pastix_float_t *ax          = NULL; /\* A times X product                                         *\/ */
    /* pastix_int_t    iparm[IPARM_SIZE];  /\* integer parameters for pastix                             *\/ */
    /* double          dparm[DPARM_SIZE];  /\* floating parameters for pastix                            *\/ */
    /* pastix_int_t   *perm        = NULL; /\* Permutation tabular                                       *\/ */
    /* pastix_int_t   *invp        = NULL; /\* Reverse permutation tabular                               *\/ */
    /* char           *type        = NULL; /\* type of the matrix                                        *\/ */
    /* char           *rhstype     = NULL; /\* type of the right hand side                               *\/ */
    /* driver_type_t  *driver_type;        /\* Matrix driver(s) requested by user                        *\/ */
    /* char          **filename;           /\* Filename(s) given by user                                 *\/ */
    /* int             nbmatrices;         /\* Number of matrices given by user                          *\/ */
    /* int             nbthread;           /\* Number of thread wanted by user                           *\/ */
    /* int             verbosemode;        /\* Level of verbose mode (0, 1, 2)                           *\/ */
    /* int             ordering;           /\* Ordering to use                                           *\/ */
    /* int             incomplete;         /\* Indicate if we want to use incomplete factorisation       *\/ */
    /* int             level_of_fill;      /\* Level of fill for incomplete factorisation                *\/ */
    /* int             amalgamation;       /\* Level of amalgamation for Kass                            *\/ */
    /* int             ooc;                /\* OOC limit (Mo/percent depending on compilation options)   *\/ */
    /* pastix_int_t    mat_type; */
    /* long            i,j; */
    /* pastix_float_t *schur; */
    /* pastix_int_t    nschur; */
    /* pastix_int_t   *schurlist = NULL; */
    /* double norme1, norme2; */


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
     * Initialize pastix
     */
    iparm[IPARM_START_TASK] = API_TASK_INIT;
    iparm[IPARM_END_TASK]   = API_TASK_INIT;
    pastix( &pastix_data, MPI_COMM_WORLD,
            -1, NULL, NULL, NULL,
            NULL, NULL, NULL, 1, iparm, dparm );

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

    iparm[IPARM_END_TASK] = API_TASK_ORDERING;
    pastix(&pastix_data, MPI_COMM_WORLD,
           spm->n, spm->colptr, spm->rowptr, spm->values,
           NULL, NULL, NULL, nrhs, iparm, dparm );

    /**
     * Scale the matrix to avoid unexpected rouding errors
     */
    normA = spmNorm( PastixFrobeniusNorm, spm );
    spmScal( 1./normA, spm );

    int baseval=spmFindBase(spm);
    pastix_data->schur_n = spm->gN/3;
    pastix_data->schur_list = (pastix_int_t*)malloc(pastix_data->schur_n*sizeof(pastix_int_t));
    for (int i=0; i<pastix_data->schur_n; i++)
        pastix_data->schur_list[i]=i+baseval;
    iparm[IPARM_SCHUR]=API_YES;

    /*******************************************/
    /* Set schur unknowns list                 */
    /*******************************************/
    iparm[IPARM_END_TASK] = API_TASK_SYMBFACT;
    pastix(&pastix_data, MPI_COMM_WORLD,
           spm->n, spm->colptr, spm->rowptr, spm->values,
           NULL, NULL, NULL, nrhs, iparm, dparm );
    orderSave(pastix_data->ordemesh,"ordername.txt");

    /*******************************************/
    /*           Call pastix                   */
    /*******************************************/
    iparm[IPARM_START_TASK]          = API_TASK_ANALYSE;
    iparm[IPARM_END_TASK]            = API_TASK_NUMFACT;
    pastix(&pastix_data, MPI_COMM_WORLD,
           spm->n, spm->colptr, spm->rowptr, spm->values,
           NULL, NULL, NULL, nrhs, iparm, dparm );


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

    /*******************************************/
    /*               Solve                     */
    /*******************************************/
    iparm[IPARM_START_TASK]          = API_TASK_SOLVE;
    iparm[IPARM_END_TASK]            = API_TASK_REFINE;

    pastix(&pastix_data, MPI_COMM_WORLD,
           spm->n, spm->colptr, spm->rowptr, spm->values,
           NULL, NULL, x, nrhs, iparm, dparm );

    if ( check ) {
        spmCheckAxb( nrhs, spm, x0, spm->n, b, spm->n, x, spm->n );
        if (x0) free(x0);
        free(x); free(b);
    }


    /*******************************************/
    /*               Clean                     */
    /*******************************************/

    iparm[IPARM_START_TASK]          = API_TASK_CLEAN;
    iparm[IPARM_END_TASK]            = API_TASK_CLEAN;

    pastix(&pastix_data, MPI_COMM_WORLD,
           spm->n, spm->colptr, spm->rowptr, spm->values,
           NULL, NULL, x, nrhs, iparm, dparm );


    //    fprintf(stdout, "Precision : ||ax -b||/||b|| = %.20lg\n", sqrt(norme1/norme2));
    /* free(ax); */
    /* free(colptr); */
    /* free(rows); */
    /* free(values); */
    /* free(perm); */
    /* free(invp); */
    /* free(rhs); */
    /* free(type); */
    /* free(rhstype); */
    /* free(schurlist); */
    /*******************************************/
    /* Do some stuff with the schur complement */
    /*******************************************/
    fprintf(stdout, "Now we have to do some stuff with our schur complement.\n");
    /* for (i = 0; i < nschur; i++) */
    /*   { */
    /*     for (j =0; j < nschur; j++) */
    /*  fprintf (stdout, "\t%lg", (double)schur[j*nschur+i]); */
    /*     fprintf(stdout, "\n"); */
    /*   } */

    /* free(schur); */


    return EXIT_SUCCESS;
}
