/**
 *
 * @file bcsc_matvec_test.c
 *
 * Tests and validate the bcsc_matvec routines.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 **/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <pastix.h>
#include "../matrix_drivers/drivers.h"
#include "common.h"
#include <csc.h>
#include <bcsc.h>
#include "sopalin_data.h"

int z_bcsc_matvec_check( int trans, const pastix_csc_t *spm, const pastix_data_t *pastix_data );
int c_bcsc_matvec_check( int trans, const pastix_csc_t *spm, const pastix_data_t *pastix_data );
int d_bcsc_matvec_check( int trans, const pastix_csc_t *spm, const pastix_data_t *pastix_data );
int s_bcsc_matvec_check( int trans, const pastix_csc_t *spm, const pastix_data_t *pastix_data );

#define PRINT_RES(_ret_)                        \
    if(_ret_) {                                 \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* mtxnames[] = { "General", "Symmetric", "Hermitian" };

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix  */
    pastix_int_t    iparm[IPARM_SIZE];  /* integer parameters for pastix                    */
    double          dparm[DPARM_SIZE];  /* floating parameters for pastix                   */
    pastix_driver_t driver;             /* Matrix driver(s) requested by user               */
    pastix_csc_t   *csc, *csc2;
    char *filename;                     /* Filename(s) given by user                        */
    int mtxtype;
    int ret = PASTIX_SUCCESS;
    int err = 0;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Get options from command line
     */
    pastix_ex_getoptions( argc, argv,
                          NULL, NULL,
                          &driver, &filename );

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
     * Run preprocessing steps required to generate the blocked csc
     */
    pastix_task_order( pastix_data, csc, NULL, NULL );
    pastix_task_symbfact( pastix_data, NULL, NULL );
    pastix_task_blend( pastix_data );

    /**
     * Generate the blocked csc
     */
    pastix_data->bcsc = malloc( sizeof(pastix_bcsc_t) );
    bcscInit( csc,
              pastix_data->ordemesh,
              pastix_data->solvmatr,
              csc->mtxtype == PastixGeneral,
              pastix_data->bcsc );

    printf(" -- BCSC MatVec Test --\n");
    printf(" Datatype: %s\n", fltnames[csc->flttype] );
    spmBase( csc, 0 );

    for( mtxtype=PastixGeneral; mtxtype<=PastixHermitian; mtxtype++ )
    {
        if ( ((mtxtype == PastixGeneral) && (csc->mtxtype != PastixGeneral)) ||
             ((mtxtype != PastixGeneral) && (csc->mtxtype == PastixGeneral)) )
        {
            continue;
        }
        if ( (mtxtype == PastixHermitian) &&
                ((csc->flttype != PastixComplex64) && (csc->flttype != PastixComplex32)) )
        {
            continue;
        }
        csc->mtxtype  = mtxtype;
        pastix_data->bcsc->mtxtype = mtxtype;

        printf("   Matrix type : %s\n", mtxnames[mtxtype - PastixGeneral] );
        printf("   -- Test Matrix * Vector : ");

        switch( csc->flttype ){
        case PastixComplex64:
            ret = z_bcsc_matvec_check( PastixNoTrans, csc, pastix_data );
            break;

        case PastixComplex32:
            ret = c_bcsc_matvec_check( PastixNoTrans, csc, pastix_data );
            break;

        case PastixFloat:
            ret = s_bcsc_matvec_check( PastixNoTrans, csc, pastix_data );
            break;

        case PastixDouble:
        default:
            ret = d_bcsc_matvec_check( PastixNoTrans, csc, pastix_data );
        }
        PRINT_RES(ret);
    }

    spmExit( csc );
    free( csc );

    pastixFinalize( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    if( err == 0 ) {
        printf(" -- All tests PASSED --\n");
        return EXIT_SUCCESS;
    }
    else
    {
        printf(" -- %d tests FAILED --\n", err);
        return EXIT_FAILURE;
    }
}
