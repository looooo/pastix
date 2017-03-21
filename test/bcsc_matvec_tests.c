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
#include <spm.h>
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
char* transnames[] = { "NoTrans", "Trans", "ConjTrans" };
char* mtxnames[] = { "General", "Symmetric", "Hermitian" };

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix  */
    pastix_int_t    iparm[IPARM_SIZE];  /* integer parameters for pastix                    */
    double          dparm[DPARM_SIZE];  /* floating parameters for pastix                   */
    pastix_driver_t driver;             /* Matrix driver(s) requested by user               */
    pastix_spm_t   *spm, *spm2;
    char *filename;                     /* Filename(s) given by user                        */
    int t;
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

    spm = malloc( sizeof( pastix_spm_t ) );
    spmReadDriver( driver, filename, spm, MPI_COMM_WORLD );
    free(filename);
    spm2 = spmCheckAndCorrect( spm );
    if ( spm2 != spm ) {
        spmExit( spm );
        free(spm);
        spm = spm2;
    }
    spmBase( spm, 0 );

    /**
     * Run preprocessing steps required to generate the blocked csc
     */
    pastix_task_analyze( pastix_data, spm );

    /**
     * Generate the blocked csc
     */
    pastix_data->bcsc = malloc( sizeof(pastix_bcsc_t) );
    bcscInit( spm,
              pastix_data->ordemesh,
              pastix_data->solvmatr,
              spm->mtxtype == PastixGeneral,
              pastix_data->bcsc );

    printf(" -- BCSC MatVec Test --\n");
    for( t=PastixNoTrans; t<=PastixConjTrans; t++ )
    {
        if ( (t == PastixConjTrans) &&
             ((spm->flttype != PastixComplex64) && (spm->flttype != PastixComplex32)) )
        {
            continue;
        }
        if ( (spm->mtxtype != PastixGeneral) && (t != PastixNoTrans) )
        {
            continue;
        }
        printf("   Case %s - %s - %s:\n",
               fltnames[spm->flttype],
               mtxnames[spm->mtxtype - PastixGeneral],
               transnames[t - PastixNoTrans] );

        switch( spm->flttype ){
        case PastixComplex64:
            ret = z_bcsc_matvec_check( t, spm, pastix_data );
            break;

        case PastixComplex32:
            ret = c_bcsc_matvec_check( t, spm, pastix_data );
            break;

        case PastixFloat:
            ret = s_bcsc_matvec_check( t, spm, pastix_data );
            break;

        case PastixDouble:
        default:
            ret = d_bcsc_matvec_check( t, spm, pastix_data );
        }
        PRINT_RES(ret);
    }

    spmExit( spm );
    free( spm );

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
