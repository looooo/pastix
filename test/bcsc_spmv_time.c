/**
 *
 * @file bcsc_spmv_time.c
 *
 * Tests performance of spmv for parallel and sequential version.
 * Validation is made in bcsc_matmvec_tests.c
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.1
 * @author Vincent Bridonneau
 * @author Esragul Korkmaz
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2023-07-21
 *
 * @precisions normal z -> c d s
 *
 **/
#include <pastix.h>
#include "common.h"

#include "z_tests.h"
#include "c_tests.h"
#include "d_tests.h"
#include "s_tests.h"

int main ( int argc, char **argv )
{
    pastix_data_t   *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t     iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    pastix_fixdbl_t  dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    spm_driver_t     driver;
    char            *filename;
    pastix_spm_t    *spm, spm2;
    int              check = 1;
    int              scatter = 0;
    int              rc, nrhs = 20;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );
    iparm[IPARM_VERBOSE]  = PastixVerboseNot;

    /**
     * Get options from command line
     */
    pastixGetOptions( argc, argv,
                      iparm, dparm,
                      &check, &scatter, &driver, &filename );

    /**
     * Initialize the PaStiX library
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Read the sparse matrix with the driver
     */
    spm = malloc( sizeof( pastix_spm_t ) );
    if ( scatter ) {
        rc = spmReadDriverDist( driver, filename, spm, MPI_COMM_WORLD );
    }
    else {
        rc = spmReadDriver( driver, filename, spm );
    }
    free( filename );
    if ( rc != SPM_SUCCESS ) {
        pastixFinalize( &pastix_data );
        return rc;
    }
    spmPrintInfo( spm, stdout );

    rc = spmCheckAndCorrect( spm, &spm2 );
    if ( rc != 0 ) {
        spmExit( spm );
        *spm = spm2;
    }

    if ( spm->flttype == SpmPattern ) {
        spmGenFakeValues( spm );
    }

    /**
     * Startup pastix to perform the analyze step
     */
    pastix_task_analyze( pastix_data, spm );
    pastix_subtask_spm2bcsc( pastix_data, spm );

    switch( spm->flttype ){
    case SpmComplex64:
        rc = z_bcsc_spmv_time( pastix_data, spm, nrhs );
        break;

    case SpmComplex32:
        rc = c_bcsc_spmv_time( pastix_data, spm, nrhs );
        break;

    case SpmFloat:
        rc = s_bcsc_spmv_time( pastix_data, spm, nrhs );
        break;

    case SpmDouble:
    default:
        rc = d_bcsc_spmv_time( pastix_data, spm, nrhs );
    }

    spmExit( spm );
    free( spm );

    pastixFinalize( &pastix_data );

    return rc;
}
