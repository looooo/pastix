/**
 *
 * @file bvec_tests.c
 *
 * Tests for preformance of parallel subroutines.
 * We test each subroutines of vectorial functions (copy, axpy, scale, ...)
 * and we compute the average time over 50 calls to each functions.
 * Size of the matrices/vectors are given by the first two dimensions of the
 * Laplacian driver.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Vincent Bridonneau
 * @author Esragul Korkmaz
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2023-07-21
 *
 **/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#include <pastix.h>
#include "refinement/z_refine_functions.h"
#include "bcsc/bcsc.h"
#include "bcsc/bcsc_z.h"
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */

#include "z_tests.h"
#include "c_tests.h"
#include "d_tests.h"
#include "s_tests.h"

char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* schednames[] = { "Sequential", "Static" };

int main ( int argc, char **argv )
{
    pastix_data_t      *pastix_data = NULL;
    pastix_int_t        iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    pastix_fixdbl_t     dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    spm_driver_t        driver = (spm_driver_t)-1;
    spmatrix_t         *spm, spm2;
    char               *filename = NULL;
    int                 check    = 1;
    int                 scatter  = 0;
    pastix_int_t        m, n, s, t;
    int rc = 0;
    pastix_int_t    dim3, dof;
    spm_coeftype_t  flttype;
    pastix_fixdbl_t alpha, beta;

    /**
     * Initialize parameters to default values
     */
    pastixInitParam( iparm, dparm );
    iparm[IPARM_VERBOSE] = 0;

    /**
     * Get options from command line
     * Prevent from failing if no arguments is given
     */
    if ( argc > 1 ) {
        pastixGetOptions( argc, argv,
                          iparm, dparm,
                          &check, &scatter, &driver, &filename );
    }

    /**
     * Initialize the PaStiX library
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    if ( driver == (spm_driver_t)-1 ) {
        driver = SpmDriverLaplacian;
        flttype = SpmDouble;
        m = 400;
        n = iparm[IPARM_GMRES_IM];
        rc = asprintf( &filename, "d:%d:%d", (int)m, (int)n );
        assert( rc != -1 );
    }
    else {
        spmParseLaplacianInfo( filename, &flttype, &m, &n, &dim3, &alpha, &beta, &dof );
    }

    /**
     * Read the sparse matrix with the driver
     */
    spm = malloc( sizeof( spmatrix_t ) );
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

    rc = spmCheckAndCorrect( spm, &spm2 );
    if ( rc != 0 ) {
        spmExit( spm );
        *spm = spm2;
    }

    /**
     * Run preprocessing steps required to generate solverMatrix
     */
    pastix_task_analyze( pastix_data, spm );

    /**
     * Generate the blocked csc
     */
    pastix_subtask_spm2bcsc( pastix_data, spm );

    if ( check ) {
        for( s=PastixSchedSequential; s<=PastixSchedStatic; s++ )
        {
            pastix_data->iparm[IPARM_SCHEDULER] = s;
            for( t=SpmFloat; t<=SpmComplex64; t++ )
            {
                if ( pastix_data->procnum == 0 ) {
                    printf("  Case %s - %s:\n",
                           fltnames[t], schednames[s] );
                }
                switch( t ){
                case SpmComplex64:
                    rc = z_bvec_check( pastix_data );
                    break;

                case SpmComplex32:
                    rc = c_bvec_check( pastix_data );
                    break;

                case SpmFloat:
                    rc = s_bvec_check( pastix_data );
                    break;

                case SpmDouble:
                default:
                    rc = d_bvec_check( pastix_data );
                }
            }
        }
    }
    else {
        switch( flttype ){
        case SpmComplex64:
            rc = z_bvec_time( pastix_data );
            break;

        case SpmComplex32:
            rc = c_bvec_time( pastix_data );
            break;

        case SpmFloat:
            rc = s_bvec_time( pastix_data );
            break;

        case SpmDouble:
        default:
            rc = d_bvec_time( pastix_data );
        }
    }

    pastixFinalize( &pastix_data );
    spmExit( spm );
    free( spm );

    return rc;
}
