/**
 * @file bvec_applyorder_tests.c
 *
 * @brief Test for applyorder with mpi distributed or replicated,
 *        in shared memory, with and without multidof.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Alycia Lisito
 * @date 2022-12-15
 *
 */
#define _GNU_SOURCE 1
#include <pastix.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#if defined(HAVE_GETOPT_H)
#include <getopt.h>
#endif  /* defined(HAVE_GETOPT_H) */

#include "common.h"
#include "blend/solver.h"
#include "bcsc/bcsc.h"

#include "z_tests.h"
#include "c_tests.h"
#include "d_tests.h"
#include "s_tests.h"

#define PRINT_RES(_ret_)                        \
    if(_ret_) {                                 \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

char* fltnames[] = { "Pattern", "", "Float", "Double", "Complex32", "Complex64" };
char* dofs[]     = { "Single DoF", "Constant Multi DoF", "Variadic Multi DoF" };
char* mpi[]      = { "Shared Memory", "Replicated", "Distributed" };

int main (int argc, char **argv)
{
    pastix_data_t  *pastix_data = NULL; /*< Pointer to the storage structure required by pastix */
    pastix_int_t    iparm[IPARM_SIZE];  /*< Integer in/out parameters for pastix                */
    double          dparm[DPARM_SIZE];  /*< Floating in/out parameters for pastix               */
    spm_driver_t    driver = (spm_driver_t)-1;
    char           *filename = NULL;
    spmatrix_t     *original, *spm, *spmdof, spmtmp;
    int             variadic, mpi_type;
    pastix_int_t    m, n, nrhs, dim3, dof;
    int             mpi_begin = 0;
    int             mpi_end   = 0;
    spm_coeftype_t  type, flttype;
    pastix_fixdbl_t alpha, beta;
    int             clustnbr, myrank;
    static int      dofmax = 4; /* Maximum degree of freedom for multi-dof cases. */
    int             dmax   = 2; /* 0: one, 1: multi-constant, 2: multi-variadic   */
    int             check  = 1;
    int             rc     = 0;
    int             err    = 0;

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
        pastixGetOptions( argc, argv, iparm, dparm, &check, &driver, &filename );
    }

    /**
     * Startup PaStiX
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Use the Laplacian driver to pass parameters to the test
     */
    if ( driver == (spm_driver_t)-1 ) {
        driver  = SpmDriverLaplacian;
        flttype = SpmPattern;
        m       = 10;
        n       = 10;
        rc      = asprintf( &filename, "p:%d:%d", (int)m, (int)n );
        assert( rc != -1 );
    }
    else {
        assert( driver == SpmDriverLaplacian );
        spmParseLaplacianInfo( filename, &flttype, &m, &n, &dim3, &alpha, &beta, &dof );
    }

    /**
     * Create the sparse matrix with the driver
     */
    original = malloc( sizeof( spmatrix_t ) );
    rc       = spmReadDriver( driver, filename, original );
    free( filename );

    if ( rc != SPM_SUCCESS ) {
        fprintf( stderr, "Failed to read the matrix\n" );
        pastixFinalize( &pastix_data );
        free( original );
        return rc;
    }

    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &clustnbr );
    if ( clustnbr > 1 ) {
        mpi_begin = 1;
        mpi_end   = 2;
        dmax      = 2; /* Variadic is not working with MPI */
    }

    for ( mpi_type = mpi_begin; mpi_type <= mpi_end; mpi_type++ ) {

        /**
         * Scatters the spm if distributed.
         */
        if ( mpi_type == 2 ) {
            rc = spmScatter( &spmtmp, -1, original, -1, NULL, 1, MPI_COMM_WORLD );
            if ( rc != SPM_SUCCESS ) {
                fprintf( stderr, "Failed to scatter the spm\n" );
                err++;
                continue;
            }
            spm = &spmtmp;
        }
        else {
            spm = original;
        }

        /**
         * Loop over dof configuration: None, Constant, Variadic
         */
        for ( dof = 0; dof < dmax; dof++ ) {

            /* Multidof. */
            if ( dof > 0 ) {
                spmdof   = malloc( sizeof(spmatrix_t) );
                variadic = dof - 1;
                rc       = spmDofExtend( spm, variadic, dofmax, spmdof );
            }
            else {
                spmdof = spm;
            }

            /* Make sure internal spm is reset, and recompute symbfact and mapping for the new matrix */
            pastix_data->csc = spmdof;
            pastix_task_analyze( pastix_data, spmdof );

            pastix_data->bcsc = calloc( 1, sizeof( pastix_bcsc_t ) );
            bcsc_init_struct( spmdof, pastix_data->solvmatr, pastix_data->bcsc );

            for ( nrhs = 1; nrhs <= 5; nrhs += 4 ) {
                for ( type = SpmFloat; type <= SpmComplex64; type ++ ) {

                    if ( myrank == 0 ) {
                        fprintf( stdout,
                                 "  Check b == (P^t (P b)) / Case %-18s - %-13s - %-9s - %ld nrhs:  ",
                                 dofs[dof], mpi[mpi_type], fltnames[type], (long)nrhs );
                    }

                    /* Checks the result. */
                    switch ( type ) {
                    case SpmFloat :
                        rc = s_bvec_applyorder_check( pastix_data, spmdof, nrhs );
                        break;
                    case SpmDouble :
                        rc = d_bvec_applyorder_check( pastix_data, spmdof, nrhs );
                        break;
                    case SpmComplex32 :
                        rc = c_bvec_applyorder_check( pastix_data, spmdof, nrhs );
                        break;
                    case SpmComplex64 :
                        rc = z_bvec_applyorder_check( pastix_data, spmdof, nrhs );
                        break;
                    default :
                        fprintf( stderr, "bvec_applyorder_test: wrong flttype. \n" );
                        break;
                    }
                    spmdof->flttype = SpmPattern;

                    if ( myrank == 0 ) {
                        PRINT_RES(rc);
                    }
                    err += rc;
                }
            }

            bcsc_exit_struct( pastix_data->bcsc );
            free( pastix_data->bcsc );
            pastix_data->bcsc = NULL;

            if ( dof > 0 ) {
                spmExit( spmdof );
                free( spmdof );
            }
        }

        if ( mpi_type == 2 ) {
            spmExit( spm );
        }
    }

    if ( pastix_data->procnum == 0 ) {
        if( err == 0 ) {
            printf(" -- All tests PASSED --\n");
        }
        else {
            printf(" -- %d tests FAILED --\n", err);
        }
    }

    pastixFinalize( &pastix_data );
    spmExit( original );
    free( original );

    (void)clustnbr;
    return rc;
}

/**
 * @endcode
 */
