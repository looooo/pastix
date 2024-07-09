/**
 * @file bvec_applyorder_tests.c
 *
 * @brief Test for applyorder with mpi distributed or replicated,
 *        in shared memory, with and without multidof.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Alycia Lisito
 * @date 2023-07-21
 *
 */
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define _GNU_SOURCE 1
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
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
    spmatrix_t     *original, *spm, *spmdof;
    int             variadic, mpi_type;
    pastix_int_t    m, n, nrhs, dim3, dof;
    spm_coeftype_t  type, flttype;
    pastix_fixdbl_t alpha, beta;
    int             clustnbr = 1;
    int             myrank   = 0;
    static int      dofmax   = 4; /* Maximum degree of freedom for multi-dof cases. */
    int             dmax     = 2; /* 0: one, 1: multi-constant, 2: multi-variadic   */
    int             check    = 1;
    int             rc       = 0;
    int             err      = 0;
    int             scatter  = 0;

    /**
     * Initializes parameters to default values.
     */
    pastixInitParam( iparm, dparm );
    iparm[IPARM_VERBOSE] = 0;

    /**
     * Gets options from command line.
     * Prevents from failing if no arguments is given.
     */
    if ( argc > 1 ) {
        pastixGetOptions( argc, argv, iparm, dparm, &check, &scatter, &driver, &filename );
    }

    /**
     * Starts-up PaStiX.
     */
    pastixInit( &pastix_data, MPI_COMM_WORLD, iparm, dparm );

    /**
     * Uses the Laplacian driver to pass parameters to the test.
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
     * Creates the sparse matrix with the driver
     */
    original = malloc( sizeof( spmatrix_t ) );
    /* Scatters the spm if scatter. */
    if ( scatter ) {
        rc = spmReadDriverDist( driver, filename, original, MPI_COMM_WORLD );
    }
    else {
        rc = spmReadDriver( driver, filename, original );
    }
    free( filename );

    if ( rc != SPM_SUCCESS ) {
        fprintf( stderr, "Failed to read the matrix\n" );
        pastixFinalize( &pastix_data );
        free( original );
        return rc;
    }

    /* Computes mpi_type. */
#if defined(PASTIX_WITH_MPI)
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &clustnbr );
    if ( scatter ) {
        mpi_type = 2;
    }
    else {
        mpi_type = 1;
    }
#else
    mpi_type = 0;
#endif

    spm = original;

    /**
     * Loops over dof configuration: None, Constant, Variadic.
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

        /* Makes sure internal spm is reset, and recomputes symbfact and mapping for the new matrix. */
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
