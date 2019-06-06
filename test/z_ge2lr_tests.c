/**
 *
 * @file z_ge2lr_tests.c
 *
 * Tests and validate the Xge2lr routine.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Gregoire Pichon
 * @date 2018-07-16
 *
 * @precisions normal z -> z c d s
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
#include "common/common.h"
#include <lapacke.h>
#include "z_tests.h"
#include "tests.h"

int main( int argc, char **argv )
{
    pastix_complex64_t *A;
    pastix_int_t n, r, lda;
    int mode, p, i, ret, rc = 0;
    double tolerance, threshold;
    double normA;
    test_param_t params;
    double eps = LAPACKE_dlamch_work('e');

    testGetOptions( argc, argv, &params, eps );

    fprintf( stdout, "%7s %4s %12s %12s %12s %12s\n",
             "Method", "Rank", "Time", "||A||_f", "||A-UVt||_f",
             "||A-UVt||_f/(||A||_f * eps)" );

    tolerance = params.tol_gen;
    threshold = params.threshold;
    for (n=params.n[0]; n<=params.n[1]; n+=params.n[2]) {
        lda = n;
        A = malloc( n * lda * sizeof(pastix_complex64_t) );

        for (p=params.prank[0]; p<=params.prank[1]; p+=params.prank[2]) {
            r = (p * n) / 100;

            for (mode=params.mode[0]; mode<=params.mode[1]; mode+=params.mode[2])
            {
                printf( "   -- Test GE2LR Tol=%e M=N=LDA=%ld R=%ld MODE=%d\n",
                        tolerance, (long)n, (long)r, mode );

                /*
                 * Generate a matrix of a given rank for the prescribed tolerance
                 */
                z_lowrank_genmat( mode, tolerance, threshold, r,
                                  n, n, A, lda, &normA );

                /* Let's test all methods we have */
                for(i=params.method[0]; i<=params.method[1]; i+=params.method[2])
                {
                    ret = z_lowrank_check_ge2lr( i, tolerance,
                                                 n, n, A, lda, normA,
                                                 ge2lrMethods[i][PastixComplex64-2] );
                    rc += (ret ? 1 : 0 );
                }
            }
        }
        free(A);
    }

    if( rc == 0 ) {
        printf( " -- All tests PASSED --\n" );
        return EXIT_SUCCESS;
    }
    else
    {
        printf( " -- %d tests FAILED --\n", rc );
        return EXIT_FAILURE;
    }
}
