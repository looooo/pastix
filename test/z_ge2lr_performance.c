/**
 *
 * @file z_ge2lr_performance.c
 *
 * Tests and validate the Xge2lr routine.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @author Esragul Korkmaz
 * @date 2018-07-16
 *
 * @precisions normal z -> z c d s
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
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
#include <cblas.h>
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"
#include "kernels/pastix_zlrcores.h"
#include "kernels/pastix_lowrank.h"
#include "flops.h"
#include "z_tests.h"
#include "tests.h"

double
z_lowrank_ge2lr_performance( FILE *f, pastix_compress_method_t method, int prank, int mode, double tol_cmp,
                             int m, int n, pastix_complex64_t *A, pastix_int_t lda,
                             double normA )
{
    fct_ge2lr_t core_ge2lr = ge2lrMethods[method][PastixComplex64-2];
    pastix_complex64_t *A2;
    pastix_lrblock_t    lrA;
    pastix_int_t minMN = pastix_imin(m, n);
    pastix_fixdbl_t flops, gflops;
    double resid, normR, timer;
    double tol_gen = tol_cmp;
    int rank, crank;

    if (m < 0) {
        fprintf(stderr, "Invalid m parameter\n");
        return -4;
    }
    if (n < 0) {
        fprintf(stderr, "Invalid n parameter\n");
        return -5;
    }
    if (lda < m) {
        fprintf(stderr, "Invalid lda parameter\n");
        return -6;
    }

    /* Backup A in A2 */
    A2 = malloc( m * n * sizeof(pastix_complex64_t));
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', m, n,
                         A, lda, A2, m );

    /* Compress A */
    clockStart(timer);
    core_ge2lr( tol_cmp, minMN, m, n, A, lda, &lrA );
    clockStop(timer);
    timer = clockVal(timer);

    crank = lrA.rk;
    rank = ( prank * n ) / 100;

    flops = FLOPS_ZGEQRF( m, rank ) +
        FLOPS_ZUNMQR( m, n - rank, rank, PastixLeft ) +
        FLOPS_ZUNGQR( m, rank, rank);
    gflops = flops * 1e-9 / timer;

    /*
     * Let's check the result
     */
    core_zlr2ge( PastixNoTrans, m, n,
                 &lrA, A2, lda );

    core_zgeadd( PastixNoTrans, m, n,
                 -1., A,  lda,
                  1., A2, lda );

    /* Frobenius norm of ||A - (U_i *V_i)|| */
    normR = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n, A2, lda, NULL );
    resid = normR / ( tol_cmp * normA );

    if ( f == stdout ) {
        fprintf( f, "%7s %5d %4d %e %6d %5d %e %e %5d %e %e %e %e %s\n",
                 compmeth_shnames[method], prank, mode, tol_gen,
                 n, rank, normA,
                 tol_cmp, crank, timer, gflops, normR, resid,
                 (resid > 10.) ? "FAILED" : "SUCCESS" );
    }
    else {
        fprintf( f, "%s;%d;%d;%e;%d;%d;%e;%e;%d;%e;%e;%e;%e;%s\n",
                 compmeth_shnames[method], prank, mode, tol_gen,
                 n, rank, normA,
                 tol_cmp, crank, timer, gflops, normR, resid,
                 (resid > 10.) ? "FAILED" : "SUCCESS" );
    }

    free(A2);
    core_zlrfree(&lrA);
    return (resid > 10.);
}

int main( int argc, char **argv )
{
    pastix_complex64_t *A;
    pastix_int_t n, r, lda;
    int mode, p, i, ret, rc = 0;
    double tolerance, tol_cmp, threshold;
    double normA;
    test_param_t params;
    double eps = LAPACKE_dlamch_work('e');

    testGetOptions( argc, argv, &params, eps );

    if ( params.output == stdout ) {
        fprintf( params.output,
                 "%7s %5s %4s %12s %6s %5s %12s %12s %5s %12s %12s %12s %12s\n",
                 "Method", "PRank", "Mode", "TolGen",
                 "N", "Rank", "NormA",
                 "TolCmp", "CRank", "Time", "GFlops", "||A-UVt||_f", "||A-UVt||_f/(||A||_f * eps)" );
    }
    else {
        fprintf( params.output,
                 "%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;Check\n",
                 "Method", "PRank", "Mode", "TolGen",
                 "N", "Rank", "NormA",
                 "TolCmp", "CRank", "Time", "GFlops", "||A-UVt||_f", "||A-UVt||_f/(||A||_f * eps)" );
    }
    tolerance = params.tol_gen;
    threshold = params.threshold;
    tol_cmp = params.tol_cmp;

    for (n=params.n[0]; n<=params.n[1]; n+=params.n[2]) {
        lda = n;
        A = malloc( n * lda * sizeof(pastix_complex64_t) );

        for (p=params.prank[0]; p<=params.prank[1]; p+=params.prank[2]) {
            r = (p * n) / 100;

            for (mode=params.mode[0]; mode<=params.mode[1]; mode+=params.mode[2])
            {
                /*
                 * Generate a matrix of a given rank for the prescribed tolerance
                 */
                z_lowrank_genmat( mode, tolerance, threshold, r,
                                  n, n, A, lda, &normA );

                /* Let's test all methods we have */
                for(i=params.method[0]; i<=params.method[1]; i+=params.method[2])
                {
                    ret = z_lowrank_ge2lr_performance( params.output, i, p, mode, tol_cmp,
                                                       n, n, A, lda, normA );
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
