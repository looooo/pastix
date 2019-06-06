/**
 *
 * @file z_ge2lr_stability.c
 *
 * Tests and validate the Xge2lr routine.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Esragul Korkmaz
 * @date 2019-03-19
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

static pastix_complex64_t mzone = -1.0;

int
z_lowrank_stability_ge2lr( int mode, pastix_compress_method_t method, double tolerance,
                           pastix_int_t m, pastix_int_t n,
                           pastix_complex64_t *A, pastix_int_t lda,
                           double normA,
                           fct_ge2lr_t core_ge2lr, pastix_int_t prank )
{
    pastix_lrblock_t    lrA;
    pastix_complex64_t *A2;
    pastix_complex64_t *u, *v;
    pastix_int_t minMN = pastix_imin(m, n);
    pastix_int_t i, ldu, ldv;
    double norm_residual;
    FILE *f;

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

    /*
     * Compress and then uncompress
     */
    core_ge2lr( -1., minMN, m, n, A, lda, &lrA );

    /* Let's check we have the maximal rank */
    assert( lrA.rk == minMN );

    /* Open the file */
    {
        char *filename;
        int rc;

        rc = asprintf( &filename, "stability_m%d_n%d_mode%d_tol%e_prank%d_%s.log",
                       (int)m, (int)n, mode, tolerance, (int)prank, compmeth_shnames[method]);
        f = fopen( filename, "w" );
        free(filename);
        (void)rc;
    }

    fprintf( f,
             "# method = %s\n"
             "# mode = %d\n"
             "# M = %d\n"
             "# N = %d\n"
             "# tol = %e\n"
             "# ||A|| = %e\n",
             compmeth_lgnames[method], mode,
             (int)m, (int)n, tolerance, normA );

    /*
     * Let's compute the frobenius norm of A - U[:,1:i] * V[:,1:i]^T for i in [1:minMN]
     */
    u = lrA.u;
    v = lrA.v;
    ldu = m;
    ldv = lrA.rkmax;
    for(i=0; i<minMN; i++) {
        norm_residual = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                             A2, m, NULL );

        fprintf( f, "%d %e\n", (int)(i+1), norm_residual/normA );

	/* Subtract the i^th outer product */
        cblas_zgerc( CblasColMajor, m, n,
                     CBLAS_SADDR(mzone),
                     u + ldu * i, 1,
                     v + i,       ldv,
                     A2,          m );
    }

    fclose(f);
    core_zlrfree(&lrA);
    free(A2);

    (void)normA;
    return 0;
}

int main( int argc, char **argv )
{
    pastix_complex64_t *A;
    pastix_int_t n, r, lda;
    int mode, p, i;
    double tolerance, threshold;
    double normA;
    test_param_t params;
    double eps = LAPACKE_dlamch_work('e');

    testGetOptions( argc, argv, &params, eps );

    tolerance = params.tol_gen;
    threshold = params.threshold;
    for (n=params.n[0]; n<=params.n[1]; n+=params.n[2]) {
        lda = n;
        A = malloc( n * lda * sizeof(pastix_complex64_t) );

        for (p=params.prank[0]; p<=params.prank[1]; p+=params.prank[2]) {
            r = (p * n) / 100;

            for (mode=params.mode[0]; mode<=params.mode[1]; mode+=params.mode[2]) {
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
                    z_lowrank_stability_ge2lr( mode, i, tolerance,
                                               n, n, A, lda, normA,
                                               ge2lrMethods[i][PastixComplex64-2], p );
                }
            }
        }
        free(A);
    }
}
