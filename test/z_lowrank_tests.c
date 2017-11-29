/**
 *
 * @file z_lowrank_tests.c
 *
 * Test functions for the low-rank kernels.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @date 2016-11-24
 *
 * @precisions normal z -> z c d s
 *
 **/
#include <stdio.h>
#include <assert.h>
#include <pastix.h>
#include "common/common.h"
#include <lapacke.h>
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"
#include "kernels/pastix_zlrcores.h"

static int ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

int
z_lowrank_genmat( int mode, double tolerance, pastix_int_t rank,
                  pastix_int_t m, pastix_int_t n, pastix_int_t lda,
                  pastix_complex64_t **Aptr,
                  pastix_lrblock_t    *lrA_svd,
                  pastix_lrblock_t    *lrA_rrqr,
                  double              *normA )
{
    pastix_complex64_t *A;
    pastix_int_t minMN    = pastix_imin(m, n);
    double       rcond    = (double) minMN;
    double       dmax     = 1.0;
    pastix_complex64_t *work;
    double *S;

    double alpha;

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
    if (rank > pastix_imin( m, n )) {
        fprintf(stderr, "Invalid rank parameter\n");
        return -3;
    }

    A      = malloc(n * lda * sizeof(pastix_complex64_t));
    *Aptr  = A;

    S    = malloc(minMN * sizeof(double));
    work = malloc(3 * pastix_imax(m, n)* sizeof(pastix_complex64_t));

    if ( (!A) || (!S) || (!work) ) {
        fprintf(stderr, "Out of Memory\n");
        free(A); free(S); free(work);
        return -2;
    }

    /*
     * Choose alpha such that alpha^rank = tolerance
     */
    alpha = exp(log(tolerance) / rank);

    if (mode == 0) {
        pastix_int_t i;
        S[0] = 1;

        if (rank == 0) {
            S[0] = 0.;
        }
        for (i=1; i<minMN; i++){
            S[i] = S[i-1] * alpha;
        }
    }

    /* Initialize A */
    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, m, n,
                         'U', ISEED,
                         'N', S, mode, rcond,
                         dmax, m, n,
                         'N', A, lda, work );

    *normA = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                  A, lda, NULL );

    /* Adjust the tolerance base on the generated matrix norm */
    if ( *normA > 0. ) {
        tolerance /= *normA;
    }

    /* Compress and then uncompress  */
    core_zge2lr_rrqr( tolerance, -1,
                      m, n,
                      A, lda,
                      lrA_rrqr );

    core_zge2lr_svd( tolerance, -1,
                      m, n,
                      A, lda,
                      lrA_svd );

    free(S); free(work);
    return 0;
}


int
z_lowrank_rradd( pastix_int_t mA, pastix_int_t nA,
                 pastix_int_t offx, pastix_int_t offy,
                 const pastix_complex64_t *A, const pastix_lrblock_t *lrA, double normA,
                 pastix_int_t mB, pastix_int_t nB,
                 const pastix_complex64_t *B, const pastix_lrblock_t *lrB, double normB,
                 pastix_lr_t *lowrank )
{
    pastix_complex64_t *Bfr, *Blr;
    pastix_complex64_t mone = -1.0;
    double norm_diff, res;
    pastix_int_t rkABmax;
    pastix_int_t rankmax  = core_get_rklimit(mB, nB);
    int          rc = 0;
    pastix_lrblock_t lrAB;

    rkABmax = lrA->rk + lrB->rk;

    if ( (lrA->rk == -1) ||
         (lrB->rk == -1) ||
         (rkABmax > rankmax) )
    {
        printf("Operation not supported\n");
        return -1;
    }

    /* Init lrAB */
    memset( &lrAB, 0, sizeof(pastix_lrblock_t));

    /* Backup B into AB */
    core_zlrcpy( lowrank, PastixNoTrans, 1.,
                 mB, nB, lrB, mB, mB, &lrAB, 0, 0 );

    /* Add A and B in their LR format */
    lowrank->core_rradd( lowrank, PastixNoTrans, &mone,
                         mA, nA, lrA,
                         mB, nB, &lrAB,
                         offx, offy );

    /*
     * Check || (A+B) - (c(A)+c(B)) || < tol * || A+B ||
     * What we really want is || (A+B) - c(A+B) || < tol * || A+B ||
     */
    Bfr = malloc( mB * nB * sizeof(pastix_complex64_t) );
    Blr = malloc( mB * nB * sizeof(pastix_complex64_t) );

    /* Dense sum */
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', mB, nB,
                        B, mB, Bfr, mB );
    core_zgeadd( PastixNoTrans, mA, nA,
                 -1.0, A,                      mA,
                  1.0, Bfr + offx + mB * offy, mB );

    /* Uncompresse the low-rank sum */
    core_zlr2ge( PastixNoTrans, mB, nB,
                 lrB, Blr, mB );

    /* Compute the diff */
    core_zgeadd( PastixNoTrans, mB, nB,
                 -1., Bfr, mB,
                  1., Blr, mB );

    norm_diff = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', mB, nB,
                                     Blr, mB, NULL );

    if ( (lrA->rk != 0) || (lrB->rk != 0) ) {
        res = norm_diff / ( lowrank->tolerance * (normA + normB) );
    }
    else {
        res = norm_diff;
    }

    free(Bfr);
    free(Blr);

    printf("A+B %d RES %e ", lrAB.rk, res);

    /*
     * Check the validity of the results
     */
    /* Check the correctness of the result */
    if ( res > 10.0 ) {
        rc += 1;
    }

    /* Check that final matrix is not full rank, we should have exited before */
    if ( lrAB.rk == -1 ) {
        rc += 2;
    }

    /* Check that final rank does not exceed the sum of the ranks */
    if ( lrAB.rk > rkABmax ) {
        rc += 4;
    }

    /* Check that final rank is larger than the minimal rank (rA, rB, abs(rA-rb)) */
    if ( lrAB.rk < pastix_imin( pastix_imin( lrA->rk, lrB->rk ), abs(lrA->rk - lrB->rk) ) )
    {
        rc += 8;
    }

    if ( rc == 0 ) {
        printf("SUCCESS\n");
    }
    else {
        printf("FAILED(%d)\n", rc);
    }

    return rc;
}

