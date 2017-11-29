/**
 *
 * @file z_rradd_tests.c
 *
 * Tests and validate the core_zrradd() routine.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Gregoire Pichon
 * @date 2016-11-24
 *
 * @precisions normal z -> c d s
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
#include <cblas.h>
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"
#include "kernels/pastix_zlrcores.h"

#define PRINT_RES(_ret_)                        \
    if(_ret_ == -1) {                           \
        printf("UNDEFINED\n");                  \
    }                                           \
    else if(_ret_ > 0) {                        \
        printf("FAILED(%d)\n", _ret_);          \
        err++;                                  \
    }                                           \
    else {                                      \
        printf("SUCCESS\n");                    \
    }

int
z_lowrank_genmat( int mode, double tolerance, pastix_int_t rank,
                  pastix_int_t m, pastix_int_t n, pastix_int_t lda,
                  pastix_complex64_t **Aptr,
                  pastix_lrblock_t    *lrA_svd,
                  pastix_lrblock_t    *lrA_rrqr,
                  double              *normA );

int
z_lowrank_rradd( pastix_int_t mA, pastix_int_t nA,
                 pastix_int_t offx, pastix_int_t offy,
                 const pastix_complex64_t *A, const pastix_lrblock_t *lrA, double normA,
                 pastix_int_t mB, pastix_int_t nB,
                 const pastix_complex64_t *B, const pastix_lrblock_t *lrB, double normB,
                 pastix_lr_t *lowrank );

int
z_rradd_test( int mode, double tolerance, pastix_int_t rankA, pastix_int_t rankB,
              pastix_int_t mA, pastix_int_t nA,
              pastix_int_t mB, pastix_int_t nB,
              pastix_int_t offx, pastix_int_t offy )
{
    pastix_complex64_t *A, *B;
    pastix_lrblock_t    lrA_svd, lrB_svd;
    pastix_lrblock_t    lrA_rrqr, lrB_rrqr;
    pastix_lr_t         lr_RRQR, lr_SVD;
    double              norm_dense_A, norm_dense_B;
    int                 rc = 0;

    /*
     * Generate matrices of rankA and rankB and thei compress SVD/RRQR versions
     */
    z_lowrank_genmat( mode, tolerance, rankA, mA, nA, mA,
                      &A, &lrA_svd, &lrA_rrqr, &norm_dense_A );
    z_lowrank_genmat( mode, tolerance, rankB, mB, nB, mB,
                      &B, &lrB_svd, &lrB_rrqr, &norm_dense_B );

    printf( " The rank of A   is: RRQR %d SVD %d rkmax %d\n",
            lrA_rrqr.rk, lrA_svd.rk, (int)core_get_rklimit( mA, nA ) );
    printf( " The rank of B   is: RRQR %d SVD %d rkmax %d\n",
            lrB_rrqr.rk, lrB_svd.rk, (int)core_get_rklimit( mB, nB ) );


    {
        lr_RRQR.compress_when       = PastixCompressWhenEnd;
        lr_RRQR.compress_method     = PastixCompressMethodRRQR;
        lr_RRQR.compress_min_width  = 0;
        lr_RRQR.compress_min_height = 0;
        lr_RRQR.tolerance  = tolerance;
        lr_RRQR.core_ge2lr = core_zge2lr_rrqr;
        lr_RRQR.core_rradd = core_zrradd_rrqr;

        printf(" RRQR: ");
        rc += z_lowrank_rradd( mA, nA, offx, offy, A, &lrA_rrqr, norm_dense_A,
                               mB, nB,             B, &lrB_rrqr, norm_dense_B,
                               &lr_RRQR );
    }

    {
        lr_SVD.compress_when       = PastixCompressWhenEnd;
        lr_SVD.compress_method     = PastixCompressMethodSVD;
        lr_SVD.compress_min_width  = 0;
        lr_SVD.compress_min_height = 0;
        lr_SVD.tolerance  = tolerance;
        lr_SVD.core_ge2lr = core_zge2lr_svd;
        lr_SVD.core_rradd = core_zrradd_svd;

        printf(" SVD:  ");
        rc += z_lowrank_rradd( mA, nA, offx, offy, A, &lrA_svd, norm_dense_A,
                               mB, nB,             B, &lrB_svd, norm_dense_B,
                               &lr_SVD );
    }

    return rc;
}

int main (int argc, char **argv)
{
    (void) argc;
    (void) argv;
    int err = 0;
    int ret;
    pastix_int_t m, r, rmax;
    double tolerance = 1e-8;

    for (m=200; m<=400; m+=100){
        rmax = core_get_rklimit( m, m );
        for (r=0; (r + (r/2)) < rmax; r += ( r + 1 ) ) {
            printf("   -- Test RRADD MA=NA=LDA=%ld MB=NB=LDB=%ld RA=%ld RB=%ld\n",
                   (long)m, (long)m, (long)r/2, (long)r);

            ret = z_rradd_test(0, tolerance, r/2, r,
                               m, m,
                               m, m,
                               0, 0);
            PRINT_RES(ret);
        }
    }


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
