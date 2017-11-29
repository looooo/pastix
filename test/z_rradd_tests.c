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

int main (int argc, char **argv)
{
    (void) argc;
    (void) argv;
    int err = 0;
    int ret, rc, mode = 0;
    pastix_int_t m, r, rmax;
    double       eps = LAPACKE_dlamch_work('e');
    double       tolerance = sqrt(eps);
    pastix_lr_t  lr_RRQR, lr_SVD;
    pastix_complex64_t *A, *B;
    pastix_lrblock_t    lrA_svd, lrB_svd;
    pastix_lrblock_t    lrA_rrqr, lrB_rrqr;
    double              norm_dense_A, norm_dense_B;

    /* Initialize the RRQR structure */
    lr_RRQR.compress_when       = PastixCompressWhenEnd;
    lr_RRQR.compress_method     = PastixCompressMethodRRQR;
    lr_RRQR.compress_min_width  = 0;
    lr_RRQR.compress_min_height = 0;
    lr_RRQR.tolerance  = tolerance;
    lr_RRQR.core_ge2lr = core_zge2lr_rrqr;
    lr_RRQR.core_rradd = core_zrradd_rrqr;

    /* Initialize the SVD structure */
    lr_SVD.compress_when       = PastixCompressWhenEnd;
    lr_SVD.compress_method     = PastixCompressMethodSVD;
    lr_SVD.compress_min_width  = 0;
    lr_SVD.compress_min_height = 0;
    lr_SVD.tolerance  = tolerance;
    lr_SVD.core_ge2lr = core_zge2lr_svd;
    lr_SVD.core_rradd = core_zrradd_svd;

    for (m=200; m<=400; m+=100){
        rmax = core_get_rklimit( m, m );
        for (r=0; (r + (r/2)) < rmax; r += ( r + 1 ) ) {
            pastix_int_t rankA = r/2;
            pastix_int_t rankB = r;
            pastix_int_t offx = 0;
            pastix_int_t offy = 0;

            printf( "  -- Test RRADD MA=NA=LDA=%ld MB=NB=LDB=%ld RA=%ld RB=%ld rkmax=%ld\n",
                    (long)m, (long)m, (long)rankA, (long)rankB,
                    (long)core_get_rklimit( m, m ) );

            /*
             * Generate matrices of rankA and rankB and thei compress SVD/RRQR versions
             */
            z_lowrank_genmat( mode, tolerance, rankA, m, m, m,
                              &A, &lrA_svd, &lrA_rrqr, &norm_dense_A );
            z_lowrank_genmat( mode, tolerance, rankB, m, m, m,
                              &B, &lrB_svd, &lrB_rrqr, &norm_dense_B );

            ret = 0;
            {
                printf("     RRQR: A %2d B %2d ", lrA_rrqr.rk, lrB_rrqr.rk );
                rc = z_lowrank_rradd( m, m, offx, offy, A, &lrA_rrqr, norm_dense_A,
                                      m, m,             B, &lrB_rrqr, norm_dense_B,
                                      &lr_RRQR );
                assert( rc >= 0 );
                ret += rc;
            }
            core_zlrfree( &lrA_rrqr );
            core_zlrfree( &lrB_rrqr );

            {
                printf("     SVD:  A %2d B %2d ", lrA_svd.rk, lrB_svd.rk );
                rc += z_lowrank_rradd( m, m, offx, offy, A, &lrA_svd, norm_dense_A,
                                       m, m,             B, &lrB_svd, norm_dense_B,
                                       &lr_SVD );
                assert( rc >= 0 );
                ret += rc;
            }
            core_zlrfree( &lrA_svd );
            core_zlrfree( &lrB_svd );

            free(A);
            free(B);
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
