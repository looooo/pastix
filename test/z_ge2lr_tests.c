/**
 *
 * @file z_ge2lr_tests.c
 *
 * Tests and validate the Xge2lr routine.
 *
 * @version 5.1.0
 * @author Gregoire Pichon
 * @date 2016-11-24
 *
 * @precisions normal z -> c d s
 *
 **/
#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <pastix.h>
#include "../common/common.h"
#include <lapacke.h>
#include <cblas.h>
#include "../blend/solver.h"
#include "../kernels/pastix_zcores.h"

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
z_ge2lr_test( double tolerance, pastix_int_t rank,
              pastix_int_t m, pastix_int_t n, pastix_int_t lda )
{

    pastix_complex64_t *A, *A_RRQR, *A_SVD;
    pastix_lrblock_t    LR_RRQR, LR_SVD;

    double norm_dense;
    double norm_diff_RRQR, norm_diff_SVD;
    double res_SVD, res_RRQR;

    pastix_int_t minMN    = pastix_imin(m, n);
    int          mode     = 0;
    double       rcond    = (double) minMN;
    double       dmax     = 1.0;
    int          ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

    pastix_complex64_t *work;
    double *S;

    double alpha;

    if (lda < m || lda < n){
        printf("Invalid lda parameter\n");
        return -3;
    }

    A      = malloc(n * lda * sizeof(pastix_complex64_t));
    A_RRQR = malloc(n * lda * sizeof(pastix_complex64_t));
    A_SVD  = malloc(n * lda * sizeof(pastix_complex64_t));

    S    = malloc(minMN * sizeof(double));
    work = malloc(3 * pastix_imax(m, n)* sizeof(pastix_complex64_t));

    if ((!A)||(!A_SVD)||(!A_RRQR)||(!S)||(!work)){
        printf("Out of Memory \n ");
        free(A); free(A_RRQR); free(A_SVD); free(S); free(work);
        return -2;
    }

    /* Chose alpha such that alpha^rank = tolerance */
    alpha = exp(log(tolerance) / rank);

    if (mode == 0) {
        pastix_int_t i;
        S[0] = 1;
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

    norm_dense = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                      A, lda, NULL );

    /* Compress and then uncompress  */
    core_zge2lr_RRQR( tolerance,
                      m, n,
                      A, lda,
                      &LR_RRQR );

    core_zge2lr_SVD( tolerance,
                      m, n,
                      A, lda,
                      &LR_SVD );

    core_zlr2ge( PastixNoTrans, m, n,
                 &LR_RRQR,
                 A_RRQR, lda );

    core_zlr2ge( PastixNoTrans, m, n,
                 &LR_SVD,
                 A_SVD, lda );

    printf(" The rank of A is: RRQR %d SVD %d\n", LR_RRQR.rk, LR_SVD.rk);

    core_zgeadd( PastixNoTrans, m, n,
                 -1., A, lda,
                  1., A_RRQR, lda );

    core_zgeadd( PastixNoTrans, m, n,
                 -1., A, lda,
                  1., A_SVD, lda );

    norm_diff_RRQR = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                          A_RRQR, lda, NULL );

    norm_diff_SVD = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', m, n,
                                         A_SVD, lda, NULL );

    res_RRQR = norm_diff_RRQR / ( tolerance * norm_dense );
    res_SVD  = norm_diff_SVD  / ( tolerance * norm_dense );

    free(A);
    free(A_SVD);
    free(A_RRQR);
    free(S);
    free(work);

    if ((res_RRQR < 10) && (res_SVD < 10) && (LR_RRQR.rk >= LR_SVD.rk || LR_RRQR.rk == -1))
        return 0;
    return 1;
}

int main (int argc, char **argv)
{
    (void) argc;
    (void) argv;
    int err = 0;
    int ret;
    pastix_int_t m, r;
    double tolerance = 0.001;

    for (m=100; m<300; m+=100){
        for (r=10; r<100; r+=10){
            printf("   -- Test GE2LR M=N=LDA=%ld R=%ld\n", (long)m, (long)r);

            ret = z_ge2lr_test(tolerance, r, m, m, m);
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
