/**
 *
 * @file z_rradd_tests.c
 *
 * Tests and validate the Xrradd routine.
 *
 * @version 5.1.0
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
z_rradd_test( int mode, double tolerance, pastix_int_t rankA, pastix_int_t rankB,
              pastix_int_t mA, pastix_int_t nA,
              pastix_int_t mB, pastix_int_t nB,
              pastix_int_t offx, pastix_int_t offy )
{
    pastix_complex64_t *A, *B, *B_tmp;
    pastix_complex64_t *C_RRQR, *C_SVD;
    pastix_lrblock_t    LR_A_SVD, LR_B_SVD;
    pastix_lrblock_t    LR_A_RRQR, LR_B_RRQR;

    double norm_dense_A, norm_dense_B;
    double norm_diff_SVD, norm_diff_RRQR;
    double res_SVD, res_RRQR;

    pastix_int_t minMN_A = pastix_imin(mA, nA);
    pastix_int_t minMN_B = pastix_imin(mB, nB);

    double rcond = (double) minMN_A;
    double dmax  = 1.0;
    int ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

    pastix_complex64_t *work;
    double *SA, *SB;
    double alphaA, alphaB;

    A      = malloc(mA * nA * sizeof(pastix_complex64_t));
    B      = malloc(mB * nB * sizeof(pastix_complex64_t));
    C_RRQR = malloc(mA * nA * sizeof(pastix_complex64_t));
    C_SVD  = malloc(mA * nA * sizeof(pastix_complex64_t));
    SA     = malloc(minMN_A * sizeof(double));
    SB     = malloc(minMN_B * sizeof(double));
    work   = malloc(3 * pastix_imax(pastix_imax(mA, nA), pastix_imax(mB, nB)) * sizeof(pastix_complex64_t));

    if ( (!A) || (!B) || (!C_SVD) || (!C_RRQR) || (!SA) || (!SB) || (!work) ) {
        printf("Out of Memory \n ");
        free(A); free(B); free(C_RRQR); free(C_SVD); free(SA); free(SB); free(work);
        return -2;
    }

    /* Chose alpha such that alpha^rank = tolerance */
    alphaA = exp(log(tolerance) / rankA);
    alphaB = exp(log(tolerance) / rankB);

    if (mode == 0){
        pastix_int_t i;
        SA[0] = 1;
        SB[0] = 1;

        if (rankA == 0)
            SA[0] = 0.;
        if (rankB == 0)
            SB[0] = 0.;

        for (i=1; i<minMN_A; i++){
            SA[i] = SA[i-1] * alphaA;
        }
        for (i=1; i<minMN_B; i++){
            SB[i] = SB[i-1] * alphaB;
        }
    }

    /* Initialize A and B */
    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, mA, nA,
                         'U', ISEED,
                         'N', SA, mode, rcond,
                         dmax, mA, nA,
                         'N', A, mA, work );

    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, mB, nB,
                         'U', ISEED,
                         'N', SB, mode, rcond,
                         dmax, mB, nB,
                         'N', B, mB, work );

    norm_dense_A = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', mA, nA,
                                        A, mA, NULL );

    norm_dense_B = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', mB, nB,
                                        B, mB, NULL );


    core_zge2lr_SVD( tolerance,
                      mA, nA,
                      A, mA,
                      &LR_A_SVD );

    core_zge2lr_SVD( tolerance,
                      mB, nB,
                      B, mB,
                      &LR_B_SVD );

    core_zge2lr_RRQR( tolerance,
                      mA, nA,
                      A, mA,
                      &LR_A_RRQR );

    core_zge2lr_RRQR( tolerance,
                      mB, nB,
                      B, mB,
                      &LR_B_RRQR );

    printf(" The rank of A is: RRQR %d SVD %d\n", LR_A_RRQR.rk, LR_A_SVD.rk);
    printf(" The rank of B is: RRQR %d SVD %d\n", LR_B_RRQR.rk, LR_B_SVD.rk);

    if (LR_A_RRQR.rk == -1 || LR_B_RRQR.rk == -1 || (LR_A_RRQR.rk + LR_B_RRQR.rk) > pastix_imin(mA, nA)){
        printf("Operation non supported\n");
        free(A); free(B); free(C_RRQR); free(C_SVD); free(SA); free(SB); free(work);
        return 0;
    }
    if (LR_A_SVD.rk == -1 || LR_B_SVD.rk == -1 || (LR_A_SVD.rk + LR_B_SVD.rk) > pastix_imin(mA, nA)){
        printf("Operation non supported\n");
        free(A); free(B); free(C_RRQR); free(C_SVD); free(SA); free(SB); free(work);
        return 0;
    }

    /* Add A and B in their LR format */
    core_zrradd_SVD( tolerance, PastixNoTrans, -1.0,
                     mA, nA, &LR_A_SVD,
                     mB, nB, &LR_B_SVD,
                     offx, offy );

    core_zrradd_RRQR( tolerance, PastixNoTrans, -1.0,
                      mA, nA, &LR_A_RRQR,
                      mB, nB, &LR_B_RRQR,
                      offx, offy );

    printf(" The rank of A+B is: RRQR %d SVD %d\n", LR_B_RRQR.rk, LR_B_SVD.rk);

    /* Build uncompressed LR+LR matrix */
    core_zlr2ge( PastixNoTrans, mB, nB,
                 &LR_B_SVD,
                 C_SVD, mB );

    core_zlr2ge( PastixNoTrans, mB, nB,
                 &LR_B_RRQR,
                 C_RRQR, mB );

    /* Compute A+B in dense */
    B_tmp = B + offx + mB * offy;
    core_zgeadd( PastixNoTrans, mA, nA,
                 -1.0, A, mA,
                 1.0, B_tmp, mB );

    /* Compute norm of dense and LR matrices */
    core_zgeadd( PastixNoTrans, mB, nB,
                 -1., B, mB,
                  1., C_SVD, mB );

    core_zgeadd( PastixNoTrans, mB, nB,
                 -1., B, mB,
                 1., C_RRQR, mB );

    norm_diff_SVD = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', mB, nB,
                                         C_SVD, mB, NULL );
    norm_diff_RRQR = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', mB, nB,
                                          C_RRQR, mB, NULL );

    if ( (rankA != 0) || (rankB != 0) ){
        res_RRQR = norm_diff_RRQR / ( tolerance * (norm_dense_A + norm_dense_B) );
        res_SVD  = norm_diff_SVD  / ( tolerance * (norm_dense_A + norm_dense_B) );
    }
    else{
        res_RRQR = norm_diff_RRQR;
        res_SVD  = norm_diff_SVD;
    }

    printf("RES SVD=%.3g RRQR=%.3g\n", res_SVD, res_RRQR);

    free(A);
    free(B);
    free(C_SVD);
    free(C_RRQR);
    free(SA);
    free(SB);
    free(work);

    if ((res_RRQR < 10) && (res_SVD < 10))
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
    double tolerance = 0.01;

    for (m=200; m<=400; m+=100){
        for (r=0; r <= (m/2); r += ( r + 1 ) ) {
            printf("   -- Test RRADD MA=NA=LDA=%ld MB=NB=LDB=%ld RA=%ld RB=%ld\n", (long)m, (long)m, (long)r, (long)(r/2));

            ret = z_rradd_test(0, tolerance, r, r/2,
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
