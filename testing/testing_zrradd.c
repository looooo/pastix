/**
 *  @file: testing_zrradd.c
 *
 *  testing LR + LR operation
 *
 * @precisions normal z -> c d s
 *
 */
#include "testing_zmain.h"

int testing_zrradd(int argc, char **argv)
{
    /* Check for number of arguments*/
    if ( argc != 9) {
        USAGE("RRADD", "tol MA NA RA MB NB RB offx offy",
              "   - tol    : tolerance for SVD compression\n"
              "   - MA     : number of rows of matrices A\n"
              "   - NA     : number of columns of matrices A\n"
              "   - RA     : rank of matrix A\n"
              "   - MB     : number of rows of matrices B\n"
              "   - NB     : number of columns of matrices B\n"
              "   - RB     : rank of matrix B\n"
              "   - offx   : first row of B with respect to A\n"
              "   - offy   : first column of B with respect to A\n");
        return -1;
    }

    double tolerance  = atof(argv[0]);
    pastix_int_t MA   = atoi(argv[1]);
    pastix_int_t NA   = atoi(argv[2]);
    pastix_int_t MB   = atoi(argv[4]);
    pastix_int_t NB   = atoi(argv[5]);
    pastix_int_t offx = atoi(argv[7]);
    pastix_int_t offy = atoi(argv[8]);

    if (MA != MB || NA != NB){
        printf("INVALID PARAMETERS\n");
        return -1;
    }

    int RA = atoi(argv[3]);
    int RB = atoi(argv[6]);

    double eps, res;
    pastix_complex64_t *A, *B, *C;
    pastix_lrblock_t    LR_A, LR_B;

    double *SA, *SB;

    pastix_int_t LDA = MA;
    pastix_int_t LDB = MB;

    pastix_complex64_t alpha = -1.0;
    pastix_complex64_t *B_tmp;

    double norm_LR, norm_dense, norm_diff;

    int mode     = 0;
    double rcond = (double) (pastix_imin(MA, NA));
    double dmax  = 1.0;

    pastix_int_t minMN = pastix_imin(MA, NA);

    if (MA + offx > MB || NA + offy > NB){
        printf("B receives a contribution from A\n");
        printf("MA + offx <= MB AND NA + offy <= NB\n");
        return -3;
    }

    MALLOC_INTERN(A, MA * NA, pastix_complex64_t);
    MALLOC_INTERN(B, MB * NB, pastix_complex64_t);
    MALLOC_INTERN(C, MB * NB, pastix_complex64_t);

    MALLOC_INTERN(SA, pastix_imin(MA, NA), double);
    MALLOC_INTERN(SB, pastix_imin(MB, NB), double);


    if ((!A)||(!B)||(!C)){
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work('e');

    printf("\n");
    printf("------ TESTS FOR PASTIX ZRRADD ROUTINE -------  \n");
    printf("            Size of the Matrix A %8ld by %8ld. Rank is %8d\n", MA, NA, RA);
    printf("            Size of the Matrix B %8ld by %8ld. Rank is %8d\n", MB, NB, RB);
    printf("\n");
    printf(" The matrix A and B are generated through zlatms routine, with singular values depending on rank and precision.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Low-rank tolerance is set to %e\n", tolerance);
    printf(" Computational tests pass if scaled residuals are less than 10.\n");

    /*----------------------------------------------------------
     *  TESTING ZRRADD
     */

    double alphaA = exp(log(tolerance) / RA);
    double alphaB = exp(log(tolerance) / RB);

    /* Initialize A, B */
    if (mode == 0){
        pastix_int_t i;
        SA[0] = 1.;
        SB[0] = 1.;
        for (i=1; i<minMN; i++){
            SA[i] = SA[i-1] * alphaA;
            SB[i] = SB[i-1] * alphaB;
        }

        /* printf("Singular values:\n"); */
        /* for (i=0; i<minMN; i++){ */
        /*     printf("INDEX %5ld: %.3g %.3g\n", i, SA[i], SB[i]); */
        /* } */
        /* printf("\n"); */
    }

    /* printf("SA_%ld = %.3g\n", RA, SA[RA]); */
    /* printf("SB_%ld = %.3g\n", RB, SB[RB]); */

    /* Initialize A */
    pastix_complex64_t *work;
    MALLOC_INTERN(work, 3*pastix_imax(MA, NA), pastix_complex64_t);

    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, MA, NA,
                         'U', ISEED,
                         'N', SA, mode, rcond,
                         dmax, MA, NA,
                         'N', A, LDA, work );

    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, MA, NA,
                         'U', ISEED,
                         'N', SB, mode, rcond,
                         dmax, MB, NB,
                         'N', B, LDB, work );

    core_zge2lr_RRQR( tolerance,
                      MA, NA,
                      A, LDA,
                      &LR_A );

    core_zge2lr_RRQR( tolerance,
                      MB, NB,
                      B, LDB,
                      &LR_B );

    printf(" The rank of A is %d B is %d\n", LR_A.rk, LR_B.rk);

    /* Add A and B in their LR format */
    core_zrradd_RRQR(tolerance, CblasNoTrans, alpha,
                     MA, NA, &LR_A,
                     MB, NB, &LR_B,
                     offx, offy);

    /* Build uncompressed LR+LR matrix */
    core_zlr2ge( MB, NB,
                 &LR_B,
                 C, LDB );

    /* Add A and B with a classical dense operation */
    B_tmp = B + offx + LDB * offy;
    core_zgeadd( CblasNoTrans, MA, NA,
                 alpha, A, LDA,
                 1.0, B_tmp, LDB );

    printf(" The rank of A+B is %d\n", LR_B.rk);

    /* Compute norm of dense and LR matrices */
    norm_LR    = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MB, NB,
                                      C, LDB, NULL );
    norm_dense = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MB, NB,
                                      B, LDB, NULL );

    core_zgeadd( PastixNoTrans, MB, NB,
                 -1., B, LDB,
                  1., C, LDB );

    norm_diff    = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MB, NB,
                                        C, LDB, NULL );

    res = norm_diff / ( tolerance * norm_dense );

    printf(" ||full(A)+full(B)|| = %e, ||comp(A)+comp(B)|| = %e\n", norm_dense, norm_LR);
    printf(" ||(full(A)+full(B)) - (comp(A)+comp(B))|| = %e\n", norm_diff);
    printf(" res = %e\n", res);

    if (res < 10){
        printf("***************************************************\n");
        printf(" ---- TESTING ZRRADD...................... PASSED !\n");
        printf("***************************************************\n");
    }
    else{
        printf("***************************************************\n");
        printf(" ---- TESTING ZRRADD.................. SUSPICIOUS !\n");
        printf("***************************************************\n");
    }

    memFree_null(A);
    memFree_null(B);
    memFree_null(C);
    return 1;
}
