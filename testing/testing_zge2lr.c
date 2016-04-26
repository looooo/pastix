/**
 *  @file: testing_zrradd.c
 *
 *  testing LR + LR operation
 *
 * @precisions normal z -> c d s
 *
 */
#include <testing_zmain.h>

int testing_zge2lr(int argc, char **argv)
{
    /* Check for number of arguments*/
    if ( argc != 4) {
        USAGE("GE2LR", "tol MA NA LDA",
              "   - tol    : tolerance for SVD compression\n"
              "   - MA     : number of rows of matrices A\n"
              "   - NA     : number of columns of matrices A\n"
              "   - LDA    : leading dimension of matrix A\n");
        return -1;
    }

    double tolerance = atof(argv[0]);
    pastix_int_t MA  = atoi(argv[1]);
    pastix_int_t NA  = atoi(argv[2]);
    pastix_int_t LDA = atoi(argv[3]);

    double eps, res;
    pastix_complex64_t *A, *C;
    pastix_lrblock_t    LR_A;

    double norm_LR, norm_dense, norm_diff;

    pastix_int_t minMN = pastix_imin(MA, NA);
    int mode           = 0;
    double rcond       = (double) minMN;
    double dmax        = 1.0;

    pastix_complex64_t *work;
    double *S;

    MALLOC_INTERN(A,  MA * LDA, pastix_complex64_t);
    MALLOC_INTERN(C,  MA * LDA, pastix_complex64_t);

    MALLOC_INTERN(S, minMN, double);
    MALLOC_INTERN(work, 3*pastix_imax(MA, NA), pastix_complex64_t);

    if ((!A)||(!C)||(!S)||(!work)){
        printf("Out of Memory \n ");
        return -2;
    }

    if (LDA < MA || LDA < NA){
        printf("Invalid LDA parameter\n");
        return -3;
    }

    eps = LAPACKE_dlamch_work('e');

    printf("\n");
    printf("------ TESTS FOR PASTIX ZGE2LR ROUTINE -------  \n");
    printf("            Size of the Matrix A %8ld by %8ld. LDA is %8ld\n", MA, NA, LDA);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Low-rank tolerance is set to %e\n", tolerance);
    printf(" Computational tests pass if scaled residuals are less than 10.\n");

    /*----------------------------------------------------------
     *  TESTING ZGE2LR
     */
    if (mode == 0){
        pastix_int_t i;
        S[0] = 1;
        for (i=1; i<minMN; i++){
            S[i] = S[i-1] / 1.1;
        }
    }

    /* Initialize A */
    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, MA, NA,
                         'U', ISEED,
                         'N', S, mode, rcond,
                         dmax, MA, NA,
                         'N', A, LDA, work );

    /* Initialize Low-Rank structure */
    LR_A.rk    = -1;
    LR_A.rkmax = -1;
    LR_A.u     = NULL;
    LR_A.v     = NULL;

    /* Compress and then uncompress  */
    core_zge2lr( tolerance,
                 MA, NA,
                 A, LDA,
                 &LR_A );

    core_zlr2ge( MA, NA,
                 &LR_A,
                 C, LDA );

    printf(" The rank of A is %d\n", LR_A.rk);

    /* Compute norm of dense and LR matrices */
    norm_LR    = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MA, NA,
                                      C, LDA, NULL );
    norm_dense = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MA, NA,
                                      A, LDA, NULL );

    core_zgeadd( PastixNoTrans, MA, NA,
                 -1., A, LDA,
                  1., C, LDA );

    norm_diff    = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'm', MA, NA,
                                        C, LDA, NULL );

    res = norm_diff / ( eps * norm_dense );

    printf(" ||full(A)|| = %e, ||comp(A)|| = %e\n", norm_dense, norm_LR);
    printf(" ||full(A) - comp(A)|| = %e\n", norm_diff);
    printf(" res = %e\n", res);

    if (res < 10){
        printf("***************************************************\n");
        printf(" ---- TESTING ZGE2LR...................... PASSED !\n");
        printf("***************************************************\n");
    }
    else{
        printf("***************************************************\n");
        printf(" ---- TESTING ZGE2LR.................. SUSPICIOUS !\n");
        printf("***************************************************\n");
    }

    memFree_null(LR_A.u);
    memFree_null(LR_A.v);
    memFree_null(A);
    memFree_null(C);
    memFree_null(S);
    memFree_null(work);
    return 1;
}
