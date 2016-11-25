/**
 *  @file: testing_zrradd.c
 *
 *  testing LR + LR operation
 *
 * @precisions normal z -> c d s
 *
 */
#include "testing_zmain.h"

static pastix_complex64_t mzone = -1.;
static pastix_complex64_t zone  =  1.;
static pastix_complex64_t zzero =  0.;

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

    double eps, res_SVD, res_RRQR;
    pastix_complex64_t *A, *B, *C;
    pastix_lrblock_t    LR_A, LR_B;

    double norm_LR_RRQR, norm_LR_SVD, norm_dense, norm_diff_RRQR, norm_diff_SVD;

    pastix_int_t minMN = pastix_imin(MA, NA);
    int mode           = 0;
    double rcond       = (double) minMN;
    double dmax        = 1.0;

    pastix_complex64_t *work;
    double *S;

    MALLOC_INTERN(A,  MA * LDA, pastix_complex64_t);
    MALLOC_INTERN(B,  MA * LDA, pastix_complex64_t);
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
        for (i=1; i<minMN/4; i++){
            S[0] *= 3.5;
        }
        for (i=1; i<minMN; i++){
            S[i] = S[i-1] / 2.1;
            /* printf("S %.3g\n", S[i]); */
        }
        for (i=1; i<minMN; i++){
            S[i] /= S[0];
        }
        S[0] = 1.;
        /* for (i=minMN/4; i<minMN; i++){ */
        /*     S[i] = 0.; */
        /* } */
    }

    /* memset( A, 0xdeadbeef, LDA * NA * sizeof(pastix_complex64_t) ); */

    /* Initialize A */
    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, MA, NA,
                         'U', ISEED,
                         'N', S, mode, rcond,
                         dmax, MA, NA,
                         'N', A, LDA, work );

    norm_dense = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MA, NA,
                                      A, LDA, NULL );

    if (0){
        pastix_int_t rank = 100; //minMN / 4;
        pastix_complex64_t *u = malloc(MA * rank * sizeof(pastix_complex64_t));
        pastix_complex64_t *v = malloc(NA * rank * sizeof(pastix_complex64_t));

        pastix_int_t i;
        for (i=0; i<MA*rank; i++){
            u[i] = (float)rand()/(float)(RAND_MAX);
        }
        for (i=0; i<NA*rank; i++){
            v[i] = (float)rand()/(float)(RAND_MAX);
        }

        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                    MA, NA, rank,
                    CBLAS_SADDR(zone),  u, MA,
                    v, rank,
                    CBLAS_SADDR(zzero), A, LDA);

        norm_dense = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MA, NA,
                                          A, LDA, NULL );

        for (i=0; i<MA*NA; i++){
            A[i] /= norm_dense;
        }

        norm_dense = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MA, NA,
                                          A, LDA, NULL );
    }

    /* Initialize Low-Rank structure */
    LR_A.rk    = -1;
    LR_A.rkmax = -1;
    LR_A.u     = NULL;
    LR_A.v     = NULL;

    LR_B.rk    = -1;
    LR_B.rkmax = -1;
    LR_B.u     = NULL;
    LR_B.v     = NULL;

    /* Compress and then uncompress  */
    core_zge2lr_RRQR( tolerance,
                      MA, NA,
                      A, LDA,
                      &LR_A );

    core_zge2lr_SVD( tolerance,
                      MA, NA,
                      A, LDA,
                      &LR_B );

    core_zlr2ge( MA, NA,
                 &LR_A,
                 C, LDA );

    core_zlr2ge( MA, NA,
                 &LR_B,
                 B, LDA );

    printf(" The rank of A is: RRQR %d SVD %d\n", LR_A.rk, LR_B.rk);

    /* Compute norm of dense and LR matrices */
    norm_LR_RRQR = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MA, NA,
                                        C, LDA, NULL );
    norm_LR_SVD  = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MA, NA,
                                        B, LDA, NULL );

    core_zgeadd( PastixNoTrans, MA, NA,
                 -1., A, LDA,
                  1., C, LDA );

    core_zgeadd( PastixNoTrans, MA, NA,
                 -1., A, LDA,
                  1., B, LDA );

    norm_diff_RRQR = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MA, NA,
                                          C, LDA, NULL );

    norm_diff_SVD = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MA, NA,
                                         B, LDA, NULL );


    printf(" ||full(A)|| = %e, ||comp(A_RRQR)|| = %e, ||comp(A_SVD)|| = %e\n",
           norm_dense, norm_LR_RRQR, norm_LR_SVD);
    printf(" ||full(A) - comp(A_RRQR)|| = %e, ||full(A) - comp(A_SVD)|| = %e\n\n",
           norm_diff_RRQR, norm_diff_SVD);
    printf(" ||full(A)|| * sqrt(tol) = %e\n",
           norm_dense * tolerance);


    res_RRQR = norm_diff_RRQR / ( tolerance * norm_dense );
    res_SVD  = norm_diff_SVD  / ( tolerance * norm_dense );
    printf(" res = %e %e\n", res_RRQR, res_SVD);

    if (res_RRQR < 100){
        printf("***************************************************\n");
        printf(" ---- TESTING ZGE2LR / RRQR...................... PASSED !\n");
        printf("***************************************************\n");
    }
    else{
        printf("***************************************************\n");
        printf(" ---- TESTING ZGE2LR / RRQR.................. SUSPICIOUS !\n");
        printf("***************************************************\n");
    }

    if (res_SVD < 100){
        printf("***************************************************\n");
        printf(" ---- TESTING ZGE2LR / SVD...................... PASSED !\n");
        printf("***************************************************\n");
    }
    else{
        printf("***************************************************\n");
        printf(" ---- TESTING ZGE2LR / SVD.................. SUSPICIOUS !\n");
        printf("***************************************************\n");
    }

    /* memFree_null(LR_A.u); */
    /* memFree_null(LR_A.v); */
    /* memFree_null(LR_B.u); */
    /* memFree_null(LR_B.v); */
    /* memFree_null(A); */
    /* memFree_null(B); */
    /* memFree_null(C); */
    /* memFree_null(S); */
    /* memFree_null(work); */
    return 1;
}
