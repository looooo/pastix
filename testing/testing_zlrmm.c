/**
 *  @file: testing_zrradd.c
 *
 *  testing LR + LR operation
 *
 * @precisions normal z -> c d s
 *
 */
#include <testing_zmain.h>

int testing_zlrmm(int argc, char **argv)
{
    /* Check for number of arguments*/
    if ( argc != 8) {
        USAGE("ZLRMM", "tol MA NB K MC NC offx offy",
              "   - tol    : tolerance for SVD compression\n"
              "   - MA     : number of rows of matrix A\n"
              "   - NB     : number of columns of matrix B\n"
              "   - K      : number of columns of matrix A and number of rows of matrix B\n"
              "   - MC     : number of rows of matrix C\n"
              "   - NC     : number of columns of matrix C\n"
              "   - offx   : first row of AB with respect to C\n"
              "   - offy   : first column of AB with respect to C\n");
        return -1;
    }

    double tolerance  = atof(argv[0]);
    pastix_int_t MA   = atoi(argv[1]);
    pastix_int_t NB   = atoi(argv[2]);
    pastix_int_t K    = atoi(argv[3]);
    pastix_int_t MC   = atoi(argv[4]);
    pastix_int_t NC   = atoi(argv[5]);
    pastix_int_t offx = atoi(argv[6]);
    pastix_int_t offy = atoi(argv[7]);

    if (MA + offx > MC || NB + offy > NC){
        printf("C receives a contribution from AB\n");
        printf("MA + offx <= MC AND NB + offy <= NC\n");
        return -3;
    }

    double eps, res;
    pastix_complex64_t *A, *B, *C, *C2;
    pastix_lrblock_t    LR_A, LR_B, LR_C;

    double norm_LR, norm_dense, norm_diff;

    int mode           = 0;
    double dmax        = 1.0;

    pastix_int_t minMN_A = pastix_imin(MA, K);
    pastix_int_t minMN_B = pastix_imin(NB, K);
    pastix_int_t minMN_C = pastix_imin(MC, NC);

    pastix_complex64_t *work;
    double *S_A, *S_B, *S_C;

    pastix_int_t LDA = MA;
    pastix_int_t LDB = NB;
    pastix_int_t LDC = MC;

    int transA = PastixNoTrans;
    int transB = PastixTrans;

    MALLOC_INTERN(A , MA * K , pastix_complex64_t);
    MALLOC_INTERN(B , NB * K , pastix_complex64_t);
    MALLOC_INTERN(C , MC * NC, pastix_complex64_t);
    MALLOC_INTERN(C2, MC * NC, pastix_complex64_t);

    MALLOC_INTERN(S_A, minMN_A, double);
    MALLOC_INTERN(S_B, minMN_B, double);
    MALLOC_INTERN(S_C, minMN_C, double);

    MALLOC_INTERN(work, 3*pastix_imax(MC, NC), pastix_complex64_t);

    if ((!A)||(!B)||(!C)||(!C2)||(!S_A)||(!S_B)||(!S_C)||(!work)){
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work('e');

    printf("\n");
    printf("------ TESTS FOR PASTIX ZZLRMM ROUTINE -------  \n");
    printf("            Size of the Matrix A %8ld by %8ld\n", MA, K);
    printf("            Size of the Matrix B %8ld by %8ld\n", K, NB);
    printf("            Size of the Matrix C %8ld by %8ld\n", MC, NC);
    printf("\n");
    printf(" The matrix A and B are randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Low-rank tolerance is set to %e\n", tolerance);
    printf(" Computational tests pass if scaled residuals are less than 10.\n");

    /*----------------------------------------------------------
     *  TESTING ZZLRMM
     */
    if (mode == 0){
        pastix_int_t i;
        S_A[0] = 1;
        S_B[0] = 1.2;
        S_C[0] = 1.125;
        for (i=1; i<minMN_A; i++){
            S_A[i] = S_A[i-1] / 1.1;
        }
        for (i=1; i<minMN_B; i++){
            S_B[i] = S_B[i-1] / 1.05;
        }
        for (i=1; i<minMN_C; i++){
            S_C[i] = S_C[i-1] / 1.12;
        }
    }

    /* Initialize A and B*/
    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, MA, K,
                         'U', ISEED,
                         'N', S_A, mode, minMN_A,
                         dmax, MA, K,
                         'N', A, LDA, work );

    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, NB, K,
                         'U', ISEED,
                         'N', S_B, mode, minMN_B,
                         dmax, NB, K,
                         'N', B, LDB, work );

    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, MC, NC,
                         'U', ISEED,
                         'N', S_C, mode, minMN_C,
                         dmax, MC, NC,
                         'N', C, LDC, work );

    /* Initialize Low-Rank structure */
    LR_A.rk    = -1;
    LR_A.rkmax = -1;
    LR_A.u     = NULL;
    LR_A.v     = NULL;
    LR_B.rk    = -1;
    LR_B.rkmax = -1;
    LR_B.u     = NULL;
    LR_B.v     = NULL;
    LR_C.rk    = -1;
    LR_C.rkmax = -1;
    LR_C.u     = NULL;
    LR_C.v     = NULL;

    /* Compress matrices */
    core_zge2lr( tolerance,
                 MA, K,
                 A, LDA,
                 &LR_A );

    core_zge2lr( tolerance,
                 NB, K,
                 B, LDB,
                 &LR_B );

    core_zge2lr( tolerance,
                 MC, NC,
                 C, LDC,
                 &LR_C );

    printf(" The rank of A is %d B is %d C is %d\n", LR_A.rk, LR_B.rk, LR_C.rk);

    /* Compute A*B + C */
    double alpha = -1.0;
    double beta  = 1.0;
    core_zlrmm( tolerance, transA, transB,
                MA, NB, K,
                MC, NC,
                offx, offy,
                alpha, &LR_A, &LR_B,
                beta , &LR_C,
                NULL, -1 );

    printf(" New rank of C is %d\n", LR_C.rk);

    core_zlr2ge( MC, NC,
                 &LR_C,
                 C2, LDC );

    pastix_complex64_t *Cptr = C + offy * LDC + offx;

    cblas_zgemm( CblasColMajor, transA, transB,
                 MA, NB, K,
                 CBLAS_SADDR(alpha), A, LDA,
                                     B, LDB,
                 CBLAS_SADDR(beta),  Cptr, LDC );


    /* Compute norm of dense and LR matrices */
    norm_LR    = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MC, NC,
                                      C2, LDC, NULL );
    norm_dense = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'f', MC, NC,
                                      C,  LDC, NULL );

    core_zgeadd( PastixNoTrans, MC, NC,
                 -1., C2, LDC,
                  1., C,  LDC );

    norm_diff    = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'm', MC, NC,
                                        C, LDC, NULL );

    res = norm_diff / ( eps * norm_dense );

    printf(" ||full(A)*full(B)+full(C)|| = %e, ||full(comp(A)*comp(B)+comp(C))|| = %e\n", norm_dense, norm_LR);
    printf(" ||full(A*B+C) - comp(A*B+C)|| = %e\n", norm_diff);
    printf(" res = %e\n", res);

    if (res < 10){
        printf("***************************************************\n");
        printf(" ---- TESTING ZZLRMM...................... PASSED !\n");
        printf("***************************************************\n");
    }
    else{
        printf("***************************************************\n");
        printf(" ---- TESTING ZZLRMM.................. SUSPICIOUS !\n");
        printf("***************************************************\n");
    }

    memFree_null(LR_A.u);
    memFree_null(LR_A.v);
    memFree_null(LR_B.u);
    memFree_null(LR_B.v);
    memFree_null(LR_C.u);
    memFree_null(LR_C.v);
    memFree_null(A);
    memFree_null(B);
    memFree_null(C);
    memFree_null(C2);
    memFree_null(S_A);
    memFree_null(S_B);
    memFree_null(S_C);
    memFree_null(work);
    return 1;
}
