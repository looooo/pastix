/**
 *  @file: testing_zrradd.c
 *
 *  testing LR + LR operation
 *
 * @precisions normal z -> c d s
 *
 */
#include <testing_zmain.h>

int testing_zgradd(int argc, char **argv)
{
    /* Check for number of arguments*/
    if ( argc != 9) {
        USAGE("GRADD", "tol MA NA RA MB NB RB offx offy",
              "   - tol    : tolerance for SVD compression\n"
              "   - MA     : number of rows of matrices A\n"
              "   - NA     : number of columns of matrices A\n"
              "   - LDA    : rank of matrix A\n"
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
    pastix_int_t LDA  = atoi(argv[3]);
    pastix_int_t MB   = atoi(argv[4]);
    pastix_int_t NB   = atoi(argv[5]);
    pastix_int_t offx = atoi(argv[7]);
    pastix_int_t offy = atoi(argv[8]);

    int RB = atoi(argv[6]);

    double eps, res;
    pastix_complex64_t *uB, *vB;
    pastix_complex64_t *A, *B, *C;
    pastix_lrblock_t    LR_B;

    pastix_int_t LDB = MB;
    pastix_complex64_t alpha = -1.0;
    pastix_complex64_t *B_tmp;

    double norm_LR, norm_dense, norm_diff;

    if (MA + offx > MB || NA + offy > NB){
        printf("B receives a contribution from A\n");
        printf("MA + offx <= MB AND NA + offy <= NB\n");
        return -3;
    }

    MALLOC_INTERN(uB, MB * MB, pastix_complex64_t);
    MALLOC_INTERN(vB, NB * NB, pastix_complex64_t);

    MALLOC_INTERN(A, MA * LDA, pastix_complex64_t);
    MALLOC_INTERN(B, MB * NB, pastix_complex64_t);
    MALLOC_INTERN(C, MB * NB, pastix_complex64_t);

    if ((!uB)||(!vB)||(!A)||(!B)||(!C)){
        printf("Out of Memory \n ");
        return -2;
    }

    eps = LAPACKE_dlamch_work('e');

    printf("\n");
    printf("------ TESTS FOR PASTIX ZGRADD ROUTINE -------  \n");
    printf("            Size of the Matrix A %8ld by %8ld. Leading dimension is %8ld\n", MA, NA, LDA);
    printf("            Size of the Matrix B %8ld by %8ld. Rank is              %8d\n", MB, NB, RB);
    printf("            Offsets are %8ld %8ld\n", offx, offy);
    printf("\n");
    printf(" The matrix A and B are randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Low-rank tolerance is set to %e\n", tolerance);
    printf(" Computational tests pass if scaled residuals are less than 10.\n");

    /*----------------------------------------------------------
     *  TESTING ZGRADD
     */

    /* Initialize A, B */
    LAPACKE_zlarnv_work(IONE, ISEED, MA*LDA, A);
    LAPACKE_zlarnv_work(IONE, ISEED, MB*RB , uB);
    LAPACKE_zlarnv_work(IONE, ISEED, NB*RB , vB);

    LR_B.rk    = RB;
    LR_B.rkmax = RB;
    LR_B.u     = uB;
    LR_B.v     = vB;

    /* Build uncompressed matrix */
    core_zlr2ge( MB, NB,
                 &LR_B,
                 B, LDB );

    /* Add dense A and LR B matrices */
    core_zgradd(tolerance, alpha,
                MA, NA, A, LDA,
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

    norm_diff    = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'm', MB, NB,
                                        C, LDB, NULL );

    res = norm_diff / ( eps * norm_dense );

    printf(" ||full(A)+full(B)|| = %e, ||comp(A)+comp(B)|| = %e\n", norm_dense, norm_LR);
    printf(" ||(full(A)+full(B)) - (comp(A)+comp(B))|| = %e\n", norm_diff);
    printf(" res = %e\n", res);

    if (res < 10){
        printf("***************************************************\n");
        printf(" ---- TESTING ZGRADD...................... PASSED !\n");
        printf("***************************************************\n");
    }
    else{
        printf("***************************************************\n");
        printf(" ---- TESTING ZGRADD.................. SUSPICIOUS !\n");
        printf("***************************************************\n");
    }

    memFree_null(LR_B.u);
    memFree_null(LR_B.v);
    memFree_null(A);
    memFree_null(B);
    memFree_null(C);
    return 1;
}
