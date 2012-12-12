/**
 *
 * @file testing_zgesvd.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Hatem Ltaief
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <plasma.h>
#include <cblas.h>
#include <lapacke.h>
#include <core_blas.h>
#include "testing_zmain.h"

#undef REAL
#define COMPLEX

static int check_orthogonality(int, int, int, PLASMA_Complex64_t*, int, double);
static int check_reduction(int, int, PLASMA_Complex64_t*, PLASMA_Complex64_t*, int, PLASMA_Complex64_t*, int, PLASMA_Complex64_t*, int, double);
static int check_solution(int, double*, double*, double);

int testing_zgesvd(int argc, char **argv)
{
    int tree = 0;

    if ( argc < 1 ){
        goto usage;
    } else {
        tree = atoi(argv[0]);
    }

    /* Check for number of arguments*/
    if ( ((tree == 0) && (argc != 4)) ||
         ((tree != 0) && (argc != 5)) ){
      usage:
        USAGE("GESVD", "MODE M N LDA [RH]",
              "   - MODE : 0: flat, 1: tree (RH needed)\n"
              "   - M    : number of rows of the matrix A\n"
              "   - N    : number of columns of the matrix A\n"
              "   - LDA  : leading dimension of the matrix A\n"
              "   - RH   : Size of each subdomains\n");
        return -1;
    }

    int M   = atoi(argv[1]);
    int N   = atoi(argv[2]);
    int LDA = atoi(argv[3]);
    int rh;
    if ( tree ) {
        rh = atoi(argv[4]);
        PLASMA_Set(PLASMA_HOUSEHOLDER_MODE, PLASMA_TREE_HOUSEHOLDER);
        PLASMA_Set(PLASMA_HOUSEHOLDER_SIZE, rh);
    }

    if (LDA < M){
        printf("LDA should be >= M !\n");
        return -1;
    }

    double eps  = LAPACKE_dlamch_work('e');
    double dmax = 1.0;
    PLASMA_enum vecu  = PlasmaNoVec;
    PLASMA_enum vecvt = PlasmaNoVec;
    int info_orthou    = 0;
    int info_orthovt   = 0;
    int info_solution  = 0;
    int info_reduction = 0;
    int minMN = min(M, N);
    int mode  = 4;
    double rcond = (double) minMN;

    PLASMA_Complex64_t *A1   = (PLASMA_Complex64_t *)malloc(LDA*N*sizeof(PLASMA_Complex64_t));
    double *S1               = (double *)            malloc(minMN*sizeof(double));
    double *S2               = (double *)            malloc(minMN*sizeof(double));
    PLASMA_Complex64_t *work = (PLASMA_Complex64_t *)malloc(3*max(M, N)* sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *A2 = NULL;
    PLASMA_Complex64_t *U  = NULL;
    PLASMA_Complex64_t *VT = NULL;
    PLASMA_desc *T;

    /* Check if unable to allocate memory */
    if ( (!A1) || (!S1) || (!S2) || (!work) ) {
        printf("Out of Memory \n ");
        return -2;
    }

    /*
    PLASMA_Disable(PLASMA_AUTOTUNING);
    PLASMA_Set(PLASMA_TILE_SIZE, 120);
    PLASMA_Set(PLASMA_INNER_BLOCK_SIZE, 40);
    */

    PLASMA_Enable(PLASMA_WARNINGS);
    PLASMA_Enable(PLASMA_ERRORS);
    PLASMA_Alloc_Workspace_zgesvd(M, N, &T);

    /*----------------------------------------------------------
    *  TESTING ZGESVD
    */
    /* Initialize A1 */
    LAPACKE_zlatms_work( LAPACK_COL_MAJOR, M, N,
                         lapack_const(PlasmaDistUniform), ISEED,
                         lapack_const(PlasmaNonsymPosv), S1, mode, rcond,
                         dmax, M, N,
                         lapack_const(PlasmaNoPacking), A1, LDA, work );
    free(work);

    /* Copy A1 for check */
    if ( (vecu == PlasmaVec) && (vecvt == PlasmaVec) ) {
        A2 = (PLASMA_Complex64_t *)malloc(LDA*N*sizeof(PLASMA_Complex64_t));
        LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, 'A', M, N, A1, LDA, A2, LDA);
    }
    if ( vecu == PlasmaVec ) {
        U = (PLASMA_Complex64_t *)malloc(M*M*sizeof(PLASMA_Complex64_t));
        LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', M, M, 0., 1., U, M);
    }
    if ( vecvt == PlasmaVec ) {
        VT = (PLASMA_Complex64_t *)malloc(N*N*sizeof(PLASMA_Complex64_t));
        LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', N, N, 0., 1., VT, N);
    }
 
    /* PLASMA ZGESVD */
    PLASMA_zgesvd(vecu, vecvt, M, N, A1, LDA, S2, U, M, VT, N, T);

    if (getenv("PLASMA_TESTING_VERBOSE"))
    {
        int i;
        printf("Singular values original\n");
        for (i = 0; i < min(N,25); i++){
            printf("%f ", S1[i]);
        }
        printf("\n");

        printf("Singular values computed\n");
        for (i = 0; i < min(N,25); i++){
            printf("%f ", S2[i]);
        }
        printf("\n");
    }

    printf("\n");
    printf("------ TESTS FOR PLASMA ZGESVD ROUTINE -------  \n");
    printf("        Size of the Matrix %d by %d\n", M, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n",eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");


    /* Check the orthogonality, reduction and the singular values */
    if ( vecu == PlasmaVec )
        info_orthou = check_orthogonality(PlasmaLeft, M, M, U, M, eps);

    if ( vecvt == PlasmaVec )
        info_orthovt = check_orthogonality(PlasmaRight, N, N, VT, N, eps);

    /* 
     * WARNING: For now, Q is associated to Band tridiagonal reduction and 
     * not to the final tridiagonal reduction, so we can not call the check
     * Need to accumulate the orthogonal transformations
     * during the bulge chasing to be able to perform the next test!
     */
    if ( (vecu == PlasmaVec) && (vecvt == PlasmaVec) && 0 )
        info_reduction = check_reduction(M, N, A2, A1, LDA, U, M, VT, N, eps);

    info_solution = check_solution(minMN, S1, S2, eps);

    if ( (info_solution == 0) & (info_orthou == 0) & 
         (info_orthovt == 0) & (info_reduction == 0) ) {
        if (M >= N) {
           printf("***************************************************\n");
           printf(" ---- TESTING ZGESVD .. M >= N ........... PASSED !\n");
           printf("***************************************************\n");
        }
        else {
           printf("***************************************************\n");
           printf(" ---- TESTING ZGESVD .. M < N ............ PASSED !\n");
           printf("***************************************************\n");
        }
    }
    else {
        if (M >= N) {
           printf("************************************************\n");
           printf(" - TESTING ZGESVD .. M >= N .. FAILED !\n");
           printf("************************************************\n");
        }
        else {
           printf("************************************************\n");
           printf(" - TESTING ZGESVD .. M < N .. FAILED !\n");
           printf("************************************************\n");
        }
    }

    if ( A2 != NULL ) free(A2); 
    if ( U  != NULL ) free(U); 
    if ( VT != NULL ) free(VT); 
    free(A1); free(S1); free(S2);
    PLASMA_Dealloc_Handle_Tile(&T);

    return 0;
}

/*-------------------------------------------------------------------
 * Check the orthogonality of Q
 */
static int check_orthogonality(int side, int M, int N, PLASMA_Complex64_t *Q, int LDQ, double eps)
{
    double  alpha =  1.0;
    double  beta  = -1.0;
    double  normQ, result;
    int     info_ortho;
    int     minMN = min(M, N);
    double *work = (double *)malloc(minMN*sizeof(double));

    /* Build the idendity matrix */
    PLASMA_Complex64_t *Id = (PLASMA_Complex64_t *) malloc(minMN*minMN*sizeof(PLASMA_Complex64_t));
    LAPACKE_zlaset_work(LAPACK_COL_MAJOR, 'A', minMN, minMN, 0., 1., Id, minMN);

    /* Perform Id - Q'Q */
    if (M >= N)
        cblas_zherk(CblasColMajor, CblasUpper, CblasConjTrans, N, M, alpha, Q, LDQ, beta, Id, N);
    else
        cblas_zherk(CblasColMajor, CblasUpper, CblasNoTrans,   M, N, alpha, Q, LDQ, beta, Id, M);

    normQ = LAPACKE_zlansy_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), 'U', minMN, Id, minMN, work);

    if (getenv("PLASMA_TESTING_VERBOSE"))
        printf( "||Q||_oo=%f\n", normQ );
    
    printf("================================\n");
    result = normQ / (M * eps);
    if (side == PlasmaLeft)
    {
        printf("Checking the orthogonality of U (Left singular vectors)\n");
        printf("||Id-U'*U||_oo / (N*eps) = %e \n", result);
    }
    else 
    {
        printf("Checking the orthogonality of VT (Right singular vectors)\n");
        printf("||Id-VT'*VT||_oo / (N*eps) = %e \n", result);
    }

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- Orthogonality is suspicious ! \n");
        info_ortho = 1;
    }
    else {
        printf("-- Orthogonality is CORRECT ! \n");
        info_ortho = 0;
    }
    
    free(work); free(Id);
    
    return info_ortho;
}

/*------------------------------------------------------------
 *  Check the bidiagonal reduction
 */
static int check_reduction(int M, int N, PLASMA_Complex64_t *A1, PLASMA_Complex64_t *A2, int LDA, 
                           PLASMA_Complex64_t *U, int LDU, PLASMA_Complex64_t *VT, int LDVT, double eps )
{
    PLASMA_Complex64_t alpha = 1.0;
    PLASMA_Complex64_t beta  = 0.0;
    double Anorm, Rnorm, result;
    int info_reduction;
    int i, j;
    int maxMN = max(M, N);

    PLASMA_Complex64_t *Aorig    = (PLASMA_Complex64_t *)malloc(M*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *Residual = (PLASMA_Complex64_t *)malloc(M*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *B        = (PLASMA_Complex64_t *)malloc(M*N*sizeof(PLASMA_Complex64_t));
    double *work = (double *)malloc(N*sizeof(double));

    memset((void*)B, 0, M*N*sizeof(PLASMA_Complex64_t));

    /* Rebuild the B */
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, PlasmaUpperLower, M, N, A2, LDA, B, M);
    memset((void*)Aorig, 0, M*N*sizeof(PLASMA_Complex64_t));

    /* Compute Aorig = U * B * VT */
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, M, CBLAS_SADDR(alpha), U,     LDU, B,  M,    CBLAS_SADDR(beta), Aorig, M);
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, N, N, CBLAS_SADDR(alpha), Aorig, M,   VT, LDVT, CBLAS_SADDR(beta), B,     M);
    
    /* Compute the Residual */
    for (i = 0; i < M; i++)
        for (j = 0 ; j < N; j++)
            Residual[j*M+i] = A1[j*LDA+i] - B[j*M+i];
    
    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, Residual, M,   work);
    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, A2,       LDA, work);
    
   if (getenv("PLASMA_TESTING_VERBOSE"))
        printf( "||A||_oo=%f\n||A - U*B*VT||_oo=%e\n", Anorm, Rnorm );

    result = Rnorm / ( Anorm * maxMN * eps);
    printf("============\n");
    printf("Checking the bidiagonal reduction \n");
    printf("-- ||A-U*B*VT||_oo/(||A||_oo.max(M,N).eps) = %e \n", result);

    if ( isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- Reduction is suspicious ! \n");
        info_reduction = 1;
    }
    else {
        printf("-- Reduction is CORRECT ! \n");
        info_reduction = 0;
    }

    free(Aorig); free(Residual); free(B); free(work);

    return info_reduction;
}

static int check_solution(int N, double *E1, double *E2, double eps)
{
    int info_solution, i;
    double resid;
    double maxtmp;
    double maxel = fabs( fabs(E1[0]) - fabs(E2[0]) );
    double maxeig = max( fabs(E1[0]), fabs(E2[0]) );
    for (i = 1; i < N; i++){
        resid   = fabs(fabs(E1[i])-fabs(E2[i]));
        maxtmp  = max(fabs(E1[i]), fabs(E2[i]));

        /* Update */
        maxeig = max(maxtmp, maxeig);
        maxel  = max(resid,  maxel );
    }

    maxel = maxel / (maxeig * N * eps);
    printf(" ======================================================\n");
    printf(" Checking the singular values of A\n");
    printf("    | D - eigcomputed | / (|D| * N * eps) : %15.3E \n",  maxel );

    if ( isnan(maxel) || isinf(maxel) || (maxel > 100) ) {
        printf("-- The singular values are suspicious ! \n");
        info_solution = 1;
    }
    else{
        printf("-- The singular values are CORRECT ! \n");
        info_solution = 0;
    }

    return info_solution;
}
