/**
 *
 * @file example_zgelqs.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @brief Example for solving undetermined linear systems
 *
 * @version 2.5.0
 * @author Bilel Hadri
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

#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif

int check_solution(int, int, int, PLASMA_Complex64_t*, int, PLASMA_Complex64_t*, PLASMA_Complex64_t*, int);

int IONE=1;
int ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

int main ()
{

    int cores = 2;
    int M     = 10;
    int N     = 15;
    int LDA   = 10;
    int NRHS  = 5;
    int LDB   = 15;

    int info;
    int info_solution;
    int i,j;
    int LDAxN = LDA*N;
    int LDBxNRHS = LDB*NRHS;

    PLASMA_Complex64_t *A1 = (PLASMA_Complex64_t *)malloc(LDA*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *A2 = (PLASMA_Complex64_t *)malloc(LDA*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *B1 = (PLASMA_Complex64_t *)malloc(LDB*NRHS*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *B2 = (PLASMA_Complex64_t *)malloc(LDB*NRHS*sizeof(PLASMA_Complex64_t));
    PLASMA_desc *T;

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)||(!B1)||(!B2)){
        printf("Out of Memory \n ");
        return EXIT_SUCCESS;
    }

    /* Plasma Initialization */
    PLASMA_Init(cores);
    printf("-- PLASMA is initialized to run on %d cores. \n",cores);

    /* Allocate T */
    PLASMA_Alloc_Workspace_zgeqrf(M, N, &T);

    /* Initialize A1 and A2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDAxN, A1);
    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            A2[LDA*j+i] = A1[LDA*j+i] ;

    /* Initialize B1 and B2 */
    LAPACKE_zlarnv_work(IONE, ISEED, LDBxNRHS, B1);
    for (i = 0; i < M; i++)
        for (j = 0; j < NRHS; j++)
             B2[LDB*j+i] = B1[LDB*j+i] ;

    /* Factorization QR of the matrix A2 */
    info = PLASMA_zgelqf(M, N, A2, LDA, T);

    /* Solve the problem */
    info = PLASMA_zgelqs(M, N, NRHS, A2, LDA, T, B2, LDB);

    /* Check the solution */
    info_solution = check_solution(M, N, NRHS, A1, LDA, B1, B2, LDB);

    if ((info_solution != 0)|(info != 0))
       printf("-- Error in ZGELQS example ! \n");
    else
       printf("-- Run of ZGELQS example successful ! \n");

    free(A1); free(A2); free(B1); free(B2); free(T);

    PLASMA_Finalize();

    return EXIT_SUCCESS;
}

/*--------------------------------------------------------------
 * Check the solution
 */

int check_solution(int M, int N, int NRHS, PLASMA_Complex64_t *A1, int LDA, PLASMA_Complex64_t *B1, PLASMA_Complex64_t *B2, int LDB)
{
    int info_solution;
    double Rnorm, Anorm, Xnorm, Bnorm;
    PLASMA_Complex64_t alpha, beta;
    double *work = (double *)malloc(max(M, N)* sizeof(double));
    double eps;

    eps = LAPACKE_dlamch_work('e');

    alpha = 1.0;
    beta  = -1.0;

    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, N, A1, LDA, work);
    Xnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, NRHS, B2, LDB, work);
    Bnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B1, LDB, work);

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, NRHS, N, CBLAS_SADDR(alpha), A1, LDA, B2, LDB, CBLAS_SADDR(beta), B1, LDB);

    if (M >= N) {
       PLASMA_Complex64_t *Residual = (PLASMA_Complex64_t *)malloc(M*NRHS*sizeof(PLASMA_Complex64_t));
       memset((void*)Residual, 0, M*NRHS*sizeof(PLASMA_Complex64_t));
       cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, NRHS, M, CBLAS_SADDR(alpha), A1, LDA, B1, LDB, CBLAS_SADDR(beta), Residual, M);
       Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), M, NRHS, Residual, M, work);
       free(Residual);
    }
    else {
       PLASMA_Complex64_t *Residual = (PLASMA_Complex64_t *)malloc(N*NRHS*sizeof(PLASMA_Complex64_t));
       memset((void*)Residual, 0, N*NRHS*sizeof(PLASMA_Complex64_t));
       cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, N, NRHS, M, CBLAS_SADDR(alpha), A1, LDA, B1, LDB, CBLAS_SADDR(beta), Residual, N);
       Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, Residual, N, work);
       free(Residual);
    }

    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||)_oo.N.eps) = %e \n",Rnorm/((Anorm*Xnorm+Bnorm)*N*eps));

    if (isnan(Rnorm / ((Anorm * Xnorm + Bnorm) * N * eps)) || (Rnorm / ((Anorm * Xnorm + Bnorm) * N * eps) > 10.0) ) {
         printf("-- The solution is suspicious ! \n");
         info_solution = 1;
    }
    else {
         printf("-- The solution is CORRECT ! \n");
         info_solution= 0 ;
    }

    free(work);

    return info_solution;
}
