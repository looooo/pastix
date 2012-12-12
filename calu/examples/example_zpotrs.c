/**
 *
 * @file example_zpotrs.c
 *
 *  PLASMA testing routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @brief Example of solving a linear system with a Cholesky factorization
 *
 * @version 2.4.6
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

int check_solution(int, int, PLASMA_Complex64_t*, int, PLASMA_Complex64_t*, PLASMA_Complex64_t*, int);

int IONE=1;
int ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

int main ()
{

    int cores = 2;
    int N     = 10 ;
    int LDA   = 10 ;
    int NRHS  = 5 ;
    int LDB   = 10 ;
    int info;
    int info_solution;

    PLASMA_Complex64_t *A1   = (PLASMA_Complex64_t *)malloc(LDA*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *A2   = (PLASMA_Complex64_t *)malloc(LDA*N*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *B1   = (PLASMA_Complex64_t *)malloc(LDB*NRHS*sizeof(PLASMA_Complex64_t));
    PLASMA_Complex64_t *B2   = (PLASMA_Complex64_t *)malloc(LDB*NRHS*sizeof(PLASMA_Complex64_t));

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)||(!B1)||(!B2)){
        printf("Out of Memory \n ");
        return EXIT_SUCCESS;
    }

    /* Plasma Initialize */
    PLASMA_Init(cores);
    printf("-- PLASMA is initialized to run on %d cores. \n",cores);

    /* Initialize A1 and A2 for Symmetric Positif Matrix (Hessenberg in the complex case) */
    PLASMA_zplghe( (double)N, N, A1, LDA, 51 );
    PLASMA_zlacpy( PlasmaUpperLower, N, N, A1, LDA, A2, LDA );

    /* Initialize B1 and B2 */
    PLASMA_zplrnt( N, NRHS, B1, LDB, 371 );
    PLASMA_zlacpy( PlasmaUpperLower, N, NRHS, B1, LDB, B2, LDB );

    /* Plasma routines */
    info = PLASMA_zpotrf(PlasmaUpper, N, A2, LDA);
    info = PLASMA_zpotrs(PlasmaUpper, N, NRHS, A2, LDA, B2, LDB);

    /* Check the solution */
    info_solution = check_solution(N, NRHS, A1, LDA, B1, B2, LDB);

    if ((info_solution != 0)|(info != 0))
       printf("-- Error in ZPOTRS example ! \n");
    else
       printf("-- Run of ZPOTRS example successful ! \n");

    free(A1); free(A2); free(B1); free(B2);

    PLASMA_Finalize();

    return EXIT_SUCCESS;
}

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution of the linear system
 */

int check_solution(int N, int NRHS, PLASMA_Complex64_t *A1, int LDA, PLASMA_Complex64_t *B1, PLASMA_Complex64_t *B2, int LDB )
{
    int info_solution;
    double Rnorm, Anorm, Xnorm, Bnorm;
    PLASMA_Complex64_t alpha, beta;
    double *work = (double *)malloc(N*sizeof(double));
    double eps;

    eps = LAPACKE_dlamch_work('e');

    /* Initialize A1 and A2 for Symmetric Positive Matrix */
    alpha = 1.0;
    beta  = -1.0;

    Xnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B2, LDB, work);
    Anorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, N, A1, LDA, work);
    Bnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B1, LDB, work);

    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, NRHS, N, CBLAS_SADDR(alpha), A1, LDA, B2, LDB, CBLAS_SADDR(beta), B1, LDB);
    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, lapack_const(PlasmaInfNorm), N, NRHS, B1, LDB, work);

    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n",Rnorm/((Anorm*Xnorm+Bnorm)*N*eps));

    if (Rnorm/((Anorm*Xnorm+Bnorm)*N*eps) > 10.0){
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
     }
    else{
        printf("-- The solution is CORRECT ! \n");
        info_solution = 0;
    }

    free(work);

    return info_solution;
}
