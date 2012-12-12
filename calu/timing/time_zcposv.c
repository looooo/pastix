/**
 *
 * @precisions mixed zc -> ds
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zposv"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_POTRF( N ) + FMULS_POTRS( N, NRHS ))
#define _FADDS (FADDS_POTRF( N ) + FADDS_POTRS( N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    int iter;

    LDA = max(LDA, N);
    
    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, PLASMA_Complex64_t, LDA, N    );
    PASTE_CODE_ALLOCATE_MATRIX( B, 1, PLASMA_Complex64_t, LDB, NRHS );

    /* Initialiaze Data */
    PLASMA_zplghe((double)N, N, A, LDA, 51 );
    PLASMA_zplrnt( N, NRHS, B, LDB, 5673 );

    PASTE_CODE_ALLOCATE_COPY( X, 1, PLASMA_Complex64_t, B, LDB, NRHS );

    /* PLASMA ZCPOSV */
    START_TIMING();
    PLASMA_zcposv(PlasmaUpper, N, NRHS, A, LDA, B, LDB, X, LDB, &iter);
    STOP_TIMING();

    /* Check the solution */
    if (check)
      {
        dparam[IPARAM_RES] = z_check_solution(N, N, NRHS, A, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
      }

    free(A); free(B); free(X); 
    return 0;
}
