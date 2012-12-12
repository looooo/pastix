/**
 *
 * @precisions mixed zc -> ds
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zcgesv"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GETRF( N, N ) + FMULS_GETRS( N, NRHS ))
#define _FADDS (FADDS_GETRF( N, N ) + FADDS_GETRS( N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    int                iter = 0;
    
    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }
    
    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, PLASMA_Complex64_t, LDA, N    );
    PASTE_CODE_ALLOCATE_MATRIX( B, 1, PLASMA_Complex64_t, LDB, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX( X, 1, PLASMA_Complex64_t, LDB, NRHS );
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, N, 1 );
    
    /* Initialiaze Data */
    PLASMA_zplrnt( N, N,    A, LDA,   51 );
    PLASMA_zplrnt( N, NRHS, B, LDB, 5673 );

    /* Save A and b  */
    PASTE_CODE_ALLOCATE_COPY( Acpy, check, PLASMA_Complex64_t, A, LDA, N    );
    PASTE_CODE_ALLOCATE_COPY( Bcpy, check, PLASMA_Complex64_t, B, LDB, NRHS );

    START_TIMING();
    PLASMA_zcgesv( N, NRHS, A, LDA, piv, B, LDB, X, LDB, &iter );
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = z_check_solution(N, N, NRHS, Acpy, LDA, Bcpy, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(Acpy); free(Bcpy);
    }

    free( piv );
    free( X );
    free( B );
    free( A );


    return 0;
}
