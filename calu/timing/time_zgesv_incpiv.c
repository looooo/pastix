/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgesv_incpiv"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GETRF( N, N ) + FMULS_GETRS( N, NRHS ))
#define _FADDS (FADDS_GETRF( N, N ) + FADDS_GETRS( N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PLASMA_desc *L;
    int *piv;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }
    
    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, PLASMA_Complex64_t, LDA, N    );
    PASTE_CODE_ALLOCATE_MATRIX( X, 1, PLASMA_Complex64_t, LDB, NRHS );
    
    /* Initialiaze Data */
    PLASMA_zplrnt( N, N,    A, LDA,   51 );
    PLASMA_zplrnt( N, NRHS, X, LDB, 5673 );

    PLASMA_Alloc_Workspace_zgesv_incpiv(N, &L, &piv);

    /* Save A and b  */
    PASTE_CODE_ALLOCATE_COPY( Acpy, check, PLASMA_Complex64_t, A, LDA, N    );
    PASTE_CODE_ALLOCATE_COPY( B,    check, PLASMA_Complex64_t, X, LDB, NRHS );

    START_TIMING();
    PLASMA_zgesv_incpiv( N, NRHS, A, N, L, piv, X, LDB );
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = z_check_solution(N, N, NRHS, Acpy, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(Acpy); free(B);
    }

    PLASMA_Dealloc_Handle_Tile( &L );
    free( piv );
    free( X );
    free( A );


    return 0;
}
