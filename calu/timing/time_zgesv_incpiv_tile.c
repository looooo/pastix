/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgesv_incpiv_Tile"
/* See Lawn 41 page 120 */
#define _FMULS (FMULS_GETRF( N, N ) + FMULS_GETRS( N, NRHS ))
#define _FADDS (FADDS_GETRF( N, N ) + FADDS_GETRS( N, NRHS ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PLASMA_desc *descL;
    int *piv;
    PASTE_CODE_IPARAM_LOCALS( iparam );
    
    if ( M != N ) {
        fprintf(stderr, "This timing works only with M == N\n");
        return -1;
    }
    
    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N    );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDB, N, NRHS );

    /* Initialize A and b */
    PLASMA_zplrnt_Tile( descA, 8796 );
    PLASMA_zplrnt_Tile( descB, 7732 );

    /* Save AT and bT in lapack layout for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, PLASMA_Complex64_t, LDA, N    );
    PASTE_TILE_TO_LAPACK( descB, b, check, PLASMA_Complex64_t, LDB, NRHS );

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_zgesv_incpiv_Tile(N, &descL, &piv);

    START_TIMING();
    PLASMA_zgesv_incpiv_Tile( descA, descL, piv, descB );
    STOP_TIMING();
    
    /* Allocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descL);

    /* Check the solution */
    if ( check )
    {
        PASTE_TILE_TO_LAPACK( descB, x, check, PLASMA_Complex64_t, LDB, NRHS );
        
        dparam[IPARAM_RES] = z_check_solution(N, N, NRHS, A, LDA, b, x, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(A); free(b); free(x);
    }

    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );
    free( piv );

    return 0;
}
