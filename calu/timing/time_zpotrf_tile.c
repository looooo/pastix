/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zpotrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_POTRF( N )
#define _FADDS FADDS_POTRF( N )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    int uplo = PlasmaUpper;

    LDA = max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1,     PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, check, PLASMA_Complex64_t, PlasmaComplexDouble, LDB, N, NRHS );
    PLASMA_zplghe_Tile( (double)N, descA, 51 );

    /* Save A for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, PLASMA_Complex64_t, LDA, N );

    /* PLASMA ZPOSV */
    START_TIMING();
    PLASMA_zpotrf_Tile(uplo, descA);
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        PLASMA_zplrnt_Tile( descB, 7672 );
        PASTE_TILE_TO_LAPACK( descB, B, check, PLASMA_Complex64_t, LDB, NRHS );

        PLASMA_zpotrs_Tile( uplo, descA, descB );

        PASTE_TILE_TO_LAPACK( descB, X, check, PLASMA_Complex64_t, LDB, NRHS );

        dparam[IPARAM_RES] = z_check_solution(N, N, NRHS, A, LDA, B, X, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));

        PASTE_CODE_FREE_MATRIX( descB );
        free( A );
        free( B );
        free( X );
      }

    PASTE_CODE_FREE_MATRIX( descA );

    return 0;
}
