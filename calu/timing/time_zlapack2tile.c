/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch

#define _NAME  "PLASMA_zlapack_to_tile"
/* See Lawn 41 page 120 */
#define _FMULS (0.0)
#define _FADDS (M * N * sizeof(_TYPE))

#include "./timing.c"

int z_check_conversion(int m, int n, int mba, int nba, int mbb, int nbb,
                      PLASMA_Complex64_t *A, PLASMA_Complex64_t *B, 
                      int (*mapA)(int, int, int, int, int, int), int (*mapB)(int, int, int, int, int, int)) {
    int i, j;

    for( j=0; j<n; j++) {
        for (i=0; i<m; i++) {
            if (A[ mapA(m, n, mba, nba, i, j) ] != B[ mapB(m, n, mbb, nbb, i, j) ] ) {
                return -1; 
            }
        }
    }
    return 0;
}

static int
RunTest(int *iparam, _PREC *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );

    LDA = M;

    dparam[IPARAM_ANORM] = (_PREC)M;
    dparam[IPARAM_BNORM] = (_PREC)_FADDS;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, PLASMA_Complex64_t, LDA, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, M, N );
    
    /* Initialize Data */
    PLASMA_zplrnt(M, N, A, LDA, 3456);

    START_TIMING();
    PLASMA_Lapack_to_Tile( (void *)A, LDA, descA);
    STOP_TIMING();
    
    /* Check the solution */
    if (check)
    {
        dparam[IPARAM_RES] = (_PREC)z_check_conversion(M, N, M, 1, MB, NB, A, descA->mat, map_CM, map_CCRB);
    }

    PASTE_CODE_FREE_MATRIX( descA );
    free( A );

    return 0;
}
