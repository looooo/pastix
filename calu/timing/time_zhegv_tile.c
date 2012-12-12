/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zhegv_Tile"
/* See Lawn 41 page 120 */
/* cholesky + 2 trsm's + trd */
#define _FMULS (((2. / 3.) * ((double)N * (double)N * (double)N)) + (N * (1.0 / 6.0 * N + 0.5) * N) + 2 * (N * N * (double)((N + 1) / 2.0) ))
#define _FADDS (((2. / 3.) * ((double)N * (double)N * (double)N)) + (N * (1.0 / 6.0 * N )      * N) + 2 * (N * N * (double)((N + 1) / 2.0) ))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    PLASMA_desc *descT;
    int itype = 1;
    int vec   = PlasmaNoVec;
    int uplo  = PlasmaUpper;

    LDA = max(LDA, N);

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descQ, (vec == PlasmaVec), PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( W, 1, double, N, 1 );

    /* Initialization */
    PLASMA_zplghe_Tile((double)0.0, descA, 51   );
    PLASMA_zplghe_Tile((double)N,   descB, 3753 );

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_zhegv(N, N, &descT);

    START_TIMING();
    PLASMA_zhegv_Tile( itype, vec, uplo, descA, descB, W, descT, descQ );
    STOP_TIMING();

    /* DeAllocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descT);

    if (vec == PlasmaVec) {
        PASTE_CODE_FREE_MATRIX( descQ );
    }
    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );
    free( W );

    return 0;
}
