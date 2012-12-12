/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgebrd_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEBRD( N, N )
#define _FADDS FADDS_GEBRD( N, N )

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    PLASMA_desc *descT;
    int vec   = PlasmaNoVec;
    
    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( W, 1, double, N, 1 );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descQ, (vec == PlasmaVec), PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );
    
    /* Initialiaze Data */
    PLASMA_zplghe_Tile((double)0.0, descA, 51 );

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_zgesvd(N, N, &descT);

    START_TIMING();
    PLASMA_zgesvd_Tile( vec, vec, descA, W, descQ, descQ, descT );
    STOP_TIMING();
    
    /* Check the solution */
    if ( check )
      {
      }

    /* DeAllocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descT);

    if (vec == PlasmaVec) {
      PASTE_CODE_FREE_MATRIX( descQ );
    }
    PASTE_CODE_FREE_MATRIX( descA );
    free( W );

    return 0;
}
