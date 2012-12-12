/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zheev_Tile"
/* See Lawn 41 page 120 */
#define _FMULS ((2. / 3.) * ((double)N * (double)N * (double)N))
#define _FADDS ((2. / 3.) * ((double)N * (double)N * (double)N))

#include "./timing.c"

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );
    PLASMA_desc *descT;
    int jobu  = PlasmaNoVec;
    int jobvt = PlasmaNoVec;

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descU, (jobu == PlasmaVec), PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descVT, (jobvt == PlasmaVec), PLASMA_Complex64_t, PlasmaComplexDouble, LDA, N, N );
    PASTE_CODE_ALLOCATE_MATRIX( S, 1, double, N, 1 );

    /* Initialiaze Data */
    PLASMA_zplrnt_Tile(descA, 51 );

    /* Save AT and bT in lapack layout for check */
    if ( check ) {
    }

    /* Allocate Workspace */
    PLASMA_Alloc_Workspace_zgesvd(N, N, &descT);

    START_TIMING();
    PLASMA_zgesvd_Tile(jobu, jobvt, descA, S, descU, descVT, descT);
    STOP_TIMING();
    
    /* Check the solution */
    if ( check )
      {
      }

    /* DeAllocate Workspace */
    PLASMA_Dealloc_Handle_Tile(&descT);

    if (jobu == PlasmaVec) {
      PASTE_CODE_FREE_MATRIX( descU );
    }
    if (jobvt == PlasmaVec) {
      PASTE_CODE_FREE_MATRIX( descVT );
    }
    PASTE_CODE_FREE_MATRIX( descA );
    free( S );

    return 0;
}
