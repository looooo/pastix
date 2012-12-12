/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgetrf_Tile"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GETRF(M, N)
#define _FADDS FADDS_GETRF(M, N)

int plasma_element_size(int type);
#include "./timing.c"
#include "../control/descriptor.h"

#define BLKLDD(A, k) ( ( (k) + (A).i/(A).mb) < (A).lm1 ? (A).mb : (A).lm%(A).mb )

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_) 
{
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, M, N );
    
    PLASMA_zplrnt_Tile(descA, 3456);

    {
        PLASMA_Complex64_t *Amat;
        int m, i, ldam;
        for(m=0; m<MT; m++) {
            ldam = BLKLDD( *descA, m );
            Amat = (PLASMA_Complex64_t*)plasma_getaddr(*descA, m, m);
            for(i=0; i<ldam; i++) {
                Amat[i*ldam+i] += max(M,N);
            }
        }
    }

    /* Save AT in lapack layout for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, PLASMA_Complex64_t, LDA, N );
    
    START_TIMING();
    PLASMA_zgetrf_nopiv_Tile( descA );
    STOP_TIMING();
    
    /* Check the solution */
    if ( check )
    {
        PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDB, N, NRHS );
        PLASMA_zplrnt_Tile( descB, 7732 );
        PASTE_TILE_TO_LAPACK( descB, b, check, PLASMA_Complex64_t, LDB, NRHS );

        PLASMA_ztrsm_Tile( PlasmaLeft, PlasmaLower, PlasmaNoTrans, PlasmaUnit,
                           1.0, descA, descB );
        PLASMA_ztrsm_Tile( PlasmaLeft, PlasmaUpper, PlasmaNoTrans, PlasmaNonUnit,
                           1.0, descA, descB );

        PASTE_TILE_TO_LAPACK( descB, x, check, PLASMA_Complex64_t, LDB, NRHS );
        dparam[IPARAM_RES] = z_check_solution(M, N, NRHS, A, LDA, b, x, LDB,
                                              &(dparam[IPARAM_ANORM]), 
                                              &(dparam[IPARAM_BNORM]), 
                                              &(dparam[IPARAM_XNORM]));
        free(A); free(b); free(x);
    }

    PASTE_CODE_FREE_MATRIX( descA );
 
    return 0;
}
