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

#include <assert.h>
#include "../control/common.h"
#include "./timing.c"

#define A(m,n)    BLKADDR(A, PLASMA_Complex64_t, m, n)

void initMatrix(PLASMA_desc A)
{
    int it, jt, i, j;
    int row, col;
    double val = 1.0;
    int *mypivot = (int *)malloc(sizeof(int)*A.m);

    int par1, par2, par3;
    int p = 5;
    assert( A.m % (p*p) == 0 );
    for( i=0; i<A.m; i++){
        par1 = i % p;
        par2 = ((i-par1)/p) % (A.m/(p*p));
        par3 = (( (i-par1)/p - par2) * p * p)/A.m;
        mypivot[i] = par2 * p * p + par3 * p + par1;
    }

    for(it=0; it<A.mt; it++) {
        int tempm = it == A.mt-1 ? A.m-it*A.mb : A.mb;
        int ldam = BLKLDD(A, it);
        for(i=0; i<tempm; i++) {
            row = mypivot[it * A.mb + i ];
            //row = it * A.mb + i;
            val = 1.0;
            for(jt=0; jt<A.nt; jt++) {
                int tempn = jt == A.nt-1 ? A.n-jt*A.nb : A.nb;
                for(j=0; j<tempn; j++) {
                    col = jt * A.nb + j;
                    if (col < row) {
                        (A(it, jt))[j*ldam+i] = col + 1.0;
                    } else {
                        (A(it, jt))[j*ldam+i] = A.m;
                    }
                }
            }
        }
    }
}

/* static void printMatrix(char *name, int M, int N, PLASMA_Complex64_t *A, int LDA) */
/* { */
/*     int i, j; */
/*     fprintf(stderr, "======== %s =========\n", name ); */

/*     for(i=0; i<M; i++) { */
/*         for(j=0; j<N; j++) { */
/*             fprintf(stderr, "%e ", A[j*LDA+i] ); */
/*         } */
/*         fprintf(stderr, "\n" ); */
/*     } */
/* } */

static double z_check2(int M, int N, PLASMA_Complex64_t *A, int LDA )
{
    double Rnorm = -1.00;
    PLASMA_Complex64_t mzone = -1.0;
    int i, j;

    PASTE_CODE_ALLOCATE_MATRIX( R, 1, PLASMA_Complex64_t, LDA, N );

    for(i=0; i<M; i++) {
        for(j=0; j<N; j++) {
            if (j < i) {
                R[j*LDA+i] = 1.0 / (M-j);
            } else {
                R[j*LDA+i] = M - i;
            }
        }
    }

    cblas_zaxpy(LDA*N, CBLAS_SADDR(mzone), A, 1, R, 1);
    Rnorm = LAPACKE_zlange_work(LAPACK_COL_MAJOR, 'm', M, N, R, LDA, NULL);

    //printMatrix("R", M, N, R, LDA );

    free( R );
    return Rnorm;
}

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

    int l = 4361;
    /* for( l = 1; l< 10000; l+=217) { */

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDA, M, N );
    PASTE_CODE_ALLOCATE_MATRIX( piv, 1, int, min(M, N), 1 );

    PLASMA_zplrnt_Tile(descA, l);
    //initMatrix(*descA);

    /* Save AT in lapack layout for check */
    PASTE_TILE_TO_LAPACK( descA, A, check, PLASMA_Complex64_t, LDA, N );

    START_TIMING();
    PLASMA_zgetrf_Tile( descA, piv );
    STOP_TIMING();

    /* Check the solution */
    if ( check )
    {
        PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, PLASMA_Complex64_t, PlasmaComplexDouble, LDB, N, NRHS );
        PLASMA_zplrnt_Tile( descB, 7732 );
        PASTE_TILE_TO_LAPACK( descB, b, check, PLASMA_Complex64_t, LDB, NRHS );

        PLASMA_zgetrs_Tile( PlasmaNoTrans, descA, piv, descB );

        PASTE_TILE_TO_LAPACK( descB, x, check, PLASMA_Complex64_t, LDB, NRHS );
        dparam[IPARAM_RES] = z_check_solution(M, N, NRHS, A, LDA, b, x, LDB,
                                              &(dparam[IPARAM_ANORM]),
                                              &(dparam[IPARAM_BNORM]),
                                              &(dparam[IPARAM_XNORM]));

        PASTE_CODE_FREE_MATRIX( descB );

        /* PASTE_TILE_TO_LAPACK( descA, A2, check, PLASMA_Complex64_t, LDA, N ); */
        /* dparam[IPARAM_BNORM] = z_check2(M, N, A2, LDA ); */

        free(A); free(b); free(x);
        /* free( A2 ); */
        //        printf( "Residual=%9.3e\n", dparam[IPARAM_RES] );
    }

    PASTE_CODE_FREE_MATRIX( descA );
    free( piv );
    /*}*/

    return 0;
}
