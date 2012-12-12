/**
 *
 * @precisions normal z -> c d s
 *
 **/
#define _TYPE  PLASMA_Complex64_t
#define _PREC  double
#define _LAMCH LAPACKE_dlamch_work

#define _NAME  "PLASMA_zgeqrf"
/* See Lawn 41 page 120 */
#define _FMULS FMULS_GEQRF(M, N)
#define _FADDS FADDS_GEQRF(M, N)

#include "./timing.c"
#include <cblas.h>

#define COMPLEX

static int
RunTest(int *iparam, double *dparam, real_Double_t *t_)
{
    int info, lwork, i;
    int *jpvt, *jpvt2;
    PLASMA_Complex64_t *tau, *tau2;
    PLASMA_Complex64_t *work;
    double *rwork;
    double tnorm;
    const PLASMA_Complex64_t mzone = -1.0;
    
    PASTE_CODE_IPARAM_LOCALS( iparam );

    if ( M != N && check ) {
        fprintf(stderr, "Check cannot be perfomed with M != N\n");
        check = 0;
    }

    /* Allocate Data */
    PASTE_CODE_ALLOCATE_MATRIX( A, 1, PLASMA_Complex64_t, LDA, N );
    jpvt  = (int*)                malloc( N   * sizeof(int)                );
    tau   = (PLASMA_Complex64_t*) malloc( 2*N   * sizeof(PLASMA_Complex64_t) );
    rwork = (double*)             malloc( 2*2*N * sizeof(double)             );
    
    #ifdef COMPLEX
    lwork = (N+1)*NB;
    #else
    lwork = (N+1)*NB + 2*N;
    #endif
    work  = (PLASMA_Complex64_t*) malloc( lwork * sizeof(PLASMA_Complex64_t) );
    
    if ( jpvt == NULL || tau == NULL || work == NULL || rwork == NULL ) {
        printf( "malloc failed\n" );
        return -1;
    }
    
    /* zero out pivots (required by LAPACK) */
    for( i = 0; i < N; ++i ) {
        jpvt[i] = 0;
    }
    
    /* Initialize Data */
    PLASMA_zplrnt(M, N, A, LDA, 123456);
    
    /* Save AT in lapack layout for check */
    PASTE_CODE_ALLOCATE_COPY( Acpy, check, PLASMA_Complex64_t, A, LDA, N );

    START_TIMING();
    info = PLASMA_zgeqp3( M, N, A, LDA, jpvt, tau, work, rwork );
    STOP_TIMING();
    
    /* Check the solution */
    if ( info != 0 ) {
        printf( "\nPLASMA_zgeqp3 returned error %d, due to a numerical instability that we detect but do not yet handle.\n", info );
    }
    else if ( check ) {
        jpvt2  = (int*)                malloc( N * sizeof(int)    );
        tau2   = (PLASMA_Complex64_t*) malloc( N * sizeof(PLASMA_Complex64_t) );
        if ( jpvt2 == NULL || tau2 == NULL ) {
            printf( "test malloc failed\n" );
        }
        
        /* zero out pivots (required by LAPACK) */
        for( i = 0; i < N; ++i ) {
            jpvt2[i] = 0;
        }
        
        dparam[IPARAM_ANORM] = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'F', N, N, Acpy, LDA, rwork );
        dparam[IPARAM_XNORM] = 1.;
        dparam[IPARAM_BNORM] = 0.;
        
        double time = cWtime();
        #ifdef COMPLEX
        info = LAPACKE_zgeqp3_work( LAPACK_COL_MAJOR, N, N, Acpy, LDA, jpvt2, tau2, work, lwork, rwork );
        #else
        info = LAPACKE_zgeqp3_work( LAPACK_COL_MAJOR, N, N, Acpy, LDA, jpvt2, tau2, work, lwork );
        #endif
        time = cWtime() - time;
        if ( info != 0 ) {
            printf( "qp3 returned error %d\n", info );
        }
        /* printf( "   %7.3f", time ); */
                
        cblas_zaxpy( LDA*N, CBLAS_SADDR(mzone), A, 1, Acpy, 1 );
        dparam[IPARAM_RES] = LAPACKE_zlange_work( LAPACK_COL_MAJOR, 'F', N, N, Acpy, LDA, rwork );
        
        cblas_zaxpy( N, CBLAS_SADDR(mzone), tau, 1, tau2, 1 );
        tnorm = cblas_dznrm2( N, tau2, 1 );
        /* printf( "|t|=%.2e\n", tnorm ); */
        
        /*
        * for( i = 0; i < N; ++i ) {
        *     if ( jpvt[i]+1 != jpvt2[i] ) {
        *         printf( "pivot mis-match jpvt[%2d]+1=%2d, jpvt2[%2d]=%2d\n", i, jpvt[i]+1, i, jpvt2[i] );
        *     }
        * }
        */
        
        free( Acpy  );
        free( jpvt2 );
        free( tau2  );
    }

    /* Free data */
    free( A     );
    free( jpvt  );
    free( tau   );
    free( work  );
    free( rwork );

    return 0;
}
