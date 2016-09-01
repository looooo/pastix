/**
 *
 * @file z_spm_matvec_test.c
 *
 * Tests and validate the spm_convert routines.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 *
 **/
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <pastix.h>
#include <common.h>
#include <spm.h>
#include "cblas.h"
#include "lapacke.h"
#include <z_spm.h>
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
int
z_spm_matvec_check( int trans, const pastix_spm_t *spm )
{
    unsigned long long int seed = 35469;
    pastix_complex64_t *A, *x, *y0, *ys, *yd;
    pastix_complex64_t alpha, beta;

    double Anorm, Xnorm, Y0norm, Ysnorm, Ydnorm, Rnorm;
    double eps, result;
    int info_solution, start = 1;

    eps = LAPACKE_dlamch_work('e');

    core_zplrnt( 1, 1, &alpha, 1, 1, start, 0, seed ); start++;
    core_zplrnt( 1, 1, &beta,  1, 1, start, 0, seed ); start++;

    x = (pastix_complex64_t*)malloc(spm->gNexp * sizeof(pastix_complex64_t));
    core_zplrnt( spm->gNexp, 1, x, spm->gNexp, 1, start, 0, seed ); start += spm->gNexp;

    y0 = (pastix_complex64_t*)malloc(spm->gNexp * sizeof(pastix_complex64_t));
    core_zplrnt( spm->gNexp, 1, y0, spm->gNexp, 1, start, 0, seed ); start += spm->gNexp;

    /* Create a dense backup of spm */
    A = z_spm2dense( spm );

    /* Allocate cs/cd */
    ys = (pastix_complex64_t*)malloc(spm->gNexp * sizeof(pastix_complex64_t));
    yd = (pastix_complex64_t*)malloc(spm->gNexp * sizeof(pastix_complex64_t));

    /* Initialize cs/cd */
    memcpy( ys, y0, spm->gNexp * sizeof(pastix_complex64_t) );
    memcpy( yd, y0, spm->gNexp * sizeof(pastix_complex64_t) );

    /* Compute the sparse matrix-vector product */
    spmMatVec( trans, &alpha, spm, x, &beta, ys );

    /* Compute the dense matrix-vector product */
    cblas_zgemm( CblasColMajor, trans, CblasNoTrans, spm->gNexp, 1, spm->gNexp,
                 CBLAS_SADDR(alpha), A, spm->gNexp,
                                     x, spm->gNexp,
                 CBLAS_SADDR(beta), yd, spm->gNexp );

    Anorm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, spm->gNexp,  A, spm->gNexp );
    Xnorm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1,        x, spm->gNexp );
    Y0norm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1,       y0, spm->gNexp );
    Ysnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1,       ys, spm->gNexp );
    Ydnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1,       yd, spm->gNexp );

    core_zgeadd(PastixNoTrans, spm->gNexp, 1,
                -1., ys, spm->gNexp,
                 1., yd, spm->gNexp);
    Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', spm->gNexp, 1, yd, spm->gNexp );

    if ( 1 ) {
        printf("  ||A||_inf = %e, ||x||_inf = %e, ||y||_inf = %e\n"
               "  ||dense(a*A*x+b*y)||_inf = %e, ||sparse(a*A*x+b*y)||_inf = %e, ||R||_m = %e\n",
               Anorm, Xnorm, Y0norm, Ydnorm, Ysnorm, Rnorm);
    }

    result = Rnorm / ((Anorm + Xnorm + Y0norm) * spm->gNexp* eps);
    if (  isinf(Ydnorm) || isinf(Ysnorm) ||
          isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    free(A); free(x); free(y0); free(ys); free(yd);

    return info_solution;
}

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
int
z_spm_norm_check( const pastix_spm_t *spm )
{
    pastix_complex64_t *A;
    double norms, normd;
    double eps, result;
    int ret = 0;

    eps = LAPACKE_dlamch_work('e');

    /* Create a dense backup of spm */
    A = z_spm2dense( spm );

    /**
     * Test Norm Max
     */
    printf(" -- Test norm Max :");
    norms = spmNorm( PastixMaxNorm, spm );
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', spm->gNexp, spm->gNexp, A, spm->gNexp );
    result = fabs(norms - normd) / (normd * eps);

    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        ret++;
    }

    printf("   Nsparse = %e, Ndense = %e\n", norms, normd );
    printf("  | Nsparse - Ndense | / Ndense = %e\n", result);

    /**
     * Test Norm Inf
     */
    printf(" -- Test norm Inf :");
    norms = spmNorm( PastixInfNorm, spm );
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, spm->gNexp, A, spm->gNexp );
    result = fabs(norms - normd) / (normd * eps);
    result = result * ((double)(spm->gnnzexp)) / ((double)(spm->gNexp));

    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        ret++;
    }

    printf("   Nsparse = %e, Ndense = %e\n", norms, normd );
    printf("  | Nsparse - Ndense | / Ndense = %e\n", result);

    /**
     * Test Norm One
     */
    printf(" -- Test norm One :");
    norms = spmNorm( PastixOneNorm, spm );
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'O', spm->gNexp, spm->gNexp, A, spm->gNexp );
    result = fabs(norms - normd) / (normd * eps);
    result = result * ((double)(spm->gnnzexp)) / ((double)(spm->gNexp));

    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        ret++;
    }

    printf("   Nsparse = %e, Ndense = %e\n", norms, normd );
    printf("  | Nsparse - Ndense | / Ndense = %e\n", result);

    /**
     * Test Norm Frobenius
     */
    printf(" -- Test norm Frb :");
    norms = spmNorm( PastixFrobeniusNorm, spm );
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'F', spm->gNexp, spm->gNexp, A, spm->gNexp );
    result = abs(norms - normd) / (normd * eps);
    result = result / ((double)spm->gnnzexp);

    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        ret++;
    }

    printf("   Nsparse = %e, Ndense = %e\n", norms, normd );
    printf("  | Nsparse - Ndense | / Ndense = %e\n", result);

    free(A);
    return ret;
}
