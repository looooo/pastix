/**
 *
 * @file z_bcsc_tests.c
 *
 * Tests and validate the bcsc routines.
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
#include <bcsc.h>
#include "lapacke.h"
#include <z_spm.h>
#include <z_bcsc.h>
#include <bcsc.h>
#include <order.h>
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
int
z_bcsc_matvec_check( int trans, const pastix_spm_t *spm, const pastix_data_t *pastix_data )
{
    unsigned long long int seed = 35469;
    pastix_complex64_t *x, *y0, *ys, *yd;
    pastix_complex64_t alpha, beta;

    double Anorm, Xnorm, Y0norm, Ysnorm, Ydnorm, Rnorm;
    double eps, result;
    int info_solution, start = 1;

    eps = LAPACKE_dlamch_work('e');

    core_zplrnt( 1, 1, &alpha, 1, 1, start, 0, seed ); start++;
    core_zplrnt( 1, 1, &beta,  1, 1, start, 0, seed ); start++;

    x = (pastix_complex64_t*)malloc(spm->gN * sizeof(pastix_complex64_t));
    core_zplrnt( spm->gN, 1, x, spm->gN, 1, start, 0, seed ); start += spm->gN;

    y0 = (pastix_complex64_t*)malloc(spm->gN * sizeof(pastix_complex64_t));
    core_zplrnt( spm->gN, 1, y0, spm->gN, 1, start, 0, seed ); start += spm->gN;

    /* Allocate cs/cd */
    ys    = (pastix_complex64_t*)malloc(spm->gN * sizeof(pastix_complex64_t));
    yd    = (pastix_complex64_t*)malloc(spm->gN * sizeof(pastix_complex64_t));

    /* Initialize cs/cd */
    memcpy( ys, y0, spm->gN * sizeof(pastix_complex64_t) );
    memcpy( yd, y0, spm->gN * sizeof(pastix_complex64_t) );

    /* Compute the spm matrix-vector product */
    spmMatVec( trans, &alpha, spm, x, &beta, ys );

    /* Compute the bcsc matrix-vector product */
    z_bcscApplyPerm( pastix_data->bcsc->gN, 1, yd, pastix_data->bcsc->gN, pastix_data->ordemesh->permtab );
    z_bcscApplyPerm( pastix_data->bcsc->gN, 1, x, pastix_data->bcsc->gN, pastix_data->ordemesh->permtab );

    bcscMatVec( trans, &alpha, pastix_data->bcsc, x, &beta, yd );

    z_bcscApplyPerm( pastix_data->bcsc->gN, 1, yd, pastix_data->bcsc->gN, pastix_data->ordemesh->peritab );
    z_bcscApplyPerm( pastix_data->bcsc->gN, 1, x, pastix_data->bcsc->gN, pastix_data->ordemesh->peritab );

    // Anorm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gN, spm->gN,  A, spm->gN );
    Anorm  = spmNorm( PastixInfNorm, spm );
    Xnorm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gN, 1,        x, spm->gN );
    Y0norm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gN, 1,       y0, spm->gN );
    Ysnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gN, 1,       ys, spm->gN );
    Ydnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gN, 1,       yd, spm->gN );

    core_zgeadd(PastixNoTrans, spm->gN, 1,
                -1., ys, spm->gN,
                 1., yd, spm->gN );
    Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', spm->gN, 1, yd, spm->gN );

    if ( 1 ) {
        printf("  ||A||_inf = %e, ||x||_inf = %e, ||y||_inf = %e\n"
               "  ||spm(a*A*x+b*y)||_inf = %e, ||bcsc(a*A*x+b*y)||_inf = %e, ||R||_m = %e\n",
               Anorm, Xnorm, Y0norm, Ysnorm, Ydnorm, Rnorm);
    }

    result = Rnorm / ((Anorm + Xnorm + Y0norm) * spm->gN* eps);
    if (  isinf(Ydnorm) || isinf(Ysnorm) ||
          isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

    free(x); free(y0); free(ys); free(yd);

    return info_solution;
}

/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
int
z_bcsc_norm_check( const pastix_spm_t *spm, const pastix_bcsc_t *bcsc )
{
    double norms, normd;
    double eps, result;
    int ret = 0;

    eps = LAPACKE_dlamch_work('e');

    /**
     * Test Norm Max
     */
    printf(" -- Test norm Max :");
    norms = spmNorm( PastixMaxNorm, spm );
    normd = z_bcscNorm( PastixMaxNorm, bcsc );
    result = fabs(norms - normd) / (norms * eps);

    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        ret++;
    }

    printf("   Nspm = %e, Nbcsc = %e\n", norms, normd );
    printf("  | Nspm - Nbcsc | / Nspm = %e\n", result);

    /**
     * Test Norm Inf
     */
    printf(" -- Test norm Inf :");
    norms = spmNorm( PastixInfNorm, spm );
    normd = z_bcscNorm( PastixInfNorm, bcsc );
    result = fabs(norms - normd) / (norms * eps);
    result = result * ((double)(spm->gN)) / ((double)(spm->gnnz));

    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        ret++;
    }

    printf("   Nspm = %e, Nbcsc = %e\n", norms, normd );
    printf("  | Nspm - Nbcsc | / Nspm = %e\n", result);

    /**
     * Test Norm One
     */
    printf(" -- Test norm One :");
    norms = spmNorm( PastixOneNorm, spm );
    normd = z_bcscNorm( PastixOneNorm, bcsc );
    result = fabs(norms - normd) / (norms * eps);
    result = result * ((double)(spm->gN)) / ((double)(spm->gnnz));

    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        ret++;
    }

    printf("   Nspm = %e, Nbcsc = %e\n", norms, normd );
    printf("  | Nspm - Nbcsc | / Nspm = %e\n", result);

    /**
     * Test Norm Frobenius
     */
    printf(" -- Test norm Frb :");
    norms = spmNorm( PastixFrobeniusNorm, spm );
    normd = z_bcscNorm( PastixFrobeniusNorm, bcsc );
    result = abs(norms - normd) / (norms * eps);
    result = result / ((double)spm->gnnz);

    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        ret++;
    }

    printf("   Nspm = %e, Nbcsc = %e\n", norms, normd );
    printf("  | Nspm - Nbcsc | / Nspm = %e\n", result);

    return ret;
}
