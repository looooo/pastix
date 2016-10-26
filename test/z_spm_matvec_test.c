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

pastix_complex64_t *
z_spm2dense( const pastix_spm_t *spm )
{
    pastix_int_t i, j, lda, baseval;
    pastix_complex64_t *A, *valptr;
    pastix_int_t *colptr, *rowptr;

    assert( spm->fmttype == PastixCSC );
    assert( spm->flttype == PastixComplex64 );

    lda = spm->gN;
    A = (pastix_complex64_t*)malloc(lda * lda * sizeof(pastix_complex64_t));
    memset( A, 0, lda * lda * sizeof(pastix_complex64_t));

    baseval = spmFindBase( spm );
    i = 0; j = 0;

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (pastix_complex64_t*)(spm->values);

    switch( spm->mtxtype ){
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PastixHermitian:
        for(i=0; i<spm->n; i++, colptr++)
        {
            for(j=colptr[0]; j<colptr[1]; j++, rowptr++, valptr++)
            {
                if( i == (*rowptr-baseval) ) {
                    /* Make sure the matrix is hermitian */
                    A[ i * lda + (*rowptr - baseval) ] = creal(*valptr) + I * 0.;
                }
                else {
                    A[ i * lda + (*rowptr - baseval) ] = *valptr;
                    A[ (*rowptr - baseval) * lda + i ] = conj(*valptr);
                }
            }
        }
        break;
#endif
    case PastixSymmetric:
        for(i=0; i<spm->n; i++, colptr++)
        {
            for(j=colptr[0]; j<colptr[1]; j++, rowptr++, valptr++)
            {
                A[ i * lda + (*rowptr - baseval) ] = *valptr;
                A[ (*rowptr - baseval) * lda + i ] = *valptr;
            }
        }
        break;
    case PastixGeneral:
    default:
        for(i=0; i<spm->n; i++, colptr++)
        {
            for(j=colptr[0]; j<colptr[1]; j++, rowptr++, valptr++)
            {
                A[ i * lda + (*rowptr - baseval) ] = *valptr;
            }
        }
    }
    return A;
}

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

    x = (pastix_complex64_t*)malloc(spm->gN * sizeof(pastix_complex64_t));
    core_zplrnt( spm->gN, 1, x,  spm->gN, 1, start, 0, seed ); start += spm->gN;

    y0 = (pastix_complex64_t*)malloc(spm->gN * sizeof(pastix_complex64_t));
    core_zplrnt( spm->gN, 1, y0, spm->gN, 1, start, 0, seed ); start += spm->gN;

    /* Create a dense backup of spm */
    A = z_spm2dense( spm );

    /* Allocate cs/cd */
    ys = (pastix_complex64_t*)malloc(spm->gN * sizeof(pastix_complex64_t));
    yd = (pastix_complex64_t*)malloc(spm->gN * sizeof(pastix_complex64_t));

    /* Initialize cs/cd */
    memcpy( ys, y0, spm->gN * sizeof(pastix_complex64_t) );
    memcpy( yd, y0, spm->gN * sizeof(pastix_complex64_t) );

    /* Compute the sparse matrix-vector product */
    spmMatVec( trans, &alpha, spm, x, &beta, ys );

    /* Compute the dense matrix-vector product */
    cblas_zgemm( CblasColMajor, trans, CblasNoTrans, spm->gN, 1, spm->gN,
                 CBLAS_SADDR(alpha), A, spm->gN,
                                     x, spm->gN,
                 CBLAS_SADDR(beta), yd, spm->gN );

    Anorm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gN, spm->gN,  A, spm->gN );
    Xnorm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gN, 1,        x, spm->gN );
    Y0norm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gN, 1,       y0, spm->gN );
    Ysnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gN, 1,       ys, spm->gN );
    Ydnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gN, 1,       yd, spm->gN );

    core_zgeadd(PastixNoTrans, spm->gN, 1,
                -1., ys, spm->gN,
                 1., yd, spm->gN);
    Rnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', spm->gN, 1, yd, spm->gN );

    if ( 1 ) {
        printf("  ||A||_inf = %e, ||x||_inf = %e, ||y||_inf = %e\n"
               "  ||dense(a*A*x+b*y)||_inf = %e, ||sparse(a*A*x+b*y)||_inf = %e, ||R||_m = %e\n",
               Anorm, Xnorm, Y0norm, Ydnorm, Ysnorm, Rnorm);
    }

    result = Rnorm / ((Anorm + Xnorm + Y0norm) * spm->gN* eps);
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
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'M', spm->gN, spm->gN,  A, spm->gN );
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
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gN, spm->gN,  A, spm->gN );
    result = fabs(norms - normd) / (normd * eps);
    result = result * ((double)(spm->gN)) / ((double)(spm->gnnz));

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
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'O', spm->gN, spm->gN,  A, spm->gN );
    result = fabs(norms - normd) / (normd * eps);
    result = result * ((double)(spm->gN)) / ((double)(spm->gnnz));

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
    normd = LAPACKE_zlange( LAPACK_COL_MAJOR, 'F', spm->gN, spm->gN,  A, spm->gN );
    result = abs(norms - normd) / (normd * eps);
    result = result / ((double)spm->gnnz);

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
