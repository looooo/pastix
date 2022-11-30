/**
 *
 * @file z_bcsc_tests.c
 *
 * Tests and validate the bcsc routines.
 *
 * @copyright 2015-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @author Vincent Bridonneau
 * @date 2022-12-05
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
#include "bcsc/bcsc.h"
#include "bcsc/bcsc_z.h"
#include <lapacke.h>
#include <spm/z_spm.h>
#include <pastix/order.h>
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"
#include "z_tests.h"

int
z_bcsc_spmv_check( spm_trans_t       trans,
                   const spmatrix_t *spm,
                   pastix_data_t    *pastix_data )
{
    unsigned long long int seed  = 35469;
    unsigned long long int seedX = 79243;
    unsigned long long int seedY = 29037;
    pastix_complex64_t *xd, *yd; /* The potentially distributed or replicated vectors    */
    pastix_complex64_t *xr, *yr; /* The always replicated vectors to validate the result */
    pastix_complex64_t alpha, beta;

    double Anorm, Xnorm, Ynorm, Yrnorm, Ydnorm, Rnorm;
    double eps, result;
    int    rc, info_solution, start = 1;
    int    nrhs = 1;

    eps = LAPACKE_dlamch_work('e');

    core_zplrnt( 1, 1, &alpha, 1, 1, start, 0, seed ); start++;
    core_zplrnt( 1, 1, &beta,  1, 1, start, 0, seed ); start++;

    /* Make sure alpha/beta are doubles */
#if defined(PRECISION_c) || defined(PRECISION_z)
    alpha = creal( alpha );
    beta  = creal( beta  );
#endif

    /* spm->nexp should always be of the correct size for both replicated and distributed vectors */
    xd = (pastix_complex64_t*)malloc( spm->nexp  * sizeof(pastix_complex64_t) );
    yd = (pastix_complex64_t*)malloc( spm->nexp  * sizeof(pastix_complex64_t) );
    yr = (pastix_complex64_t*)malloc( spm->gNexp * sizeof(pastix_complex64_t) );

    if ( spm->loc2glob == NULL ) {
        /* The vectors are replicated */
        core_zplrnt( spm->nexp, nrhs, xd, spm->nexp, spm->gNexp, 0, 0, seedX );
        core_zplrnt( spm->nexp, nrhs, yd, spm->nexp, spm->gNexp, 0, 0, seedY );
        xr = xd;
    }
    else {
        assert(0);
        /* The vectors are distributed (TODO: Fix initialization based on spm instead of bcsc) */
        z_init( pastix_data, seedX, nrhs, xd, spm->nexp );
        z_init( pastix_data, seedY, nrhs, yd, spm->nexp );

        xr = (pastix_complex64_t*)malloc( spm->gNexp * sizeof(pastix_complex64_t) );
        core_zplrnt( spm->gNexp, nrhs, xr, spm->gNexp, spm->gNexp, 0, 0, seedX );
    }

    core_zplrnt( spm->gNexp, nrhs, yr, spm->gNexp, spm->gNexp, 0, 0, seedY );

    if ( trans == SpmNoTrans ) {
        Anorm = spmNorm( SpmInfNorm, spm );
    }
    else {
        Anorm = spmNorm( SpmOneNorm, spm );
    }
    Xnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1, xr, spm->gNexp );
    Ynorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, 1, yr, spm->gNexp );

    /* Compute the spm matrix-vector product */
    {
        rc = spmMatVec( trans, alpha, spm, xr, beta, yr );
        if ( rc != SPM_SUCCESS ) {
            info_solution = 1;
            goto end;
        }
    }

    /* Compute the bcsc matrix-vector product */
    {
        pastix_rhs_t Pxd;
        pastix_rhs_t Pyd;

        pastixRhsInit( &Pxd );
        pastixRhsInit( &Pyd );
        bvec_zlapmr( pastix_data, PastixDirForward, spm->nexp, nrhs, xd, spm->nexp, Pxd );
        bvec_zlapmr( pastix_data, PastixDirForward, spm->nexp, nrhs, yd, spm->nexp, Pyd );

        bcsc_zspmv( pastix_data, (pastix_trans_t)trans, alpha, Pxd->b, beta, Pyd->b );

        bvec_zlapmr( pastix_data, PastixDirBackward, spm->nexp, nrhs, xd, spm->nexp, Pxd );
        bvec_zlapmr( pastix_data, PastixDirBackward, spm->nexp, nrhs, yd, spm->nexp, Pyd );
        pastixRhsFinalize( Pxd );
        pastixRhsFinalize( Pyd );
    }

    info_solution = 0;
    Yrnorm = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->gNexp, nrhs, yr, spm->gNexp ); /* Spm computation  */
    Ydnorm = spmNormMat( SpmInfNorm, spm, nrhs, yd, spm->nexp );                        /* bcsc computation */

    /* Compute the max norm of the difference of the two computations */
    {
        pastix_complex64_t *yref = yr;
        pastix_complex64_t *ynew = yd;
        pastix_int_t        i, j;
        double              val;

        Rnorm = 0.;

        if ( spm->loc2glob == NULL ) {
            for( j=0; j<nrhs; j++ ) {
                for( i=0; i<spm->gNexp; i++, ynew++, yref++ ) {
                    val = cabs(*yref - *ynew);

                    if ( val > Rnorm ) {
                        Rnorm = val;
                    }
                }
            }
        }
        else {
            pastix_int_t *loc2glob = spm->loc2glob;
            pastix_int_t  ig, baseval = spm->baseval;

            for( j=0; j<nrhs; j++ ) {
                for( i=0; i<spm->nexp; i++, ynew++, loc2glob++ ) {
                    ig  = *loc2glob - baseval;
                    val = cabs(yref[ig] - *ynew);

                    if ( val > Rnorm ) {
                        Rnorm = val;
                    }
                }
            }

#if defined(PASTIX_WITH_MPI)
            MPI_Allreduce( MPI_IN_PLACE, &Rnorm, 1, MPI_DOUBLE,
                           MPI_MAX, pastix_data->inter_node_comm );
#endif
        }
    }

    if ( 1 ) {
        printf("[%2d]  ||A||_inf = %e, ||x||_inf = %e, ||y||_inf = %e\n"
               "[%2d]  ||spm(a*A*x+b*y)||_inf = %e, ||bcsc(a*A*x+b*y)||_inf = %e, ||R||_m = %e\n",
               pastix_data->procnum, Anorm, Xnorm, Ynorm,
               pastix_data->procnum, Yrnorm, Ydnorm, Rnorm );
    }

    /*
     * Use the average degree as a limit, the exact max degree should be computed to be correct
     */
    result = Rnorm / ((Anorm + Xnorm + Ynorm) * ((double)spm->gnnzexp / (double)spm->gNexp) * eps);
    if (  isinf(Ydnorm) || isinf(Yrnorm) ||
          isnan(result) || isinf(result) || (result > 10.0) ) {
        info_solution = 1;
    }
    else {
        info_solution = 0;
    }

  end:
#if defined(PASTIX_WITH_MPI)
    MPI_Allreduce( MPI_IN_PLACE, &info_solution, 1, MPI_INT,
                   MPI_SUM, pastix_data->inter_node_comm );
#endif

    if ( spm->loc2glob != NULL ) {
        free( xr );
    }
    free(xd); free(yr); free(yd);

    return info_solution;
}

int
z_bcsc_norm_check( const spmatrix_t   *spm, const pastix_bcsc_t *bcsc )
{
    double norms, normd;
    double eps, result;
    int ret = 0;

    eps = LAPACKE_dlamch_work('e');

    /**
     * Test Norm Max
     */
    printf(" -- Test norm Max :");
    norms = spmNorm( SpmMaxNorm, spm );
    normd = bcsc_znorm( PastixMaxNorm, bcsc );
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
    norms = spmNorm( SpmInfNorm, spm );
    normd = bcsc_znorm( PastixInfNorm, bcsc );
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
    norms = spmNorm( SpmOneNorm, spm );
    normd = bcsc_znorm( PastixOneNorm, bcsc );
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
    norms = spmNorm( SpmFrobeniusNorm, spm );
    normd = bcsc_znorm( PastixFrobeniusNorm, bcsc );
    result = fabs(norms - normd) / (norms * eps);
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
