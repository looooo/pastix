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
#include <csc.h>
#include <bcsc.h>
#include <lapacke.h>
#include <z_spm.h>
#include <z_bcsc.h>


/*------------------------------------------------------------------------
 *  Check the accuracy of the solution
 */
int
z_bcsc_norm_check( const pastix_csc_t *spm, const pastix_bcsc_t *bcsc )
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
    result = fabs(norms - normd) / (normd * eps);

    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        ret++;
    }

    printf("   Nspm = %e, Nbcsc = %e\n", norms, normd );
    printf("  | Nspm - Nbcsc | / Nbcsc = %e\n", result);

    /**
     * Test Norm Inf
     */
    printf(" -- Test norm Inf :");
    norms = spmNorm( PastixInfNorm, spm );
    normd = z_bcscNorm( PastixInfNorm, bcsc );
    result = fabs(norms - normd) / (normd * eps);
    result = result * ((double)(spm->gnnz)) / ((double)(spm->gN));

    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        ret++;
    }

    printf("   Nspm = %e, Nbcsc = %e\n", norms, normd );
    printf("  | Nspm - Nbcsc | / Nbcsc = %e\n", result);

    /**
     * Test Norm One
     */
    printf(" -- Test norm One :");
    norms = spmNorm( PastixOneNorm, spm );
    normd = z_bcscNorm( PastixOneNorm, bcsc );
    result = fabs(norms - normd) / (normd * eps);
    result = result * ((double)(spm->gnnz)) / ((double)(spm->gN));

    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        ret++;
    }

    printf("   Nspm = %e, Nbcsc = %e\n", norms, normd );
    printf("  | Nspm - Nbcsc | / Nbcsc = %e\n", result);

    /**
     * Test Norm Frobenius
     */
    printf(" -- Test norm Frb :");
    norms = spmNorm( PastixFrobeniusNorm, spm );
    normd = z_bcscNorm( PastixFrobeniusNorm, bcsc );
    result = abs(norms - normd) / (normd * eps);
    result = result / ((double)spm->gnnz);

    if ( (result >= 0.) && (result < 1.) ) {
        printf("SUCCESS !\n");
    } else {
        printf("FAILED !\n");
        ret++;
    }

    printf("   Nspm = %e, Nbcsc = %e\n", norms, normd );
    printf("  | Nspm - Nbcsc | / Nbcsc = %e\n", result);

    return ret;
}