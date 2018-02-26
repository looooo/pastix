/**
 *
 * @file z_spm_genrhs.c
 *
 * SParse Matrix package right hand side generators.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 * @precisions normal z -> c s d
 **/
#include "cblas.h"
#include "lapacke.h"
#include "common.h"
#include "solver.h"
#include "spm.h"
#include "z_spm.h"
#include "kernels/pastix_zcores.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define Rnd64_A  6364136223846793005ULL
#define Rnd64_C  1ULL
#define RndF_Mul 5.4210108624275222e-20f
#define RndD_Mul 5.4210108624275222e-20


static pastix_complex64_t mzone = (pastix_complex64_t)-1.;
static pastix_complex64_t zone  = (pastix_complex64_t) 1.;
static pastix_complex64_t zzero = (pastix_complex64_t) 0.;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Random generator from the HPL library
 *
 *******************************************************************************
 *
 * @param[in] n
 *         Number of element to jump over in the generator cycle.
 *
 * @param[in] seed
 *
 *******************************************************************************
 *
 * @retval a random integer value
 *
 *******************************************************************************/
static inline unsigned long long int
Rnd64_jump(unsigned long long int n, unsigned long long int seed ) {
  unsigned long long int a_k, c_k, ran;
  int i;

  a_k = Rnd64_A;
  c_k = Rnd64_C;

  ran = seed;
  for (i = 0; n; n >>= 1, ++i) {
    if (n & 1)
      ran = a_k * ran + c_k;
    c_k *= (a_k + 1);
    a_k *= a_k;
  }

  return ran;
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#if defined(PRECISION_z) || defined(PRECISION_c)
#define NBELEM   2
#else
#define NBELEM   1
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate a vector of random values.
 *
 *******************************************************************************
 *
 * @param[in] scale
 *         Scaling factor for each randomized values.
 *
 * @param[in] m
 *         The number of rows of the tile A. m >= 0.
 *
 * @param[in] n
 *         The number of columns of the tile A. n >= 0.
 *
 * @param[inout] A
 *         On entry, the m-by-n tile to be initialized.
 *         On exit, the tile initialized in the mtxtype format.
 *
 * @param[in] lda
 *         The leading dimension of the tile A. lda >= max(1,m).
 *
 * @param[in] gM
 *         The global number of rows of the full matrix, A is belonging to. gM >= (m0+M).
 *
 * @param[in] m0
 *         The index of the first row of tile A in the full matrix. m0 >= 0.
 *
 * @param[in] n0
 *         The index of the first column of tile A in the full matrix. n0 >= 0.
 *
 * @param[in] seed
 *         The seed used for random generation. Must be the same for
 *         all tiles initialized with this routine.
 *
 ******************************************************************************/
void
z_spmRndVect( double scale, int m, int n, pastix_complex64_t *A, int lda,
              int gM, int m0, int n0, unsigned long long int seed )
{
    pastix_complex64_t *tmp = A;
    int64_t i, j;
    unsigned long long int ran, jump;

    jump = (unsigned long long int)m0 + (unsigned long long int)n0 * (unsigned long long int)gM;

    for (j=0; j<n; ++j ) {
        ran = Rnd64_jump( NBELEM*jump, seed );
        for (i = 0; i < m; ++i) {
            *tmp = (0.5f - ran * RndF_Mul) * scale;
            ran  = Rnd64_A * ran + Rnd64_C;
#ifdef COMPLEX
            *tmp += (I*(0.5f - ran * RndF_Mul)) * scale;
            ran   = Rnd64_A * ran + Rnd64_C;
#endif
            tmp++;
        }
        tmp  += lda-i;
        jump += gM;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Generate nrhs right hand side vectors associated to a given
 * matrix to test a problem with a solver.
 *
 *******************************************************************************
 *
 * @param[in] type
 *          Defines how to compute the vector b.
 *          - PastixRhsOne:  b is computed such that x = 1 [ + I ]
 *          - PastixRhsI:    b is computed such that x = i [ + i * I ]
 *          - PastixRhsRndX: b is computed by matrix-vector product, such that
 *            is a random vector in the range [-0.5, 0.5]
 *          - PastixRhsRndB: b is computed randomly and x is not computed.
 *
 * @param[in] nrhs
 *          Defines the number of right hand side that must be generated.
 *
 * @param[in] spm
 *          The sparse matrix uses to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[out] x
 *          On exit, if x != NULL, then the x vector(s) generated to compute b
 *          is returned. Must be of size at least ldx * spm->n.
 *
 * @param[in] ldx
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= spm->n.
 *
 * @param[inout] b
 *          b must be an allocated matrix of size at least ldb * nrhs.
 *          On exit, b is initialized as defined by the type parameter.
 *
 * @param[in] ldb
 *          Defines the leading dimension of b when multiple right hand sides
 *          are available. ldb >= spm->n.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the b vector has been computed succesfully,
 * @retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_spmGenRHS( pastix_rhstype_t type, int nrhs,
             const pastix_spm_t  *spm,
             void                *x, int ldx,
             void                *b, int ldb )
{
    pastix_complex64_t *xptr = (pastix_complex64_t*)x;
    pastix_complex64_t *bptr = (pastix_complex64_t*)b;
    pastix_int_t i, j;
    int rc;

    if (( spm == NULL ) ||
        ( spm->values == NULL )) {
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Other format not supported for now */
    if( spm->fmttype != PastixCSC ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    if( spm->gN <= 0 ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    if( nrhs <= 0 ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    if( (nrhs > 1) && (ldx < spm->n) ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    if( (nrhs > 1) && (ldb < spm->n) ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    if( spm->dof != 1 ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    if (nrhs == 1) {
        ldb = spm->n;
        ldx = spm->n;
    }

    /* We don't handle distributed spm for now */
    assert( spm->n == spm->gN );

    /* If random b, we do it and exit */
    if ( type == PastixRhsRndB ) {
        /* Compute the spm norm to scale the b vector */
        double norm = z_spmNorm( PastixFrobeniusNorm, spm );

        z_spmRndVect( norm, spm->n, nrhs, bptr, ldb,
                      spm->gN, 0, 0, 24356 );
        return PASTIX_SUCCESS;
    }

    if ( (type == PastixRhsOne  ) ||
         (type == PastixRhsI    ) ||
         (type == PastixRhsRndX ) )
    {
        if ( xptr == NULL ) {
            MALLOC_INTERN( xptr, ldx * nrhs, pastix_complex64_t );
        }

        switch( type ) {
        case PastixRhsOne:
            for( j=0; j<nrhs; j++ )
            {
                for( i=0; i<spm->n; i++, xptr++ )
                {
#if defined(PRECISION_z) || defined(PRECISION_c)
                    *xptr = (pastix_complex64_t)(1.+1.*I);
#else
                    *xptr = (pastix_complex64_t)1.;
#endif
                }
                xptr += ldx-i;
            }
            xptr -= nrhs * ldx;
            break;

        case PastixRhsI:
            for( j=0; j<nrhs; j++ )
            {
                for( i=0; i<spm->n; i++, xptr++ )
                {
#if defined(PRECISION_z) || defined(PRECISION_c)
                    *xptr = (pastix_complex64_t)(i + i * I);
#else
                    *xptr = (pastix_complex64_t)i;
#endif
                }
                xptr += ldx-i;
            }
            xptr -= nrhs * ldx;
            break;

        case PastixRhsRndX:
        default:
            z_spmRndVect( 1., spm->n, nrhs, xptr, ldx,
                          spm->gN, 0, 0, 24356 );
        }

        /* Compute B */
        rc = z_spmCSCMatMat( PastixNoTrans, nrhs, &zone, spm, xptr, ldx, &zzero, bptr, ldb );

        if ( x == NULL ) {
            memFree_null(xptr);
        }
        return rc;
    }

    fprintf(stderr, "z_spmGenRHS: Generator not implemented yet\n");

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_rhs
 *
 * @brief Check the backward error, and the forward error if x0 is provided.
 *
 *******************************************************************************
 *
 * @param[in] nrhs
 *          Defines the number of right hand side that must be generated.
 *
 * @param[in] spm
 *          The sparse matrix uses to generate the right hand side, and the
 *          solution of the full problem.
 *
 * @param[inout] x0
 *          If x0 != NULL, the forward error is computed.
 *          On exit, x0 stores (x0-x)
 *
 * @param[in] ldx0
 *          Defines the leading dimension of x0 when multiple right hand sides
 *          are available. ldx0 >= spm->n.
 *
 * @param[inout] b
 *          b is a matrix of size at least ldb * nrhs.
 *          On exit, b stores Ax-b.
 *
 * @param[in] ldb
 *          Defines the leading dimension of b when multiple right hand sides
 *          are available. ldb >= spm->n.
 *
 * @param[in] x
 *          Contains the solution computed by the solver.
 *
 * @param[in] ldx
 *          Defines the leading dimension of x when multiple right hand sides
 *          are available. ldx >= spm->n.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the tests are succesfull
 * @retval 1, if one of the test failed
 *
 *******************************************************************************/
int
z_spmCheckAxb( int nrhs,
               const pastix_spm_t  *spm,
                     void *x0, int ldx0,
                     void *b,  int ldb,
               const void *x,  int ldx )
{
    const pastix_complex64_t *zx  = (const pastix_complex64_t *)x;
    pastix_complex64_t       *zx0 = (pastix_complex64_t *)x0;
    pastix_complex64_t       *zb  = (pastix_complex64_t *)b;
    double normA, normB, normX, normX0, normR;
    double backward, forward, eps;
    int failure = 0;
    int i;

    assert( spm->nexp == spm->n );
    assert( spm->dof == 1 );

    eps = LAPACKE_dlamch('e');

    /**
     * Compute the starting norms
     */
    normA = spmNorm( PastixOneNorm, spm );

    normB = 0.;
    normX = 0.;
    for( i=0; i<nrhs; i++ ) {
        double norm;

        norm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->n, 1, zb + i * ldb, ldb );
        normB = (norm > normB ) ? norm : normB;
        norm  = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->n, 1, zx + i * ldx, ldx );
        normX = (norm > normX ) ? norm : normX;
    }
    printf( "   || A ||_1                                               %e\n"
            "   max(|| b_i ||_oo)                                       %e\n"
            "   max(|| x_i ||_oo)                                       %e\n",
            normA, normB, normX );

    /**
     * Compute r = b - A * x
     */
    spmMatMat( PastixNoTrans, nrhs, &mzone, spm, x, ldx, &zone, b, ldb );

    normR    = 0.;
    backward = 0.;
    failure  = 0;

    for( i=0; i<nrhs; i++ ) {
        double nx   = cblas_dzasum( spm->n, zx + i * ldx, 1 );
        double nr   = cblas_dzasum( spm->n, zb + i * ldb, 1 );
        double back =  ((nr / normA) / nx) / eps;
        int fail = 0;

        normR    = (nr   > normR   ) ? nr   : normR;
        backward = (back > backward) ? back : backward;

        fail = isnan(nr) || isinf(nr) || isnan(back) || isinf(back) || (back > 1.e3);
        if ( fail ) {
            printf( "   || b_%d - A x_%d ||_1                                 %e\n"
                    "   || b_%d - A x_%d ||_1 / (||A||_1 * ||x_%d||_oo * eps) %e (%s)\n",
                    i, i, nr,
                    i, i, i, back,
                    fail ? "FAILED" : "SUCCESS" );
        }

        failure = failure || fail;
    }

    printf( "   max(|| b_i - A x_i ||_1)                                %e\n"
            "   max(|| b_i - A x_i ||_1 / (||A||_1 * ||x_i||_oo * eps)) %e (%s)\n",
            normR, backward,
            failure ? "FAILED" : "SUCCESS" );

    /**
     * Compute r = x0 - x
     */
    if ( x0 != NULL ) {
        double forw, nr, nx;
        int fail;

        forward = 0.;
        normR   = 0.;
        normX0  = 0.;
        failure = 0;

        for( i=0; i<nrhs; i++, zx += ldx, zx0 += ldx0 ) {

            nx = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->n, 1, zx0, ldx0 );

            cblas_zaxpy( spm->n, CBLAS_SADDR(mzone),
                         zx, 1, zx0, 1);

            nr = LAPACKE_zlange( LAPACK_COL_MAJOR, 'I', spm->n, 1, zx0, ldx0 );

            forw = (nr / nx) / eps;

            normX0  = ( nx   > normX0  ) ? nx   : normX0;
            normR   = ( nr   > normR   ) ? nr   : normR;
            forward = ( forw > forward ) ? forw : forward;

            fail = isnan(nx) || isinf(nx) || isnan(forw) || isinf(forw) || (forw > 1.e3);
            if ( fail ) {
                printf( "   || x0_%d ||_oo                                 %e\n"
                        "   || x0_%d - x_%d ||_oo / (||x0_%d||_oo * eps)   %e (%s)\n",
                        i, nr,
                        i, i, i, forw,
                        fail ? "FAILED" : "SUCCESS" );
            }

            failure = failure || fail;
        }

        printf( "   max(|| x0_i ||_oo)                                      %e\n"
                "   max(|| x0_i - x_i ||_oo)                                %e\n"
                "   max(|| x0_i - x_i ||_oo / || x0_i ||_oo)                %e (%s)\n",
                normX0, normR, forward,
                failure ? "FAILED" : "SUCCESS" );
    }

    return - failure;
}
