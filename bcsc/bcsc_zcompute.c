/**
 *
 * @file bcsc_zcompute.c
 *
 *  Functions computing operations on the BCSC.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Théophile terraz
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include <math.h>
#include "lapacke.h"
#include "bcsc.h"
#include "z_bcsc.h"
#include "frobeniusupdate.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * @brief Compute the Frobenius norm of a vector.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The size of the vector x.
 *
 * @param[in] x
 *          The vector x of size n.
 *
 *******************************************************************************
 *
 * @retval the Frobenius norm of x.
 *
 *******************************************************************************/
double
bvec_znrm2( pastix_int_t              n,
            const pastix_complex64_t *x )
{
    double scale = 0.;
    double sum = 1.;
    double norm;
    double *valptr = (double*)x;
    pastix_int_t i;

    for( i=0; i < n; i++, valptr++ )
    {
        frobenius_update( 1, &scale, &sum, valptr );
#if defined(PRECISION_z) || defined(PRECISION_c)
        valptr++;
        frobenius_update( 1, &scale, &sum, valptr );
#endif
    }

    norm = scale*sqrt(sum);

    return norm;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * @brief Scale a vector by the scalar alpha.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The size of the vector x.
 *
 * @param[in] alpha
 *          The scalar to sclae the vector x.
 *
 * @param[inout] x
 *          The vector x to scale.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the x vector has been computed succesfully,
 * @retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
bvec_zscal( pastix_int_t        n,
            pastix_complex64_t  alpha,
            pastix_complex64_t *x )
{
    pastix_complex64_t *xptr = x;
    pastix_int_t i;

    if( x == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    if( alpha == (pastix_complex64_t)0.0 )
    {
        memset( x, 0., n * sizeof(pastix_complex64_t) );
        return PASTIX_SUCCESS;
    }

    for( i=0; i<n; i++, xptr++ )
    {
        *xptr *= alpha;
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * @brief Compute y <- alpha * x + y.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The size of the vectors.
 *
 * @param[in] alpha
 *          A scalar.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[inout] y
 *          The vector y.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the y vector has been computed succesfully,
 * @retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
bvec_zaxpy( pastix_int_t              n,
            pastix_complex64_t        alpha,
            const pastix_complex64_t *x,
            pastix_complex64_t       *y)
{
    const pastix_complex64_t *xptr = x;
    pastix_complex64_t       *yptr = y;
    pastix_int_t i;

    if( (y == NULL) || (x == NULL) ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    if( alpha == (pastix_complex64_t)0.0 ) {
        return PASTIX_SUCCESS;
    }

    for(i = 0; i < n; i++)
    {
        *yptr = *yptr + alpha * (*xptr);
        yptr++;
        xptr++;
    }

    return PASTIX_SUCCESS;
}

#if defined(PRECISION_z) || defined(PRECISION_c)
/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * @brief Compute the scalar product x.conj(y).
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The size of the vectors.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] y
 *          The vector y.
 *
 *******************************************************************************
 *
 * @retval the scalar product of x and conj(y).
 *
 *******************************************************************************/
pastix_complex64_t
bvec_zdotc( pastix_int_t              n,
            const pastix_complex64_t *x,
            const pastix_complex64_t *y )
{
    int i;
    pastix_complex64_t *xptr = (pastix_complex64_t*)x;
    pastix_complex64_t *yptr = (pastix_complex64_t*)y;
    pastix_complex64_t r = 0.0;

    for (i=0; i<n; i++, xptr++, yptr++)
    {
        r = r + *xptr * conj(*yptr);
    }

    return r;
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * @brief Compute the scalar product x.y.
 *
 *******************************************************************************
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] y
 *          The vector y.
 *
 * @param[in] n
 *          The size of the vectors.
 *
 *******************************************************************************
 *
 * @retval the scalar product of x and y.
 *
 *******************************************************************************/
pastix_complex64_t
bvec_zdotu( pastix_int_t              n,
            const pastix_complex64_t *x,
            const pastix_complex64_t *y )
{
    int i;
    const pastix_complex64_t *xptr = x;
    const pastix_complex64_t *yptr = y;
    pastix_complex64_t r = 0.0;

    for (i=0; i<n; i++, xptr++, yptr++)
    {
        r = r + *xptr * (*yptr);
    }

    return r;
}

int
bvec_zswap( pastix_int_t        m,
            pastix_int_t        n,
            pastix_complex64_t *A,
            pastix_int_t        lda,
            pastix_int_t       *perm )
{
    pastix_complex64_t tmp;
    pastix_int_t i, j, k, jj;

    for(k=0; k<m; k++) {
        i = k;
        j = perm[i];

        /* Cycle already seen */
        if ( j < 0 ) {
            continue;
        }

        /* Mark the i^th element as being seen */
        perm[i] = -j-1;

        while( j != k ) {

            for(jj=0; jj<n; jj++) {
                tmp             = A[j + jj * lda];
                A[j + jj * lda] = A[k + jj * lda];
                A[k + jj * lda] = tmp;
            }

            i = j;
            j = perm[i];
            perm[i] = -j-1;

            assert( (j != i) && (j >= 0) );
        }
    }

    for(k=0; k<m; k++) {
        assert(perm[k] < 0);
        perm[k] = - perm[k] - 1;
    }

    return PASTIX_SUCCESS;
}
