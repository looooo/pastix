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
 * z_bcscNormErr - Computes the norm 2 of r and the norm 2 of b
 *                 and return the quotient of these two values:
 *
 *                         || r ||_2 / ||b||_2
 *
 *******************************************************************************
 *
 * @param[in] r
 *          The vector r.
 *
 * @param[in] b
 *          The vector b.
 *
 * @param[in] n
 *          The size of the vectors.
 *
 *******************************************************************************
 *
 * @retval the quotient.
 *
 *******************************************************************************/
double
z_bcscNormErr( pastix_complex64_t *r,
               pastix_complex64_t *b,
               pastix_int_t        n )
{
    double norm1, norm2;

    norm1 = z_vectFrobeniusNorm( n, r );
    norm2 = z_vectFrobeniusNorm( n, b );

    return norm1 / norm2;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscBerr - Compute the operation $$ berr= max_{i}(\\frac{|r1_{i}|}{|r2_{i}|}) $$.
 *
 *******************************************************************************
 *
 * @param[in] r1
 *          The vector r1.
 *
 * @param[in] r2
 *          The vector r2.
 *
 * @param[in] n
 *          The size of the vectors.
 *
 *******************************************************************************
 *
 * @retval the error.
 *
 *******************************************************************************/
double
z_bcscBerr( pastix_complex64_t *r1,
            pastix_complex64_t *r2,
            pastix_int_t        n )
{
    pastix_complex64_t *r1ptr = (pastix_complex64_t*)r1;
    pastix_complex64_t *r2ptr = (pastix_complex64_t*)r2;
    double berr = 0.;
    pastix_int_t i;

    for( i = 0; i < n; i++)
    {
        double module1 = cabs(r1ptr[i]);
        double module2 = cabs(r2ptr[i]);
        if( (module2 > 0.) && ((module1 / module2) > berr) )
            berr = module1 / module2;
    }

    return berr;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscScal - Multiply a vector by a scalaire x <- alpha*x.
 *
 *******************************************************************************
 *
 * @param[inout] x
 *          The vector x.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] n
 *          The size of the vectors.
 *
 * @param[in] smxnbr
 *          The number of vectors (multi-right-hand-side method).
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the x vector has been computed succesfully,
 * @retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_bcscScal( pastix_complex64_t *x,
            pastix_complex64_t  alpha,
            pastix_int_t        n,
            pastix_int_t        smxnbr )
{
    pastix_complex64_t *xptr = (pastix_complex64_t*)x;
    pastix_int_t i;

    if( x == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    if( alpha == (pastix_complex64_t)0.0 )
    {
        memset(xptr,0.0,smxnbr*n*sizeof(pastix_complex64_t));
        return PASTIX_SUCCESS;
    }

    for( i = 0; i < n*smxnbr; i++)
    {
        *xptr *=alpha;
        xptr ++;
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscAxpy - Compute y <- alpha * x + y.
 * TODO: Look at zgeadd
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The size of the vectors.
 *
 * @param[in] smxnbr
 *          The number of vectors (multi-right-hand-side method).
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
z_bcscAxpy( pastix_int_t              n,
            pastix_int_t              smxnbr,
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

    for(i = 0; i < n*smxnbr; i++)
    {
        *yptr = *yptr + alpha * (*xptr);
        yptr++;
        xptr++;
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscAxpb - compute r = |A||x| + |b|
 * TODO: Check that we needs it and if yes, rename it
 *
 *******************************************************************************
 *
 * @param[in] A
 *          The Pastix bcsc.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] b
 *          The vector b.
 *
 * @param[out] r
 *          The result.
 *******************************************************************************/
void
z_bcscAxpb( pastix_trans_t       trans,
            const pastix_bcsc_t *bcsc,
            pastix_complex64_t  *x,
            pastix_complex64_t  *b,
            pastix_complex64_t  *r )
{
    pastix_complex64_t *Lvalptr = NULL;
    pastix_complex64_t *xptr    = (pastix_complex64_t*)x;
    pastix_complex64_t *bptr    = (pastix_complex64_t*)b;
    pastix_complex64_t *rptr    = (pastix_complex64_t*)r;
    pastix_int_t        bloc, col, i, j, n;

    Lvalptr = (pastix_complex64_t*)bcsc->Lvalues;
    n = bcsc->n;

    switch (trans) {
#if defined(PRECISION_c) || defined(PRECISION_z)
    case PastixConjTrans:
        col = 0;
        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
            {
                for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                {
                    rptr[col] += cabs( conj( Lvalptr[i] ) ) * cabs( xptr[bcsc->rowtab[i]] );
                }
                col += 1;
            }
        }
    break;
#endif
    case PastixTrans:
        col = 0;
        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
            {
                for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                {
                    rptr[col] += cabs( Lvalptr[i] ) * cabs( xptr[bcsc->rowtab[i]] );
                }
                col += 1;
            }
        }
    break;

    case PastixNoTrans:
    default:
        col = 0;
        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
            {
                for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                {
                    rptr[bcsc->rowtab[i]] += cabs( Lvalptr[i] ) * cabs( xptr[col] );
                }
                col += 1;
            }
        }
    }

    for( i=0; i<n; i++, rptr++, bptr++)
        *rptr += cabs( *bptr );
}

#if defined(PRECISION_z) || defined(PRECISION_c)
/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscDotc - compute the scalar product x.conj(y).
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
 * @retval the scalar product of x and conj(y).
 *
 *******************************************************************************/
pastix_complex64_t
z_bcscDotc( pastix_int_t              n,
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
 * z_bcscDotu - compute the scalar product x.y.
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
z_bcscDotu( pastix_int_t               n,
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

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_vectFrobeniusNorm - compute the Frobenius norm of a vector.
 *
 *******************************************************************************
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] n
 *          The size of the vector x.
 *
 *******************************************************************************
 *
 * @retval the Frobenius norm of x.
 *
 *******************************************************************************/
double
z_vectFrobeniusNorm( pastix_int_t              n,
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

int z_bcscApplyPerm( pastix_int_t        m,
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
