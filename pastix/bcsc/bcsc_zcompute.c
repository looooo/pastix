/**
 * 
 * @file bcsc_zcompute.c
 * 
 *  Functions computing operations on the BCSC.
 * 
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
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
#include <lapacke.h>
#include "bcsc.h"
#include "z_bcsc.h"
#include "frobeniusupdate.h"

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
 * @param[out] berr
 *          The returned result.
 *
 *******************************************************************************
 *
 * @return
 *      \retval the error.
 *
 *******************************************************************************/
double
z_bcscBerr( void         *r1,
            void         *r2,
            pastix_int_t  n )
{
    pastix_complex64_t *r1ptr = (pastix_complex64_t*)r1;
    pastix_complex64_t *r2ptr = (pastix_complex64_t*)r2;
    double module1, module2;
    double berr = 0.;
    pastix_int_t i;

    if(r1==NULL || r1== NULL)
        return PASTIX_ERR_BADPARAMETER;

    for( i = 0; i < n; i++)
    {
        module1 = cabs(r1ptr[i]);
        module2 = cabs(r2ptr[i]);
        if( module2 > 0.)
            if( module1 / module2 > berr )
                berr = module1 / module2;
    }

    return berr;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscNormErr - Computes the norm 2 of r1 and the norm 2 of r2
 *                 and return the quotient of these two vectors.
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
 * @return
 *      \retval The error
 *
 *******************************************************************************/
double
z_bcscNormErr( void         *r1,
               void         *r2,
               pastix_int_t  n )
{
    double norm2r1;
    double norm2r2;
    double scale = 0.;
    double sum = 1.;
    double* valptr = (double*)r1;
    int i;

    if(r1==NULL || r2== NULL)
        return PASTIX_ERR_BADPARAMETER;

    for( i = 0; i < n; i++, valptr++ )
    {
        frobenius_update( 1, &scale, &sum, valptr);
#if defined(PRECISION_z) || defined(PRECISION_c)
        valptr++;
        frobenius_update( 1, &scale, &sum, valptr);
#endif
    }
    norm2r1 = scale*sqrt(sum);

    scale = 0.;
    sum=1.;
    valptr = (double*)r2;
    for( i = 0; i < n; i++ )
    {
        frobenius_update( 1, &scale, &sum, valptr);
#if defined(PRECISION_z) || defined(PRECISION_c)
        valptr++;
        frobenius_update( 1, &scale, &sum, valptr);
#endif
    }
    norm2r2 = scale*sqrt(sum);

    return norm2r1/norm2r2;
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
 * @param[in,out] x
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
 * @return
 *      \retval PASTIX_SUCCESS if the x vector has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_bcscScal( void               *x,
            pastix_complex64_t  alpha,
            pastix_int_t        n,
            pastix_int_t        smxnbr )
{
    pastix_complex64_t *xptr = (pastix_complex64_t*)x;
    pastix_int_t i;

    if(x==NULL)
        return PASTIX_ERR_BADPARAMETER;

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
 * z_bcscAxpy - compute y<-alpha*x+y.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          A scalar.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] n
 *          The size of the vectors.
 *
 * @param[in,out] y
 *          The vector y.
 *
 * @param[in] smxnbr
 *          The number of vectors (multi-right-hand-side method).
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the y vector has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_bcscAxpy(pastix_complex64_t  alpha,
           void               *x,
           pastix_int_t        n, /* size of the vectors */
           void               *y,
           pastix_int_t        smxnbr )
{
    pastix_complex64_t *xptr = (pastix_complex64_t*)x;
    pastix_complex64_t *yptr = (pastix_complex64_t*)y;
    pastix_int_t i;

    if(y==NULL || x== NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }
    if( alpha == (pastix_complex64_t)0.0 )
    {
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
 * z_bcscAxpb - compute |A||x| + |b|
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
 *******************************************************************************
 *
 * @return 
 *      \retval The solution
 *
 *******************************************************************************/
double
z_bcscAxpb( const pastix_bcsc_t *bcsc,
            void                *x,
            void                *b )
{
    double normA, normX, normB;

    normX = LAPACKE_zlange( LAPACK_COL_MAJOR, 'F', bcsc->gN, 1, x, bcsc->gN );
    normB = LAPACKE_zlange( LAPACK_COL_MAJOR, 'F', bcsc->gN, 1, b, bcsc->gN );
    normA = z_bcscNorm( PastixFrobeniusNorm, bcsc );

    return (normA * normX + normB);
}


/* Pivot : */
/* r=b-ax */
/* tmp_berr =  max_i(|r_i|/|r'_i|)*/
/* rberror = ||r||/||b|| */
/* Copy */
/* GEMM */

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscDotc - compute the scalar product x.y
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
 * @return
 *      \retval the scalar product of x and y.
 *
 *******************************************************************************/
pastix_complex64_t
z_bcscDotc( void                *x,
            void                *y,
            pastix_int_t         n )
{
    int i;
    pastix_complex64_t *xptr = (pastix_complex64_t*)x;
    pastix_complex64_t *yptr = (pastix_complex64_t*)y;
    pastix_complex64_t r = 0.0;

    for (i=0; i<n; i++, xptr++, yptr++)
    {
#if defined(PRECISION_z) || defined(PRECISION_c)
        r = r + *xptr * conj(*yptr);
#else
        r = r + *xptr * (*yptr);
#endif
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
 * @return
 *      \retval the Frobenius norm of x.
 *
 *******************************************************************************/
double
z_vectFrobeniusNorm( void        *x,
                     pastix_int_t n )
{
    double scale = 0.;
    double sum = 1.;
    double norm;
    double *valptr = (double*)x;
    pastix_int_t i;

    for( i=0; i < n; i++, valptr++ )
    {
        frobenius_update( 1, &scale, &sum, valptr);
#if defined(PRECISION_z) || defined(PRECISION_c)
        valptr++;
        frobenius_update( 1, &scale, &sum, valptr);
#endif
    }

    norm = scale*sqrt(sum);

    return norm;
}

int z_bcscApplyPerm( pastix_int_t m,
                     pastix_int_t n,
                     pastix_complex64_t *A,
                     pastix_int_t lda,
                     pastix_int_t *perm )
{
    pastix_int_t i, j, k;
    (void)lda;

    if ( n == 1 ) {
        pastix_complex64_t tmp;
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

                tmp = A[j];
                A[j] = A[k];
                A[k] = tmp;

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
    }

    return PASTIX_SUCCESS;
}
