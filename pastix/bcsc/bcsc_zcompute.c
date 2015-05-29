/**
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
/*
  File: bcsc_zcompute.c

  Functions computing operations on the BCSC.

*/

#include "common.h"
#include "bcsc.h"
#include "math.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscGemv - compute the matrix-vector product y=alpha*A**trans*x+beta*y.
 * A is a PastixGeneral bcsc,
 * trans specifies the operation to be performed as follows:
 *              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y
 *              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y
 *              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y
 * x and y are two vectors of size n,
 * alpha and beta are scalars.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          The operation to be performed.
 *
 * @param[in] n
 *          Number of columns of the matrix.
 *
 * @param[in] alpha
 *          A scalar.
 *
 * @param[in] csc
 *          The PastixGeneral csc.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          A scalar.
 *
 * @param[in,out] y
 *          The vector y.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the y vector has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_bcscGemv(int                 trans,
           pastix_int_t        n, /* size of the matrix */
           pastix_complex64_t  alpha,
           pastix_bcsc_t      *bcsc,
           void               *x,
           pastix_complex64_t  beta,
           void               *y )
{
    pastix_complex64_t *Lvalptr = NULL;
//     pastix_complex64_t *Uvalptr = NULL;
    pastix_complex64_t *yptr    = (pastix_complex64_t*)y;
    pastix_complex64_t *xptr    = (pastix_complex64_t*)x;
    pastix_int_t        bloc, col, i, j;

    if(bcsc==NULL || y==NULL || x== NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }
    Lvalptr = bcsc->Lvalues;
//     Uvalptr = bcsc->Uvalues;

    /* first, y = beta*y */
    if( beta != (pastix_complex64_t)0.0 )
    {
        for( col = 0; col < n; col++ )
        {
            yptr[col] *= beta;
        }
    }
    else if( beta == (pastix_complex64_t)0.0 )
    {
        memset(yptr,0.0,n*sizeof(pastix_complex64_t));
    }

    switch (trans) {
    case PastixNoTrans:
        col = 0;
        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
            {
                for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                {
                    yptr[bcsc->rowtab[i]] += alpha * Lvalptr[i] * xptr[col];
                }
                col += 1;
            }
        }
    break;
    case PastixTrans:
        col = 0;
        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
            {
                for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                {
//                     yptr[bcsc->rowtab[i]] += alpha * Uvalptr[i] * xptr[col];
                    yptr[col] += alpha * Lvalptr[i] * xptr[bcsc->rowtab[i]];
                }
                col += 1;
            }
        }
    break;
#if defined(PRECISION_c) || defined(PRECISION_z)
    case PastixConjTrans:
        col = 0;
        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
            {
                for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                {
//                     yptr[bcsc->rowtab[i]] += alpha * conj( Uvalptr[i] ) * xptr[col];
#if defined(PRECISION_z)
                    yptr[col] += alpha * conj( Lvalptr[i] ) * xptr[bcsc->rowtab[i]];
#elseif defined(PRECISION_c)
                    yptr[col] += alpha * conjf( Lvalptr[i] ) * xptr[bcsc->rowtab[i]];
#endif
                }
                col += 1;
            }
        }
    break;
#endif
    default:
        return PASTIX_ERR_BADPARAMETER;
    }

    return PASTIX_SUCCESS;

}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscNormMax - compute the max norm of a general matrix.
 *
 *******************************************************************************
 *
 * @param[in] values
 *          The values array of the matrix.
 *
 * @param[in] n
 *          The number of elements in the array A.
 *
 * @param[out] norm
 *          The norm of the matrix A.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the norm has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_bcscNormMax( void         *values,
               pastix_int_t  n,
               double       *norm )
{
    double temp1;
#if defined(PRECISION_c) || defined(PRECISION_z)
    double temp2;
#endif
    pastix_complex64_t *valptr = values;
    pastix_int_t i;

    *norm = 0.;
    for( i=0; i < n; i++ )
    {
#if defined(PRECISION_c) || defined(PRECISION_z)
        temp1 = pow(creal(valptr[i]),2.);
        temp2 = pow(cimag(valptr[i]),2.);
        temp1 = sqrt(temp1 + temp2);
#else
        temp1 = fabs((double)valptr[i]);
#endif
        if(*norm < temp1)
        {
            *norm = temp1;
        }
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscNorm2 - compute the frobenius norm of a general matrix.
 *
 *******************************************************************************
 *
 * @param[in] values
 *          The values array of the matrix.
 *
 * @param[in] n
 *          The number of elements in the array A.
 *
 * @param[out] norm
 *          The norm of the matrix A.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the norm has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_bcscNorm2( void         *values,
             pastix_int_t  n,
             double       *norm )
{
    double scale = 0.;
    double sum = 1.;
    double temp;
    pastix_complex64_t *valptr = values;
    pastix_int_t i;

    for( i=0; i < n; i++ )
    {
#if defined(PRECISION_c) || defined(PRECISION_z)
        temp = fabs(creal(valptr[i]));
#else
        temp = fabs((double)valptr[i]);
#endif
        if(temp != 0.)
        {
            if(scale < temp)
            {
                sum = 1. + sum*pow((scale / temp), 2.);
                scale = temp;
            }else{
                sum = sum + pow((double)(temp / scale), 2.);
            }
        }
#if defined(PRECISION_c) || defined(PRECISION_z)
        temp = fabs(cimag(valptr[i]));
        if(temp != 0.)
        {
            if(scale < temp)
            {
                sum = 1. + sum*pow((double)(scale / temp), 2.);
                scale = temp;
            }else{
                sum = sum + pow((double)(temp / scale), 2.);
            }
        }
#endif
    }
    *norm = scale*sqrt(sum);

    return PASTIX_SUCCESS;
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
 * @param[out] berr
 *          The returned result.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the berr has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_bcscBerr( void         *r1,
            void         *r2,
            pastix_int_t  n,
            double       *berr )
{
    pastix_complex64_t *r1ptr = (pastix_complex64_t*)r1;
    pastix_complex64_t *r2ptr = (pastix_complex64_t*)r2;
#if defined(PRECISION_c) || defined(PRECISION_z)
    pastix_complex64_t Re1, Re2, Im1, Im3, module1, module2;
#endif
    pastix_int_t i;

    if(r1==NULL || r1== NULL)
        return PASTIX_ERR_BADPARAMETER;

    *berr = 0.;

    for( i = 0; i < n; i++)
    {
#if defined(PRECISION_c) || defined(PRECISION_z)
        Re1 = creal(r1ptr[i]);
        Re2 = creal(r1ptr[i]);
        Im1 = cimag(r2ptr[i]);
        Im2 = cimag(r2ptr[i]);
        module1 = sqrt(Re1*Re1 + Im1*Im1);
        module2 = sqrt(Re2*Re2 + Im2*Im2);
        if( module2 > 0)
            if( module1 / module2 > *berr )
                *berr = module1 / module2;
#else
        if(r2ptr != 0)
            if( fabs(r1ptr[i] / r2ptr[i]) > *berr )
                *berr = fabs(r1ptr[i] / r2ptr[i]);
#endif
    }

    return PASTIX_SUCCESS;
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
 * @param[out] err
 *          The returned result.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the err has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_bcscNormErr( void         *r1,
               void         *r2,
               pastix_int_t  n,
               double       *err )
{
    double              norm2r1;
    double              norm2r2;

    if(r1==NULL || r1== NULL)
        return PASTIX_ERR_BADPARAMETER;

    if( z_bcscNorm2(r1, n, &norm2r1) != PASTIX_SUCCESS )
        return PASTIX_ERR_BADPARAMETER;

    if( z_bcscNorm2(r2, n, &norm2r2) != PASTIX_SUCCESS )
        return PASTIX_ERR_BADPARAMETER;

    *err = norm2r1/norm2r2;

    return PASTIX_SUCCESS;
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
           pastix_int_t        n, /* size of the matrix */
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
        *yptr = *yptr + alpha * *xptr;
        yptr++;
        xptr++;
    }

    return PASTIX_SUCCESS;
}


/* Pivot : */
/* r=b-ax */
/* r'=|A||x|+|b| */
/* tmp_berr =  max_i(|r_i|/|r'_i|)*/
/* rberror = ||r||/||b|| */
/* Copy */
/* GEMM */
