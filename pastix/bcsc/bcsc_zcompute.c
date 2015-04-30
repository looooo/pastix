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
 * x and y are two vectors of size csc->gN,
 * alpha and beta are scalars.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          The operation to be performed.
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
z_bcscGemv(char                trans,
           pastix_complex64_t  alpha,
           pastix_bcsc_t      *bcsc,
           pastix_complex64_t *x,
           pastix_complex64_t  beta,
           pastix_complex64_t *y )
{
    pastix_complex64_t *valptr  = (pastix_complex64_t*)bcsc->Lvalues;
    pastix_complex64_t *yptr    = y;
    pastix_complex64_t *xptr    = x;
    pastix_int_t        bloc, col, row, i, j, baseval;

    if(bcsc==NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    if(bcsc->mtxtype!=PastixGeneral || x==NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    baseval = bcsc->colptr[0];

    /* first, y = beta*y */
    for( i=0; i < csc->gN; i++ )
    {
        yptr[i] *= beta;
    }

    if( trans == 'n' || trans == 'N' )
    {
        col = 0;
        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc]->colnbr; j++ )
            {
                col += 1;
                for( i=bcsc->coltab[col]; i<csc->coltab[col+1]; i++ )
                {
                    row = bcsc->rowtab[i];
                    yptr[row] += alpha * valptr[i] * xptr[col];
                }
            }
        }
    }
    else if( trans == 't' || trans == 'T' )
    {
        col = 0;
        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc]->colnbr; j++ )
            {
                col += 1;
                for( i=bcsc->coltab[col]; i<csc->coltab[col+1]; i++ )
                {
                    row = bcsc->rowtab[i];
                    yptr[col] += alpha * valptr[i] * xptr[row];
                }
            }
        }
    }
#if defined(PRECISION_c) || defined(PRECISION_z)
    else if( trans == 'c' || trans == 'C' )
    {
        col = 0;
        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc]->colnbr; j++ )
            {
                col += 1;
                for( i=bcsc->coltab[col]; i<csc->coltab[col+1]; i++ )
                {
                    row = bcsc->rowtab[i];
                    yptr[col] += alpha * conj( valptr[i-baseval] ) * xptr[row];
                }
            }
        }
    }
#endif
    else
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    return PASTIX_SUCCESS;

}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscNorm1 - compute the norme norm of a general matrix.
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
 *      \retval PASTIX_SUCCESS if the y vector has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_bcscNorm1( void         *values,
             pastix_int_t  n,
             double        norm )
{
    double scale = 0.;
    double sum = 1.;
    double temp;
    pastix_complex64_t valptr = values;
    pastix_int_t i, j, col, bloc, n;
    
    col = 0;
    for( i=0, i < n, i++ )
    {
#if defined(PRECISION_c) || defined(PRECISION_z)
        temp = fabs(creal(valptr[i]));
#else
        temp = fabs((double)valptr[i]);
#endif
        if(scale < temp)
        {
            sum = 1. + sum*(scale / temp);
            scale = temp;
        }else{
            sum = sum + (temp / scale);
        }
#if defined(PRECISION_c) || defined(PRECISION_z)
        temp = fabs(cimag(valptr[i]));
        if(scale < temp)
        {
            sum = 1. + sum*(scale / temp);
            scale = temp;
        }else{
            sum = sum + (temp / scale);
        }
#endif
        norm = scale*sum;
    }
    
    return PASTIX_SUCCESS
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
 *      \retval PASTIX_SUCCESS if the y vector has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_bcscNorm2( void         *values,
             pastix_int_t  n,
             double        norm )
{
    double scale = 0.;
    double sum = 1.;
    double temp;
    pastix_complex64_t valptr = values;
    pastix_int_t i, j, col, bloc, n;
    
    col = 0;
    for( i=0, i < n, i++ )
    {
#if defined(PRECISION_c) || defined(PRECISION_z)
        temp = fabs(creal(valptr[i]));
#else
        temp = fabs((double)valptr[i]);
#endif
        if(scale < temp)
        {
            sum = 1. + sum*pow((scale / temp), 2.);
            scale = temp;
        }else{
            sum = sum + pow((double)(temp / scale), 2.);
        }
#if defined(PRECISION_c) || defined(PRECISION_z)
        temp = fabs(cimag(valptr[i]));
        if(scale < temp)
        {
            sum = 1. + sum*pow((double)(scale / temp), 2.);
            scale = temp;
        }else{
            sum = sum + pow((double)(temp / scale), 2.);
        }
#endif
        norm = scale*sqrt(sum);
    }
    
    return PASTIX_SUCCESS
}






/* Pivot : */
/* r=b-ax */
/* r'=|A||x|+|b| */
/* tmp_berr =  max_i(|r_i|/|r'_i|)*/
/* rberror = ||r||/||b|| */
/* Copy */
/* GEMM */