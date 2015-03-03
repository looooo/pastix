/**
 *
 * @file z_spm_matrixvector.c
 *
 *  PaStiX csc routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 **/
#include "common.h"
#include "csc.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * z_spmGeCSCv - compute y=alpha*A*x+beta*y.
 * A is a PastixGeneral csc, 
 * x and y are two vectors of size csc->gN,
 * alpha and beta are scalars.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          A scalar.
 * 
 * @param[in] csc
 *          The PastixGeneral csc.
 * 
 * @param[in] beta
 *          A scalar.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in,out] y
 *          The vector y, can be unallocated at enter.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the y vector has been computed succesfully,
 *      \retval PASTIX_ERR_MATRIX if the matrix is not a PastixGeneral csc.
 *
 *******************************************************************************/
int
z_spmGeCSCv(pastix_complex64_t  alpha,
            pastix_csc_t       *csc  ,  
            pastix_complex64_t  beta ,
            pastix_complex64_t *x    ,
            pastix_complex64_t *y     )
{
    pastix_complex64_t *valptr  = (pastix_complex64_t*)csc->avals;
    pastix_complex64_t *bptr    = y;
    pastix_complex64_t *xptr    = x;
    pastix_int_t        col, row, i, baseval;

    if(csc==NULL)
    {
        return PASTIX_ERR_MATRIX;
    }

    if(csc->mtxtype!=PastixGeneral || x==NULL)
    {
        return PASTIX_ERR_MATRIX;
    }

    baseval = pastix_imin( *(csc->colptr), *(csc->rows) );

    for( i=0; i < csc->gN; i++ )
    {
        bptr[i] *= beta;
    }
    
    for( col=0; col < csc->gN; col++ )
    {
        for( i=csc->colptr[col]; i<csc->colptr[col+1]; i++ )
        {
            row = csc->rows[i-baseval]-baseval;
            bptr[row] += alpha * valptr[i-baseval] * xptr[col];
        }
    }

    return PASTIX_SUCCESS;

}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * z_spmSyCSCv - compute y=alpha*A*x+beta*y.
 * A is a PastixSymmetric csc, 
 * x and y are two vectors of size csc->gN,
 * alpha and beta are scalars.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          A scalar.
 * 
 * @param[in] csc
 *          The PastixSymmetric csc.
 * 
 * @param[in] beta
 *          A scalar.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in,out] y
 *          The vector y, can be unallocated at enter.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the y vector has been computed succesfully,
 *      \retval PASTIX_ERR_MATRIX if the matrix is not a PastixSymmetric csc.
 *
 *******************************************************************************/
int
z_spmSyCSCv(pastix_complex64_t  alpha,
            pastix_csc_t       *csc  ,  
            pastix_complex64_t  beta ,
            pastix_complex64_t *x    ,
            pastix_complex64_t *y     )
{
    pastix_complex64_t *valptr  = (pastix_complex64_t*)csc->avals;
    pastix_complex64_t *bptr    = y;
    pastix_complex64_t *xptr    = x;
    pastix_int_t        col, row, i, baseval;

    if(csc==NULL)
    {
        return PASTIX_ERR_MATRIX;
    }

    if(csc->mtxtype!=PastixSymmetric || x==NULL)
    {
        return PASTIX_ERR_MATRIX;
    }

    baseval = pastix_imin( *(csc->colptr), *(csc->rows) );

    for( i=0; i < csc->gN; i++ )
    {
        bptr[i] *= beta;
    }

    for( col=0; col < csc->gN; col++ )
    {
        for( i=csc->colptr[col]; i < csc->colptr[col+1]; i++ )
        {
            row = csc->rows[i-baseval]-baseval;
            bptr[row] += alpha * valptr[i-baseval] * xptr[col];
            if( col != row )
            {
                bptr[col] += alpha * valptr[i-baseval] * xptr[row];
            }
        }
    }

    return PASTIX_SUCCESS;

}

#if !defined(PRECISION_s) && !defined(PRECISION_d)
/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * z_spmGeCSCv - compute y=alpha*A*x+beta*y.
 * A is a PastixHermitian csc, 
 * x and y are two vectors of size csc->gN,
 * alpha and beta are scalars.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          A scalar.
 * 
 * @param[in] csc
 *          The PastixHermitian csc.
 * 
 * @param[in] beta
 *          A scalar.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in,out] y
 *          The vector y, can be unallocated at enter.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the y vector has been computed succesfully,
 *      \retval PASTIX_ERR_MATRIX if the matrix is not a PastixHermitian csc.
 *
 *******************************************************************************/
int
z_spmHeCSCv(pastix_complex64_t  alpha,
            pastix_csc_t       *csc  ,  
            pastix_complex64_t  beta ,
            pastix_complex64_t *x    ,
            pastix_complex64_t *y     )
{
    pastix_complex64_t *valptr  = (pastix_complex64_t*)csc->avals;
    pastix_complex64_t *bptr    = y;
    pastix_complex64_t *xptr    = x;
    pastix_int_t        col, row, i, baseval;

    if(csc==NULL)
    {
        return PASTIX_ERR_MATRIX;
    }

    if(csc->mtxtype!=PastixHermitian || x==NULL)
    {
        return PASTIX_ERR_MATRIX;
    }

    for( i=0; i < csc->gN; i++ )
    {
        bptr[i] *= beta;
    }

    baseval = pastix_imin( *(csc->colptr), *(csc->rows) );

    for( col=0; col < csc->gN; col++ )
    {
        for( i=csc->colptr[col]; i < csc->colptr[col+1]; i++ )
        {
            row=csc->rows[i-baseval]-baseval;
            bptr[row] += alpha*valptr[i-baseval]*xptr[col];
            if( col != row )
                bptr[col] += alpha*( creal( valptr[i-baseval] ) - cimag(valptr[i-baseval])*I )*xptr[row];
        }
    }

    return PASTIX_SUCCESS;
  
}
#endif