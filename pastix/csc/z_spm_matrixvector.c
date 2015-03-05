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
 * z_spmGeCSCv - compute the matrix-vector product y=alpha*A**trans*x+beta*y.
 * A is a PastixGeneral csc, 
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
z_spmGeCSCv(char                trans,
            pastix_complex64_t  alpha,
            pastix_csc_t       *csc,
            pastix_complex64_t *x,
            pastix_complex64_t  beta,
            pastix_complex64_t *y )
{
    pastix_complex64_t *valptr  = (pastix_complex64_t*)csc->avals;
    pastix_complex64_t *yptr    = y;
    pastix_complex64_t *xptr    = x;
    pastix_int_t        col, row, i, baseval;

    if(csc==NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    if(csc->mtxtype!=PastixGeneral || x==NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    baseval = pastix_imin( *(csc->colptr), *(csc->rows) );

    /* first, y = beta*y */
    for( i=0; i < csc->gN; i++ )
    {
        yptr[i] *= beta;
    }

    if( trans == 'n' || trans == 'N' )
    {
        for( col=0; col < csc->gN; col++ )
        {
            for( i=csc->colptr[col]; i<csc->colptr[col+1]; i++ )
            {
                row = csc->rows[i-baseval]-baseval;
                yptr[row] += alpha * valptr[i-baseval] * xptr[col];
            }
        }
    }
    else if( trans == 't' || trans == 'T' )
    {
        for( col=0; col < csc->gN; col++ )
        {
            for( i=csc->colptr[col]; i<csc->colptr[col+1]; i++ )
            {
                row = csc->rows[i-baseval]-baseval;
                yptr[col] += alpha * valptr[i-baseval] * xptr[row];
            }
        }
    }
#if defined(PRECISION_c) || defined(PRECISION_z)
    else if( trans == 'c' || trans == 'C' )
    {
        for( col=0; col < csc->gN; col++ )
        {
            for( i=csc->colptr[col]; i<csc->colptr[col+1]; i++ )
            {
                row = csc->rows[i-baseval]-baseval;
                yptr[col] += alpha * conj( valptr[i-baseval] ) * xptr[row];
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
 * @ingroup pastix_csc
 *
 * z_spmSyCSCv - compute the matrix-vector product y=alpha*A**trans*x+beta*y.
 * A is a PastixSymetric csc, 
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
z_spmSyCSCv(pastix_complex64_t  alpha,
            pastix_csc_t       *csc,
            pastix_complex64_t *x,
            pastix_complex64_t  beta,
            pastix_complex64_t *y )
{
    pastix_complex64_t *valptr  = (pastix_complex64_t*)csc->avals;
    pastix_complex64_t *yptr    = y;
    pastix_complex64_t *xptr    = x;
    pastix_int_t        col, row, i, baseval;

    if(csc==NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    if(csc->mtxtype!=PastixSymmetric || x==NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    baseval = pastix_imin( *(csc->colptr), *(csc->rows) );

    for( i=0; i < csc->gN; i++ )
    {
        yptr[i] *= beta;
    }

    for( col=0; col < csc->gN; col++ )
    {
        for( i=csc->colptr[col]; i < csc->colptr[col+1]; i++ )
        {
            row = csc->rows[i-baseval]-baseval;
            yptr[row] += alpha * valptr[i-baseval] * xptr[col];
            if( col != row )
            {
                yptr[col] += alpha * valptr[i-baseval] * xptr[row];
            }
        }
    }

    return PASTIX_SUCCESS;

}

#if defined(PRECISION_c) || defined(PRECISION_z)
/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * z_spmHeCSCv - compute the matrix-vector product y=alpha*A**trans*x+beta*y.
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
z_spmHeCSCv(pastix_complex64_t  alpha,
            pastix_csc_t       *csc,
            pastix_complex64_t *x,
            pastix_complex64_t  beta,
            pastix_complex64_t *y )
{
    pastix_complex64_t *valptr  = (pastix_complex64_t*)csc->avals;
    pastix_complex64_t *yptr    = y;
    pastix_complex64_t *xptr    = x;
    pastix_int_t        col, row, i, baseval;

    if(csc==NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    if(csc->mtxtype!=PastixHermitian || x==NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    for( i=0; i < csc->gN; i++ )
    {
        yptr[i] *= beta;
    }

    baseval = pastix_imin( *(csc->colptr), *(csc->rows) );

    for( col=0; col < csc->gN; col++ )
    {
        for( i=csc->colptr[col]; i < csc->colptr[col+1]; i++ )
        {
            row=csc->rows[i-baseval]-baseval;
            yptr[row] += alpha * valptr[i-baseval] * xptr[col];
            if( col != row )
                yptr[col] += alpha * conj( valptr[i-baseval] ) * xptr[row];
        }
    }

    return PASTIX_SUCCESS;
  
}
#endif