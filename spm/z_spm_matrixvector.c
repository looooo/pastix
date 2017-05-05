/**
 *
 * @file z_spm_matrixvector.c
 *
 * SParse Matrix package matrix-vector multiplication routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 **/
#include "common.h"
#include "spm.h"
#include "z_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_matvec
 *
 * @brief compute the matrix-vector product:
 *          y = alpha * op( A ) * x + beta * y
 *
 * A is a PastixGeneral spm, where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = conjg( X' )
 *
 *  alpha and beta are scalars, and x and y are vectors.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix spm is transposed, not transposed or
 *          conjugate transposed:
 *          = PastixNoTrans:   A is not transposed;
 *          = PastixTrans:     A is transposed;
 *          = PastixConjTrans: A is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] spm
 *          The PastixGeneral spm.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
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
z_spmGeCSCv(const pastix_trans_t      trans,
                  pastix_complex64_t  alpha,
            const pastix_spm_t       *spm,
            const pastix_complex64_t *x,
                  pastix_complex64_t  beta,
                  pastix_complex64_t *y )
{
    const pastix_complex64_t *valptr = (pastix_complex64_t*)(spm->values);
    const pastix_complex64_t *xptr   = (const pastix_complex64_t*)x;
    pastix_complex64_t *yptr = (pastix_complex64_t*)y;
    pastix_int_t col, row, i, baseval;

    if ( (spm == NULL) || (x == NULL) || (y == NULL ) )
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    if( spm->mtxtype != PastixGeneral )
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    baseval = spmFindBase( spm );

    /* first, y = beta*y */
    if( beta == 0. ) {
        memset( yptr, 0, spm->gN * sizeof(pastix_complex64_t) );
    }
    else {
        for( i=0; i<spm->gN; i++, yptr++ ) {
            (*yptr) *= beta;
        }
        yptr = y;
    }

    if( alpha != 0. ) {
        /**
         * PastixNoTrans
         */
        if( trans == PastixNoTrans )
        {
            for( col=0; col < spm->gN; col++ )
            {
                for( i=spm->colptr[col]; i<spm->colptr[col+1]; i++ )
                {
                    row = spm->rowptr[i-baseval]-baseval;
                    yptr[row] += alpha * valptr[i-baseval] * xptr[col];
                }
            }
        }
        /**
         * PastixTrans
         */
        else if( trans == PastixTrans )
        {
            for( col=0; col < spm->gN; col++ )
            {
                for( i=spm->colptr[col]; i<spm->colptr[col+1]; i++ )
                {
                    row = spm->rowptr[i-baseval]-baseval;
                    yptr[col] += alpha * valptr[i-baseval] * xptr[row];
                }
            }
        }
#if defined(PRECISION_c) || defined(PRECISION_z)
        else if( trans == PastixConjTrans )
        {
            for( col=0; col < spm->gN; col++ )
            {
                for( i=spm->colptr[col]; i<spm->colptr[col+1]; i++ )
                {
                    row = spm->rowptr[i-baseval]-baseval;
                    yptr[col] += alpha * conj( valptr[i-baseval] ) * xptr[row];
                }
            }
        }
#endif
        else
        {
            return PASTIX_ERR_BADPARAMETER;
        }
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_matvec
 *
 * @brief compute the matrix-vector product:
 *          y = alpha * A + beta * y
 *
 * A is a PastixSymmetric spm, alpha and beta are scalars, and x and y are
 * vectors, and A a symm.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] spm
 *          The PastixSymmetric spm.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
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
z_spmSyCSCv(      pastix_complex64_t  alpha,
            const pastix_spm_t       *spm,
            const pastix_complex64_t *x,
                  pastix_complex64_t  beta,
                  pastix_complex64_t *y )
{
    const pastix_complex64_t *valptr = (pastix_complex64_t*)spm->values;
    const pastix_complex64_t *xptr   = x;
    pastix_complex64_t *yptr = y;
    pastix_int_t col, row, i, baseval;

    if ( (spm == NULL) || (x == NULL) || (y == NULL ) )
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    if( spm->mtxtype != PastixSymmetric )
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    baseval = spmFindBase( spm );

    /* First, y = beta*y */
    if( beta == 0. ) {
        memset( yptr, 0, spm->gN * sizeof(pastix_complex64_t) );
    }
    else {
        for( i=0; i<spm->gN; i++, yptr++ ) {
            (*yptr) *= beta;
        }
        yptr = y;
    }

    if( alpha != 0. ) {
        for( col=0; col < spm->gN; col++ )
        {
            for( i=spm->colptr[col]; i < spm->colptr[col+1]; i++ )
            {
                row = spm->rowptr[i-baseval]-baseval;
                yptr[row] += alpha * valptr[i-baseval] * xptr[col];
                {
                    yptr[col] += alpha * valptr[i-baseval] * xptr[row];
                }
            }
        }
    }

    return PASTIX_SUCCESS;
}

#if defined(PRECISION_c) || defined(PRECISION_z)
/**
 *******************************************************************************
 *
 * @ingroup spm_dev_matvec
 *
 * @brief compute the matrix-vector product:
 *          y = alpha * A + beta * y
 *
 * A is a PastixHermitian spm, alpha and beta are scalars, and x and y are
 * vectors, and A a symm.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] spm
 *          The PastixHermitian spm.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
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
z_spmHeCSCv(      pastix_complex64_t  alpha,
            const pastix_spm_t       *spm,
            const pastix_complex64_t *x,
                  pastix_complex64_t  beta,
                  pastix_complex64_t *y )
{
    const pastix_complex64_t *valptr = (pastix_complex64_t*)spm->values;
    const pastix_complex64_t *xptr   = x;
    pastix_complex64_t *yptr = y;
    pastix_int_t col, row, i, baseval;

    if ( (spm == NULL) || (x == NULL) || (y == NULL ) )
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    if( spm->mtxtype != PastixHermitian )
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    /* First, y = beta*y */
    if( beta == 0. ) {
        memset( yptr, 0, spm->gN * sizeof(pastix_complex64_t) );
    }
    else {
        for( i=0; i<spm->gN; i++, yptr++ ) {
            (*yptr) *= beta;
        }
        yptr = y;
    }

    baseval = spmFindBase( spm );

    if( alpha != 0. ) {
        for( col=0; col < spm->gN; col++ )
        {
            for( i=spm->colptr[col]; i < spm->colptr[col+1]; i++ )
            {
                row=spm->rowptr[i-baseval]-baseval;
                if( col != row ) {
                    yptr[row] += alpha * valptr[i-baseval] * xptr[col];
                    yptr[col] += alpha * conj( valptr[i-baseval] ) * xptr[row];
                }
                else {
                    yptr[row] += alpha * creal(valptr[i-baseval]) * xptr[col];
                }
            }
        }
    }

    return PASTIX_SUCCESS;
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_matvec
 *
 * @brief compute the matrix-vector product:
 *          y = alpha * A + beta * y
 *
 * A is a PastixHermitian spm, alpha and beta are scalars, and x and y are
 * vectors, and A a symm.
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          TODO
 *
 * @param[in] alphaptr
 *          alpha specifies the scalar alpha
 *
 * @param[in] spm
 *          The PastixHermitian spm.
 *
 * @param[in] xptr
 *          The vector x.
 *
 * @param[in] betaptr
 *          beta specifies the scalar beta
 *
 * @param[inout] yptr
 *          The vector y.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the y vector has been computed succesfully,
 * @retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_spmCSCMatVec(const pastix_trans_t  trans,
               const void           *alphaptr,
               const pastix_spm_t   *spm,
               const void           *xptr,
               const void           *betaptr,
                     void           *yptr )
{
    const pastix_complex64_t *x = (const pastix_complex64_t*)xptr;
    pastix_complex64_t *y       = (pastix_complex64_t*)yptr;
    pastix_complex64_t alpha, beta;

    alpha = *((const pastix_complex64_t *)alphaptr);
    beta  = *((const pastix_complex64_t *)betaptr);

    switch (spm->mtxtype) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PastixHermitian:
        return z_spmHeCSCv( alpha, spm, x, beta, y );
#endif
    case PastixSymmetric:
        return z_spmSyCSCv( alpha, spm, x, beta, y );
    case PastixGeneral:
    default:
        return z_spmGeCSCv( trans, alpha, spm, x, beta, y );
    }
}
