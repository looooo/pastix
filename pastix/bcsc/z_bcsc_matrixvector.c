/**
 *
 * @file z_bcsc_matrixvector.c
 *
 * Functions computing matrix-vector products for the BCSC
 *
 * PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 * LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Théophile terraz
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s
 *
 **/
#include <math.h>
#include "common.h"
#include "bcsc.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscGemv - compute the matrix-vector product y=alpha*A**trans*x+beta*y.
 * A is a PastixGeneral bcsc,
 * trans specifies the operation to be performed as follows:
 *              trans = PastixNoTrans   y := alpha*A*x + beta*y
 *              trans = PastixTrans     y := alpha*A**T*x + beta*y
 *              trans = PastixConjTrans y := alpha*A**H*x + beta*y
 * x and y are two vectors of size n,
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
 * @param[in] bcsc
 *          The PastixGeneral bcsc.
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
z_bcscGemv(pastix_trans_t      trans,
           pastix_complex64_t  alpha,
           pastix_bcsc_t      *bcsc,
           void               *x,
           pastix_complex64_t  beta,
           void               *y )
{
    pastix_complex64_t *Lvalptr = NULL;
    pastix_complex64_t *yptr    = (pastix_complex64_t*)y;
    pastix_complex64_t *xptr    = (pastix_complex64_t*)x;
    pastix_int_t        bloc, col, i, j, n;

    if(bcsc==NULL || y==NULL || x== NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }
    Lvalptr = (pastix_complex64_t*)bcsc->Lvalues;
    n = bcsc->n;

    /* first, y = beta*y */
    if( beta != (pastix_complex64_t)0.0 )
    {
        for( col = 0; col < n; col++, yptr++ )
        {
            (*yptr) *= beta;
        }
    }
    else if( beta == (pastix_complex64_t)0.0 )
    {
        memset(yptr, 0, n * sizeof(pastix_complex64_t));
    }

    yptr = (pastix_complex64_t*)y;
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
                    yptr[col] += alpha * conj( Lvalptr[i] ) * xptr[bcsc->rowtab[i]];
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
                    yptr[col] += alpha * Lvalptr[i] * xptr[bcsc->rowtab[i]];
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
                    yptr[bcsc->rowtab[i]] += alpha * Lvalptr[i] * xptr[col];
                }
                col += 1;
            }
        }
    }

    return PASTIX_SUCCESS;

}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscSymv - compute the matrix-vector product y=alpha*A*x+beta*y.
 * A is a PastixSymmetric bcsc,
 * x and y are two vectors of size n,
 * alpha and beta are scalars.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          A scalar.
 *
 * @param[in] bcsc
 *          The PastixSymmetric bcsc.
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
z_bcscSymv(pastix_complex64_t  alpha,
           pastix_bcsc_t      *bcsc,
           void               *x,
           pastix_complex64_t  beta,
           void               *y )
{
    pastix_complex64_t *Lvalptr = NULL;
    pastix_complex64_t *yptr    = (pastix_complex64_t*)y;
    pastix_complex64_t *xptr    = (pastix_complex64_t*)x;
    pastix_int_t        bloc, col, i, j, n;

    if(bcsc==NULL || y==NULL || x== NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }
    if(bcsc->mtxtype != PastixSymmetric)
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    Lvalptr = (pastix_complex64_t*)bcsc->Lvalues;
    n = bcsc->n;

    /* first, y = beta*y */
    if( beta != (pastix_complex64_t)0.0 )
    {
        for( col = 0; col < n; col++, yptr++ )
        {
            (*yptr) *= beta;
        }
    }
    else if( beta == (pastix_complex64_t)0.0 )
    {
        memset(yptr,0,n*sizeof(pastix_complex64_t));
    }

    yptr = (pastix_complex64_t*)y;
    col = 0;
    for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
    {
        for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
        {
            for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
            {
                yptr[bcsc->rowtab[i]] += alpha * Lvalptr[i] * xptr[col];
                if( bcsc->rowtab[i] != col )
                    yptr[col] += alpha * Lvalptr[i] * xptr[bcsc->rowtab[i]];
            }
            col += 1;
        }
    }

    return PASTIX_SUCCESS;

}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscHemv - compute the matrix-vector product y=alpha*A*x+beta*y.
 * A is a PastixHermitian bcsc,
 * x and y are two vectors of size n,
 * alpha and beta are scalars.
 *
 *******************************************************************************
 *
 * @param[in] alpha
 *          A scalar.
 *
 * @param[in] bcsc
 *          The PastixHermitian bcsc.
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
z_bcscHemv(pastix_complex64_t  alpha,
           pastix_bcsc_t      *bcsc,
           void               *x,
           pastix_complex64_t  beta,
           void               *y )
{
    pastix_complex64_t *Lvalptr = NULL;
    pastix_complex64_t *yptr    = (pastix_complex64_t*)y;
    pastix_complex64_t *xptr    = (pastix_complex64_t*)x;
    pastix_int_t        bloc, col, i, j, n;

    if(bcsc==NULL || y==NULL || x== NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }
    if(bcsc->mtxtype != PastixHermitian)
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    Lvalptr = (pastix_complex64_t*)bcsc->Lvalues;
    n = bcsc->n;

    /* first, y = beta*y */
    if( beta != (pastix_complex64_t)0.0 )
    {
        for( col = 0; col < n; col++, yptr++ )
        {
            (*yptr) *= beta;
        }
    }
    else if( beta == (pastix_complex64_t)0.0 )
    {
        memset(yptr,0,n*sizeof(pastix_complex64_t));
    }

    yptr = (pastix_complex64_t*)y;
    col = 0;
    for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
    {
        for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
        {
            for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
            {
                yptr[bcsc->rowtab[i]] += alpha * Lvalptr[i] * xptr[col];
                if( bcsc->rowtab[i] != col )
                {
                    yptr[col] += alpha * conj( Lvalptr[i] ) * xptr[bcsc->rowtab[i]];
                }
            }
            col += 1;
        }
    }

    return PASTIX_SUCCESS;
}
