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
#include "common.h"
#include <math.h>
#include "bcsc.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscMatVec - Compute the matrix-vector product
 *         y = alpha * op(A) * x + beta * y,
 * where A is given in the bcsc format, x and y are two vectors of size n, and
 * alpha and beta are two scalars.
 * The op function is specified by the trans parameter and performs the
 * operation as follows:
 *              trans = PastixNoTrans   y := alpha*A       *x + beta*y
 *              trans = PastixTrans     y := alpha*A'      *x + beta*y
 *              trans = PastixConjTrans y := alpha*conj(A')*x + beta*y
 *
 *******************************************************************************
 *
 * @param[in] trans
 *          Specifies whether the matrix A from the bcsc is transposed, not
 *          transposed or conjugate transposed:
 *            = PastixNoTrans:   A is not transposed;
 *            = PastixTrans:     A is transposed;
 *            = PastixConjTrans: A is conjugate transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] bcsc
 *          The bcsc structure describing the matrix A.
 *
 * @param[in] x
 *          The vector x.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
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
    pastix_complex64_t *valptr = NULL;
    pastix_complex64_t *yptr, *xptr;
    pastix_int_t        bloc, i, j, n;

    if(bcsc==NULL || y==NULL || x== NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }
    n = bcsc->n;

    yptr = (pastix_complex64_t*)y;

    /* first, y = beta*y */
    if( beta != (pastix_complex64_t)0.0 )
    {
        for( j=0; j<n; j++, yptr++ )
        {
            (*yptr) *= beta;
        }
    }
    else if( beta == (pastix_complex64_t)0.0 )
    {
        memset(yptr, 0, n * sizeof(pastix_complex64_t));
    }

    /**
     * There are three cases:
     *    We can use the Lvalues pointer directly:
     *          - The matrix is general and we use A^t
     *          - the matrix is symmetric or hermitian
     *    We can use the Uvalues pointer directly
     *          - The matrix is general and we use A
     *    We have to use Lvalues per row (instead of column)
     *          - The matrix A is general and Uvalues is unavailable
     *
     * To this, we have to add the conj call if ConjTrans or Hermitian
     *
     *     Mtxtype   | trans asked | algo applied
     *     ++++++++++++++++++++++++++++++++++++
     +     General   | NoTrans     | U if possible, otherwise indirect L
     +     General   | Trans       | L
     +     General   | ConjTrans   | conj(L)
     +     Symmetric | NoTrans     | L
     +     Symmetric | Trans       | L
     +     Symmetric | ConjTrans   | conj(L)
     +     Hermitian | NoTrans     | conj(L)
     +     Hermitian | Trans       | conj(L)
     +     Hermitian | ConjTrans   | L
     */
    yptr = (pastix_complex64_t*)y;
    xptr = (pastix_complex64_t*)x;

    if (bcsc->mtxtype == PastixGeneral && trans == PastixNoTrans )
    {
        /* U */
        if ( bcsc->Uvalues != NULL ) {
            valptr = (pastix_complex64_t*)bcsc->Uvalues;

            for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
            {
                for( j=0; j < bcsc->cscftab[bloc].colnbr; j++, yptr++ )
                {
                    for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                    {
                        *yptr += alpha * valptr[i] * xptr[ bcsc->rowtab[i] ];
                    }
                }
            }
        }
        /* Indirect L */
        else {
            valptr = (pastix_complex64_t*)bcsc->Lvalues;
            for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
            {
                for( j=0; j < bcsc->cscftab[bloc].colnbr; j++, xptr++ )
                {
                    for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                    {
                        yptr[ bcsc->rowtab[i] ] += alpha * valptr[i] * (*xptr);
                    }
                }
            }
        }
    }
#if defined(PRECISION_z) || defined(PRECISION_c)
    /* Conj(L) */
    else if ( (bcsc->mtxtype == PastixGeneral   && trans == PastixConjTrans ) ||
              (bcsc->mtxtype == PastixSymmetric && trans == PastixConjTrans ) ||
              (bcsc->mtxtype == PastixHermitian && trans != PastixConjTrans ) )
    {
        valptr = (pastix_complex64_t*)bcsc->Lvalues;

        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc].colnbr; j++, yptr++ )
            {
                for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                {
                    *yptr += alpha * conj( valptr[i] ) * xptr[ bcsc->rowtab[i] ];
                }
            }
        }
    }
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */
    /* L */
    else {
        valptr = (pastix_complex64_t*)bcsc->Lvalues;

        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc].colnbr; j++, yptr++ )
            {
                for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                {
                    *yptr += alpha * valptr[i] * xptr[ bcsc->rowtab[i] ];
                }
            }
        }
    }

    return PASTIX_SUCCESS;
}
