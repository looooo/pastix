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
 * z_spmGeCSCv - compute b=alpha*A*x+beta*b.
 * A is a PastixGeneral csc, 
 * x and b are two vectors of size csc->gN,
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
 * @param[in,out] b
 *          The vector b, can be unallocated at enter.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the b vector has been computed succesfully,
 *      \retval PASTIX_ERR_MATRIX if the matrix is not a PastixGeneral csc.
 *
 *******************************************************************************/
int
z_spmGeCSCv(pastix_complex64_t   alpha,
            pastix_csc_t        *csc  ,  
            pastix_complex64_t   beta ,
            void               **x    ,
            void               **b     )
{
    pastix_complex64_t *temp_Ax = NULL;
    pastix_complex64_t *valptr  = csc->avals;
    pastix_complex64_t *bptr    = *b;
    pastix_complex64_t *xptr    = *x;
    pastix_int_t        col, row, i, baseval;

    if(csc->mtxtype!=PastixGeneral || *x==NULL)
    {
        return PASTIX_ERR_MATRIX;
    }

    baseval = pastix_imin( *(csc->colptr), *(csc->rows) );

    temp_Ax = calloc(csc->gN,sizeof(pastix_complex64_t));
    assert( temp_Ax );
    if(*b==NULL)
    {
        *b = calloc(csc->gN,sizeof(pastix_complex64_t));
        bptr = *b;
    }

    for(col=0;col<csc->gN;col++)
    {
        for(i=csc->colptr[col];i<csc->colptr[col+1];i++)
        {
            row=csc->rows[i-baseval]-baseval;
            temp_Ax[row]+=alpha*valptr[i-baseval]*xptr[col];
        }
    }

    for(i=0;i<csc->gN;i++)
    {
        bptr[i] = temp_Ax[i]+bptr[i]*beta;
    }

    memFree_null(temp_Ax);

    return PASTIX_SUCCESS;

}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * z_spmSyCSCv - compute b=alpha*A*x+beta*b.
 * A is a PastixSymmetric csc, 
 * x and b are two vectors of size csc->gN,
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
 * @param[in,out] b
 *          The vector b, can be unallocated at enter.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the b vector has been computed succesfully,
 *      \retval PASTIX_ERR_MATRIX if the matrix is not a PastixSymmetric csc.
 *
 *******************************************************************************/
int
z_spmSyCSCv(pastix_complex64_t   alpha,
            pastix_csc_t        *csc  ,  
            pastix_complex64_t   beta ,
            void               **x    ,
            void               **b     )
{
    pastix_complex64_t *temp_Ax = NULL;
    pastix_complex64_t *valptr  = csc->avals;
    pastix_complex64_t *bptr    = *b;
    pastix_complex64_t *xptr    = *x;
    pastix_int_t        col, row, i, baseval;

    if(csc->mtxtype!=PastixSymmetric || *x==NULL)
    {
        return PASTIX_ERR_MATRIX;
    }

    baseval = pastix_imin( *(csc->colptr), *(csc->rows) );

    temp_Ax = calloc(csc->gN,sizeof(pastix_complex64_t));
    assert( temp_Ax );
    if(*b==NULL)
    {
        *b = calloc(csc->gN,sizeof(pastix_complex64_t));
        bptr = *b;
    }

    for(col=0;col<csc->gN;col++)
    {
        for(i=csc->colptr[col];i<csc->colptr[col+1];i++)
        {
            row=csc->rows[i-baseval]-baseval;
            temp_Ax[row]+=alpha*valptr[i-baseval]*xptr[col];
            if(col!=row)
            {
                temp_Ax[col]+=alpha*valptr[i-baseval]*xptr[row];
            }
        }
    }

    for(i=0;i<csc->gN;i++)
    {
        bptr[i] = temp_Ax[i]+bptr[i]*beta;
    }

    memFree_null(temp_Ax);

    return PASTIX_SUCCESS;

}

#if !defined(PRECISION_s) && !defined(PRECISION_d)
/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * z_spmGeCSCv - compute b=alpha*A*x+beta*b.
 * A is a PastixHermitian csc, 
 * x and b are two vectors of size csc->gN,
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
 * @param[in,out] b
 *          The vector b, can be unallocated at enter.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the b vector has been computed succesfully,
 *      \retval PASTIX_ERR_MATRIX if the matrix is not a PastixHermitian csc.
 *
 *******************************************************************************/
int
z_spmHeCSCv(pastix_complex64_t   alpha,
            pastix_csc_t        *csc  ,  
            pastix_complex64_t   beta ,
            void               **x    ,
            void               **b     )
{
    pastix_complex64_t *temp_Ax = NULL;
    pastix_complex64_t *valptr  = csc->avals;
    pastix_complex64_t *bptr    = *b;
    pastix_complex64_t *xptr    = *x;
    pastix_int_t        col, row, i, baseval;

    if(csc->mtxtype!=PastixHermitian || *x==NULL)
    {
        return PASTIX_ERR_MATRIX;
    }

    baseval = pastix_imin( *(csc->colptr), *(csc->rows) );

    temp_Ax = calloc(csc->gN,sizeof(pastix_complex64_t));
    assert( temp_Ax );
    if(*b==NULL)
    {
        *b = calloc(csc->gN,sizeof(pastix_complex64_t));
        bptr = *b;
    }

    for(col=0;col<csc->gN;col++)
    {
        for(i=csc->colptr[col];i<csc->colptr[col+1];i++)
        {
            row=csc->rows[i-baseval]-baseval;
            temp_Ax[row]+=alpha*valptr[i-baseval]*xptr[col];
            if(col!=row)
                temp_Ax[col]+=alpha*(creal(valptr[i-baseval]) - cimag(valptr[i-baseval])*I)*xptr[row];
        }
    }

    for(i=0;i<csc->gN;i++)
    {
        bptr[i] = temp_Ax[i]+bptr[i]*beta;
    }

    memFree_null(temp_Ax);

    return PASTIX_SUCCESS;
  
}
#endif