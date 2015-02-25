/**
 *
 * @file z_spm_genrhs.c
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
 * @precisions normal z -> c s d
 **/
#include "common.h"
#include "csc.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * z_spm_genRHS - generate a RHS such as the solution of Ax=rhs
 * is a vector of size csc->gN with all its values equal to 1+I
 * in complex cases, 1 in real cases.
 *
 *******************************************************************************
 *
 * @param[in] csc
 *          The PastixGeneral csc.
 * 
 * @param[in,out] rhs
 *          The generated rhight hand side member, 
 *          reallocated if allocated at enter.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the b vector has been computed succesfully,
 *      \retval PASTIX_ERR_MATRIX if the csc matrix is not correct.
 *
 *******************************************************************************/
int
z_spm_genRHS(pastix_csc_t  *csc,
             void         **rhs )
{
    void *x = NULL;

    if(csc->avals==NULL)
        return PASTIX_ERR_MATRIX;

    if(csc->fmttype!=PastixCSC)
        return PASTIX_ERR_MATRIX;

    if(csc->gN<=0)
        return PASTIX_ERR_MATRIX;

    x=malloc(csc->gN*sizeof(pastix_complex64_t));
    if(csc->flttype!=PastixComplex64 && csc->flttype!=PastixComplex32)
    {
        memset(x,1,csc->gN*sizeof(pastix_complex64_t));
    }else{
        memset(x,1+I,csc->gN*sizeof(pastix_complex64_t));
    }
        
    if(*rhs != NULL)
        memFree_null(*rhs);

    if(csc->mtxtype==PastixGeneral)
    {
        if(z_spmGeCSCv((pastix_complex64_t)1.,csc,(pastix_complex64_t)0.,&x,rhs)!=PASTIX_SUCCESS)
            return PASTIX_ERR_MATRIX;
    }
    else if(csc->mtxtype==PastixSymmetric)
    {
        if(z_spmSyCSCv((pastix_complex64_t)1.,csc,(pastix_complex64_t)0.,&x,rhs)!=PASTIX_SUCCESS)
            return PASTIX_ERR_MATRIX;
    }
    else if(csc->mtxtype==PastixHermitian)
    {
        if(z_spmHeCSCv((pastix_complex64_t)1.,csc,(pastix_complex64_t)0.,&x,rhs)!=PASTIX_SUCCESS)
            return PASTIX_ERR_MATRIX;
    }

    memFree_null(x);

    return PASTIX_SUCCESS;
}