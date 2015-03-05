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

static int (*CSCv[3])(char, pastix_complex64_t, pastix_csc_t*, pastix_complex64_t*, pastix_complex64_t, pastix_complex64_t*) =
{
    z_spmGeCSCv,
    z_spmSyCSCv,
#if defined(PRECISION_z) || defined(PRECISION_c)
    z_spmHeCSCv
#else
    NULL
#endif
};

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
 *      \retval PASTIX_ERR_BADPARAMETER otherwise.
 *
 *******************************************************************************/
int
z_spm_genRHS(pastix_csc_t  *csc,
             void         **rhs )
{
    void *x = NULL;
    pastix_complex64_t *xptr;
    char n  = 'n';
    pastix_int_t i;

    if(csc->avals==NULL)
        return PASTIX_ERR_BADPARAMETER;

    if(csc->fmttype!=PastixCSC)
        return PASTIX_ERR_BADPARAMETER;

    if(csc->gN<=0)
        return PASTIX_ERR_BADPARAMETER;

    x=malloc(csc->gN*sizeof(pastix_complex64_t));
    xptr = x;

#if defined(PRECISION_z) || defined(PRECISION_c)
    for(i=0; i < csc->gN; i++, xptr++)
    {
        *xptr = (pastix_complex64_t)(1.+1.*I);
    }
#else
    for(i=0; i < csc->gN; i++, xptr++)
    {
        *xptr = (pastix_complex64_t)1.;
    }
#endif

    if(*rhs == NULL)
    {
        *rhs=malloc(csc->gN*sizeof(pastix_complex64_t));
        memset(*rhs,0,csc->gN*sizeof(pastix_complex64_t));
    }

    if(CSCv[csc->mtxtype-PastixGeneral])
    {
        if( CSCv[csc->mtxtype-PastixGeneral](n, 1., csc, x, 1., *rhs ) != PASTIX_SUCCESS )
        {
            return PASTIX_ERR_BADPARAMETER;
        }
    }

    memFree_null(x);

    return PASTIX_SUCCESS;
}