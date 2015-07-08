/**
 *
 * @file genrhs.c
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

#include "z_spm.h"
#include "c_spm.h"
#include "d_spm.h"
#include "s_spm.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * genRHS - generate a RHS such as the solution of Ax=rhs
 * is a vector of size csc->gN with all its values equal to 1+I
 * in complex cases, 1 in real cases.
 *
 *******************************************************************************
 *
 * @param[in] csc
 *          The PastixGeneral csc.
 *
 * @param[in,out] rhs
 *          The generated right hand side member, reallocated if allocated on
 *          entry.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the rhs vector has been computed succesfully,
 *      \retval PASTIX_ERR_BADPARAMETER if the csc matrix is not correct.
 *
 *******************************************************************************/
int
genRHS(pastix_csc_t  *csc,
       void         **rhs )
{
    if(csc->flttype==PastixFloat)
    {
        return s_spm_genRHS(csc, rhs);
    }
    else if(csc->flttype==PastixDouble)
    {
        return d_spm_genRHS(csc, rhs);
    }
    else if(csc->flttype==PastixComplex32)
    {
        return c_spm_genRHS(csc, rhs);
    }
    else if(csc->flttype==PastixComplex64)
    {
        return z_spm_genRHS(csc, rhs);
    }
    else
    {
        return PASTIX_ERR_BADPARAMETER;
    }
}
