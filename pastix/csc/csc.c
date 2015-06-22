/**
 *
 * @file csc.c
 *
 *  PaStiX csc routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "csc.h"

#include "z_spm.h"
#include "c_spm.h"
#include "d_spm.h"
#include "s_spm.h"
#include "p_spm.h"

static int (*conversionTable[3][3][6])(pastix_csc_t*) = {
    /* From CSC */
    {{ NULL, NULL, NULL, NULL, NULL, NULL },
     { p_spmConvertCSC2CSR,
       NULL,
       s_spmConvertCSC2CSR,
       d_spmConvertCSC2CSR,
       c_spmConvertCSC2CSR,
       z_spmConvertCSC2CSR },
     { p_spmConvertCSC2IJV,
       NULL,
       s_spmConvertCSC2IJV,
       d_spmConvertCSC2IJV,
       c_spmConvertCSC2IJV,
       z_spmConvertCSC2IJV }},
    /* From CSR */
    {{ p_spmConvertCSR2CSC,
       NULL,
       s_spmConvertCSR2CSC,
       d_spmConvertCSR2CSC,
       c_spmConvertCSR2CSC,
       z_spmConvertCSR2CSC },
     { NULL, NULL, NULL, NULL, NULL, NULL },
     { p_spmConvertCSR2IJV,
       NULL,
       s_spmConvertCSR2IJV,
       d_spmConvertCSR2IJV,
       c_spmConvertCSR2IJV,
       z_spmConvertCSR2IJV }},
    /* From IJV */
    {{ p_spmConvertIJV2CSC,
       NULL,
       s_spmConvertIJV2CSC,
       d_spmConvertIJV2CSC,
       c_spmConvertIJV2CSC,
       z_spmConvertIJV2CSC },
     { p_spmConvertIJV2CSR,
       NULL,
       s_spmConvertIJV2CSR,
       d_spmConvertIJV2CSR,
       c_spmConvertIJV2CSR,
       z_spmConvertIJV2CSR },
     { NULL, NULL, NULL, NULL, NULL, NULL }}
};



/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * spmConvert - Convert the storage format of the given sparse matrix from any
 * of the following format: PastixCSC, PastixCSR, or PastixIJV to one of these.
 *
 *******************************************************************************
 *
 * @param[in] ofmttype
 *          The output format of the sparse matrix. It might be PastixCSC,
 *          PastixCSR, or PastixIJV.
 *
 * @param[in,out] spm
 *          The sparse matrix structure to convert.
 *
 ********************************************************************************
 *
 * @return
 *        \retval PASTIX_SUCCESS if the conversion happened successfuly
 *        \retval PASTIX_ERR_BADPARAMETER if one the parameter is incorrect.
 *
 *******************************************************************************/
int
spmConvert( int ofmttype, pastix_csc_t *ospm )
{
    if ( conversionTable[ospm->fmttype][ofmttype][ospm->flttype] ) {
        return conversionTable[ospm->fmttype][ofmttype][ospm->flttype]( ospm );
    }
    else {
        return PASTIX_SUCCESS;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * spmFindBase - Search the base used in the spm structure given as parameter.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 ********************************************************************************
 *
 * @return  The baseval used in the given sparse matrix structure.
 *
 *******************************************************************************/
pastix_int_t
spmFindBase( const pastix_csc_t *spm )
{

    pastix_int_t i, *tmp, baseval;

    /*
     * Check the baseval, we consider that arrays are sorted by columns or rows
     */
    baseval = pastix_imin( *(spm->colptr), *(spm->rowptr) );
    /*
     * if not:
     */
    if ( ( baseval != 0 ) &&
         ( baseval != 1 ) )
    {
        baseval = spm->n;
        tmp = spm->colptr;
        for(i=0; i<spm->nnz; i++, tmp++){
            baseval = pastix_imin( *tmp, baseval );
        }
    }

    return baseval;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * spmNorm - Return the ntype norm of the sparse matrix spm.
 *
 *     spmNorm = ( max(abs(spm(i,j))), NORM = PastixMaxNorm
 *               (
 *               ( norm1(spm),         NORM = PastixOneNorm
 *               (
 *               ( normI(spm),         NORM = PastixInfNorm
 *               (
 *               ( normF(spm),         NORM = PastixFrobeniusNorm
 *
 *  where norm1 denotes the one norm of a matrix (maximum column sum),
 *  normI denotes the infinity norm of a matrix (maximum row sum) and
 *  normF denotes the Frobenius norm of a matrix (square root of sum
 *  of squares). Note that max(abs(spm(i,j))) is not a consistent matrix
 *  norm.
 *
 *******************************************************************************
 *
 * @param[in] ntype
 *          = PastixMaxNorm: Max norm
 *          = PastixOneNorm: One norm
 *          = PastixInfNorm: Infinity norm
 *          = PastixFrobeniusNorm: Frobenius norm
 *
 * @param[in] spm
 *          The sparse matrix structure.
 *
 ********************************************************************************
 *
 * @return
 *          \retval the norm described above. Note that for simplicity, even if
 *          the norm of single real or single complex matrix is computed with
 *          single precision, the returned norm is stored in double precision
 *          number.
 *          \retval -1., if the floating point of the sparse matrix is
 *          undefined.
 *
 *******************************************************************************/
double
spmNorm( int ntype,
         const pastix_csc_t *csc )
{
    double tmp;

    switch (csc->flttype) {
    case PastixFloat:
        tmp = (double)s_spmNorm( ntype, csc );
        return tmp;

    case PastixDouble:
        return d_spmNorm( ntype, csc );

    case PastixComplex32:
        tmp = (double)c_spmNorm( ntype, csc );
        return tmp;

    case PastixComplex64:
        return z_spmNorm( ntype, csc );

    default:
        return -1.;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * spmExit - Free the spm structure given as parameter
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The sparse matrix to free.
 *
 *******************************************************************************/
void
spmExit( pastix_csc_t *spm )
{
    if(spm->colptr != NULL)
        memFree_null(spm->colptr);
    if(spm->rowptr != NULL)
        memFree_null(spm->rowptr);
    if(spm->loc2glob != NULL)
        memFree_null(spm->loc2glob);
    if(spm->values != NULL)
        memFree_null(spm->values);
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * spmBase - Rebase the spm to the given value.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The sparse matrix to rebase.
 *
 * @param[in] baseval
 *          The base value to use in the graph (0 or 1).
 *
 *******************************************************************************/
void spmBase( pastix_csc_t *spm,
              int           baseval )
{
    pastix_int_t baseadj;
    pastix_int_t i, n, nnz;

    /* Parameter checks */
    if ( spm == NULL ) {
        errorPrint("spmBase: spm pointer is NULL");
        return;
    }
    if ( (spm->colptr == NULL) ||
         (spm->rowptr == NULL) )
    {
        errorPrint("spmBase: spm pointer is not correctly initialized");
        return;
    }
    if ( (baseval != 0) &&
         (baseval != 1) )
    {
        errorPrint("spmBase: baseval is incorrect, must be 0 or 1");
        return;
    }

    baseadj = baseval - spmFindBase( spm );
    if (baseadj == 0)
	return;

    n   = spm->n;
    nnz = spm->colptr[n] - spm->colptr[0];

    for (i = 0; i <= n; i++) {
        spm->colptr[i] += baseadj;
    }
    for (i = 0; i < nnz; i++) {
        spm->rowptr[i] += baseadj;
    }

    if (spm->loc2glob != NULL) {
        for (i = 0; i < n; i++) {
            spm->loc2glob[i] += baseadj;
        }
    }
    return;
}


/**
 * TODO: Maybe we should move down the cast of the parameters to the lowest
 * functions, and simplify this one to have identical calls to all subfunction
 */
int
spmMatVec(      int           trans,
          const void         *alpha,
          const pastix_csc_t *csc,
          const void         *x,
          const void         *beta,
                void         *y )
{
    switch (csc->mtxtype) {
    case PastixHermitian:
        switch (csc->flttype) {
        case PastixFloat:
            return s_spmSyCSCv( *((const float*)alpha), csc, (const float*)x, *((const float*)beta), (float*)y );
        case PastixComplex32:
            return c_spmHeCSCv( *((const pastix_complex32_t*)alpha), csc, (const pastix_complex32_t*)x, *((const pastix_complex32_t*)beta), (pastix_complex32_t*)y );
        case PastixComplex64:
            return z_spmHeCSCv( *((const pastix_complex64_t*)alpha), csc, (const pastix_complex64_t*)x, *((const pastix_complex64_t*)beta), (pastix_complex64_t*)y );
        case PastixDouble:
        default:
            return d_spmSyCSCv( *((const double*)alpha), csc, (const double*)x, *((const double*)beta), (double*)y );
        }
    case PastixSymmetric:
        switch (csc->flttype) {
        case PastixFloat:
            return s_spmSyCSCv( *((const float*)alpha), csc, (const float*)x, *((const float*)beta), (float*)y );
        case PastixComplex32:
            return c_spmSyCSCv( *((const pastix_complex32_t*)alpha), csc, (const pastix_complex32_t*)x, *((const pastix_complex32_t*)beta), (pastix_complex32_t*)y );
        case PastixComplex64:
            return z_spmSyCSCv( *((const pastix_complex64_t*)alpha), csc, (const pastix_complex64_t*)x, *((const pastix_complex64_t*)beta), (pastix_complex64_t*)y );
        case PastixDouble:
        default:
            return d_spmSyCSCv( *((const double*)alpha), csc, (const double*)x, *((const double*)beta), (double*)y );
        }
    case PastixGeneral:
    default:
        switch (csc->flttype) {
        case PastixFloat:
            return s_spmGeCSCv( trans, *((const float*)alpha), csc, (const float*)x, *((const float*)beta), (float*)y );
        case PastixComplex32:
            return c_spmGeCSCv( trans, *((const pastix_complex32_t*)alpha), csc, (const pastix_complex32_t*)x, *((const pastix_complex32_t*)beta), (pastix_complex32_t*)y );
        case PastixComplex64:
            return z_spmGeCSCv( trans, *((const pastix_complex64_t*)alpha), csc, (const pastix_complex64_t*)x, *((const pastix_complex64_t*)beta), (pastix_complex64_t*)y );
        case PastixDouble:
        default:
            return d_spmGeCSCv( trans, *((const double*)alpha), csc, (const double*)x, *((const double*)beta), (double*)y );
        }
    }
}

