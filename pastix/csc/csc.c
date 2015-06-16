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
