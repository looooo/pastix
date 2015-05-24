/**
 * @file z_spm_norm.c
 *
 *  PaSTiX SPM management routines.
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "csc.h"

static inline void
frobenius_update( double *scale, double *sumsq, double *value )
{
    double ratio;
    if (*value != 0. ){
        if ( *scale < *value ) {
            ratio = *scale / * value;
            *sumsq = 1. + (*sumsq) * ratio * ratio;
            *scale = *value;
        } else {
            ratio = *value / *scale;
            *sumsq = *sumsq + ratio * ratio;
        }
    }
}

double
z_spmMaxNorm( const pastix_csc_t *csc )
{
    pastix_int_t i;
    pastix_complex64_t *valptr = (pastix_complex64_t*)csc->values;
    double tmp, norm = 0.;

    for(i=0; i <csc->nnz; i++, valptr++) {
        tmp = cabs( *valptr );
        norm = norm > tmp ? norm : tmp;
    }

    // TODO:
    // Thread reduction, and distributed reduction
    return norm;
}

double
z_spmInfNorm( const pastix_csc_t *csc )
{
    pastix_int_t i;
    pastix_complex64_t *valptr = (pastix_complex64_t*)csc->values;
    double tmp, norm = 0.;

    for(i=0; i <csc->nnz; i++, valptr++) {
        tmp = cabs( *valptr );
        norm = norm > tmp ? norm : tmp;
    }

    // TODO: Fix to really compute the Inf norm
    // Thread reduction, and distributed reduction
    return norm;
}

double
z_spmOneNorm( const pastix_csc_t *csc )
{
    pastix_int_t i;
    pastix_complex64_t *valptr = (pastix_complex64_t*)csc->values;
    double tmp, norm = 0.;

    for(i=0; i <csc->nnz; i++, valptr++) {
        tmp = cabs( *valptr );
        norm = norm > tmp ? norm : tmp;
    }

    // TODO: Fix to really compute the One norm
    // Thread reduction, and distributed reduction
    return norm;
}

double
z_spmFrobeniusNorm( const pastix_csc_t *csc )
{
    pastix_int_t i;
    double *valptr = (double*)csc->values;
    double scale = 1.;
    double sumsq = 0.;

    for(i=0; i <csc->nnz; i++, valptr++) {
        frobenius_update( &scale, &sumsq, valptr );

#if defined(PRECISION_z) || defined(PRECISION_c)
        valptr++;
        frobenius_update( &scale, &sumsq, valptr );
#endif
    }

    // TODO:
    // Thread reduction, and distributed reduction

    return scale * sqrt( sumsq );
}

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
double
z_spmNorm( int ntype,
           const pastix_csc_t *csc )
{
    double norm = 0.;

    if(csc == NULL)
    {
        return PASTIX_ERR_BADPARAMETER;
    }

    switch( ntype ) {
    case PastixMaxNorm:
        norm = z_spmMaxNorm( csc );
        break;

    case PastixInfNorm:
        norm = z_spmInfNorm( csc );
        break;

    case PastixOneNorm:
        norm = z_spmOneNorm( csc );
        break;

    case PastixFrobeniusNorm:
        norm = z_spmFrobeniusNorm( csc );
        break;

    default:
        fprintf(stderr, "z_spmNorm: invalid norm type\n");
        return 0.;
    }

    return norm;
}
