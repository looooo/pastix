 /**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Théophile terraz
 * @date 2015-01-01
 * @precisions normal z -> c d s
 *
 **/
/*
  File: z_bcsc_norm.c

  Functions computing norms on the BCSC.

*/

#include "common.h"
#include "bcsc.h"
#include <math.h>
#include "frobeniusupdate.h"


/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscMaxNorm - compute the max norm of a bcsc matrix.
 *
 *******************************************************************************
 *
 * @param[in] values
 *          The values array of the matrix.
 *
 *******************************************************************************
 *
 * @return
 *      \retval The norm of the matrix.
 *
 *******************************************************************************/
double
z_bcscMaxNorm( const pastix_bcsc_t *bcsc )
{
    double temp;
    double norm = 0.;
    pastix_complex64_t *valptr = (pastix_complex64_t*)bcsc->Lvalues;
    pastix_int_t i, j, bloc;

    for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
    {
        for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
        {
            for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
            {
                temp = cabs(valptr[i]);
                if(norm < temp)
                {
                    norm = temp;
                }
            }
        }
    }

    return norm;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscInfNorm - compute the infinity norm of a bcsc matrix.
 * The infinity norm is equal to the maximum value of the sum of the
 * Absolute values of the elements of each rows.
 *
 *******************************************************************************
 *
 * @param[in] bcsc
 *          The Pastix bcsc.
 *
 *******************************************************************************
 *
 * @return
 *      \retval The norm of the matrix.
 *
 *******************************************************************************/
double
z_bcscInfNorm( const pastix_bcsc_t *bcsc )
{
    double *summcol;
    double norm = 0.;
    pastix_complex64_t *valptr = (pastix_complex64_t*)bcsc->Lvalues;
    int n,i,j,bloc;

    n = bcsc->n;
    MALLOC_INTERN( summcol, n, double);
    memset( summcol, 0, n * sizeof(double) );
    for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
    {
        for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
        {
            for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
            {
                summcol[bcsc->rowtab[i]] += cabs(valptr[i]);
            }
        }
    }

    if(bcsc->mtxtype != PastixGeneral)
    {
// TODO check if this is good
        valptr = (pastix_complex64_t*)bcsc->Uvalues;
        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
            {
                for( i = bcsc->cscftab[bloc].coltab[j]+1; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                {
                    summcol[bcsc->rowtab[i]] += cabs(valptr[i]);
                }
            }
        }
    }

    for( i=0; i<n; i++)
    {
        if(norm < summcol[i])
        {
            norm = summcol[i];
        }
    }
    memFree_null( summcol );

    return norm;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscOneNorm - compute the norm 1 of a bcsc matrix.
 * Norm 1 is equal to the maximum value of the sum of the
 * Absolute values of the elements of each columns.
 *
 *******************************************************************************
 *
 * @param[in] bcsc
 *          The Pastix bcsc.
 *
 *******************************************************************************
 *
 * @return
 *      \retval The norm of the matrix.
 *
 *******************************************************************************/
double
z_bcscOneNorm( const pastix_bcsc_t *bcsc )
{
    double *summrow;
    double norm = 0.;
    pastix_complex64_t *valptr = (pastix_complex64_t*)bcsc->Lvalues;
    int n,i,j,bloc,col;

    n = bcsc->n;
    MALLOC_INTERN( summrow, n, double);
    memset( summrow, 0, n * sizeof(double) );
    col = 0;
    for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
    {
        for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
        {
            for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
            {
                summrow[col] += cabs(valptr[i]);
            }
            col += 1;
        }
    }

    if(bcsc->mtxtype != PastixGeneral)
    {
        valptr = (pastix_complex64_t*)bcsc->Uvalues;
        col = 0;
        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
            {
// TODO check if this is good
                for( i = bcsc->cscftab[bloc].coltab[j]+1; i < bcsc->cscftab[bloc].coltab[j+1]; i++ )
                {
                    summrow[col] += cabs(valptr[i]);
                }
                col += 1;
            }
        }
    }

    for( i=0; i<n; i++)
    {
        if(norm < summrow[i])
        {
            norm = summrow[i];
        }
    }
    memFree_null( summrow );

    return norm;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscFrobeniusNorm - compute the frobenius norm of a bcsc matrix.
 *
 *******************************************************************************
 *
 * @param[in] bcsc
 *          The Pastix bcsc.
 *
 *******************************************************************************
 *
 * @return
 *          The norm of the matrix
 *
 *******************************************************************************/
double
z_bcscFrobeniusNorm( const pastix_bcsc_t *bcsc)
{
    double scale = 0.;
    double sum = 1.;
    double norm;
    double *valptr = (double*)bcsc->Lvalues;
    pastix_int_t i, j, bloc;

    for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
    {
        for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
        {
            for( i = bcsc->cscftab[bloc].coltab[j]; i < bcsc->cscftab[bloc].coltab[j+1]; i++, valptr++ )
            {
                frobenius_update( 1, &scale, &sum, valptr);
#if defined(PRECISION_z) || defined(PRECISION_c)
                valptr++;
                frobenius_update( 1, &scale, &sum, valptr);
#endif
            }
        }
    }

    if(bcsc->mtxtype != PastixGeneral)
    {
// TODO check if this is good
        valptr = (double*)bcsc->Uvalues;

        for( bloc=0; bloc < bcsc->cscfnbr; bloc++ )
        {
            for( j=0; j < bcsc->cscftab[bloc].colnbr; j++ )
            {
                for( i = bcsc->cscftab[bloc].coltab[j]+1; i < bcsc->cscftab[bloc].coltab[j+1]; i++, valptr++ )
                {
                frobenius_update( 1, &scale, &sum, valptr);
#if defined(PRECISION_z) || defined(PRECISION_c)
                valptr++;
                frobenius_update( 1, &scale, &sum, valptr);
#endif
                }
            }
        }
    }

    norm = scale*sqrt(sum);

    return norm;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_bcsc
 *
 * z_bcscNorm - Compute the norm of an bcsc matrix
 *
 *******************************************************************************
 *
 * @param[in] type
 *          = PastixMaxNorm: Max norm
 *          = PastixOneNorm: One norm
 *          = PastixInfNorm: Infinity norm
 *          = PastixFrobeniusNorm: Frobenius norm
 *
 * @param[in] bcsc
 *          The bcsc structure describing the matrix.
 *
 *******************************************************************************
 *
 * @return
 *      \retval The norm of the matrix.
 *
 *******************************************************************************/
double
z_bcscNorm( pastix_normtype_t ntype,
            const pastix_bcsc_t *bcsc )
{
    double norm = 0.;

    if(bcsc == NULL)
    {
        return -1.;
    }

    switch( ntype ) {
    case PastixMaxNorm:
        norm = z_bcscMaxNorm( bcsc );
        break;

    case PastixInfNorm:
        norm = z_bcscInfNorm( bcsc );
        break;

    case PastixOneNorm:
        norm = z_bcscOneNorm( bcsc );
        break;

    case PastixFrobeniusNorm:
        norm = z_bcscFrobeniusNorm( bcsc );
        break;

    default:
        fprintf(stderr, "z_spmNorm: invalid norm type\n");
        return -1.;
    }

    return norm;
}