/**
 * @file spm_gen_fake_values.c
 *
 * SParse Matrix generic laplacian value generator routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 **/
#include "common.h"
#include "spm.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 *
 * @brief Compute the degree of each vertex.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          At start, an allocated spm structure.
 *          Contains the size of the laplacian in spm->n.
 *          At exit, contains the matrix in csc format.
 *
 * @param[inout] degrees
 *          Array of size spm->n allocated on entry. On exit, contains the
 *          degree of each vertex in the spm matrix.
 *
 *******************************************************************************/
static inline void
spm_compute_degrees( const pastix_spm_t *spm,
                     pastix_int_t *degrees )
{
    pastix_int_t i, j, k;
    pastix_int_t *colptr = spm->colptr;
    pastix_int_t *rowptr = spm->rowptr;
    pastix_int_t baseval;

    baseval = spmFindBase( spm );
    memset( degrees, 0, spm->n * sizeof(pastix_int_t) );

    switch(spm->fmttype)
    {
    case PastixCSR:
        /* Swap pointers to call CSC */
        colptr = spm->rowptr;
        rowptr = spm->colptr;

    case PastixCSC:
        for(j=0; j<spm->n; j++, colptr++) {
            for(k=colptr[0]; k<colptr[1]; k++, rowptr++) {
                i = *rowptr - baseval;

                if ( i != j ) {
                    degrees[j] += 1;
                    if ( spm->mtxtype != PastixGeneral ) {
                        degrees[i] += 1;
                    }
                }
            }
        }
        break;
    case PastixIJV:
        for(k=0; k<spm->nnz; k++, rowptr++, colptr++)
        {
            i = *rowptr - baseval;
            j = *colptr - baseval;

            if ( i != j ) {
                degrees[j] += 1;

                if ( spm->mtxtype != PastixGeneral ) {
                    degrees[i] += 1;
                }
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 *
 * @brief Generate the fake values array such that \[ M =  \alpha * D - \beta * A \]
 *
 * D is the degree matrix, and A the adjacency matrix.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The spm structure for which the values array must be generated.
 *
 * @param[in] degrees
 *          Array of size spm->n that containes the degree of each vertex in the
 *          spm structure.
 *
 *******************************************************************************/
static inline void
spm_generate_fake_values( pastix_spm_t *spm,
                          const pastix_int_t *degrees,
                          double alpha, double beta )
{
    double *values;
    pastix_int_t i, j, k;
    pastix_int_t *colptr = spm->colptr;
    pastix_int_t *rowptr = spm->rowptr;
    pastix_int_t baseval;

    baseval = spmFindBase( spm );

    spm->values = malloc( spm->nnzexp * sizeof(double) );
    values = spm->values;

    switch(spm->fmttype)
    {
    case PastixCSR:
        /* Swap pointers to call CSC */
        colptr = spm->rowptr;
        rowptr = spm->colptr;

    case PastixCSC:
        for(j=0; j<spm->n; j++, colptr++) {
            for(k=colptr[0]; k<colptr[1]; k++, rowptr++, values++) {
                i = *rowptr - baseval;

                if ( i == j ) {
                    *values = alpha * degrees[j];
                }
                else {
                    *values = - beta;
                }
            }
        }
        break;
    case PastixIJV:
        for(k=0; k<spm->nnz; k++, rowptr++, colptr++, values++)
        {
            i = *rowptr - baseval;
            j = *colptr - baseval;

            if ( i == j ) {
                *values = alpha * degrees[j];
            }
            else {
                *values = - beta;
            }
        }
    }

    spm->flttype = PastixDouble;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Generate the fake values array such that \[ M =  \alpha * D - \beta * A \]
 *
 * D is the degree matrix, and A the adjacency matrix. The resulting matrix uses
 * real double.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          The spm structure for which the values array must be generated.
 *
 *******************************************************************************/
void
spmGenFakeValues( pastix_spm_t *spm )
{
    pastix_int_t *degrees;
    double alpha = 10.;
    double beta = 1.;

    assert( spm->flttype == PastixPattern );
    assert( spm->values == NULL );
    assert( spm->dof == 1 );

    /*
     * Read environment values for alpha/beta
     */
    {
        char *str = pastix_getenv( "PASTIX_FAKE_ALPHA" );
        double value;

        if ( str != NULL ) {
            value = strtod( str, NULL );
            if ( (value != HUGE_VAL) && (value != 0.) &&
                 !isnan(value) && !isinf(value) )
            {
                alpha = value;
            }
            pastix_cleanenv( str );
        }

        str = pastix_getenv( "PASTIX_FAKE_BETA" );
        if ( str != NULL ) {
            value = strtod( str, NULL );
            if ( (value != HUGE_VAL) && (value != 0.) &&
                 !isnan(value) && !isinf(value) )
            {
                beta = value;
            }
            pastix_cleanenv( str );
        }
    }

    degrees = malloc( spm->n * sizeof(pastix_int_t));
    spm_compute_degrees( spm, degrees );
    spm_generate_fake_values( spm, degrees, alpha, beta );

    return;
}
