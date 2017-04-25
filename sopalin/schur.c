/**
 *
 * @file pastix.c
 *
 * PaStiX schur interface functions
 *
 * @copyright 2011-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathias Hastaran
 * @date 2011-11-11
 *
 * @addtogroup pastix_schur
 * @{
 *
 **/
#include "common.h"
#include "spm.h"

/**
 *******************************************************************************
 *
 * @brief Set the list of unknowns that belongs to the schur complement.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix data structure of the solver to store the list of Schur
 *          unknowns.
 *
 * @param[in] n
 *          The number of unknowns in the Schur complement.
 *
 * @param[in] list
 *          Array of integer of size n.
 *          The list of unknowns belonging to the Schur complement with the same
 *          baseval as the associated spm.
 *
 *******************************************************************************/
void
pastix_setSchurUnknownList( pastix_data_t      *pastix_data,
                            pastix_int_t        n,
                            const pastix_int_t *list)
{
    pastix_data->schur_n    = n;
    pastix_data->schur_list = (pastix_int_t*)malloc(n * sizeof(pastix_int_t));
    memcpy( pastix_data->schur_list, list, n * sizeof(pastix_int_t) );
}

/**
 *******************************************************************************
 *
 * @brief Return the Schur complement.
 *
 * The Schur complement is returned in the column major layout used by the
 * classic linear algebra libraries such as Blas or Lapack.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix data structure of the problem solved.
 *
 * @param[inout] S
 *          Array of size spm->n -by- lds of arithmetic spm->flttype, where spm
 *          is the spm of the original problem.
 *          On exit, the array contains the Schur complement of the factorized
 *          matrix.
 *
 * @param[in] lds
 *          The leading dimension of the S array.
 *
 *******************************************************************************/
void
pastix_getSchur( const pastix_data_t *pastix_data,
                 void                *S,
                 pastix_int_t         lds )
{
    return;
}

