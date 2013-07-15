/**
 *
 * @file csc_base.c
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

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * cscBase - Rebase the csc to the given value.
 *
 *******************************************************************************
 *
 * @param[in,out] csc
 *          The csc matrix to rebase.
 *
 * @param[in] baseval
 *          The base value to use in the CSC.
 *
 *******************************************************************************/
void cscBase( pastix_csc_t *csc,
              int           baseval )
{
    pastix_int_t baseadj;
    pastix_int_t i, n, nnz;

    baseadj = baseval - csc->colptr[0];
    if (baseadj == 0)
	return;

    n   = csc->n;
    nnz = csc->colptr[n] - csc->colptr[0];

    if (csc->colptr != NULL) {
        for (i = 0; i <= n; i++) {
            csc->colptr[i]   += baseadj;
        }
    }
    if (csc->rows != NULL) {
        for (i = 0; i < nnz; i++)
            csc->rows[i] += baseadj;
    }

    if (csc->loc2glob != NULL) {
        for (i = 0; i < n; i++) {
            csc->loc2glob[i] += baseadj;
        }
    }
    return;
}
