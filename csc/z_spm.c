/**
 *
 * @file z_spm.c
 *
 *  PaStiX spm routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> c d s p
 *
 **/
#include "common.h"
#include "csc.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * z_spmSort - This routine sorts the subarray of edges of each vertex in a
 * centralized spm stored in CSC or CSR format. Nothing is performed if IJV
 * format is used.
 *
 * WARNING: This function should NOT be called if dof is greater than 1.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the same sparse matrix with subarrays of edges sorted by
 *          ascending order.
 *
 *******************************************************************************/
void
z_spmSort( pastix_csc_t *csc )
{
    pastix_int_t       *colptr = csc->colptr;
    pastix_int_t       *rowptr = csc->rowptr;
    pastix_complex64_t *values = csc->values;
    void *sortptr[2];
    pastix_int_t n = csc->n;
    pastix_int_t i, size;
    (void)sortptr;

    if (csc->dof > 1){
        fprintf(stderr, "z_spmSort: Number of dof (%d) greater than one not supported\n", (int)csc->dof);
        exit(1);
    }

    /* Sort in place each subset */
    if ( csc->fmttype == PastixCSC ) {
        for (i=0; i<n; i++, colptr++)
        {
            size = colptr[1] - colptr[0];

#if defined(PRECISION_p)
            intSort1asc1( rowptr, size );
#else
            sortptr[0] = rowptr;
            sortptr[1] = values;
            z_qsortIntFloatAsc( sortptr, size );
#endif
            rowptr += size;
            values += size;
        }
    }
    else if ( csc->fmttype == PastixCSR ) {
        for (i=0; i<n; i++, rowptr++)
        {
            size = rowptr[1] - rowptr[0];

#if defined(PRECISION_p)
            intSort1asc1( colptr, size );
#else
            sortptr[0] = colptr;
            sortptr[1] = values;
            z_qsortIntFloatAsc( sortptr, size );
#endif
            colptr += size;
            values += size;
        }
    }
}

