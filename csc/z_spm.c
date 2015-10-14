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

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 * @ingroup pastix_internal
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
pastix_int_t
z_spmMergeDuplicate( pastix_csc_t *csc )
{
    pastix_int_t       *colptr = csc->colptr;
    pastix_int_t       *oldrow = csc->rowptr;
    pastix_int_t       *newrow = csc->rowptr;
    pastix_complex64_t *newval = csc->values;
    pastix_complex64_t *oldval = csc->values;
    pastix_int_t n       = csc->n;
    pastix_int_t baseval = csc->colptr[0];
    pastix_int_t dof2    = csc->dof * csc->dof;
    pastix_int_t idx, i, j, d, size;
    pastix_int_t merge = 0;
    (void)d;

    if ( csc->fmttype == PastixCSC ) {
        idx = 0;
        for (i=0; i<n; i++, colptr++)
        {
            size = colptr[1] - colptr[0];
            for (j = 0; j < size;
                 j++, oldrow++, oldval+=dof2, newrow++, newval+=dof2, idx++)
            {
                if ( newrow != oldrow ) {
                    newrow[0] = oldrow[0];
#if !defined(PRECISION_p)
                    memcpy( newval, oldval, dof2 * sizeof(pastix_complex64_t) );
#endif
                }

                while( ((j+1) < size) && (newrow[0] == oldrow[1]) ) {
                    j++; oldrow++; oldval+=dof2;
#if !defined(PRECISION_p)
                    /* Merge the two sets of values */
                    for (d = 0; d < dof2; d++) {
                        newval[d] += oldval[d];
                    }
#endif
                    merge++;
                }
            }
            assert( ((merge == 0) && ( colptr[1] == idx+baseval)) ||
                    ((merge != 0) && ( colptr[1] >  idx+baseval)) );

            colptr[1] = idx + baseval;
        }
        assert( ((merge == 0) && (csc->nnz         == idx)) ||
                ((merge != 0) && (csc->nnz - merge == idx)) );

        if (merge > 0) {
            printf("spmMergeDuplicate: %ld elements were merged\n", (int64_t)merge);

            csc->nnz = csc->nnz - merge;

            newrow = malloc( csc->nnz * sizeof(pastix_int_t) );
            memcpy( newrow, csc->rowptr, csc->nnz * sizeof(pastix_int_t) );
            free(csc->rowptr);
            csc->rowptr = newrow;

#if !defined(PRECISION_p)
            newval = malloc( csc->nnz * dof2 * sizeof(pastix_int_t) );
            memcpy( newval, csc->values, csc->nnz * dof2 * sizeof(pastix_complex64_t) );
            free(csc->values);
            csc->values = newval;
#endif
        }
    }

    return merge;
}

