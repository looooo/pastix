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
 * z_spmMergeDuplicate - This routine merge the multiple entries in a sparse
 * matrix by suming their values together. The sparse matrix needs to be sorted
 * first (see z_spmSort()).
 *
 * WARNING: This function should NOT be called if dof is greater than 1.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          On entry, the pointer to the sparse matrix structure.
 *          On exit, the reduced sparse matrix of multiple entries were present
 *          in it. The multiple values for a same vertex are sum up together.
 *
 ********************************************************************************
 *
 * @return
 *          \retval The number of vertices that were merged
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
z_spmSymmetrize( pastix_csc_t *csc )
{
    pastix_int_t *colptr, *coltmp;
    pastix_int_t *rowptr, *rowtmp;
    pastix_int_t *newol, *toaddtab, toaddcnt, taddsze;
    pastix_int_t  n       = csc->n;
    pastix_int_t  dof2    = csc->dof * csc->dof;
    pastix_int_t  i, j, k, r, size;
    pastix_int_t  baseval;

    if ( (csc->fmttype == PastixCSC) || (csc->fmttype == PastixCSR) ) {
        if (csc->fmttype == PastixCSC) {
            colptr = csc->colptr;
            coltmp = csc->colptr;
            rowptr = csc->rowptr;
            rowtmp = csc->rowptr;
        }
        else {
            colptr = csc->rowptr;
            coltmp = csc->rowptr;
            rowptr = csc->colptr;
            rowtmp = csc->colptr;
        }

        baseval  = colptr[0];
        toaddcnt = 0;
        toaddsze = 0;
        for (i=0; i<n; i++, coltmp++)
        {
            size = coltmp[1] - coltmp[0];
            for (r=0; r<size; r++, rowtmp++ )
            {
                j = rowtmp[0]-baseval;
                if ( i != j ) {
                    /* Look for the element (j, i) */
                    pastix_int_t frow = colptr[ j ];
                    pastix_int_t lrow = colptr[ j+1 ];
                    int found = 0;

                    for (k = frow; (k < lrow); k++)
                    {
                        if (i == (rowptr[k]-baseval))
                        {
                            found = 1;
                            break;
                        }
                        else if ( i < (rowptr[k]-baseval))
                        {
                            /* The csc is sorted so we will not find it later */
                            break;
                        }
                    }

                    if ( !found ) {
                        if ( newcol == NULL ) {
                            newcol = malloc( (csc->n+1) * sizeof(pastix_int_t) );
                            for (k=0; k<csc->n; k++) {
                                newcol[k] = colptr[k+1] - colptr[k];
                            }

                            /* Let's start with a buffer that will contain 5% of extra elements */
                            toaddsze = pastix_imax( 0.05 * (double)csc->nnz, 1 );
                            MALLOC_INTERN(toaddtab, 2*toaddsze, pastix_int_t);
                        }

                        if (toaddcnt >= toaddsze) {
                            toaddsze *= 2;
                            toaddtab = (pastix_int_t*)memRealloc(toadd, 2*toaddsze*sizeof(pastix_int_t));
                        }

                        toaddtab[ toaddcnt * 2     ] = j;
                        toaddtab[ toaddcnt * 2 + 1 ] = i;
                        toaddcnt++;
                    }
                }
            }
        }
    }

    return merge;
}

