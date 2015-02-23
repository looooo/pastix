/**
 * @file z_spm_laplacian.c
 *
 * $COPYRIGHTS$
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2015-01-01
 * @precisions normal z -> c d s p
 *
 **/

#include "common.h"
#include "csc.h"
#include "laplacian.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * z_spmLaplacian1D - Generate a 1D laplacian matrix and an associated right
 * hand side member.
 *
 * Example:
 * >  2 -1  0  0
 * > -1  2 -1  0
 * >  0 -1  2 -1
 * >  0  0 -1  2
 *
 *******************************************************************************
 *
 * @param[in,out] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
 *
 * @param[out] rhs
 *          At exit, contains the right hand side member.
 *
 * @param[in] dim1
 *          contains the dimension of the 1D laplacian.
 *
 *******************************************************************************/
void
z_spmLaplacian1D( pastix_csc_t  *csc,
                  void         **rhs,
                  pastix_int_t   dim1 )
{
    pastix_complex64_t *rhsptr;
    pastix_complex64_t *valptr;
    pastix_int_t *colptr, *rowptr;
    pastix_int_t i, j;
    pastix_int_t nnz = 2*(csc->gN) - 1;
    (void)rhsptr;

    csc->mtxtype  = PastixSymmetric;
    csc->fmttype  = PastixCSC;
    csc->gnnz     = nnz;
    csc->nnz      = nnz;
    csc->dof      = 1;
    csc->colptr   = NULL;
    csc->rows     = NULL;
    csc->loc2glob = NULL;
    csc->avals    = NULL;

    assert( csc->gN == dim1 );

    /* Allocating */
    csc->colptr = malloc((csc->n+1)*sizeof(pastix_int_t));
    csc->rows   = malloc(nnz       *sizeof(pastix_int_t));
    assert( csc->colptr );
    assert( csc->rows   );

#if !defined(PRECISION_p)
    csc->avals  = malloc(nnz       *sizeof(pastix_complex64_t));
    assert( csc->avals  );
#endif

    /* Building ia, ja and avals*/
    colptr = csc->colptr;
    rowptr = csc->rows;
    valptr = (pastix_complex64_t*)(csc->avals);

    j = 0;
    *colptr = 1; colptr++; /* baseval */
    for (i=0; i<csc->gN; i++, colptr++)
    {
        *rowptr = i+1;
#if !defined(PRECISION_p)
        *valptr = (pastix_complex64_t)2.;
#endif

        j++; valptr++; rowptr++;

        if (i < csc->gN-1) {
            *rowptr = i+2;

#if defined(PRECISION_z) || defined(PRECISION_c)
            *valptr = -1. + 2. * _Complex_I;
#elif defined(PRECISION_d) || defined(PRECISION_s)
            *valptr = -1.;
#endif
            j++; rowptr++; valptr++;
        }

        *colptr = j+1;
    }

    assert( (csc->colptr[ csc->gN ] - csc->colptr[0]) == nnz );

#if !defined(PRECISION_p)
    /* Initialize RHS */
    *rhs = calloc(csc->n, sizeof(pastix_complex64_t));
    assert( *rhs );
    rhsptr = (pastix_complex64_t*)(*rhs);
    rhsptr[0]         = (pastix_complex64_t)1.;
    rhsptr[csc->gN-1] = (pastix_complex64_t)1.;
#endif
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * z_spmLaplacian2D - Generate a 2D laplacian matrix and an associated right
 * hand side member.
 *
 * Example:
 * >  4 -1  0 -1  0  0
 * > -1  4 -1  0 -1  0
 * >  0 -1  4  0  0 -1
 * > -1  0  0  4 -1  0
 * >  0 -1  0 -1  4 -1
 * >  0  0 -1  0 -1  4
 *
 *******************************************************************************
 *
 * @param[in,out] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
 *
 * @param[out] rhs
 *          At exit, contains the right hand side member.
 *
 * @param[in] dim1
 *          contains the first dimension of the 2D grid of the laplacian.
 *
 * @param[in] dim2
 *          contains the second dimension of the 2D grid of the laplacian.
 *
 *******************************************************************************/
void
z_spmLaplacian2D( pastix_csc_t  *csc,
                  void         **rhs,
                  pastix_int_t   dim1,
                  pastix_int_t   dim2 )
{
    pastix_complex64_t *rhsptr;
    pastix_complex64_t *valptr;
    pastix_int_t i, j, k;
    pastix_int_t nnz = (2*(dim1)-1)*dim2 + (dim2-1)*dim1;
    (void)rhsptr;
    
    csc->mtxtype  = PastixSymmetric;
    csc->fmttype  = PastixCSC;
    csc->gnnz     = nnz;
    csc->nnz      = nnz;
    csc->dof      = 1;
    csc->colptr  = NULL;
    csc->rows    = NULL;
    csc->avals   = NULL;
    *rhs         = NULL;

    assert( csc->gN == dim1*dim2 );

    /* Allocating */
    csc->colptr = malloc((csc->n+1)*sizeof(pastix_int_t));
    csc->rows   = malloc(nnz       *sizeof(pastix_int_t));
    assert( csc->colptr );
    assert( csc->rows   );

#if !defined(PRECISION_p)
    csc->avals  = malloc(nnz       *sizeof(pastix_complex64_t));
    assert( csc->avals  );
#endif

    /* Building ia, ja and avals */
    csc->colptr[0] = 1;
    valptr = csc->avals;
    k = 0;

    for(i=0; i<dim2; i++)
    {
        for(j=1; j<=dim1; j++)
        {
            // column k = i*dim1+j of the matrix
            k += 1;
            if(j!=dim1 && i!=dim2-1)
            {
                csc->colptr[k] = csc->colptr[k-1]+3;
                csc->rows[csc->colptr[k-1]-1]=k;
                csc->rows[csc->colptr[k-1]  ]=k+1;
                csc->rows[csc->colptr[k-1]+1]=k+dim1;
#if !defined(PRECISION_p)
                *valptr = (pastix_complex64_t)4.;
                *(valptr+1) = -1.;
                *(valptr+2) = -1.;
                valptr += 3;
#endif
            }
            else if(j==dim1 && i!=dim2-1)
            {
                csc->colptr[k] = csc->colptr[k-1]+2;
                csc->rows[csc->colptr[k-1]-1]=k;
                csc->rows[csc->colptr[k-1]  ]=k+dim1;
#if !defined(PRECISION_p)
                *valptr = (pastix_complex64_t)4.;
                *(valptr+1) = -1.;
                valptr += 2;
#endif
            }
            else if(j!=dim1 && i==dim2-1)
            {
                csc->colptr[k] = csc->colptr[k-1]+2;
                csc->rows[csc->colptr[k-1]-1]=k;
                csc->rows[csc->colptr[k-1]  ]=k+1;
#if !defined(PRECISION_p)
                *valptr = (pastix_complex64_t)4.;
                *(valptr+1) = -1;
                valptr += 2;
#endif
            }
            else /* if(j==dim1 && i==dim2-1) */
            {
                csc->colptr[k] = csc->colptr[k-1]+1;
                csc->rows[csc->colptr[k-1]-1]=k;
#if !defined(PRECISION_p)
                *valptr = (pastix_complex64_t)4.;
                valptr += 1;
#endif
            }
        }
    }
    
    assert( (csc->colptr[ csc->gN ] - csc->colptr[0]) == nnz );

#if !defined(PRECISION_p)
    /* Initialize RHS */
    *rhs = calloc(csc->n, sizeof(pastix_complex64_t));
    assert( *rhs );
    rhsptr = (pastix_complex64_t*)(*rhs);
    rhsptr[0]         = (pastix_complex64_t)1.;
    rhsptr[csc->gN-1] = (pastix_complex64_t)1.;
#endif
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * z_spmLaplacian3D - Generate a 3D laplacian matrix and an associated right
 * hand side member.
 *
 * Example:
 * >  6 -1 -1  0 -1  0  0  0
 * > -1  6  0 -1  0 -1  0  0
 * > -1  0  6 -1  0  0 -1  0
 * >  0 -1 -1  6  0  0  0 -1
 * > -1  0  0  0  6 -1 -1  0
 * >  0 -1  0  0 -1  6  0 -1
 * >  0  0 -1  0 -1  0  6 -1
 * >  0  0  0 -1  0 -1 -1  6
 *
 *******************************************************************************
 *
 * @param[in,out] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
 *
 * @param[out] rhs
 *          At exit, contains the right hand side member.
 *
 * @param[in] dim1
 *          contains the first dimension of the 3D grid of the laplacian.
 *
 * @param[in] dim2
 *          contains the second dimension of the 3D grid of the laplacian.
 *
 * @param[in] dim3
 *          contains the third dimension of the 3D grid of the laplacian.
 *
 *******************************************************************************/
void
z_spmLaplacian3D( pastix_csc_t  *csc,
                  void         **rhs,
                  pastix_int_t   dim1,
                  pastix_int_t   dim2,
                  pastix_int_t   dim3 )
{

    pastix_complex64_t *rhsptr;
    pastix_complex64_t *valptr;
    pastix_int_t i, j, k, l;
    pastix_int_t nnz = (2*(dim1)-1)*dim2*dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1);
    (void)rhsptr;
    
    csc->mtxtype  = PastixSymmetric;
    csc->fmttype  = PastixCSC;
    csc->gnnz     = nnz;
    csc->nnz      = nnz;
    csc->dof      = 1;
    csc->colptr  = NULL;
    csc->rows    = NULL;
    csc->avals   = NULL;
    *rhs         = NULL;

    assert( csc->gN == dim1*dim2*dim3 );

    /* Allocating */
    csc->colptr = malloc((csc->n+1)*sizeof(pastix_int_t));
    csc->rows   = malloc(nnz       *sizeof(pastix_int_t));
    assert( csc->colptr );
    assert( csc->rows   );

#if !defined(PRECISION_p)
    csc->avals  = malloc(nnz       *sizeof(pastix_complex64_t));
    assert( csc->avals  );
#endif

    /* Building ia, ja and avals */
    csc->colptr[0] = 1;
    valptr = csc->avals;
    l = 0;

    for(i=0; i<dim3; i++)
    {
        for(j=0; j<dim2; j++)
        {
            for(k=1; k<=dim1; k++)
            {
                // column l = i*dim1*dim2+j*dim1+k of the matrix
                l += 1;
                if(k!=dim1 && j!=dim2-1 && i!=dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+4;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+1;
                    csc->rows[csc->colptr[l-1]+1]=l+dim1;
                    csc->rows[csc->colptr[l-1]+2]=l+dim1*dim2;
#if !defined(PRECISION_p)
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    *(valptr+2) = -1.;
                    *(valptr+3) = -1.;
                    valptr += 4;
#endif
                }
                else if(k==dim1 && j!=dim2-1 && i!=dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+3;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+dim1;
                    csc->rows[csc->colptr[l-1]+1]=l+dim1*dim2;
#if !defined(PRECISION_p)
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    *(valptr+2) = -1.;
                    valptr += 3;
#endif
                }
                else if(k!=dim1 && j==dim2-1 && i!=dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+3;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+1;
                    csc->rows[csc->colptr[l-1]+1]=l+dim1*dim2;
#if !defined(PRECISION_p)
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    *(valptr+2) = -1.;
                    valptr += 3;
#endif
                }
                else if(k!=dim1 && j!=dim2-1 && i==dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+3;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+1;
                    csc->rows[csc->colptr[l-1]+1]=l+dim1;
#if !defined(PRECISION_p)
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    *(valptr+2) = -1.;
                    valptr += 3;
#endif
                }
                else if(k==dim1 && j==dim2-1 && i!=dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+2;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+dim1*dim2;
#if !defined(PRECISION_p)
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    valptr += 2;
#endif
                }
                else if(k!=dim1 && j==dim2-1 && i==dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+2;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+1;
#if !defined(PRECISION_p)
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    valptr += 2;
#endif
                }
                else if(k==dim1 && j!=dim2-1 && i==dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+2;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+dim1;
#if !defined(PRECISION_p)
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    valptr += 2;
#endif
                }
                else if(k==dim1 && j==dim2-1 && i==dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+1;
                    csc->rows[csc->colptr[l-1]-1]=l;
#if !defined(PRECISION_p)
                    *valptr = 6.;
                    valptr += 1;
#endif
                }
            }
        }
    }
    
    assert( (csc->colptr[ csc->gN ] - csc->colptr[0]) == nnz );

#if !defined(PRECISION_p)
    /* Initialize RHS */
    *rhs = calloc(csc->n, sizeof(pastix_complex64_t));
    assert( *rhs );
    rhsptr = (pastix_complex64_t*)(*rhs);
    rhsptr[0]         = (pastix_complex64_t)1.;
    rhsptr[csc->gN-1] = (pastix_complex64_t)1.;
#endif
}
