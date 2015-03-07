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
 * >  1 -1  0  0
 * > -1  2 -1  0
 * >  0 -1  2 -1
 * >  0  0 -1  1
 *
 *******************************************************************************
 *
 * @param[in,out] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
 *
 * @param[in] dim1
 *          contains the dimension of the 1D laplacian.
 *
 *******************************************************************************/
void
z_spmLaplacian1D( pastix_csc_t  *csc,
                  pastix_int_t   dim1 )
{
    pastix_complex64_t *valptr;
    pastix_int_t *colptr, *rowptr;
    pastix_int_t i, j;
    pastix_int_t nnz = 2*(csc->gN) - 1;

    csc->mtxtype  = PastixSymmetric;
    csc->flttype  = PastixComplex64;
    csc->fmttype  = PastixCSC;
    csc->gnnz     = nnz;
    csc->nnz      = nnz;
    csc->dof      = 1;

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
    for (i=0; i<dim1; i++, colptr++)
    {
        *rowptr = i+1;
#if !defined(PRECISION_p)
        if ( (i == 0) || (i == dim1-1) ) {
            *valptr = (pastix_complex64_t)1.;
        }
        else {
            *valptr = (pastix_complex64_t)2.;
        }
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

    assert( (csc->colptr[ dim1 ] - csc->colptr[0]) == nnz );
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
 * >  2 -1  0 -1  0  0
 * > -1  3 -1  0 -1  0
 * >  0 -1  2  0  0 -1
 * > -1  0  0  2 -1  0
 * >  0 -1  0 -1  3 -1
 * >  0  0 -1  0 -1  2
 *
 *******************************************************************************
 *
 * @param[in,out] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
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
                  pastix_int_t   dim1,
                  pastix_int_t   dim2 )
{
    pastix_complex64_t *valptr;
    pastix_int_t *colptr, *rowptr;
    pastix_int_t i, j, k;
    pastix_int_t nnz = (2*(dim1)-1)*dim2 + (dim2-1)*dim1;

    csc->mtxtype  = PastixSymmetric;
    csc->flttype  = PastixComplex64;
    csc->fmttype  = PastixCSC;
    csc->gnnz     = nnz;
    csc->nnz      = nnz;
    csc->dof      = 1;

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

    /* Building ia, ja and avals*/
    colptr = csc->colptr;
    rowptr = csc->rows;
    valptr = (pastix_complex64_t*)(csc->avals);

    /* Building ia, ja and avals */
    *colptr = 1;
    k = 1; /* Column index in the matrix ((i-1) * dim1 + j-1) */
    for(i=1; i<=dim2; i++)
    {
        for(j=1; j<=dim1; j++)
        {
            colptr[1] = colptr[0];

            /* Diagonal value */
#if !defined(PRECISION_p)
            *rowptr = k;
            *valptr = (pastix_complex64_t) 4.;
            if (j == dim1 || j == 1)
                *valptr -= (pastix_complex64_t) 1.;
            if (i == dim2 || i == 1)
                *valptr -= (pastix_complex64_t) 1.;
#endif
            valptr++; rowptr++; colptr[1]++;

            /* Connexion along dimension 1 */
            if (j < dim1) {
#if !defined(PRECISION_p)
                *rowptr = k+1;
                *valptr = (pastix_complex64_t)-1.;
#endif
                valptr++; rowptr++; colptr[1]++;
            }

            /* Connexion along dimension 2 */
            if (i < dim2) {
#if !defined(PRECISION_p)
                *rowptr = k+dim1;
                *valptr = (pastix_complex64_t)-1.;
#endif
                valptr++; rowptr++; colptr[1]++;
            }

            colptr++; k++;
        }
    }

    assert( (csc->colptr[ csc->gN ] - csc->colptr[0]) == nnz );
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
 * >  3 -1 -1  0 -1  0  0  0
 * > -1  3  0 -1  0 -1  0  0
 * > -1  0  3 -1  0  0 -1  0
 * >  0 -1 -1  3  0  0  0 -1
 * > -1  0  0  0  3 -1 -1  0
 * >  0 -1  0  0 -1  3  0 -1
 * >  0  0 -1  0 -1  0  3 -1
 * >  0  0  0 -1  0 -1 -1  3
 *
 *******************************************************************************
 *
 * @param[in,out] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
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
                  pastix_int_t   dim1,
                  pastix_int_t   dim2,
                  pastix_int_t   dim3 )
{

    pastix_complex64_t *valptr;
    pastix_int_t *colptr, *rowptr;
    pastix_int_t i, j, k, l;
    pastix_int_t nnz = (2*(dim1)-1)*dim2*dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1);

    csc->mtxtype  = PastixSymmetric;
    csc->flttype  = PastixComplex64;
    csc->fmttype  = PastixCSC;
    csc->gnnz     = nnz;
    csc->nnz      = nnz;
    csc->dof      = 1;

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

    /* Building ia, ja and avals*/
    colptr = csc->colptr;
    rowptr = csc->rows;
    valptr = (pastix_complex64_t*)(csc->avals);

    /* Building ia, ja and avals */
    *colptr = 1;
    l = 1; /* Column index in the matrix ((i-1) * dim1 * dim2 + (j-1) * dim1 + k-1) */
    for(i=1; i<=dim3; i++)
    {
        for(j=1; j<=dim2; j++)
        {
            for(k=1; k<=dim1; k++)
            {

                colptr[1] = colptr[0];

                /* Diagonal value */
#if !defined(PRECISION_p)
                *rowptr = l;
                *valptr = (pastix_complex64_t) 6.;
                if (k == dim1 || k == 1)
                    *valptr -= (pastix_complex64_t) 1.;
                if (j == dim2 || j == 1)
                    *valptr -= (pastix_complex64_t) 1.;
                if (i == dim3 || i == 1)
                    *valptr -= (pastix_complex64_t) 1.;
#endif
                valptr++; rowptr++; colptr[1]++;

                /* Connexion along dimension 1 */
                if (k < dim1) {
#if !defined(PRECISION_p)
                    *rowptr = l+1;
                    *valptr = (pastix_complex64_t)-1.;
#endif
                    valptr++; rowptr++; colptr[1]++;
                }

                /* Connexion along dimension 2 */
                if (j < dim2) {
#if !defined(PRECISION_p)
                    *rowptr = l+dim1;
                    *valptr = (pastix_complex64_t)-1.;
#endif
                    valptr++; rowptr++; colptr[1]++;
                }

                /* Connexion along dimension 3 */
                if (i < dim3) {
#if !defined(PRECISION_p)
                    *rowptr = l+dim1*dim2;
                    *valptr = (pastix_complex64_t)-1.;
#endif
                    valptr++; rowptr++; colptr[1]++;
                }

                colptr++; l++;
            }
        }
    }

    assert( (csc->colptr[ csc->gN ] - csc->colptr[0]) == nnz );
}
