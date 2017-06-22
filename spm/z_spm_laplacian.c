/**
 * @file z_spm_laplacian.c
 *
 * SParse Matrix package laplacian generator routines.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
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
#include "spm.h"
#include "laplacian.h"

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 *
 * @brief Generate a laplacian matrix for a 7-points stencil
 * \[ M = \alpha * D - \beta * A \]
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
 * @param[inout] spm
 *          At start, an allocated spm structure.
 *          Contains the size of the laplacian in spm->n.
 *          At exit, contains the matrix in csc format.
 *
 * @param[in] dim1
 *          contains the first dimension of the grid of the laplacian.
 *
 * @param[in] dim2
 *          contains the second dimension of the grid of the laplacian.
 *
 * @param[in] dim3
 *          contains the third dimension of the grid of the laplacian.
 *
 * @param[in] alpha
 *          The alpha coefficient for the degree matrix
 *
 * @param[in] beta
 *          The beta coefficient for the adjacency matrix
 *
 *******************************************************************************/
void
z_spmLaplacian_7points( pastix_spm_t   *spm,
                        pastix_int_t    dim1,
                        pastix_int_t    dim2,
                        pastix_int_t    dim3,
                        pastix_fixdbl_t alpha,
                        pastix_fixdbl_t beta )
{

    pastix_complex64_t *valptr;
    pastix_int_t *colptr, *rowptr;
    pastix_int_t i, j, k, l;
    pastix_int_t nnz = (2*(dim1)-1)*dim2*dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1);

    spm->mtxtype  = PastixHermitian;
    spm->flttype  = PastixComplex64;
    spm->fmttype  = PastixCSC;
    spm->gnnz     = nnz;
    spm->nnz      = nnz;
    spm->dof      = 1;

    assert( spm->gN == dim1*dim2*dim3 );

    /* Allocating */
    spm->colptr = malloc((spm->n+1)*sizeof(pastix_int_t));
    spm->rowptr = malloc(nnz       *sizeof(pastix_int_t));
    assert( spm->colptr );
    assert( spm->rowptr );

#if !defined(PRECISION_p)
    spm->values = malloc(nnz       *sizeof(pastix_complex64_t));
    assert( spm->values );
#endif

    /* Building ia, ja and values*/
    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (pastix_complex64_t*)(spm->values);

    /* Building ia, ja and values*/
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
                *rowptr = l;
#if !defined(PRECISION_p)
                *valptr = 6. * alpha;
                if (k == 1)
                    *valptr -= alpha;
                if (k == dim1)
                    *valptr -= alpha;
                if (j == 1)
                    *valptr -= alpha;
                if (j == dim2)
                    *valptr -= alpha;
                if (i == 1)
                    *valptr -= alpha;
                if (i == dim3)
                    *valptr -= alpha;
#endif
                valptr++; rowptr++; colptr[1]++;

                /* Connexion along dimension 1 */
                if (k < dim1) {
                    *rowptr = l+1;
#if !defined(PRECISION_p)
                    *valptr = -beta;
#endif
                    valptr++; rowptr++; colptr[1]++;
                }

                /* Connexion along dimension 2 */
                if (j < dim2) {
                    *rowptr = l+dim1;
#if !defined(PRECISION_p)
                    *valptr = -beta;
#endif
                    valptr++; rowptr++; colptr[1]++;
                }

                /* Connexion along dimension 3 */
                if (i < dim3) {
                    *rowptr = l+dim1*dim2;
#if !defined(PRECISION_p)
                    *valptr = -beta;
#endif
                    valptr++; rowptr++; colptr[1]++;
                }

                colptr++; l++;
            }
        }
    }

    assert( (spm->colptr[ spm->gN ] - spm->colptr[0]) == nnz );
    (void)alpha; (void)beta;
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 *
 * z_spmExtendedLaplacian2D - Generate a 2D extended laplacian matrix.
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          At start, an allocated spm structure.
 *          Contains the size of the laplacian in spm->n.
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
z_spmExtendedLaplacian2D( pastix_spm_t  *spm,
                          pastix_int_t   dim1,
                          pastix_int_t   dim2 )
{
    pastix_complex64_t *valptr;
    pastix_int_t *colptr, *rowptr;
    pastix_int_t i, j, k;
    pastix_int_t nnz = (2*(dim1)-1)*dim2 + (dim2-1)*(3*dim1-2);

    spm->mtxtype  = PastixSymmetric;
    spm->flttype  = PastixComplex64;
    spm->fmttype  = PastixCSC;
    spm->gnnz     = nnz;
    spm->nnz      = nnz;
    spm->dof      = 1;

    assert( spm->gN == dim1*dim2 );

    /* Allocating */
    spm->colptr = malloc((spm->n+1)*sizeof(pastix_int_t));
    spm->rowptr = malloc(nnz       *sizeof(pastix_int_t));
    assert( spm->colptr );
    assert( spm->rowptr );

#if !defined(PRECISION_p)
    spm->values = malloc(nnz       *sizeof(pastix_complex64_t));
    assert( spm->values );
#endif

    /* Building ia, ja and values*/
    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (pastix_complex64_t*)(spm->values);

    /* Building ia, ja and values*/
    *colptr = 1;
    k = 1; /* Column index in the matrix ((i-1) * dim1 + j-1) */
    for(i=1; i<=dim2; i++)
    {
        for(j=1; j<=dim1; j++)
        {
            colptr[1] = colptr[0];

            /* Diagonal value */
            *rowptr = k;
#if !defined(PRECISION_p)
            if ( (j == dim1 || j == 1) && (i == dim2 || i == 1) )
                *valptr = (pastix_complex64_t) 2.5;
            else if (j == dim1 || j == 1 || i == dim2 || i == 1)
                *valptr = (pastix_complex64_t) 4.;
            else
                *valptr = (pastix_complex64_t) 6.;
#endif
            valptr++; rowptr++; colptr[1]++;

            /* Connexion along dimension 1 */
            if (j < dim1) {
                *rowptr = k+1;
#if !defined(PRECISION_p)
                *valptr = (pastix_complex64_t)-1.;
#endif
                valptr++; rowptr++; colptr[1]++;
            }

            /* Connexion along dimension 2 */
            if (i < dim2)
            {
                if (j > 1)
                {
                    *rowptr = k+dim1-1;
#if !defined(PRECISION_p)
                    *valptr = (pastix_complex64_t)-0.5;
#endif
                    valptr++; rowptr++; colptr[1]++;

                }

                *rowptr = k+dim1;
#if !defined(PRECISION_p)
                *valptr = (pastix_complex64_t)-1.;
#endif
                valptr++; rowptr++; colptr[1]++;

                if (j < dim1)
                {
                    *rowptr = k+dim1+1;
#if !defined(PRECISION_p)
                    *valptr = (pastix_complex64_t)-0.5;
#endif
                    valptr++; rowptr++; colptr[1]++;

                }
            }

            colptr++; k++;
        }
    }

    assert( (spm->colptr[ spm->gN ] - spm->colptr[0]) == nnz );
}

/**
 *******************************************************************************
 *
 * @ingroup spm_dev_driver
 *
 * @brief Generate an extended laplacian matrix for a 27-points stencil with
 * \[ M = \alpha * D - \beta * A \]
 *
 *******************************************************************************
 *
 * @param[inout] spm
 *          At start, an allocated spm structure.
 *          Contains the size of the laplacian in spm->n.
 *          At exit, contains the matrix in csc format.
 *
 * @param[in] dim1
 *          contains the first dimension of the grid of the laplacian.
 *
 * @param[in] dim2
 *          contains the second dimension of the grid of the laplacian.
 *
 * @param[in] dim3
 *          contains the third dimension of the grid of the laplacian.
 *
 *******************************************************************************/
void
z_spmExtendedLaplacian3D( pastix_spm_t   *spm,
                          pastix_int_t    dim1,
                          pastix_int_t    dim2,
                          pastix_int_t    dim3 )
{

    pastix_complex64_t *valptr;
    pastix_int_t *colptr, *rowptr;
    pastix_int_t i, j, k, l;
    pastix_int_t nnz = (2*dim1-1) * dim2     * dim3
        +              (3*dim1-2) * (dim2-1) * dim3
        +             ((3*dim1-2) * dim2 + 2 * (3*dim1-2) *(dim2-1)) * (dim3-1);

    spm->mtxtype  = PastixSymmetric;
    spm->flttype  = PastixComplex64;
    spm->fmttype  = PastixCSC;
    spm->gnnz     = nnz;
    spm->nnz      = nnz;
    spm->dof      = 1;

    assert( spm->gN == dim1*dim2*dim3 );

    /* Allocating */
    spm->colptr = malloc((spm->n+1)*sizeof(pastix_int_t));
    spm->rowptr = malloc(nnz       *sizeof(pastix_int_t));
    assert( spm->colptr );
    assert( spm->rowptr );

#if !defined(PRECISION_p)
    spm->values = malloc(nnz       *sizeof(pastix_complex64_t));
    assert( spm->values );
#endif

    /* Building ia, ja and values*/
    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (pastix_complex64_t*)(spm->values);

    /* Building ia, ja and values*/
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
                *rowptr = l;
#if !defined(PRECISION_p)
                if ( (j == dim2 || j == 1) && (i == dim3 || i == 1) && (k == dim1 || i == 1) )
                    *valptr = (pastix_complex64_t) 4.75;
                else if ( (j != dim2 || j != 1) && (i == dim3 || i == 1) && (k == dim1 || i == 1) )
                    *valptr = (pastix_complex64_t) 10.;
                else if ( (j == dim2 || j == 1) && (i != dim3 || i != 1) && (k == dim1 || i == 1) )
                    *valptr = (pastix_complex64_t) 10.;
                else if ( (j == dim2 || j == 1) && (i == dim3 || i == 1) && (k != dim1 || i != 1) )
                    *valptr = (pastix_complex64_t) 10.;
                else if ( (j != dim2 || j != 1) && (i != dim3 || i != 1) && (k == dim1 || i == 1) )
                    *valptr = (pastix_complex64_t) 7.;
                else if ( (j == dim2 || j == 1) && (i != dim3 || i != 1) && (k != dim1 || i != 1) )
                    *valptr = (pastix_complex64_t) 7.;
                else if ( (j != dim2 || j != 1) && (i == dim3 || i == 1) && (k != dim1 || i != 1) )
                    *valptr = (pastix_complex64_t) 7.;
                else
                    *valptr = (pastix_complex64_t) 14.;
#endif
                valptr++; rowptr++; colptr[1]++;

                /* Connexion along dimension 1 */
                if (k < dim1) {
                    *rowptr = l+1;
#if !defined(PRECISION_p)
                    *valptr = (pastix_complex64_t)-1.;
#endif
                    valptr++; rowptr++; colptr[1]++;
                }

                /* Connexion along dimension 2 */
                if (j < dim2)
                {
                    if (k > 1)
                    {
                        *rowptr = l+dim1-1;
#if !defined(PRECISION_p)
                        *valptr = (pastix_complex64_t)-0.5;
#endif
                        valptr++; rowptr++; colptr[1]++;

                    }

                    *rowptr = l+dim1;
#if !defined(PRECISION_p)
                    *valptr = (pastix_complex64_t)-1.;
#endif
                    valptr++; rowptr++; colptr[1]++;

                    if (k < dim1)
                    {
                        *rowptr = l+dim1+1;
#if !defined(PRECISION_p)
                        *valptr = (pastix_complex64_t)-0.5;
#endif
                        valptr++; rowptr++; colptr[1]++;

                    }
                }

                /* Connexion along dimension 3 */
                if (i < dim3) {
                    if( j > 1 )
                    {
                        if (k > 1)
                        {
                            *rowptr = l+dim1*dim2-dim1-1;
#if !defined(PRECISION_p)
                            *valptr = (pastix_complex64_t)-0.25;
#endif
                            valptr++; rowptr++; colptr[1]++;

                        }

                        *rowptr = l+dim1*dim2-dim1;
#if !defined(PRECISION_p)
                        *valptr = (pastix_complex64_t)-0.5;
#endif
                        valptr++; rowptr++; colptr[1]++;

                        if (k < dim1)
                        {
                            *rowptr = l+dim1*dim2-dim1+1;
#if !defined(PRECISION_p)
                            *valptr = (pastix_complex64_t)-0.25;
#endif
                            valptr++; rowptr++; colptr[1]++;

                        }
                    }
                    if (k > 1)
                    {
                        *rowptr = l+dim1*dim2-1;
#if !defined(PRECISION_p)
                        *valptr = (pastix_complex64_t)-0.5;
#endif
                        valptr++; rowptr++; colptr[1]++;

                    }

                    *rowptr = l+dim1*dim2;
#if !defined(PRECISION_p)
                    *valptr = (pastix_complex64_t)-1.;
#endif
                    valptr++; rowptr++; colptr[1]++;

                    if (k < dim1)
                    {
                        *rowptr = l+dim1*dim2+1;
#if !defined(PRECISION_p)
                        *valptr = (pastix_complex64_t)-0.5;
#endif
                        valptr++; rowptr++; colptr[1]++;

                    }

                    if( j < dim2 )
                    {
                        if (k > 1)
                        {
                            *rowptr = l+dim1*dim2+dim1-1;
#if !defined(PRECISION_p)
                            *valptr = (pastix_complex64_t)-0.25;
#endif
                            valptr++; rowptr++; colptr[1]++;

                        }

                        *rowptr = l+dim1*dim2+dim1;
#if !defined(PRECISION_p)
                        *valptr = (pastix_complex64_t)-0.5;
#endif
                        valptr++; rowptr++; colptr[1]++;

                        if (k < dim1)
                        {
                            *rowptr = l+dim1*dim2+dim1+1;
#if !defined(PRECISION_p)
                            *valptr = (pastix_complex64_t)-0.25;
#endif
                            valptr++; rowptr++; colptr[1]++;

                        }
                    }
                }

                colptr++; l++;
            }
        }
    }

    assert( (spm->colptr[ spm->gN ] - spm->colptr[0]) == nnz );
}
