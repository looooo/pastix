/**
 * @file laplacian.c
 *
 *  $COPYRIGHTS$
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "common.h"
#include "drivers.h"
#include "pastix.h"

static inline void
laplacian_usage(void)
{
    fprintf(stderr,
            "Usage: genLaplacian( \"[<type>:]<dim1>[:<dim2>[:<dim3>]]\" )\n"
            "   <type> p = pattern only\n"
            "          s = real simple\n"
            "          d = real double [default]\n"
            "          c = complex simple\n"
            "          z = complex double\n"
            "   <dim1> size of the first dimension of the 1D|2D|3D laplacian\n"
            "   <dim2> size of the second dimension of the 2D|3D laplacian\n"
            "   <dim3> size of the third dimension of the 3D laplacian\n"
            "   Example:\n"
            "     genLaplacian( \"z:10:20\" ) generates a 2D complex laplacian matrix of size 200.\n"
            );
}

static inline int
laplacian_parse_info( const char   *filename,
                      pastix_csc_t *csc,
                      pastix_int_t *dim1,
                      pastix_int_t *dim2,
                      pastix_int_t *dim3 )
{
    long tmp1, tmp2, tmp3;

    /* Look for the datatype */
    {
        char flt;
        char *tmpf = strdup( filename );

        if ( sscanf( filename, "%c:%s", &flt, tmpf ) == 2 ) {
            filename += 2;
            switch( flt ){
            case 'Z':
            case 'z':
                csc->flttype = PastixComplex64;
                break;

            case 'C':
            case 'c':
                csc->flttype = PastixComplex32;
                break;

            case 'D':
            case 'd':
                csc->flttype = PastixDouble;
                break;

            case 'S':
            case 's':
                csc->flttype = PastixFloat;
                break;

            case 'P':
            case 'p':
                csc->flttype = PastixPattern;
                break;

            default:
                laplacian_usage();
                return PASTIX_ERR_BADPARAMETER;
            }
        }
        else {
            csc->flttype = PastixDouble;
        }

        free(tmpf);
    }

    /* Scan the dimensions */
    *dim1 = *dim2 = *dim3 = 0;

    if ( sscanf( filename, "%ld:%ld:%ld", &tmp1, &tmp2, &tmp3 ) == 3 ) {
        *dim1 = (pastix_int_t)tmp1;
        *dim2 = (pastix_int_t)tmp2;
        *dim3 = (pastix_int_t)tmp3;
        csc->gN = (*dim1)*(*dim2)*(*dim3);
    }
    else if ( sscanf( filename, "%ld:%ld", &tmp1, &tmp2 ) == 2 ) {
        *dim1 = (pastix_int_t)tmp1;
        *dim2 = (pastix_int_t)tmp2;
        csc->gN = (*dim1)*(*dim2);
    }
    else if ( sscanf( filename, "%ld", &tmp1 ) == 1 ) {
        *dim1 = (pastix_int_t)tmp1;
        csc->gN = *dim1;
    }
    else {
        laplacian_usage();
        return PASTIX_ERR_BADPARAMETER;
    }

    /* One of the dimension was set to 0 */
    if ( csc->gN == 0 ) {
        laplacian_usage();
        return PASTIX_ERR_BADPARAMETER;
    }

    csc->n = csc->gN;
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * File: laplacian.c
 *
 * z_gen1Dlaplacian - Generate a 1D laplacian
 *
 * Example :
 * >  2 -1  0  0
 * > -1  2 -1  0
 * >  0 -1  2 -1
 * >  0  0 -1  2
 *
 *******************************************************************************
 *
 * @param[out] csc
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
static inline void
z_gen1Dlaplacian( pastix_csc_t  *csc,
                  void         **rhs,
                  long           dim1 )
{

    pastix_int_t i;
    pastix_int_t j;
    pastix_complex64_t * valptr;
    pastix_complex64_t * rhsptr;
    pastix_int_t nnzeros = 2*(csc->gN) - 1;
    csc->n       = csc->gN;
    csc->colptr  = NULL;
    csc->rows    = NULL;
    csc->avals   = NULL;

    fprintf(stderr, "Laplacien 1D, n = %d\n",csc->n);

    assert( csc->gN == dim1 );

    /* Allocating */
    if ((NULL == (csc->colptr = (pastix_int_t *)      malloc(((csc->n)+1) *sizeof(pastix_int_t))  ))     ||
        (NULL == (csc->rows   = (pastix_int_t *)      malloc((nnzeros)    *sizeof(pastix_int_t))  ))     ||
        (NULL == (csc->avals  = (pastix_complex64_t *)malloc((nnzeros)    *sizeof(pastix_complex64_t)))) ||
        (NULL == (*rhs        = (pastix_complex64_t *)malloc((csc->n)     *sizeof(pastix_complex64_t))))/* ||
                                                                                                         (NULL == (*type       = (char *)          malloc(4         *sizeof(char))      )) ||
                                                                                                         (NULL == (*rhstype    = (char *)          malloc(4         *sizeof(char))          ))*/)
    {
        fprintf(stderr, "Laplacien 1D, error in CSC allocation\n");
        //     if (*type != NULL)
        //     {
        //       free(*type);
        //       *type = NULL;
        //     }
        if (*rhs != NULL)
        {
            free(*rhs);
            *rhs = NULL;
        }
        if (csc->avals != NULL)
        {
            free(csc->avals);
            csc->avals = NULL;
        }
        if (csc->rows != NULL)
        {
            free(csc->rows);
            csc->rows = NULL;
        }
        if (csc->colptr != NULL)
        {
            free(csc->colptr);
            csc->colptr = NULL;
        }
        return;
    }

    /* Building ia, ja and avals and rhs*/
    j=0;
    valptr = csc->avals;
    rhsptr = *rhs;
    for (i = 0; i < (csc->gN); i++,rhsptr++)
    {
        (csc->colptr)[i] = j+1;
        /* ONLY triangular inferior matrix */
        /*       if (i != 0) */
        /*  { */
        /*    (csc->rows)[j]    = i; */
        /*    (csc->avals)[j] = -1; */
        /*    j++; */
        /*  } */
        (csc->rows)[j]    = i+1;
        *valptr = 2;
        j++;
        valptr++;
        if (i != (csc->gN)-1)
        {
            (csc->rows)[j]    = i+2;

            if (csc->flttype==PastixComplex64 || csc->flttype==PastixComplex32)
            {
                *valptr = - 1 +  2* _Complex_I;
            }else{
                *valptr = -1;
            }
            j++;
            valptr++;
        }
        *rhsptr = 0;
        if (csc->flttype==PastixComplex64 || csc->flttype==PastixComplex32)
        {
            *rhsptr = -4 + 4 * _Complex_I;
        }
    }
    (csc->colptr)[i] = j+1;
    if (csc->flttype==PastixComplex64 || csc->flttype==PastixComplex32)
    {
        rhsptr = *rhs;
        *rhsptr = -1 + 3 *  _Complex_I;
        rhsptr += csc->gN -1;
        *rhsptr = -1 + 3 *  _Complex_I;
    }else{
        rhsptr = *rhs;
        *rhsptr = 1;
        rhsptr += csc->gN -1;
        *rhsptr = 1;
    }
    /* type and rhstype */
    // #if (defined TYPE_COMPLEX && !defined SYMMETRIC_LAPLACIAN)
    //   sprintf (*type, "RHA");
    // #else
    //   sprintf (*type, "RSA");
    // #endif
    //   sprintf (*rhstype,"???");
    csc->mtxtype = PastixSymmetric;
    csc->fmttype = PastixCSC;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * File: laplacian.c
 *
 * z_gen2Dlaplacian - Generate a 2D laplacian
 *
 * Example :
 * >  4 -1  0 -1  0  0
 * > -1  4 -1  0 -1  0
 * >  0 -1  4  0  0 -1
 * > -1  0  0  4 -1  0
 * >  0 -1  0 -1  4 -1
 * >  0  0 -1  0 -1  4
 *
 *******************************************************************************
 *
 * @param[out] csc
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
static inline void
z_gen2Dlaplacian( pastix_csc_t  *csc,
                  void         **rhs,
                  long           dim1,
                  long           dim2 )
{

    pastix_int_t i;
    pastix_int_t j;
    pastix_int_t k;
    pastix_complex64_t * valptr;
    pastix_complex64_t * rhsptr;
    pastix_int_t nnzeros = (2*(dim1)-1)*dim2 + (dim2-1)*dim1;
    csc->n       = csc->gN;
    csc->colptr  = NULL;
    csc->rows    = NULL;
    csc->avals   = NULL;
    *rhs         = NULL;

    fprintf(stderr, "Laplacien 2D, n = %d\n",csc->n);

    assert( csc->gN == dim1*dim2 );

    /* Allocating */
    if ((NULL == (csc->colptr = (pastix_int_t *)      malloc(((csc->gN)+1)*sizeof(pastix_int_t))  ))     ||
        (NULL == (csc->rows   = (pastix_int_t *)      malloc((nnzeros)    *sizeof(pastix_int_t))  ))     ||
        (NULL == (csc->avals  = (pastix_complex64_t *)malloc((nnzeros)    *sizeof(pastix_complex64_t)))) ||
        (NULL == (*rhs        = (pastix_complex64_t *)calloc((csc->gN)    ,sizeof(pastix_complex64_t)))) )
    {
        fprintf(stderr, "Laplacien 2D, error in CSC allocation\n");
        if (*rhs != NULL)
        {
            free(*rhs);
            *rhs = NULL;
        }
        if (csc->avals != NULL)
        {
            free(csc->avals);
            csc->avals = NULL;
        }
        if (csc->rows != NULL)
        {
            free(csc->rows);
            csc->rows = NULL;
        }
        if (csc->colptr != NULL)
        {
            free(csc->colptr);
            csc->colptr = NULL;
        }
        return;
    }

    /* Building ia, ja and avals and rhs*/

    csc->colptr[0] = 1;
    valptr = csc->avals;
    rhsptr = *rhs;
    *rhsptr = 1.;
    *(rhsptr+csc->gN-1) = 1.;
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
                *valptr = 4.;
                *(valptr+1) = -1.;
                *(valptr+2) = -1.;
                valptr += 3;
            }
            else if(j==dim1 && i!=dim2-1)
            {
                csc->colptr[k] = csc->colptr[k-1]+2;
                csc->rows[csc->colptr[k-1]-1]=k;
                csc->rows[csc->colptr[k-1]  ]=k+dim1;
                *valptr = 4.;
                *(valptr+1) = -1.;
                valptr += 2;
            }
            else if(j!=dim1 && i==dim2-1)
            {
                csc->colptr[k] = csc->colptr[k-1]+2;
                csc->rows[csc->colptr[k-1]-1]=k;
                csc->rows[csc->colptr[k-1]  ]=k+1;
                *valptr = 4.;
                *(valptr+1) = -1;
                valptr += 2;
            }
            else /* if(j==dim1 && i==dim2-1) */
            {
                csc->colptr[k] = csc->colptr[k-1]+1;
                csc->rows[csc->colptr[k-1]-1]=k;
                *valptr = 4.;
                valptr += 1;
            }
        }
    }

    csc->mtxtype = PastixSymmetric;
    csc->fmttype = PastixCSC;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * File: laplacian.c
 *
 * z_gen3Dlaplacian - Generate a 3D laplacian
 *
 * Example :
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
 * @param[out] csc
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
static inline void
z_gen3Dlaplacian( pastix_csc_t  *csc,
                  void         **rhs,
                  long           dim1,
                  long           dim2,
                  long           dim3 )
{

    pastix_int_t i;
    pastix_int_t j;
    pastix_int_t k;
    pastix_int_t l;
    pastix_complex64_t * valptr;
    pastix_complex64_t * rhsptr;
    pastix_int_t nnzeros = (2*(dim1)-1)*dim2*dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1);
    csc->n       = csc->gN;
    csc->colptr  = NULL;
    csc->rows    = NULL;
    csc->avals   = NULL;
    *rhs         = NULL;

    assert( csc->gN == dim1*dim2*dim3 );

    /* Allocating */
    if ((NULL == (csc->colptr = (pastix_int_t *)      malloc(((csc->gN)+1)*sizeof(pastix_int_t))  ))     ||
        (NULL == (csc->rows   = (pastix_int_t *)      malloc((nnzeros)    *sizeof(pastix_int_t))  ))     ||
        (NULL == (csc->avals  = (pastix_complex64_t *)malloc((nnzeros)    *sizeof(pastix_complex64_t)))) ||
        (NULL == (*rhs        = (pastix_complex64_t *)calloc((csc->gN)    ,sizeof(pastix_complex64_t)))) )
    {
        fprintf(stderr, "Laplacien 3D, error in CSC allocation\n");
        if (*rhs != NULL)
        {
            free(*rhs);
            *rhs = NULL;
        }
        if (csc->avals != NULL)
        {
            free(csc->avals);
            csc->avals = NULL;
        }
        if (csc->rows != NULL)
        {
            free(csc->rows);
            csc->rows = NULL;
        }
        if (csc->colptr != NULL)
        {
            free(csc->colptr);
            csc->colptr = NULL;
        }
        return;
    }

    /* Building ia, ja and avals and rhs*/

    csc->colptr[0] = 1;
    valptr = csc->avals;
    rhsptr = *rhs;
    *rhsptr = 1.;
    *(rhsptr+csc->gN-1) = 1.;
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
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    *(valptr+2) = -1.;
                    *(valptr+3) = -1.;
                    valptr += 4;
                }
                else if(k==dim1 && j!=dim2-1 && i!=dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+3;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+dim1;
                    csc->rows[csc->colptr[l-1]+1]=l+dim1*dim2;
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    *(valptr+2) = -1.;
                    valptr += 3;
                }
                else if(k!=dim1 && j==dim2-1 && i!=dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+3;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+1;
                    csc->rows[csc->colptr[l-1]+1]=l+dim1*dim2;
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    *(valptr+2) = -1.;
                    valptr += 3;
                }
                else if(k!=dim1 && j!=dim2-1 && i==dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+3;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+1;
                    csc->rows[csc->colptr[l-1]+1]=l+dim1;
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    *(valptr+2) = -1.;
                    valptr += 3;
                }
                else if(k==dim1 && j==dim2-1 && i!=dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+2;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+dim1*dim2;
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    valptr += 2;
                }
                else if(k!=dim1 && j==dim2-1 && i==dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+2;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+1;
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    valptr += 2;
                }
                else if(k==dim1 && j!=dim2-1 && i==dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+2;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    csc->rows[csc->colptr[l-1]  ]=l+dim1;
                    *valptr = 6.;
                    *(valptr+1) = -1.;
                    valptr += 2;
                }
                else if(k==dim1 && j==dim2-1 && i==dim3-1)
                {
                    csc->colptr[l] = csc->colptr[l-1]+1;
                    csc->rows[csc->colptr[l-1]-1]=l;
                    *valptr = 6.;
                    valptr += 1;
                }
            }
        }
    }

    csc->mtxtype = PastixSymmetric;
    csc->fmttype = PastixCSC;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * genLaplacian - Generate a Laplacian of size csc->n
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          Path to the directory containing matrix.
 *
 * @param[out] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
 *
 * @param[out] rhs
 *          At exit, contains the right hand side member.
 *
 *******************************************************************************/
int
genLaplacian( const char    *filename,
              pastix_csc_t  *csc,
              void         **rhs )
{
    pastix_int_t dim1, dim2, dim3;
    int rc;

    rc = laplacian_parse_info(filename, csc, &dim1, &dim2, &dim3);
    if (rc != PASTIX_SUCCESS)
        return rc;

    if( dim3 > 0 ) {
        z_gen3Dlaplacian(csc, rhs, dim1, dim2, dim3);
    }
    else if (dim2 != 0) {
        z_gen2Dlaplacian(csc, rhs, dim1, dim2);
    }
    else {
        z_gen1Dlaplacian(csc, rhs, dim1);
    }

    return PASTIX_SUCCESS;
}
