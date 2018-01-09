/**
 *
 * @file laplacian.c
 *
 * @copyright 2011-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @date 2011-11-11
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "common.h"
#include "spm_drivers.h"
#include "drivers/laplacian.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_driver
 *
 * laplacian_usage - Print the usage information to generate correct Laplacian
 * matrices.
 *
 *******************************************************************************/
static inline void
laplacian_usage(void)
{
    fprintf(stderr,
            "Usage: genLaplacian( \"[<type>:]<dim1>[:<dim2>[:<dim3>[:<alpha>[:<beta>]]]]\" )\n"
            "   Generate a Laplacian matrix M, of the form alpha * D - beta * A,\n"
            "   where D is the degree matrix, and A the adjacency matrix.\n"
            "   <type> p = pattern only\n"
            "          s = real simple\n"
            "          d = real double [default]\n"
            "          c = complex simple\n"
            "          z = complex double\n"
            "   <dim1> size of the first dimension of the laplacian\n"
            "   <dim2> size of the second dimension of the laplacian\n"
            "   <dim3> size of the third dimension of the laplacian\n"
            "   Example:\n"
            "     genLaplacian( \"z:10:20\" )        generates a 2D complex double laplacian matrix of size 200.\n"
            "     genLaplacian( \"10:1:10:2.:0.5\" ) generates a 2D real double laplacian matrix of size 100 where M = 2. * D - 0.5 * A.\n"
            "     genLaplacian( \"s:10\" )           generates a 1D real single laplacian matrix of size 10.\n"
            );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_driver
 *
 * @brief Parse information given through the filename string to configure the
 * laplacian matrix to generate.
 *
 * The laplacian will be of size dim1 * dim2 * dim3, and will be equal to
 *     \[ M = \alpha * D - \beta * A \]
 *
 * where D is the degree matrix, and A the adjacency matrix.
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          Configuration string of the Laplacian. See laplacian_usage() for
 *          more information.
 *
 * @param[inout] spm
 *          At start, an allocated spm structure that will store the Lapalcian
 *          matrix.
 *          At exit, the fields of the spm are initialized and especially the
 *          type, symmetry and number of unknows are setup.
 *
 * @param[out] dim1
 *          The first dimension of the laplacian
 *
 * @param[out] dim2
 *          The second dimension of the laplacian
 *
 * @param[out] dim3
 *          The third dimension of the laplacian
 *
 * @param[out] alpha
 *          The alpha coefficient for the degree matrix
 *
 * @param[out] beta
 *          The beta coefficient for the adjacency matrix
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the matrix has been generated successfully
 * @retval PASTIX_ERR_BADPARAMETER if the configuration string is incorrect
 *
 *******************************************************************************/
int
laplacian_parse_info( const char   *filename,
                      pastix_spm_t *spm,
                      pastix_int_t *dim1,
                      pastix_int_t *dim2,
                      pastix_int_t *dim3,
                      double       *alpha,
                      double       *beta )
{
    double val1, val2;
    long tmp1, tmp2, tmp3;
    spm->colptr   = NULL;
    spm->rowptr   = NULL;
    spm->values   = NULL;
    spm->loc2glob = NULL;

    *alpha = 1.;
    *beta = 1.;

    /* Look for the datatype */
    {
        char flt;
        char *tmpf = strndup( filename, 256 );

        if ( sscanf( filename, "%c:%254s", &flt, tmpf ) == 2 ) {
            filename += 2;
            switch( flt ){
            case 'Z':
            case 'z':
                spm->flttype = PastixComplex64;
                break;

            case 'C':
            case 'c':
                spm->flttype = PastixComplex32;
                break;

            case 'D':
            case 'd':
                spm->flttype = PastixDouble;
                break;

            case 'S':
            case 's':
                spm->flttype = PastixFloat;
                break;

            case 'P':
            case 'p':
                spm->flttype = PastixPattern;
                break;

            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                spm->flttype = PastixDouble;
                /*
                 * The first dimension is only one character long so we come
                 * back to the beginning of the string
                 */
                filename -= 2;
                break;

            default:
                laplacian_usage();
                return PASTIX_ERR_BADPARAMETER;
            }
        }
        else {
            spm->flttype = PastixDouble;
        }

        free(tmpf);
    }

    /* Scan the dimensions */
    *dim1 = *dim2 = *dim3 = 1;

    if ( sscanf( filename, "%ld:%ld:%ld:%lf:%lf", &tmp1, &tmp2, &tmp3, &val1, &val2 ) == 5 ) {
        *dim1 = (pastix_int_t)tmp1;
        *dim2 = (pastix_int_t)tmp2;
        *dim3 = (pastix_int_t)tmp3;
        *alpha = val1;
        *beta  = val2;
        spm->gN = (*dim1)*(*dim2)*(*dim3);
    }
    else if ( sscanf( filename, "%ld:%ld:%ld:%lf", &tmp1, &tmp2, &tmp3, &val1 ) == 4 ) {
        *dim1 = (pastix_int_t)tmp1;
        *dim2 = (pastix_int_t)tmp2;
        *dim3 = (pastix_int_t)tmp3;
        *alpha = val1;
        spm->gN = (*dim1)*(*dim2)*(*dim3);
    }
    else if ( sscanf( filename, "%ld:%ld:%ld", &tmp1, &tmp2, &tmp3 ) == 3 ) {
        *dim1 = (pastix_int_t)tmp1;
        *dim2 = (pastix_int_t)tmp2;
        *dim3 = (pastix_int_t)tmp3;
        spm->gN = (*dim1)*(*dim2)*(*dim3);
    }
    else if ( sscanf( filename, "%ld:%ld", &tmp1, &tmp2 ) == 2 ) {
        *dim1 = (pastix_int_t)tmp1;
        *dim2 = (pastix_int_t)tmp2;
        spm->gN = (*dim1)*(*dim2);
    }
    else if ( sscanf( filename, "%ld", &tmp1 ) == 1 ) {
        *dim1 = (pastix_int_t)tmp1;
        spm->gN = *dim1;
    }
    else {
        laplacian_usage();
        return PASTIX_ERR_BADPARAMETER;
    }

    /* One of the dimension was set to 0 */
    if ( spm->gN == 0 ) {
        laplacian_usage();
        return PASTIX_ERR_BADPARAMETER;
    }

    spm->n = spm->gN;
    return PASTIX_SUCCESS;
}

static void (*laplacian_7points[6])(pastix_spm_t *, pastix_int_t, pastix_int_t, pastix_int_t, pastix_fixdbl_t, pastix_fixdbl_t) =
{
    p_spmLaplacian_7points,
    NULL,
    s_spmLaplacian_7points,
    d_spmLaplacian_7points,
    c_spmLaplacian_7points,
    z_spmLaplacian_7points
};

static void (*extended_laplacian_table2D[6])(pastix_spm_t *, pastix_int_t, pastix_int_t) =
{
    p_spmExtendedLaplacian2D,
    NULL,
    s_spmExtendedLaplacian2D,
    d_spmExtendedLaplacian2D,
    c_spmExtendedLaplacian2D,
    z_spmExtendedLaplacian2D
};

static void (*extended_laplacian_table3D[6])(pastix_spm_t *, pastix_int_t, pastix_int_t, pastix_int_t) =
{
    p_spmExtendedLaplacian3D,
    NULL,
    s_spmExtendedLaplacian3D,
    d_spmExtendedLaplacian3D,
    c_spmExtendedLaplacian3D,
    z_spmExtendedLaplacian3D
};

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_driver
 *
 * genLaplacian - Generate a Laplacian of size spm->n
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          Configuration string of the Laplacian.
 *          [<type>:]<dim1>[:<dim2>[:<dim3>]]
 *             <type> p = pattern only\n"
 *                    s = real simple\n"
 *                    d = real double [default]\n"
 *                    c = complex simple\n"
 *                    z = complex double\n"
 *             <dim1> size of the first dimension of the 1D|2D|3D laplacian\n"
 *             <dim2> size of the second dimension of the 2D|3D laplacian\n"
 *             <dim3> size of the third dimension of the 3D laplacian\n"
 *
 * @param[inout] spm
 *          At start, an allocated spm structure.
 *          At exit, contains a laplacian matrix in the spm format.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the matrix has been generated successfully
 *      \retval PASTIX_ERR_BADPARAMETER if the configuration string is incorrect
 *
 *******************************************************************************/
int
genLaplacian( const char    *filename,
              pastix_spm_t  *spm )
{
    pastix_int_t dim1, dim2, dim3;
    double alpha = 1.;
    double beta = 1.;
    int rc;

    rc = laplacian_parse_info(filename, spm, &dim1, &dim2, &dim3, &alpha, &beta );
    if (rc != PASTIX_SUCCESS)
        return rc;

    laplacian_7points[spm->flttype](spm, dim1, dim2, dim3, alpha, beta);

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_driver
 *
 * genExtendedLaplacian - Generate a extended Laplacian of size spm->n
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          Configuration string of the Laplacian.
 *          [<type>:]<dim1>[:<dim2>[:<dim3>]]
 *             <type> p = pattern only
 *                    s = real simple
 *                    d = real double [default]
 *                    c = complex simple
 *                    z = complex double
 *             <dim1> size of the first dimension of the 1D|2D|3D laplacian
 *             <dim2> size of the second dimension of the 2D|3D laplacian
 *             <dim3> size of the third dimension of the 3D laplacian
 *
 * @param[inout] spm
 *          At start, an allocated spm structure.
 *          At exit, contains a laplacian matrix in the spm format.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the matrix has been generated successfully
 *      \retval PASTIX_ERR_BADPARAMETER if the configuration string is incorrect
 *
 *******************************************************************************/
int
genExtendedLaplacian( const char    *filename,
                      pastix_spm_t  *spm )
{
    pastix_int_t dim1, dim2, dim3;
    double alpha = 1.;
    double beta = 1.;
    int rc;

    rc = laplacian_parse_info(filename, spm, &dim1, &dim2, &dim3, &alpha, &beta);
    if (rc != PASTIX_SUCCESS)
        return rc;

    if( dim3 > 0 ) {
        extended_laplacian_table3D[spm->flttype](spm, dim1, dim2, dim3);
    }
    else if (dim2 > 0) {
        extended_laplacian_table2D[spm->flttype](spm, dim1, dim2);
    }

    return PASTIX_SUCCESS;
}
