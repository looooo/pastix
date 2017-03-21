/**
 * @file laplacian.c
 *
 *  This file contains the user routine to generate 1, 2 and 3 dimensional Laplacian
 *  matrices.
 *  $COPYRIGHTS$
 *
 * @version 1.0.0
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
 * @ingroup pastix_csc_driver
 *
 * laplacian_usage - Print the usage information to generate correct Laplacian
 * matrices.
 *
 *******************************************************************************/
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

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * laplacian_parse_info - Parse information given through the filename string to
 * configure the laplacian matrix to generate.
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          Configuration string of the Laplacian. See laplacian_usage() for
 *          more information.
 *
 * @param[inout] csc
 *          At start, an allocated csc structure that will store the Lapalcian
 *          matrix.
 *          At exit, the fields of the csc are initialized and especially the
 *          type, symmetry and number of unknows are setup.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the matrix has been generated successfully
 *      \retval PASTIX_ERR_BADPARAMETER if the configuration string is incorrect
 *
 *******************************************************************************/
static inline int
laplacian_parse_info( const char   *filename,
                      pastix_csc_t *csc,
                      pastix_int_t *dim1,
                      pastix_int_t *dim2,
                      pastix_int_t *dim3 )
{
    long tmp1, tmp2, tmp3;
    csc->colptr   = NULL;
    csc->rowptr   = NULL;
    csc->values   = NULL;
    csc->loc2glob = NULL;

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

            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
                csc->flttype = PastixDouble;
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

static void (*laplacian_table1D[6])(pastix_csc_t*, pastix_int_t) =
{
    p_spmLaplacian1D,
    NULL,
    s_spmLaplacian1D,
    d_spmLaplacian1D,
    c_spmLaplacian1D,
    z_spmLaplacian1D
};

static void (*laplacian_table2D[6])(pastix_csc_t*, pastix_int_t, pastix_int_t) =
{
    p_spmLaplacian2D,
    NULL,
    s_spmLaplacian2D,
    d_spmLaplacian2D,
    c_spmLaplacian2D,
    z_spmLaplacian2D
};

static void (*laplacian_table3D[6])(pastix_csc_t*, pastix_int_t, pastix_int_t, pastix_int_t) =
{
    p_spmLaplacian3D,
    NULL,
    s_spmLaplacian3D,
    d_spmLaplacian3D,
    c_spmLaplacian3D,
    z_spmLaplacian3D
};

static void (*extended_laplacian_table2D[6])(pastix_csc_t*, pastix_int_t, pastix_int_t) =
{
    p_spmLaplacian2D,
    NULL,
    s_spmExtendedLaplacian2D,
    d_spmExtendedLaplacian2D,
    c_spmExtendedLaplacian2D,
    z_spmExtendedLaplacian2D
};

static void (*extended_laplacian_table3D[6])(pastix_csc_t*, pastix_int_t, pastix_int_t, pastix_int_t) =
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
 * @ingroup pastix_csc_driver
 *
 * genLaplacian - Generate a Laplacian of size csc->n
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
 * @param[inout] csc
 *          At start, an allocated csc structure.
 *          At exit, contains a laplacian matrix in the csc format.
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
              pastix_csc_t  *csc )
{
    pastix_int_t dim1, dim2, dim3;
    int rc;

    rc = laplacian_parse_info(filename, csc, &dim1, &dim2, &dim3);
    if (rc != PASTIX_SUCCESS)
        return rc;

    if( dim3 > 0 ) {
        laplacian_table3D[csc->flttype](csc, dim1, dim2, dim3);
    }
    else if (dim2 > 0) {
        laplacian_table2D[csc->flttype](csc, dim1, dim2);
    }
    else {
        laplacian_table1D[csc->flttype](csc, dim1);
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * genExtendedLaplacian - Generate a extended Laplacian of size csc->n
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
 * @param[inout] csc
 *          At start, an allocated csc structure.
 *          At exit, contains a laplacian matrix in the csc format.
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
                      pastix_csc_t  *csc )
{
    pastix_int_t dim1, dim2, dim3;
    int rc;

    rc = laplacian_parse_info(filename, csc, &dim1, &dim2, &dim3);
    if (rc != PASTIX_SUCCESS)
        return rc;

    if( dim3 > 0 ) {
        extended_laplacian_table3D[csc->flttype](csc, dim1, dim2, dim3);
    }
    else if (dim2 > 0) {
        extended_laplacian_table2D[csc->flttype](csc, dim1, dim2);
    }
    else {
        laplacian_table1D[csc->flttype](csc, dim1);
    }

    return PASTIX_SUCCESS;
}
