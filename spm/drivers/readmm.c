/**
 *
 * @file readmm.c
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
#include "common.h"
#include "spm_drivers.h"
#include "drivers/mmio.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_driver
 *
 * @brief Read the data part of a complex matrix in Matrix Market file.
 *
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] file
 *          The file opened in readMM which contains the matrix stored in Matrix
 *          Market format.
 *
 * @param[inout] spm
 *          At exit, the data of the matrix are stored in the spm structure.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the matrix has been read successfully
 * @retval PASTIX_ERR_IO if a problem occured in the RSA driver
 *
 *******************************************************************************/
int
z_readMM( FILE *file,
          pastix_spm_t *spm )
{
    pastix_complex64_t *valptr;
    pastix_int_t *colptr;
    pastix_int_t *rowptr;
    pastix_int_t i;
    long row, col;
    double re, im;

    spm->values = malloc( spm->nnz * sizeof(pastix_complex64_t) );

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (pastix_complex64_t*)(spm->values);

    for (i=0; i<spm->nnz; i++, colptr++, rowptr++, valptr++)
    {
        if (4 != fscanf(file,"%ld %ld %lg %lg\n", &row, &col, &re, &im))
        {
            fprintf(stderr, "readmm: erro while reading matrix file (line %ld)\n", (long)i);
            return PASTIX_ERR_IO;
        }

        *rowptr = (pastix_int_t)row;
        *colptr = (pastix_int_t)col;
        *valptr = (pastix_complex64_t)(re + im * I);
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_driver
 *
 * @brief Read the data part of a real matrix in Matrix Market file.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] file
 *          The file opened in readMM which contains the matrix stored in Matrix
 *          Market format.
 *
 * @param[inout] spm
 *          At exit, the data of the matrix are stored in the spm structure.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the matrix has been read successfully
 * @retval PASTIX_ERR_IO if a problem occured in the RSA driver
 *
 *******************************************************************************/
int
d_readMM( FILE *file,
          pastix_spm_t *spm )
{
    double       *valptr;
    pastix_int_t *colptr;
    pastix_int_t *rowptr;
    pastix_int_t i;
    long row, col;
    double re;

    spm->values = malloc( spm->nnz * sizeof(double) );

    colptr = spm->colptr;
    rowptr = spm->rowptr;
    valptr = (double*)(spm->values);

    for (i=0; i<spm->nnz; i++, colptr++, rowptr++, valptr++)
    {
        if (3 != fscanf(file,"%ld %ld %lg\n", &row, &col, &re))
        {
            fprintf(stderr, "readmm: erro while reading matrix file (line %ld)\n", (long)i);
            return PASTIX_ERR_IO;
        }

        *rowptr = (pastix_int_t)row;
        *colptr = (pastix_int_t)col;
        *valptr = re;
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_driver
 *
 * @brief Read the data part of a pattern matrix in Matrix Market file.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] file
 *          The file opened in readMM which contains the matrix stored in Matrix
 *          Market format.
 *
 * @param[inout] spm
 *          At exit, the data of the matrix are stored in the spm structure.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the matrix has been read successfully
 * @retval PASTIX_ERR_IO if a problem occured in the RSA driver
 *
 *******************************************************************************/
int
p_readMM( FILE *file,
          pastix_spm_t *spm )
{
    pastix_int_t *colptr;
    pastix_int_t *rowptr;
    pastix_int_t i;
    long row, col;

    spm->values = NULL;

    colptr = spm->colptr;
    rowptr = spm->rowptr;

    for (i=0; i<spm->nnz; i++, colptr++, rowptr++)
    {
        if (2 != fscanf(file,"%ld %ld\n", &row, &col))
        {
            fprintf(stderr, "readmm: erro while reading matrix file (line %ld)\n", (long)i);
            return PASTIX_ERR_IO;
        }

        *rowptr = (pastix_int_t)row;
        *colptr = (pastix_int_t)col;
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_driver
 *
 * @brief Read a matrix in Matrix Market fill. This corresponds to
 * IJV format with (%d %d[ %lf[ %lf]]) format per line.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The filename that contains the matrix stored in Matrix Market format.
 *
 * @param[in] spm
 *          At exit, contains the matrix in spm format.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if the matrix has been read successfully
 * @retval PASTIX_ERR_IO if a problem occured in the RSA driver
 * @retval PASTIX_ERR_BADPARAMETER if the matrix is no in a supported format
 *
 *******************************************************************************/
int
readMM( const char   *filename,
        pastix_spm_t *spm )
{
    MM_typecode matcode;
    FILE *file;
    int rc;

    file = fopen(filename,"r");
    if (file == NULL)
    {
        fprintf(stderr,"readmm: Cannot open the file (%s)\n", filename);
        return PASTIX_ERR_BADPARAMETER;
    }

    if (mm_read_banner(file, &matcode) != 0)
    {
        fprintf(stderr,"readmm: Could not process Matrix Market banner.\n");
        return PASTIX_ERR_IO;
    }

    /* Float values type */

    if (mm_is_complex(matcode)) {
        spm->flttype = PastixComplex64;
    }
    else if (mm_is_real(matcode)) {
        spm->flttype = PastixDouble;
    }
    else if (mm_is_pattern(matcode)) {
        spm->flttype = PastixPattern;
    }
    else {
        fprintf(stderr,"readmm: Unsupported type of matrix.\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Matrix structure */
    if (mm_is_general(matcode)) {
        spm->mtxtype = PastixGeneral;
    }
    else if (mm_is_symmetric(matcode)) {
        spm->mtxtype = PastixSymmetric;
    }
    else if (mm_is_hermitian(matcode)) {
        spm->mtxtype = PastixHermitian;
    }
    else {
        fprintf(stderr,"readmm: Unsupported type of matrix.\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    spm->fmttype = PastixIJV;
    spm->dof     = 1;
    spm->loc2glob= NULL;

    /* Read the size */
    {
        int m, n, nnz;
        if (mm_read_mtx_crd_size(file, &m, &n, &nnz) != 0) {
            fprintf(stderr, "readmm: error while reading matrix sizes\n");
            return PASTIX_ERR_IO;
        }

        spm->gN   = n;
        spm->n    = n;
        spm->gnnz = nnz;
        spm->nnz  = nnz;
    }

    spm->colptr = (pastix_int_t*)malloc(spm->nnz * sizeof(pastix_int_t));
    spm->rowptr = (pastix_int_t*)malloc(spm->nnz * sizeof(pastix_int_t));

    switch( spm->flttype ) {
    case PastixComplex64:
        rc = z_readMM(file, spm);
        break;

    case PastixDouble:
        rc = d_readMM(file, spm);
        break;

    case PastixPattern:
    default:
        rc = p_readMM(file, spm);
    }

    fclose(file);
    return rc;
}
