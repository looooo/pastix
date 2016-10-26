/**
 * @file readmm.c
 *
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
#include "common.h"
#include "spm_drivers.h"
#include "drivers/mmio.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * z_readMM - Read the data part of a complex matrix in Matrix Market file.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] file
 *          The file opened in readMM which contains the matrix stored in Matrix
 *          Market format.
 *
 * @param[in,out] csc
 *          At exit, the data of the matrix are stored in the csc structure.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the matrix has been read successfully
 *      \retval PASTIX_ERR_IO if a problem occured in the RSA driver
 *
 *******************************************************************************/
int
z_readMM( FILE *file,
          pastix_csc_t *csc )
{
    pastix_complex64_t *valptr;
    pastix_int_t *colptr;
    pastix_int_t *rowptr;
    pastix_int_t i;
    long row, col;
    double re, im;

    csc->values = malloc( csc->nnz * sizeof(pastix_complex64_t) );

    colptr = csc->colptr;
    rowptr = csc->rowptr;
    valptr = (pastix_complex64_t*)(csc->values);

    for (i=0; i<csc->nnz; i++, colptr++, rowptr++, valptr++)
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
 * @ingroup pastix_csc_driver
 *
 * d_readMM - Read the data part of a real matrix in Matrix Market file.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] file
 *          The file opened in readMM which contains the matrix stored in Matrix
 *          Market format.
 *
 * @param[in,out] csc
 *          At exit, the data of the matrix are stored in the csc structure.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the matrix has been read successfully
 *      \retval PASTIX_ERR_IO if a problem occured in the RSA driver
 *
 *******************************************************************************/
int
d_readMM( FILE *file,
          pastix_csc_t *csc )
{
    double       *valptr;
    pastix_int_t *colptr;
    pastix_int_t *rowptr;
    pastix_int_t i;
    long row, col;
    double re;

    csc->values = malloc( csc->nnz * sizeof(double) );

    colptr = csc->colptr;
    rowptr = csc->rowptr;
    valptr = (double*)(csc->values);

    for (i=0; i<csc->nnz; i++, colptr++, rowptr++, valptr++)
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
 * @ingroup pastix_csc_driver
 *
 * p_readMM - Read the data part of a pattern matrix in Matrix Market file.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] file
 *          The file opened in readMM which contains the matrix stored in Matrix
 *          Market format.
 *
 * @param[in,out] csc
 *          At exit, the data of the matrix are stored in the csc structure.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the matrix has been read successfully
 *      \retval PASTIX_ERR_IO if a problem occured in the RSA driver
 *
 *******************************************************************************/
int
p_readMM( FILE *file,
          pastix_csc_t *csc )
{
    pastix_int_t *colptr;
    pastix_int_t *rowptr;
    pastix_int_t i;
    long row, col;

    csc->values = NULL;

    colptr = csc->colptr;
    rowptr = csc->rowptr;

    for (i=0; i<csc->nnz; i++, colptr++, rowptr++)
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
 * @ingroup pastix_csc_driver
 *
 * readMM - Read a matrix in Matrix Market fill. This corresponds to
 * IJV format with (%d %d[ %lf[ %lf]]) format per line.
 * For more information about matrix market format see mmio.c/mmio.h
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The filename that contains the matrix stored in Matrix Market format.
 *
 * @param[in] csc
 *          At exit, contains the matrix in csc format.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the matrix has been read successfully
 *      \retval PASTIX_ERR_IO if a problem occured in the RSA driver
 *      \retval PASTIX_ERR_BADPARAMETER if the matrix is no in a supported format
 *
 *******************************************************************************/
int
readMM( const char   *filename,
        pastix_csc_t *csc )
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
        csc->flttype = PastixComplex64;
    }
    else if (mm_is_real(matcode)) {
        csc->flttype = PastixDouble;
    }
    else if (mm_is_pattern(matcode)) {
        csc->flttype = PastixPattern;
    }
    else {
        fprintf(stderr,"readmm: Unsupported type of matrix.\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Matrix structure */
    if (mm_is_general(matcode)) {
        csc->mtxtype = PastixGeneral;
    }
    else if (mm_is_symmetric(matcode)) {
        csc->mtxtype = PastixSymmetric;
    }
    else if (mm_is_hermitian(matcode)) {
        csc->mtxtype = PastixHermitian;
    }
    else {
        fprintf(stderr,"readmm: Unsupported type of matrix.\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    csc->fmttype = PastixIJV;
    csc->dof     = 1;
    csc->loc2glob= NULL;

    /* Read the size */
    {
        int m, n, nnz;
        if (mm_read_mtx_crd_size(file, &m, &n, &nnz) != 0) {
            fprintf(stderr, "readmm: error while reading matrix sizes\n");
            return PASTIX_ERR_IO;
        }

        csc->gN   = n;
        csc->n    = n;
        csc->gnnz = nnz;
        csc->nnz  = nnz;
    }

    csc->colptr = (pastix_int_t*)malloc(csc->nnz * sizeof(pastix_int_t));
    csc->rowptr = (pastix_int_t*)malloc(csc->nnz * sizeof(pastix_int_t));

    switch( csc->flttype ) {
    case PastixComplex64:
        rc = z_readMM(file,csc);
        break;

    case PastixDouble:
        rc = d_readMM(file,csc);
        break;

    case PastixPattern:
    default:
        rc = p_readMM(file,csc);
    }

    fclose(file);
    return rc;
}
