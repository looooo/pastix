/**
 * @file readijv.c
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
#include <stdio.h>
#include <stdlib.h>
#include "common.h"
#include "drivers.h"
#include "pastix.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * threeFilesReadHeader - Read header from three file IJV format.
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The file containing the matrix.
 *
 * @param[out] Nrow
 *          At exit, contains the number of rows of the matrix.
 *
 * @param[out] Ncol
 *          At exit, contains the number of columns of the matrix.
 *
 * @param[out] Nnzero
 *          At exit, contains the number of non zero entries of the matrix.
 *
 *******************************************************************************/
int
threeFilesReadHeader(FILE         *infile,
                     pastix_int_t *Nrow,
                     pastix_int_t *Ncol,
                     pastix_int_t *Nnzero)
{
    long temp1,temp2,temp3;

    /* ncol nrow nnzero */
    if (fscanf(infile, "%ld %ld %ld\n", &temp1, &temp2, &temp3) != 3) {
        Nrow = Ncol = Nnzero = 0;
        fprintf(stderr, "readijv: Wrong format in header file\n");
        return PASTIX_ERR_BADPARAMETER;
    }
    *Nrow   = (pastix_int_t)temp1;
    *Ncol   = (pastix_int_t)temp2;
    *Nnzero = (pastix_int_t)temp3;

    return PASTIX_SUCCESS;
}

/**
 * ******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * readIJV - Read matrix from three files IJV
 *
 * header file is "filename"/header
 * columns file is "filename"/ia_threeFiles
 * rows file is "filename"/ja_threeFiles
 * values file is "filename"/ra_threeFiles
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          The file containing the matrix.
 *
 * @param[out] csc
 *          At exit, contains the matrix in csc format.
 *
 *******************************************************************************/
int
readIJV( const char   *dirname,
         pastix_csc_t *csc )
{

    FILE *iafile, *jafile, *rafile;
    FILE *hdrfile;
    char *filename;
    pastix_int_t *tempcol;
    pastix_int_t *temprow;
    double       *tempval;
    pastix_int_t  i, Nrow, Ncol, Nnzero;

    filename = malloc(strlen(dirname)+10);

    csc->flttype = PastixDouble;
    csc->mtxtype = PastixGeneral;
    csc->fmttype = PastixIJV;
    csc->dof     = 1;
    csc->loc2glob= NULL;

    /* Read the header information */
    {
        sprintf(filename,"%s/header",dirname);
        hdrfile = fopen (filename,"r");
        if (hdrfile == NULL)
        {
            fprintf(stderr,"readijv: Cannot open the header file (%s)\n", filename);
            return PASTIX_ERR_BADPARAMETER;
        }
        threeFilesReadHeader(hdrfile, &Nrow, &Ncol, &Nnzero);
        fclose(hdrfile);
    }

    csc->gN      = Ncol;
    csc->n       = Ncol;
    csc->gnnz    = Nnzero;
    csc->nnz     = Nnzero;
    csc->colptr = (pastix_int_t *) malloc(Nnzero*sizeof(pastix_int_t));
    csc->rows   = (pastix_int_t *) malloc(Nnzero*sizeof(pastix_int_t));
    csc->avals  = (double *)       malloc(Nnzero*sizeof(double));

    /* Open the 3 files */
    sprintf(filename,"%s/ia_threeFiles",dirname);
    iafile = fopen(filename,"r");
    if (iafile == NULL)
    {
        fprintf(stderr,"readijv: Cannot open the ia file (%s)\n", filename);
        return PASTIX_ERR_BADPARAMETER;
    }

    sprintf(filename,"%s/ja_threeFiles",dirname);
    jafile = fopen(filename,"r");
    if (jafile == NULL)
    {
        fprintf(stderr,"readijv: Cannot open the ja file (%s)\n", filename);
        fclose(iafile);
        return PASTIX_ERR_BADPARAMETER;
    }

    sprintf(filename,"%s/ra_threeFiles",dirname);
    rafile = fopen(filename,"r");
    if (rafile == NULL)
    {
        fprintf(stderr,"readijv: Cannot open the ra file (%s)\n", filename);
        fclose(iafile);
        fclose(jafile);
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Read the files */
    tempcol = csc->colptr;
    temprow = csc->rows;
    tempval = csc->avals;

    for (i=0; i<Nnzero; i++, tempcol++, temprow++, tempval++)
    {
        long temp1, temp2;
        double temp3;

        if (( 1 != fscanf(iafile,"%ld\n", &temp1)) ||
            ( 1 != fscanf(jafile,"%ld\n", &temp2)) ||
            ( 1 != fscanf(rafile,"%le\n", &temp3)) )
        {
            fprintf(stderr, "ERROR: reading matrix\n");
            return PASTIX_ERR_IO;
        }
        *temprow = (pastix_int_t)temp1;
        *tempcol = (pastix_int_t)temp2;
        *tempval = temp3;
    }
    fclose(iafile);
    fclose(jafile);
    fclose(rafile);

    return PASTIX_SUCCESS;
}
