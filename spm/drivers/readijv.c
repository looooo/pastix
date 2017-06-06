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
#include "spm_drivers.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_driver
 *
 * threeFilesReadHeader - Read header from three file IJV format.
 *
 *******************************************************************************
 *
 * @param[in] infile
 *          The opened header file
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
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the information has been read successfully
 *      \retval PASTIX_ERR_BADPARAMETER if the header has a wrong format
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
 * @ingroup pastix_spm_driver
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
 * @param[in] dirname
 *          Directory that contains the files.
 *
 * @param[out] spm
 *          At exit, contains the matrix in ijv format.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS if the matrix has been read successfully
 *      \retval PASTIX_ERR_IO if a problem occurs while reading the files
 *      \retval PASTIX_ERR_BADPARAMETER if a problem occurs while opening the
 *              files
 *
 *******************************************************************************/
int
readIJV( const char   *dirname,
         pastix_spm_t *spm )
{

    FILE *iafile, *jafile, *rafile;
    FILE *hdrfile;
    char *filename;
    pastix_int_t *tempcol;
    pastix_int_t *temprow;
    double       *tempval;
    pastix_int_t  i, Nrow, Ncol, Nnzero;

    filename = malloc(strlen(dirname)+10);

    spm->flttype = PastixDouble;
    spm->mtxtype = PastixGeneral;
    spm->fmttype = PastixIJV;
    spm->dof     = 1;
    spm->loc2glob= NULL;

    /* Read the header information */
    {
        sprintf(filename,"%s/header",dirname);
        hdrfile = fopen (filename,"r");
        if (hdrfile == NULL)
        {
            fprintf(stderr,"readijv: Cannot open the header file (%s)\n", filename);
            free(filename);
            return PASTIX_ERR_BADPARAMETER;
        }
        threeFilesReadHeader(hdrfile, &Nrow, &Ncol, &Nnzero);
        fclose(hdrfile);
        free(filename);
    }

    spm->gN      = Ncol;
    spm->n       = Ncol;
    spm->gnnz    = Nnzero;
    spm->nnz     = Nnzero;
    spm->colptr = (pastix_int_t *) malloc(Nnzero*sizeof(pastix_int_t));
    spm->rowptr = (pastix_int_t *) malloc(Nnzero*sizeof(pastix_int_t));
    spm->values = (double *)       malloc(Nnzero*sizeof(double));

    /* Open the 3 files */
    sprintf(filename,"%s/ia_threeFiles",dirname);
    iafile = fopen(filename,"r");
    if (iafile == NULL)
    {
        fprintf(stderr,"readijv: Cannot open the ia file (%s)\n", filename);
        free(filename);
        return PASTIX_ERR_BADPARAMETER;
    }
    free(filename);

    sprintf(filename,"%s/ja_threeFiles",dirname);
    jafile = fopen(filename,"r");
    if (jafile == NULL)
    {
        fprintf(stderr,"readijv: Cannot open the ja file (%s)\n", filename);
        fclose(iafile);
        free(filename);
        return PASTIX_ERR_BADPARAMETER;
    }
    free(filename);

    sprintf(filename,"%s/ra_threeFiles",dirname);
    rafile = fopen(filename,"r");
    if (rafile == NULL)
    {
        fprintf(stderr,"readijv: Cannot open the ra file (%s)\n", filename);
        fclose(iafile);
        fclose(jafile);
        free(filename);
        return PASTIX_ERR_BADPARAMETER;
    }
    free(filename);

    /* Read the files */
    tempcol = spm->colptr;
    temprow = spm->rowptr;
    tempval = spm->values;

    for (i=0; i<Nnzero; i++, tempcol++, temprow++, tempval++)
    {
        long temp1, temp2;
        double temp3;

        if (( 1 != fscanf(iafile,"%ld\n", &temp1)) ||
            ( 1 != fscanf(jafile,"%ld\n", &temp2)) ||
            ( 1 != fscanf(rafile,"%le\n", &temp3)) )
        {
            fprintf(stderr, "ERROR: reading matrix\n");
            fclose(iafile);
            fclose(jafile);
            fclose(rafile);
            free(filename);
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
