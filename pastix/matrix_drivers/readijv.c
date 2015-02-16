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
 * ******************************************************************************
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
void
threeFilesReadHeader(FILE         *infile,
                          pastix_int_t *Nrow,
                          pastix_int_t *Ncol,
                          pastix_int_t *Nnzero)
{
  char line[BUFSIZ];
  long temp1,temp2,temp3;

  /* ncol nrow nnzero */
  fgets(line,BUFSIZ,infile);

  sscanf(line, "%ld %ld %ld", &temp1,&temp2,&temp3);
  if (temp1!=temp2)
    {
      temp2=temp1;
      fgets(line,BUFSIZ,infile);
      sscanf(line, "%ld", &temp3);
    }

  *Nrow = (pastix_int_t)temp1;
  *Ncol = (pastix_int_t)temp2;
  *Nnzero = (pastix_int_t)temp3;
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
void
readIJV( const char   *dirname,
						  pastix_csc_t *csc )
{

  FILE * iaFile;
  FILE * jaFile;
  FILE * raFile;
  FILE * headerFile;
  char * filename;
  pastix_int_t * tempcol;
  pastix_int_t * temprow;
  pastix_int_t Nrow;
  pastix_int_t Ncol;
  pastix_int_t Nnzero;
  double * tempval;
  double * values;
  pastix_int_t iter,baseval;
  pastix_int_t tmp,total,pos,limit;

  filename = malloc(strlen(dirname)+10);

	csc->flttype = PastixDouble;
	csc->mtxtype = PastixGeneral;
	
#ifdef TYPE_COMPLEX
  fprintf(stderr, "\nWARNING: This drivers reads non complex matrices, imaginary part will be 0\n\n");
#endif

  sprintf(filename,"%s/header",dirname);
  headerFile = fopen (filename,"r");
  if (headerFile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      exit(-1);
    }
  threeFilesReadHeader(headerFile,&Nrow,&Ncol,&Nnzero);
  fclose (headerFile);

  sprintf(filename,"%s/ia_threeFiles",dirname);
  iaFile = fopen(filename,"r");
  if (iaFile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      exit(-1);
    }

  sprintf(filename,"%s/ja_threeFiles",dirname);
  jaFile = fopen(filename,"r");
  if (jaFile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      exit(-1);
    }
  sprintf(filename,"%s/ra_threeFiles",dirname);
  raFile = fopen(filename,"r");
  if (raFile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      exit(-1);
    }

  /* Allocation memoire */
  tempcol = (pastix_int_t *) malloc(Nnzero*sizeof(pastix_int_t));
  temprow = (pastix_int_t *) malloc(Nnzero*sizeof(pastix_int_t));
  tempval = (double *) malloc(Nnzero*sizeof(double));

  if ((tempcol==NULL) || (temprow == NULL) || (tempval == NULL))
    {
      fprintf(stderr, "threeFilesRead : Not enough memory (Nnzero %ld)\n",(long)Nnzero);
      exit(-1);
    }

  /* Remplissage */
  for (iter=0; iter<Nnzero; iter++)
    {
      long temp1,temp2;
      double tempv;
      if ( 1 != fscanf(iaFile,"%ld\n", &temp1))
        {
          fprintf(stderr, "ERROR: reading matrix\n");
          exit(1);
        }
      temprow[iter]=(pastix_int_t)temp1;
      if (1 != fscanf(jaFile,"%ld\n", &temp2))
        {
          fprintf(stderr, "ERROR: reading matrix\n");
          exit(1);
        }
      tempcol[iter]=(pastix_int_t)temp2;
      if (1 != fscanf(raFile,"%le\n", &tempv))
        {
          fprintf(stderr, "ERROR: reading matrix\n");
          exit(1);
        }
      tempval[iter]= (double)tempv;
    }

  fclose (iaFile);
  fclose (jaFile);
  fclose (raFile);

  csc->colptr = (pastix_int_t *) malloc((Nrow+1)*sizeof(pastix_int_t));
  memset(csc->colptr,0,(Nrow+1)*sizeof(pastix_int_t));
  csc->rows = (pastix_int_t *) malloc(Nnzero*sizeof(pastix_int_t));
  memset(csc->rows,0,Nnzero*sizeof(pastix_int_t));
  csc->avals = (double *) malloc(Nnzero*sizeof(double));
  values = (double *) malloc(Nnzero*sizeof(double));
  if ((csc->colptr==NULL) || (csc->rows == NULL) || (csc->avals == NULL))
    {
      fprintf(stderr, "threeFilesRead : Not enough memory (Nnzero %ld)\n",(long)Nnzero);
      exit(-1);
    }

  for (iter = 0; iter < Nnzero; iter ++)
    {
      csc->colptr[tempcol[iter]-1]++;
    }

  baseval=1; /* Attention on base a 1 */
  total = baseval;

  for (iter = 0; iter < Ncol+1; iter ++)
    {
      tmp = csc->colptr[iter];
      csc->colptr[iter]=total;
      total+=tmp;
    }

  for (iter = 0; iter < Nnzero; iter ++)
    {

      pos = csc->colptr[tempcol[iter]-1]-1;
      limit = csc->colptr[tempcol[iter]]-1;
      while(csc->rows[pos] != 0 && pos < limit)
        {
          pos++;
        }
      if (pos == limit)
        fprintf(stderr, "Erreur de lecture\n");

      csc->rows[pos] = temprow[iter];
      values[pos] = tempval[iter];
    }
	memcpy(csc->avals, values, csc->dof*sizeof(double));

  csc->n=(pastix_int_t)Nrow;
  csc->gN=csc->n;
  csc->dof=(pastix_int_t)Nnzero;
  memFree_null(values);
  memFree_null(tempval);
  memFree_null(temprow);
  memFree_null(tempcol);
}
