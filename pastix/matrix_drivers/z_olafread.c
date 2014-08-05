/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
/*
  File: olafread.c

  Driver for the olaf matrix format.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <stdint.h>

#include "z_pastix.h"
#include "common_drivers.h"
#include "z_olafread.h"

/*
  Function: z_olafReadHeader

  Reads header from file *infile*

  header format is :
  > Nrow
  > Nnzero

  Nrow is equal to Ncol

  Parameters:
    infile - File to read from
    Nrow   - Number of rows
    Ncol   - Number of columns
    Nnzero - Number of non zeros
    Type   - Type of the matrix
 */
void z_olafReadHeader(FILE         *infile,
                    pastix_int_t *Nrow,
                    pastix_int_t *Ncol,
                    pastix_int_t *Nnzero,
                    char         *Type)
{
  long temp1;
  if (1 != fscanf(infile, "%ld\n", &temp1))
    {
      fprintf(stderr, "ERROR: Reading matrix header\n");
      exit(1);
    }
  *Nrow = (pastix_int_t)temp1;
  if (1 != fscanf(infile, "%ld\n", &temp1))
    {
      fprintf(stderr, "ERROR: Reading matrix header\n");
      exit(1);
    }
  *Nnzero = (pastix_int_t)temp1;
  *Ncol = *Nrow;
  Type[0] = 'R';
  Type[1] = 'S';
  Type[2] = 'A';
  Type[3] = '\0';
}

/*
  Function: z_olafRead

  Reads a matrix in olaf format.

  Header format is described in <z_olafReadHeader>,
  Olaf files contains :
  colptr, row and avals in the CSC format
  > colptr[0]
  > colptr[1]
  >....
  > row[0]
  > row[1]
  > ...
  > avals[0]
  > avals[1]
  > ...


  Parameters:
    filename - Path to the directory containing hfile, ifile, jfile and afile
    Nrow     - Number of rows
    Ncol     - Number of columns
    Nnzero   - Number of non zeros
    col      - Index of first element of each column in *row* and *val*
    row      -	Row of eah element
    val      -	Value of each element
    Type     -	Type of the matrix
    RhsType  -	Type of the right hand side.
 */
void z_olafRead(char const      *filename,
              pastix_int_t    *Nrow,
              pastix_int_t    *Ncol,
              pastix_int_t    *Nnzero,
              pastix_int_t   **col,
              pastix_int_t   **row,
              pastix_complex64_t **val,
              char           **Type,
              char           **RhsType,
              pastix_complex64_t **rhs)
{
  FILE *infile;
  pastix_int_t iter,size;
  long temp1;
  double temp2;

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));
  (*RhsType)[0] = 'A';
  (*RhsType)[1] = 'A';
  (*RhsType)[2] = 'A';

#ifdef TYPE_COMPLEX
  fprintf(stderr, "\nWARNING: This drivers reads non complex matrices, imaginary part will be 0\n\n");
#endif

  infile = fopen(filename, "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "olafcsr");
      exit(-1);
    }
  z_olafReadHeader(infile, Nrow, Ncol, Nnzero, *Type);

  printf("Nrow %ld Ncol %ld Nnzero %ld\n", (long)*Nrow, (long)*Ncol, (long)*Nnzero);

  (*col) = (pastix_int_t *) malloc((*Ncol+1)*sizeof(pastix_int_t));
  (*row) = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  (*val) = (pastix_complex64_t *) malloc((*Nnzero)*sizeof(pastix_complex64_t));
  (*rhs) = (pastix_complex64_t *) malloc((*Ncol)*sizeof(pastix_complex64_t));

  if (((*col) == NULL) || ((*row) == NULL) || ((*val) == NULL) || ((*rhs) == NULL))
    fprintf(stderr, "z_olafRead : Not enough memory for \n");

  for (iter=0; iter<(*Ncol+1); iter++)
    {
      if (1 != fscanf(infile, "%ld", &temp1))
        {
          fprintf(stderr, "ERROR: Reading matrix header\n");
          exit(1);
        }

      (*col)[iter] = (pastix_int_t)temp1;
    }

  size=*Nnzero;

  for (iter=0; iter<size; iter++)
    {
      if (1 != fscanf(infile, "%ld", &temp1))
        {
          fprintf(stderr, "ERROR: Reading matrix header\n");
          exit(1);
        }

      (*row)[iter] = (pastix_int_t)temp1;
    }

  for (iter=0; iter<size; iter++)
    {
      if (1 != fscanf(infile, "%lf", &temp2))
        {
          fprintf(stderr, "ERROR: Reading matrix header\n");
          exit(1);
        }

      (*val)[iter] = (pastix_complex64_t)temp2;
    }

  fclose(infile);

  infile = fopen("olafrhs", "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "olafrhs");
      exit(-1);
    }

  for (iter=0; iter<(*Ncol); iter++)
    {
      if (1 != fscanf(infile, "%lf", &temp2))
        {
          fprintf(stderr, "ERROR: Reading matrix header\n");
          exit(1);
        }

      (*rhs)[iter] = (pastix_complex64_t)temp2;
    }

  fclose(infile);
}
