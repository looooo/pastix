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
          "\n    Laplacian usage:\n"
          "  -9 <type>:<dim1>[:<dim2>[:<dim3>]]\n"
          "  <type> (mandatory): s = real simple\n"
          "                      d = real double\n"
          "                      c = complex simple\n"
          "                      z = complex double\n"
          "  <dim1> (mandatory): size of the first dimension of the 1D|2D|3D laplacian\n"
          "  <dim2> (optionnal): size of the second dimension of the 2D|3D laplacian\n"
          "  <dim3> (optionnal): size of the third dimension of the 3D laplacian\n"
          "  Example:\n"
          "  -9 z:10:20\n"
          "  generate a 2D complex laplacian matrix of size 200.\n"
          "\n"
  );
}

static inline void
get_size_from_name( const char    *filename,
                    pastix_csc_t  *csc,
                    long          *dim1,
                    long          *dim2,
                    long          *dim3 )
{
  char *end1 = NULL;
  char *end2 = NULL;
  char *end3 = NULL;

  if(filename[0] == 's' || filename[0] == 'S')
  {
    csc->flttype=PastixFloat;
  }
  else if(filename[0] == 'd' || filename[0] == 'D')
  {
    csc->flttype=PastixDouble;
  }
  else if(filename[0] == 'c' || filename[0] == 'C')
  {
    csc->flttype=PastixComplex32;
  }
  else if(filename[0] == 'z' || filename[0] == 'Z')
  {
    csc->flttype=PastixComplex64;
  }
  else{
    /* TODO error message */
    laplacian_usage(); exit(0);
  }
  if(filename[1] != ':')
  {
    /* TODO error message */
    laplacian_usage(); exit(0);
  }
  *dim1 = strtol( &filename[2], &end1, 10 );
  if(*dim1 == 0)
  {
    /* TODO error message */
    laplacian_usage(); exit(0);
  }
  if(*end1 != ':')
  {
    fprintf(stderr, "1D laplacian, dim1 = %ld\n",*dim1);
    csc->gN = *dim1;
    return;
  }
  *dim2 = strtol( end1+1, &end2, 10 );
  if(*dim2 == 0)
  {
    /* TODO error message */
    laplacian_usage(); exit(0);
  }
  if(*end2 != ':')
  {
    fprintf(stderr, "2D laplacian, dim1 = %ld, dim2 = %ld\n",*dim1,*dim2);
    csc->gN = *dim1*(*dim2);
    return;
  }
  *dim3 = strtol( end2+1, &end3, 10 );
  if(*dim3 == 0)
  {
    /* TODO error message */
    laplacian_usage(); exit(0);
  }
  if(*end3 != '\0')
  {
    /* TODO error message */
    laplacian_usage(); exit(0);
  }
  fprintf(stderr, "3D laplacian, dim1 = %ld, dim2 = %ld, dim2 = %ld\n",*dim1,*dim2,*dim3);
  csc->gN = *dim1*(*dim2)*(*dim3);
}

static inline void
gen1Dlaplacian( const char    *filename,
                pastix_csc_t  *csc,
                void         **rhs,
                long           dim1 )
{

  pastix_int_t i;
  pastix_int_t j;
  pastix_complex64_t * values;
  pastix_complex64_t * rhs_temp;
  pastix_int_t nnzeros = 2*(csc->gN) - 1;
  csc->n   = csc->gN;
  csc->colptr      = NULL;
  csc->rows      = NULL;
  csc->avals   = NULL;
  *rhs     = NULL;
//   *type    = NULL;
//   *rhstype = NULL;
    fprintf(stderr, "Laplacien, n = %d\n",csc->n);

  /* Allocating */
  if ((NULL == (csc->colptr = (pastix_int_t *)  malloc(((csc->n)+1)     *sizeof(pastix_int_t))  )) ||
      (NULL == (csc->rows   = (pastix_int_t *)  malloc((nnzeros)*sizeof(pastix_int_t))  )) ||
      (NULL == (csc->avals  = (pastix_complex64_t *)malloc((nnzeros)*sizeof(pastix_complex64_t)))) ||
      (NULL == (values      = (pastix_complex64_t *)malloc((nnzeros)*sizeof(pastix_complex64_t)))) ||
      (NULL == (*rhs        = (pastix_complex64_t *)malloc((csc->n)         *sizeof(pastix_complex64_t)))) ||
      (NULL == (rhs_temp    = (pastix_complex64_t *)malloc((csc->n)         *sizeof(pastix_complex64_t))))/* ||
      (NULL == (*type       = (char *)          malloc(4         *sizeof(char))      )) ||
      (NULL == (*rhstype    = (char *)          malloc(4         *sizeof(char))          ))*/)
  {
    fprintf(stderr, "Error in CSC allocation\n");
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
  for (i = 0; i < (csc->n); i++)
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
    values[j] = 2;
    j++;
    if (i != (csc->n)-1)
    {
      (csc->rows)[j]    = i+2;
#ifdef TYPE_COMPLEX
      values[j] = - 1 +  2* _Complex_I;
#else
      values[j] = -1;
#endif
      j++;
    }
    (rhs_temp)[i] = 0;
#if (defined TYPE_COMPLEX && defined SYMMETRIC_LAPLACIAN)
    (rhs_temp)[i] = -4 + 4 * _Complex_I;
#endif
  }
  (csc->colptr)[i] = j+1;
#ifdef TYPE_COMPLEX
  (rhs_temp)[0] = 3 - _Complex_I;
#ifdef SYMMETRIC_LAPLACIAN
  (rhs_temp)[0] =  -1 + 3 *  _Complex_I;
#endif
  (rhs_temp)[(csc->n)-1] =  -1 + 3 *  _Complex_I;
#else
  (rhs_temp)[0] = 1;
  (rhs_temp)[(csc->n)-1] = 1;
#endif
  /* type and rhstype */
// #if (defined TYPE_COMPLEX && !defined SYMMETRIC_LAPLACIAN)
//   sprintf (*type, "RHA");
// #else
//   sprintf (*type, "RSA");
// #endif
//   sprintf (*rhstype,"???");
  csc->mtxtype = PastixSymmetric;
  csc->fmttype = PastixCSC;
	memcpy(*rhs, rhs_temp, csc->n*sizeof(pastix_complex64_t));
  memFree_null(rhs_temp);
	memcpy(csc->avals, values, nnzeros*sizeof(pastix_complex64_t));
  memFree_null(values);
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * File: laplacian.c
 *
 * genlaplacian - Generate a laplacian of size csc->n
 *
 * Example :
 * >  2 -1  0  0
 * > -1  2 -1  0
 * >  0 -1  2 -1
 * >  0  0 -1  2
 *
 *******************************************************************************
 *
 * @param[in] filename
 *          Path to the directory containing matrix.
 * 
 * @param[in] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
 * 
 * @param[in] rhs
 *          At exit, contains the right hand side member.
 *
 *******************************************************************************/
void
genlaplacian( const char    *filename,
                pastix_csc_t  *csc,
                void         **rhs )
{
  long  dim1 = 0;
  long  dim2 = 0;
  long  dim3 = 0;

  get_size_from_name(filename, csc, &dim1, &dim2, &dim3);
  if(dim3 == 0)
  {
    if(dim2 == 0)
    {
      gen1Dlaplacian(filename, csc, rhs, dim1);
    }else{
//       gen2Dlaplacian(filename, csc, rhs, dim1, dim2)
      gen1Dlaplacian(filename, csc, rhs, dim1);
    }
  }else{
//     gen3Dlaplacian(filename, csc, rhs, dim1, dim2, dim3)
    gen1Dlaplacian(filename, csc, rhs, dim1);
  }
}
