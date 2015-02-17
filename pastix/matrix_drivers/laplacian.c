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
#include "common.h"
#include "drivers.h"
#include "pastix.h"

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
 * @param[in] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
 *
 *******************************************************************************/
void
genlaplacian(pastix_csc_t *csc
                 , void **rhs )
{

  pastix_int_t i;
  pastix_int_t j;
  pastix_complex64_t * values;
  pastix_complex64_t * rhs_temp;

  pastix_int_t nnzeros = 2*(csc->n) - 1;
  csc->gN   = csc->n;
  csc->colptr      = NULL;
  csc->rows      = NULL;
  csc->avals   = NULL;
  *rhs     = NULL;
//   *type    = NULL;
//   *rhstype = NULL;
    fprintf(stderr, "Laplacien, n=%d\n",csc->n);

  /* Allocating */
  if ((NULL == (csc->colptr = (pastix_int_t *)  malloc(((csc->n)+1)     *sizeof(pastix_int_t))  )) ||
      (NULL == (csc->rows   = (pastix_int_t *)  malloc((nnzeros)*sizeof(pastix_int_t))  )) ||
      (NULL == (csc->avals  = (pastix_complex64_t *)malloc((nnzeros)*sizeof(pastix_complex64_t)))) ||
      (NULL == (values      = (pastix_complex64_t *)malloc((nnzeros)*sizeof(pastix_complex64_t)))) ||
      (NULL == (*rhs        = (pastix_complex64_t *)malloc((csc->n)         *sizeof(pastix_complex64_t)))) ||
      (NULL == (rhs_temp   = (pastix_complex64_t *)malloc((csc->n)         *sizeof(pastix_complex64_t))))/* ||
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
	memcpy(*rhs, rhs_temp, csc->n*sizeof(pastix_complex64_t));
  memFree_null(rhs_temp);
	memcpy(csc->avals, values, csc->dof*sizeof(pastix_complex64_t));
  memFree_null(values);
}
