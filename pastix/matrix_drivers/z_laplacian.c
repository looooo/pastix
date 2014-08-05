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
 * File: laplacian.c
 *
 * Generate a laplacian
 *
 * Example :
 * >  2 -1  0  0
 * > -1  2 -1  0
 * >  0 -1  2 -1
 * >  0  0 -1  2
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <stdint.h>

#include "z_pastix.h"
#include "common_drivers.h"
#include "z_laplacian.h"
/*
 * Function: genlaplacien
 *
 * Generate a laplacien of size *n*
 *
 * Parameters:
 *   n       - Size of the wanted matrix
 *   nnzeros - Number of non zeros in the produced matrice
 *   ia      - Index of first element of each column in *row* and *val*
 *   ja      - Row of eah element
 *   avals   - Value of each element
 *   rhs     - Right-hand-side member
 *   type    - Type of the matrix
 *   rhstype - Type of the right hand side.
 */
int z_genlaplacian(pastix_int_t     n,
                 pastix_int_t    *nnzeros,
                 pastix_int_t   **ia,
                 pastix_int_t   **ja,
                 pastix_complex64_t **avals,
                 pastix_complex64_t **rhs,
                 char           **type,
                 char           **rhstype)
{

  pastix_int_t i;
  pastix_int_t j;

  *nnzeros = 2*n - 1;
  *ia      = NULL;
  *ja      = NULL;
  *avals   = NULL;
  *rhs     = NULL;
  *type    = NULL;
  *rhstype = NULL;

  /* Allocating */
  if ((NULL == (*ia       = (pastix_int_t *)  malloc((n+1)     *sizeof(pastix_int_t))  )) ||
      (NULL == (*ja       = (pastix_int_t *)  malloc((*nnzeros)*sizeof(pastix_int_t))  )) ||
      (NULL == (*avals    = (pastix_complex64_t *)malloc((*nnzeros)*sizeof(pastix_complex64_t)))) ||
      (NULL == (*rhs      = (pastix_complex64_t *)malloc(n         *sizeof(pastix_complex64_t)))) ||
      (NULL == (*type     = (char *)          malloc(4         *sizeof(char))      )) ||
      (NULL == (*rhstype  = (char *)          malloc(4         *sizeof(char))          )))
  {
    fprintf(stderr, "Error in CSC allocation\n");
    if (*type != NULL)
    {
      free(*type);
      *type = NULL;
    }
    if (*rhs != NULL)
    {
      free(*rhs);
      *rhs = NULL;
    }
    if (*avals != NULL)
    {
      free(*avals);
      *avals = NULL;
    }
    if (*ja != NULL)
    {
      free(*ja);
      *ja = NULL;
    }
    if (*ia != NULL)
    {
      free(*ia);
      *ia = NULL;
    }
    return EXIT_FAILURE;
  }

  /* Building ia, ja and avals and rhs*/
  j=0;
  for (i = 0; i < n; i++)
    {
      (*ia)[i] = j+1;
      /* ONLY triangular inferior matrix */
      /*       if (i != 0) */
      /*  { */
      /*    (*ja)[j]    = i; */
      /*    (*avals)[j] = -1; */
      /*    j++; */
      /*  } */
      (*ja)[j]    = i+1;
      (*avals)[j] = 2;
      j++;
      if (i != n-1)
      {
        (*ja)[j]    = i+2;
#ifdef TYPE_COMPLEX
        (*avals)[j] = - 1 +  2* _Complex_I;
#else
        (*avals)[j] = -1;
#endif
        j++;
      }
      (*rhs)[i] = 0;
#if (defined TYPE_COMPLEX && defined SYMMETRIC_LAPLACIAN)
      (*rhs)[i] = -4 + 4 * _Complex_I;
#endif
    }
  (*ia)[i] = j+1;
#ifdef TYPE_COMPLEX
  (*rhs)[0] = 3 - _Complex_I;
#ifdef SYMMETRIC_LAPLACIAN
  (*rhs)[0] =  -1 + 3 *  _Complex_I;
#endif
  (*rhs)[n-1] =  -1 + 3 *  _Complex_I;
#else
  (*rhs)[0] = 1;
  (*rhs)[n-1] = 1;
#endif
  /* type and rhstype */
#if (defined TYPE_COMPLEX && !defined SYMMETRIC_LAPLACIAN)
  sprintf (*type, "RHA");
#else
  sprintf (*type, "RSA");
#endif
  sprintf (*rhstype,"???");
  return EXIT_SUCCESS;
}
