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
  fprintf(stderr, "3D laplacian, dim1 = %ld, dim2 = %ld, dim3 = %ld\n",*dim1,*dim2,*dim3);
  csc->gN = *dim1*(*dim2)*(*dim3);
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * File: laplacian.c
 *
 * z_gen1Dlaplacian - Generate a 1D laplacian
 *
 * Example :
 * >  2 -1  0  0
 * > -1  2 -1  0
 * >  0 -1  2 -1
 * >  0  0 -1  2
 *
 *******************************************************************************
 * 
 * @param[out] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
 * 
 * @param[out] rhs
 *          At exit, contains the right hand side member.
 * 
 * @param[in] dim1
 *          contains the dimension of the 1D laplacian.
 *
 *******************************************************************************/
static inline void
z_gen1Dlaplacian( pastix_csc_t  *csc,
                  void         **rhs,
                  long           dim1 )
{

  pastix_int_t i;
  pastix_int_t j;
  pastix_complex64_t * valptr;
  pastix_complex64_t * rhsptr;
  pastix_int_t nnzeros = 2*(csc->gN) - 1;
  csc->n       = csc->gN;
  csc->colptr  = NULL;
  csc->rows    = NULL;
  csc->avals   = NULL;
  
  fprintf(stderr, "Laplacien 1D, n = %d\n",csc->n);
  
  assert( csc->gN == dim1 );

  /* Allocating */
  if ((NULL == (csc->colptr = (pastix_int_t *)      malloc(((csc->n)+1) *sizeof(pastix_int_t))  ))     ||
      (NULL == (csc->rows   = (pastix_int_t *)      malloc((nnzeros)    *sizeof(pastix_int_t))  ))     ||
      (NULL == (csc->avals  = (pastix_complex64_t *)malloc((nnzeros)    *sizeof(pastix_complex64_t)))) ||
      (NULL == (*rhs        = (pastix_complex64_t *)malloc((csc->n)     *sizeof(pastix_complex64_t))))/* ||
      (NULL == (*type       = (char *)          malloc(4         *sizeof(char))      )) ||
      (NULL == (*rhstype    = (char *)          malloc(4         *sizeof(char))          ))*/)
  {
    fprintf(stderr, "Laplacien 1D, error in CSC allocation\n");
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
  valptr = csc->avals;
  rhsptr = *rhs;
  for (i = 0; i < (csc->gN); i++,rhsptr++)
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
    *valptr = 2;
    j++;
    valptr++;
    if (i != (csc->gN)-1)
    {
      (csc->rows)[j]    = i+2;
      
      if (csc->flttype==PastixComplex64 || csc->flttype==PastixComplex32)
      {
        *valptr = - 1 +  2* _Complex_I;
      }else{
        *valptr = -1;
      }
      j++;
      valptr++;
    }
    *rhsptr = 0;
    if (csc->flttype==PastixComplex64 || csc->flttype==PastixComplex32)
    {
      *rhsptr = -4 + 4 * _Complex_I;
    }
  }
  (csc->colptr)[i] = j+1;
  if (csc->flttype==PastixComplex64 || csc->flttype==PastixComplex32)
  {
    rhsptr = *rhs;
    *rhsptr = -1 + 3 *  _Complex_I;
    rhsptr += csc->gN -1;
    *rhsptr = -1 + 3 *  _Complex_I;
  }else{
    rhsptr = *rhs;
    *rhsptr = 1;
    rhsptr += csc->gN -1;
    *rhsptr = 1;
  }
  /* type and rhstype */
// #if (defined TYPE_COMPLEX && !defined SYMMETRIC_LAPLACIAN)
//   sprintf (*type, "RHA");
// #else
//   sprintf (*type, "RSA");
// #endif
//   sprintf (*rhstype,"???");
  csc->mtxtype = PastixSymmetric;
  csc->fmttype = PastixCSC;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * File: laplacian.c
 *
 * z_gen2Dlaplacian - Generate a 2D laplacian
 *
 * Example :
 * >  4 -1  0 -1  0  0
 * > -1  4 -1  0 -1  0
 * >  0 -1  4  0  0 -1
 * > -1  0  0  4 -1  0
 * >  0 -1  0 -1  4 -1
 * >  0  0 -1  0 -1  4
 *
 *******************************************************************************
 * 
 * @param[out] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
 * 
 * @param[out] rhs
 *          At exit, contains the right hand side member.
 * 
 * @param[in] dim1
 *          contains the first dimension of the 2D grid of the laplacian.
 * 
 * @param[in] dim2
 *          contains the second dimension of the 2D grid of the laplacian.
 *
 *******************************************************************************/
static inline void
z_gen2Dlaplacian( pastix_csc_t  *csc,
                  void         **rhs,
                  long           dim1,
                  long           dim2 )
{

  pastix_int_t i;
  pastix_int_t j;
  pastix_int_t k;
  pastix_complex64_t * valptr;
  pastix_complex64_t * rhsptr;
  pastix_int_t nnzeros = (2*(dim1)-1)*dim2 + (dim2-1)*dim1;
  csc->n       = csc->gN;
  csc->colptr  = NULL;
  csc->rows    = NULL;
  csc->avals   = NULL;
  *rhs         = NULL;
  
  fprintf(stderr, "Laplacien 2D, n = %d\n",csc->n);
  
  assert( csc->gN == dim1*dim2 );

  /* Allocating */
  if ((NULL == (csc->colptr = (pastix_int_t *)      malloc(((csc->gN)+1)*sizeof(pastix_int_t))  ))     ||
      (NULL == (csc->rows   = (pastix_int_t *)      malloc((nnzeros)    *sizeof(pastix_int_t))  ))     ||
      (NULL == (csc->avals  = (pastix_complex64_t *)malloc((nnzeros)    *sizeof(pastix_complex64_t)))) ||
      (NULL == (*rhs        = (pastix_complex64_t *)calloc((csc->gN)    ,sizeof(pastix_complex64_t)))) )
  {
    fprintf(stderr, "Laplacien 2D, error in CSC allocation\n");
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
  
  csc->colptr[0] = 1;
  valptr = csc->avals;
  rhsptr = *rhs;
  *rhsptr = 1.;
  *(rhsptr+csc->gN-1) = 1.;
  k = 0;
  
  for(i=0; i<dim2; i++)
  {
    for(j=1; j<=dim1; j++)
    {
      // column k = i*dim1+j of the matrix
      k += 1;
      if(j!=dim1 && i!=dim2-1)
      {
        csc->colptr[k] = csc->colptr[k-1]+3;
        csc->rows[csc->colptr[k-1]-1]=k;
        csc->rows[csc->colptr[k-1]  ]=k+1;
        csc->rows[csc->colptr[k-1]+1]=k+dim1;
        *valptr = 4.;
        *(valptr+1) = -1.;
        *(valptr+2) = -1.;
        valptr += 3;
      }
      else if(j==dim1 && i!=dim2-1)
      {
        csc->colptr[k] = csc->colptr[k-1]+2;
        csc->rows[csc->colptr[k-1]-1]=k;
        csc->rows[csc->colptr[k-1]  ]=k+dim1;
        *valptr = 4.;
        *(valptr+1) = -1.;
        valptr += 2;
      }
      else if(j!=dim1 && i==dim2-1)
      {
        csc->colptr[k] = csc->colptr[k-1]+2;
        csc->rows[csc->colptr[k-1]-1]=k;
        csc->rows[csc->colptr[k-1]  ]=k+1;
        *valptr = 4.;
        *(valptr+1) = -1;
        valptr += 2;
      }
      else /* if(j==dim1 && i==dim2-1) */
      {
        csc->colptr[k] = csc->colptr[k-1]+1;
        csc->rows[csc->colptr[k-1]-1]=k;
        *valptr = 4.;
        valptr += 1;
      }
    }
  }
  
  csc->mtxtype = PastixSymmetric;
  csc->fmttype = PastixCSC;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc_driver
 *
 * File: laplacian.c
 *
 * z_gen3Dlaplacian - Generate a 3D laplacian
 *
 * Example :
 * >  6 -1 -1  0 -1  0  0  0
 * > -1  6  0 -1  0 -1  0  0
 * > -1  0  6 -1  0  0 -1  0
 * >  0 -1 -1  6  0  0  0 -1
 * > -1  0  0  0  6 -1 -1  0
 * >  0 -1  0  0 -1  6  0 -1
 * >  0  0 -1  0 -1  0  6 -1
 * >  0  0  0 -1  0 -1 -1  6
 *
 *******************************************************************************
 * 
 * @param[out] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
 * 
 * @param[out] rhs
 *          At exit, contains the right hand side member.
 * 
 * @param[in] dim1
 *          contains the first dimension of the 3D grid of the laplacian.
 * 
 * @param[in] dim2
 *          contains the second dimension of the 3D grid of the laplacian.
 * 
 * @param[in] dim3
 *          contains the third dimension of the 3D grid of the laplacian.
 *
 *******************************************************************************/
static inline void
z_gen3Dlaplacian( pastix_csc_t  *csc,
                  void         **rhs,
                  long           dim1,
                  long           dim2,
                  long           dim3 )
{

  pastix_int_t i;
  pastix_int_t j;
  pastix_int_t k;
  pastix_int_t l;
  pastix_complex64_t * valptr;
  pastix_complex64_t * rhsptr;
  pastix_int_t nnzeros = (2*(dim1)-1)*dim2*dim3 + (dim2-1)*dim1*dim3 + dim2*dim1*(dim3-1);
  csc->n       = csc->gN;
  csc->colptr  = NULL;
  csc->rows    = NULL;
  csc->avals   = NULL;
  *rhs         = NULL;
  
  assert( csc->gN == dim1*dim2 );

  /* Allocating */
  if ((NULL == (csc->colptr = (pastix_int_t *)      malloc(((csc->gN)+1)*sizeof(pastix_int_t))  ))     ||
      (NULL == (csc->rows   = (pastix_int_t *)      malloc((nnzeros)    *sizeof(pastix_int_t))  ))     ||
      (NULL == (csc->avals  = (pastix_complex64_t *)malloc((nnzeros)    *sizeof(pastix_complex64_t)))) ||
      (NULL == (*rhs        = (pastix_complex64_t *)calloc((csc->gN)    ,sizeof(pastix_complex64_t)))) )
  {
    fprintf(stderr, "Laplacien 3D, error in CSC allocation\n");
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
  
  csc->colptr[0] = 1;
  valptr = csc->avals;
  rhsptr = *rhs;
  *rhsptr = 1.;
  *(rhsptr+csc->gN-1) = 1.;
  l = 0;
  
  for(i=0; i<dim3; i++)
  {
    for(j=0; j<dim2; j++)
    {
      for(k=1; k<=dim1; k++)
      {
        // column l = i*dim1*dim2+j*dim1+k of the matrix
        l += 1;
        if(k!=dim1 && j!=dim2-1 && i!=dim3-1)
        {
          csc->colptr[l] = csc->colptr[l-1]+4;
          csc->rows[csc->colptr[l-1]-1]=l;
          csc->rows[csc->colptr[l-1]  ]=l+1;
          csc->rows[csc->colptr[l-1]+1]=l+dim1;
          csc->rows[csc->colptr[l-1]+2]=l+dim1*dim2;
          *valptr = 6.;
          *(valptr+1) = -1.;
          *(valptr+2) = -1.;
          *(valptr+3) = -1.;
          valptr += 4;
        }
        else if(k==dim1 && j!=dim2-1 && i!=dim3-1)
        {
          csc->colptr[l] = csc->colptr[l-1]+3;
          csc->rows[csc->colptr[l-1]-1]=l;
          csc->rows[csc->colptr[l-1]  ]=l+dim1;
          csc->rows[csc->colptr[l-1]+1]=l+dim1*dim2;
          *valptr = 6.;
          *(valptr+1) = -1.;
          *(valptr+2) = -1.;
          valptr += 3;
        }
        else if(k!=dim1 && j==dim2-1 && i!=dim3-1)
        {
          csc->colptr[l] = csc->colptr[l-1]+3;
          csc->rows[csc->colptr[l-1]-1]=l;
          csc->rows[csc->colptr[l-1]  ]=l+1;
          csc->rows[csc->colptr[l-1]+1]=l+dim1*dim2;
          *valptr = 6.;
          *(valptr+1) = -1.;
          *(valptr+2) = -1.;
          valptr += 3;
        }
        else if(k!=dim1 && j!=dim2-1 && i==dim3-1)
        {
          csc->colptr[l] = csc->colptr[l-1]+3;
          csc->rows[csc->colptr[l-1]-1]=l;
          csc->rows[csc->colptr[l-1]  ]=l+1;
          csc->rows[csc->colptr[l-1]+1]=l+dim1;
          *valptr = 6.;
          *(valptr+1) = -1.;
          *(valptr+2) = -1.;
          valptr += 3;
        }
        else if(k==dim1 && j==dim2-1 && i!=dim3-1)
        {
          csc->colptr[l] = csc->colptr[l-1]+2;
          csc->rows[csc->colptr[l-1]-1]=l;
          csc->rows[csc->colptr[l-1]  ]=l+dim1*dim2;
          *valptr = 6.;
          *(valptr+1) = -1.;
          valptr += 2;
        }
        else if(k!=dim1 && j==dim2-1 && i==dim3-1)
        {
          csc->colptr[l] = csc->colptr[l-1]+2;
          csc->rows[csc->colptr[l-1]-1]=l;
          csc->rows[csc->colptr[l-1]  ]=l+1;
          *valptr = 6.;
          *(valptr+1) = -1.;
          valptr += 2;
        }
        else if(k==dim1 && j!=dim2-1 && i==dim3-1)
        {
          csc->colptr[l] = csc->colptr[l-1]+2;
          csc->rows[csc->colptr[l-1]-1]=l;
          csc->rows[csc->colptr[l-1]  ]=l+dim1;
          *valptr = 6.;
          *(valptr+1) = -1.;
          valptr += 2;
        }
        else if(k==dim1 && j==dim2-1 && i==dim3-1)
        {
          csc->colptr[l] = csc->colptr[l-1]+1;
          csc->rows[csc->colptr[l-1]-1]=l;
          *valptr = 6.;
          valptr += 1;
        }
      }
    }
  }
  
  csc->mtxtype = PastixSymmetric;
  csc->fmttype = PastixCSC;
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
 *******************************************************************************
 *
 * @param[in] filename
 *          Path to the directory containing matrix.
 * 
 * @param[out] csc
 *          At start, contains the size of the laplacian in csc->n.
 *          At exit, contains the matrix in csc format.
 * 
 * @param[out] rhs
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
      z_gen1Dlaplacian(csc, rhs, dim1);
    }else{
      z_gen2Dlaplacian(csc, rhs, dim1, dim2);
    }
  }else{
    z_gen3Dlaplacian(csc, rhs, dim1, dim2, dim3);
  }
}
