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
  File: tools.c

  Some tools used in sopalin.

 */
#include "common.h"
#include "z_sopalin_compute.h"
#include "z_tools.h"

/*
   Function: z_dim_dgeam


   Computes b = alpha * a + b.

     - if transa and transb equals 'N'.
        > b := alpha * a + b
     - if transa = 'T' and transb ='N'.
        > b := alpha * trans(a) + b
     - if transa = 'N' and transb ='T'.
        > trans(b) := alpha * a + trans(b)
     - if transa = 'T' and transb ='T'.
        > trans(b) := alpha * trans(a) + trans(b)

   Parameters:
     transa - indicates if a needs to be transposed.
     transb - indicates if b needs to be transposed.
     m      - number of row in a and b.
     n      - number of colonnes in a and b.
     alpha  - scalar.
     a      - Matrix a.
     lda    - Stride between 2 columns of a.
     b      - Matrix b.
     ldb    - Stride between 2 columns of b.
*/
void z_dim_dgeam(char *transa,char *transb,pastix_int_t m,pastix_int_t n,pastix_complex64_t alpha,pastix_complex64_t *a,
               pastix_int_t lda,pastix_complex64_t *b,pastix_int_t ldb)
{
  pastix_int_t i;

  if (*transa=='N')
    {
      if (*transb=='N')
        for (i=0;i<n;i++)
          {
#ifdef TYPE_COMPLEX
          {SOPALIN_AXPY(m,alpha,(pastix_complex64_t *) &(a[i*lda]),1, (pastix_complex64_t *) &(b[i*ldb]),1);}
#else
          {SOPALIN_AXPY(m,alpha,&(a[i*lda]),1,&(b[i*ldb]),1);}
#endif
          }
      else
        for (i=0;i<n;i++)
          {
#ifdef TYPE_COMPLEX
          {SOPALIN_AXPY(m,alpha,(pastix_complex64_t *) &(a[i*lda]),1, (pastix_complex64_t *) &(b[i]),ldb);}
#else
          {SOPALIN_AXPY(m,alpha,&(a[i*lda]),1,&(b[i]),ldb);}
#endif
          }
    }
  else
    {
      if (*transb=='N')
        for (i=0;i<n;i++)
          {
#ifdef TYPE_COMPLEX
          {SOPALIN_AXPY(m,alpha,(pastix_complex64_t *) &(a[i]),lda, (pastix_complex64_t *) &(b[i*ldb]),1);}
#else
          {SOPALIN_AXPY(m,alpha,&(a[i]),lda,&(b[i*ldb]),1);}
#endif
          }
      else
        for (i=0;i<n;i++)
          {
#ifdef TYPE_COMPLEX
          {SOPALIN_AXPY(m,alpha,(pastix_complex64_t *) &(a[i]),lda, (pastix_complex64_t *) &(b[i]),ldb);}
#else
          {SOPALIN_AXPY(m,alpha,&(a[i]),lda,&(b[i]),ldb);}
#endif
          }
    }
}

/*
   Function: z_GetMpiType

   Construct a MPI type to store complex values.

   Parameters:
     none

   Return:

   the constructed MPI type for complex floating values.
 */
MPI_Datatype z_GetMpiType(void)
{
  static int first=1;
  static MPI_Datatype floattype;

  if (first==1)
    {
#ifdef FORCE_NOMPI
#ifdef PREC_DOUBLE
      floattype=MPI_DOUBLE_COMPLEX;
#else
      floattype=MPI_COMPLEX;
#endif
#else
#ifdef PREC_DOUBLE
#ifdef MPI_DOUBLE_COMPLEX
      floattype=MPI_DOUBLE_COMPLEX;
#else
      MPI_Type_contiguous(2,MPI_DOUBLE,&floattype);
      MPI_Type_commit(&floattype);
#endif
#else
#ifdef MPI_COMPLEX
      floattype=MPI_COMPLEX;
#else
      MPI_Type_contiguous(2,MPI_FLOAT,&floattype);
      MPI_Type_commit(&floattype);
#endif
#endif
#endif
      first=0;
    }

  return floattype;
}

/*
   Function: z_FreeMpiType

   Free the MPI type to store complex values.

   Parameters:
     none

 */
void z_FreeMpiType(void)
{

#ifndef FORCE_NOMPI
#if ((defined PREC_DOUBLE && (!defined MPI_DOUBLE_COMPLEX)) || (!defined MPI_COMPLEX))
  MPI_Datatype floattype = z_GetMpiType();
  MPI_Type_free(&(floattype));
#endif
#endif
}

/*
  Function: mysum

  Computes the sum of *in* and *inout* vectors
  and stores it in *inout*.


  Parameters:
    in    - First vector.
    inout - Second vector wich will store the sum.
    len   - Size of each vector.
    dptr  - MPI datatype
 */
void z_mysum(void *in, void *inout, int *len, MPI_Datatype *dptr)
{
  pastix_int_t i;
  pastix_complex64_t *a = (pastix_complex64_t *) in;
  pastix_complex64_t *b = (pastix_complex64_t *) inout;
  (void)dptr;

  for (i=0; i<*len; ++i)
    {
      b[i] += a[i];
    }
}

/*
  Function: z_GetMpiSum

  Creates a new MPI sum operation.

  Parameters:
    none

  Return: the new operation created.

 */
MPI_Op z_GetMpiSum(void)
{
  static MPI_Op myop;
#if ((defined PREC_DOUBLE && (!defined MPI_DOUBLE_COMPLEX)) || (!defined MPI_COMPLEX))
  static int first=1;
  if (first==1)
    {
      MPI_Op_create(&mysum,1,&myop);
      first=0;
    }
#else
  myop = MPI_SUM;
#endif
  return myop;
}

/*
  Function: z_FreeMpiSum

  Free the MPI sum operation.

  Parameters:
    none

 */
void z_FreeMpiSum(void)
{
#if ((defined PREC_DOUBLE && (!defined MPI_DOUBLE_COMPLEX)) || (!defined MPI_COMPLEX))
  MPI_Op myop = z_GetMpiSum();
  MPI_Op_free(&(myop));
#endif
}
