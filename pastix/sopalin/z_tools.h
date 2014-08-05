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
  File: z_tools.h

  headers for some tools used in sopalin.

 */
#ifndef Z_TOOLS_H
#define Z_TOOLS_H

/* 
   Function: z_dim_dgeam


   Computes b = alph * a + b.
 
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
void         z_dim_dgeam(char *transa, char *transb, pastix_int_t m, pastix_int_t n, pastix_complex64_t alpha,
		       pastix_complex64_t *a, pastix_int_t lda, pastix_complex64_t *b, pastix_int_t ldb);

/* 
   Function: z_GetMpiType

   Construct a MPI type to store complex values.

   Parameters:
     none

   Return:
   
   the constructed MPI type for complex floating values.
 */
MPI_Datatype z_GetMpiType(void);

/* 
   Function: z_FreeMpiType

   Free the MPI type to store complex values.

   Parameters:
     none

 */
void z_FreeMpiType(void);

/*
  Function: z_mysum
  
  Computes the sum of *in* and *inout* vectors 
  and stores it in *inout*.
  
  
  Parameters: 
    in    - First vector.
    inout - Second vector wich will store the sum.
    len   - Size of each vector.
    dptr  - MPI datatype
 */
void         z_mysum(void *in, void *inout, int *len, MPI_Datatype *dptr);
/*
  Function: z_GetMpiSum
  
  Creates a new MPI sum operation.

  Parameters:
    none
    
  Return: the new operation created.
  
 */
MPI_Op       z_GetMpiSum(void);

/*
  Function: z_FreeMpiSum
  
  Free the MPI sum operation.

  Parameters:
    none
    
 */
void z_FreeMpiSum(void);

#endif
