/*
  File: tools.h

  headers for some tools used in sopalin.

 */
#ifndef TOOLS_H
#define TOOLS_H

/* 
   Function: dim_dgeam


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
void         dim_dgeam(char *transa, char *transb, PASTIX_INT m, PASTIX_INT n, PASTIX_FLOAT alpha,
		       PASTIX_FLOAT *a, PASTIX_INT lda, PASTIX_FLOAT *b, PASTIX_INT ldb);

/* 
   Function: GetMpiType

   Construct a MPI type to store complex values.

   Parameters:
     none

   Return:
   
   the constructed MPI type for complex floating values.
 */
MPI_Datatype GetMpiType(void);

/* 
   Function: FreeMpiType

   Free the MPI type to store complex values.

   Parameters:
     none

 */
void FreeMpiType(void);

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
void         mysum(void *in, void *inout, int *len, MPI_Datatype *dptr);
/*
  Function: GetMpiSum
  
  Creates a new MPI sum operation.

  Parameters:
    none
    
  Return: the new operation created.
  
 */
MPI_Op       GetMpiSum(void);

/*
  Function: FreeMpiSum
  
  Free the MPI sum operation.

  Parameters:
    none
    
 */
void FreeMpiSum(void);

#endif
