/* 
   File: cscdread.h
   
   Reads file in cscd format.

 */
#ifndef CSCDREAD_H
#define CSCDREAD_H

/*
  Function: cscdRead

  Reads a matrix in cscd format
  Nrow is equal to Ncol.

  File format is like :
  ...

  
  Parameters:
    dirname     - Path to the directory containing matrix
    colptr      - Index of first element of each column in *row* and *val*
    row         - Row of eah element
    loc2glb     - Correspondance between local and global numbering
    avals       - Value of each element
    rhs         - Right Hand Side
    colnbr      - Number of columns
    nnz         - Number of non-zeros
    pastix_comm - MPI communicator
 */
void  
cscdRead(char const      *dirname, 
	 pastix_int_t   **colptr, 
	 pastix_int_t   **row, 
	 pastix_int_t   **loc2glb, 
	 pastix_float_t **avals, 
	 pastix_float_t **rhs, 
	 pastix_int_t    *colnbr, 
	 pastix_int_t    *nnz, 
	 MPI_Comm         pastix_comm);


#endif /* CSCDREAD_H */
