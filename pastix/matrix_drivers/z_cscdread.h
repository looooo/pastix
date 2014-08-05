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
   File: z_cscdread.h
   
   Reads file in cscd format.

 */
#ifndef Z_CSCDREAD_H
#define Z_CSCDREAD_H

/*
  Function: z_cscdRead

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
z_cscdRead(char const      *dirname, 
	 pastix_int_t   **colptr, 
	 pastix_int_t   **row, 
	 pastix_int_t   **loc2glb, 
	 pastix_complex64_t **avals, 
	 pastix_complex64_t **rhs, 
	 pastix_int_t    *colnbr, 
	 pastix_int_t    *nnz, 
	 MPI_Comm         pastix_comm);


#endif /* Z_CSCDREAD_H */
