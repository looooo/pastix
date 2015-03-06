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
  File: olafread.h
  
  Driver for the olaf matrix format.
 */

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
		    char         *Type);

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
	      pastix_complex64_t **rhs);
