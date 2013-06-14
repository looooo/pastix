/*
  File: olafread.h
  
  Driver for the olaf matrix format.
 */

/*
  Function: olafReadHeader
  
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
void olafReadHeader(FILE         *infile,
		    pastix_int_t *Nrow, 
		    pastix_int_t *Ncol, 
		    pastix_int_t *Nnzero, 
		    char         *Type);

/*
  Function: olafRead

  Reads a matrix in olaf format.

  Header format is described in <olafReadHeader>, 
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
void olafRead(char const      *filename, 
	      pastix_int_t    *Nrow, 
	      pastix_int_t    *Ncol, 
	      pastix_int_t    *Nnzero, 
	      pastix_int_t   **col, 
	      pastix_int_t   **row, 
	      pastix_float_t **val, 
	      char           **Type, 
	      char           **RhsType, 
	      pastix_float_t **rhs);
