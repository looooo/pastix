/* File: peerread.h

   Reads a matrix in PEER format.
 */

/*
  Function: peerRead
  
  Reads a matrix in PEER format.
  
  first file contain :
  > NumberOfFiles
  > file1
  > file2
  > ...

  each file contains:
  > %ld%ld          (globaln localn)
  > %ld             (local nnzeros)
  four elements from col by line, nnzeros local elements in total,
  four elements from row by line, nnzeros local elements in total,
  four elements from val by line, nnzeros local elements in total,
  four elements from rhs by line, nnzeros local elements in total,

  for each part, last line can be with 1, 2 ,3 or 4 elements.

  Parameters:
    filename - Path to file to read from
    Nrow     - Number of rows						 
    Ncol     - Number of columns					 
    Nnzero,  - Number of non zeros					 
    col      - Index of first element of each column in *row* and *val* 
    row      - Row of eah element				       	 
    val      - Value of each element				       	 
    Type     - Type of the matrix				       	 
    RhsType  - Type of the right-hand-side.			         
    rhs      - right-hand-side term(s)
*/
void peerRead(char const      *filename, 
	      pastix_int_t    *Nrow, 
	      pastix_int_t    *Ncol, 
	      pastix_int_t    *Nnzero, 
	      pastix_int_t   **col, 
	      pastix_int_t   **row, 
	      pastix_float_t **val, 
	      char           **Type, 
	      char           **RhsType, 
	      pastix_float_t **rhs);


/*
  Function: peerRead2
  
  Reads a matrix in PEER format.
  
  first file contain :
  > NumberOfFiles
  > file1
  > file2
  > ...

  each file contains:
  > %ld %ld %ld  (globaln localn  localnnzeros)
  six elements from col by line, localnnzeros elements in total,
  six elements from row by line, localnnzeros elements in total,
  six elements from val by line, localnnzeros elements in total,
  six elements from rhs by line, localnnzeros elements in total,

  for each part, last line can be with 1, 2 ,3 , 4, 5 or 6 elements.

  Parameters:
    filename - Path to file to read from
    Nrow     - Number of rows						 
    Ncol     - Number of columns					 
    Nnzero,  - Number of non zeros					 
    col      - Index of first element of each column in *row* and *val* 
    row      - Row of eah element				       	 
    val      - Value of each element				       	 
    Type     - Type of the matrix				       	 
    RhsType  - Type of the right-hand-side.			         
    rhs      - right-hand-side term(s)
*/
void peerRead2(char const      *filename, 
	       pastix_int_t    *Nrow, 
	       pastix_int_t    *Ncol, 
	       pastix_int_t    *Nnzero, 
	       pastix_int_t   **col, 
	       pastix_int_t   **row, 
	       pastix_float_t **val, 
	       char           **Type, 
	       char           **RhsType, 
	       pastix_float_t **rhs);
