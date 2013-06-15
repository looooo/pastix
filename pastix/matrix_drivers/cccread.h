#ifndef _CCCREAD_H_
#define _CCCREAD_H_

/*
   File: cccread.h

   Reads file in ccc format.

 */

/* set to extract only symmetric part in CCC format */
#define SYMPART

/*
  Function: cccReadHeader

  Reads header from a file in ccc matrix format.
  Nrow is equal to Ncol.
  Type is "CSA" if SYMPART is defined (default) "CUA" otherwise.

  File format is like :
  > Ncol Nnzero


  Parameters:
    infile - File to read from.
    Nrow   - Number of rows
    Ncol   - Number of columns
    Nnzero - Number of non zeros
    Type   - Type of the matrix
 */
void cccReadHeader(FILE         *infile,
                   pastix_int_t *Nrow,
                   pastix_int_t *Ncol,
                   pastix_int_t *Nnzero,
                   char         *Type);


/*
   Function: cccRead

   Reads Matrix in ccc format.

   Header format is described in <cccReadHeader>,
   "filename"/hfile contains columns
   Enf of the matrix is in three files in CSC format :
   "filename"/ifile contains Ncol columns,
   "filename"/jfile contains Nnzeros rows,
   "filename"/afile contains Nnzeros values.

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
void cccRead(char const      *filename,
             pastix_int_t    *Nrow,
             pastix_int_t    *Ncol,
             pastix_int_t    *Nnzero,
             pastix_int_t   **col,
             pastix_int_t   **row,
             pastix_float_t **val,
             char           **Type,
             char           **RhsType);

#endif
