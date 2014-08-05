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
#ifndef Z__CCCREAD_H_
#define Z_CCCREAD_H_

/*
   File: z_cccread.h

   Reads file in ccc format.

 */

/* set to extract only symmetric part in CCC format */
#define SYMPART

/*
  Function: z_cccReadHeader

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
void z_cccReadHeader(FILE         *infile,
                   pastix_int_t *Nrow,
                   pastix_int_t *Ncol,
                   pastix_int_t *Nnzero,
                   char         *Type);


/*
   Function: z_cccRead

   Reads Matrix in ccc format.

   Header format is described in <z_cccReadHeader>,
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
void z_cccRead(char const      *filename,
             pastix_int_t    *Nrow,
             pastix_int_t    *Ncol,
             pastix_int_t    *Nnzero,
             pastix_int_t   **col,
             pastix_int_t   **row,
             pastix_complex64_t **val,
             char           **Type,
             char           **RhsType);

#endif /* Z_CCCREAD_H_ */
