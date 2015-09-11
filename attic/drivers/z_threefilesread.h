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
  File: threetilesread.h

  Reads matrix from three files in IJV separated format.

 */

/*
   Function: z_threeFilesReadHeader

   Read header from three file IJV format.

   Header contains:
   > Nrow Ncol Nnzero
   or
   > Ncol
   > Nnzero

   Parameters:
     infile - file to read header from
     Nrow   - Number of rows
     Ncol   - Number of columns
     Nnzero - Number of non zeros
     Type   - Type of the matrix (always "RUA")
 */
void z_threeFilesReadHeader(FILE         *infile,
                          pastix_int_t *Nrow,
                          pastix_int_t *Ncol,
                          pastix_int_t *Nnzero,
                          char         *Type);

/*
  Function: z_threeFilesRead

  Read matrix from three files IJV

  header file is "filename"/header
  columns file is "filename"/ia_threeFiles
  rows file is "filename"/ja_threeFiles
  values file is "filename"/ra_threeFiles

  header is describde in <z_threeFilesReadHeader>
  each other file contain one element by line.

  Parameters:
    dirname - Path to the directory containing matrix
    Ncol    - Number of columns
    Nrow    - Number of rows
    Nnzero  - Number of non zeros
    col     - Index of first element of each column in *row* and *val*
    row     - Row of eah element
    val     - Value of each element
    Type    - Type of the matrix
    RhsType - Type of the right-hand-side.
 */
void z_threeFilesRead(char const      *dirname,
                    pastix_int_t    *Ncol,
                    pastix_int_t    *Nrow,
                    pastix_int_t    *Nnzero,
                    pastix_int_t   **col,
                    pastix_int_t   **row,
                    pastix_complex64_t **val,
                    char           **Type,
                    char           **RhsType);
