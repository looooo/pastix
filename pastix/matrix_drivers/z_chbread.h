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
  File: z_chbread.h

  Read a matrix in chb format.
 */



/*
  Function: z_chbReadHeader

  Reads header from a chb matrix file.

  header format is:
  > title 73-80 Key
  > 1-14 totcrd 15-28 ptrcrd 29-42 indcrd 43-56 valcrd 57-70 rhscrd
  > 1-3 mxtype 15-28 nrow 29-42 ncol 43-56 nnzero 57-70 neltvl
  > 1-16 ptrfmt 17-32 indfmt 33-52 valfmt 53-72 rhsfmt
  > 1  2 rhstyp 3  15-28 nrhs 29-42 nrhsix

  Parameters
    infile  - File to read from
    Type    - Type of the matrix
    Nrow    - Number of row in the matrix
    Ncol    - Number of columns in the matrix
    Nnzero  - Number of non zeros in the matrix
    Nrhs    - Number of right-hand-side terms
    Ptrfmt  -
    Indfmt  -
    Valfmt  -
    Rhsfmt  -
    Ptrcrd  -
    Indcrd  -
    Valcrd  -
    Rhscrd  -
    RhsType - Type of right-hand-side term(s)

 */
void z_chbReadHeader(FILE         *infile,
                   char         *Type,
                   pastix_int_t *Nrow,
                   pastix_int_t *Ncol,
                   pastix_int_t *Nnzero,
                   pastix_int_t *Nrhs,
                   char         *Ptrfmt,
                   char         *Indfmt,
                   char         *Valfmt,
                   char         *Rhsfmt,
                   pastix_int_t *Ptrcrd,
                   pastix_int_t *Indcrd,
                   pastix_int_t *Valcrd,
                   pastix_int_t *Rhscrd,
                   char         *RhsType);
/*
  Function: z_chbRead

  Reads a matrix in chb format.

  Header is described in <z_chbReadHeader>
  Formats are sicribed in <chbParseRfmt> and <chbParseIfmt>

  In our file we have
  header
  valuesFormat
  rowFormat
  columnFormat
  (rhsFormat)
  then the columns,
  the rows,
  the values,
  (the rhs)


  Parameters:
    filename - Path to the file to read from
    Nrow     - Number of rows
    Ncol     - Number of columns
    Nnzero   - Number of non zeros
    col      - Index of first element of each column in *row* and *val*
    row      - Row of eah element
    val      - Value of each element
    Type     - Type of the matrix
    RhsType  - Type of the right-hand-side terms.
    rhs      - right-hand-side term(s)


 */
void z_chbRead(char const      *filename,
             pastix_int_t    *Nrow,
             pastix_int_t    *Ncol,
             pastix_int_t    *Nnzero,
             pastix_int_t   **col,
             pastix_int_t   **row,
             pastix_complex64_t **val,
             char           **Type,
             char           **RhsType,
             pastix_complex64_t **rhs);


/*
  Function: hbParseRfmt

  CHB float format parser

  Format is like :
  > (3(1P,E25.16))
  or
  > (1P,3E25.16)
  or
  > (1P3E25.16)
  or
  > (3E25.16)
  or
  > (3E25)
  for perline = 3, format = E, width = 25 and prec = 16


  Parameters:
    fmt      - format to parse
    perline  - number of element per line
    width    -
    prec     - Precision
    flag     -
*/
void chbParseRfmt(char *fmt, pastix_int_t *perline, pastix_int_t *width, pastix_int_t *prec, char *flag);

/*
  Function:chbParseIfmt

  CHB integer format parser

  format is :
  > (perlineIwidth) or (X,perlineIwidth)
  Parameters:
    fmt      - format to parse
    perline  - number of element per line
    width    -
    flag     -
*/
void chbParseIfmt(char *fmt, pastix_int_t *perline, pastix_int_t *width, pastix_int_t *flag);
