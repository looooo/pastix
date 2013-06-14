/*
 * File: petscread.h
 *
 * Read Matrix written by PETSc in binary format.
 *
 */

/*
 *  Function: PETScRead
 *
 *  Reads a matrix in a binary PETScformat/
 *
 * Parameters:
 *   dirname - Path to the directory containing matrix
 *   Ncol    - Number of columns
 *   Nrow    - Number of rows
 *   Nnzero  - Number of non zeros
 *   col     - Index of first element of each column in *row* and *val*
 *   row     - Row of eah element
 *   val     - Value of each element
 *   Type    - Type of the matrix
 *   RhsType - Type of the right-hand-side.
 *
 */
void PETScRead(char const      *filename,
               pastix_int_t    *Ncol,
               pastix_int_t    *Nrow,
               pastix_int_t    *Nnzero,
               pastix_int_t   **col,
               pastix_int_t   **row,
               pastix_float_t **val,
               char           **Type,
               char           **RhsType);
