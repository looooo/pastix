/*
 *  File: mmdread.h
 *
 * Distributed Matrix Market driver.
 */
#ifndef MMDREAD_H
#define MMDREAD_H
/*
 *  Function: DistributedMatrixMarketRead
 *
 *  Reads a matrix in distributed matrix market format
 *
 *  For more information about matrix market format see mmio.c/mmio.h
 *
 * Parameters:
 *   filename - root file of the MMD matrix.
 *   Ncol     - Number of columns
 *   Nrow     - Number of rows
 *   Nnzero   - Number of non zeros
 *   col      - Index of first element of each column in *row* and *val*
 *   row      - Row of eah element
 *   val      - Value of each element
 *   Type     - Type of the matrix
 *   RhsType  - Type of the right-hand-side.
 *
 */
void DistributedMatrixMarketRead(char const      *filename,
                                 pastix_int_t    *Ncol,
                                 pastix_int_t    *Nrow,
                                 pastix_int_t    *Nnzero,
                                 pastix_int_t   **col,
                                 pastix_int_t   **row,
                                 pastix_float_t **val,
                                 pastix_int_t   **l2g,
                                 char           **Type,
                                 char           **RhsType);
#endif /* not MMDREAD_H */
