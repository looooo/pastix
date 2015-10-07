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
#ifndef Z_CSC_H
#define Z_CSC_H
#include "common.h"

/* Section: Macros */
/*
  Macro: CSC_FNBR

  Accessor to the number of column block.

  Parameters:
    a - Pointer to the CSC.

  Returns:
    Number of column block.
 */
#define CSC_FNBR(a)     (a)->cscfnbr /* cblk nbr */
/*
  Macro: CSC_FTAB

  Accessor to the array of column blocks.

  Parameters:
    a - Pointer to the CSC.

  Returns:
    Address of the array of column blocks.
*/
#define CSC_FTAB(a)     (a)->cscftab
/*
  Macro: CSC_COLNBR

  Accessor to the number of column in a block column.

  Parameters:
    a - Pointer to the CSC matrix.
    b - Column block index.

  Returns:
    The number of column in the block column
*/
#define CSC_COLNBR(a,b) (a)->cscftab[b].colnbr
/*
  Macro: CSC_COLTAB

  Accessor to the array of start for each column
  in the rows and values arrays.

  Parameters:
    a - Pointer to the CSC matrix.
    b - Column block index.

  Reurns:
    Address of the array of indexes of start for each column
    in the rows and values arrays.

*/
#define CSC_COLTAB(a,b) (a)->cscftab[b].coltab
/*
  Macro: CSC_COL

  Accessor to the index of first element of a column in rows
  and values.

  Parameters:
    a - Pointer to the CSC matrix.
    b - Column block index.
    c - Column index.
*/
#define CSC_COL(a,b,c)  (a)->cscftab[b].coltab[c]
/*
   Macro: CSC_ROWTAB

   Accessor to the array of rows.

   Parameters:
     a - Pointer to the CSC matrix.
*/
#define CSC_ROWTAB(a)   (a)->rowtab
/*
   Macro: CSC_ROW

   Accessor to a row in the CSC.

   Parameters:
     a - Pointer to the CSC matrix.
     b - Index of the row.
*/
#define CSC_ROW(a,b)    (a)->rowtab[b]
/*
   Macro: CSC_VALTAB

   Accessor to the array of values.

   Parameters:
     a - Pointer to the CSC matrix.
*/
#define CSC_VALTAB(a)   (a)->valtab
/*
   Macro: CSC_VAL

   Accessor to a value in the CSC.

   Parameters:
     a - Pointer to the CSC matrix.
     b - Index of the value.
*/
#define CSC_VAL(a,b)    (a)->valtab[b]
/*
   Macro: CSC_FROW

   Accessor to the first row of the column $c$ in
   the column block $b$.

   Parameters:
     a - Pointer to the CSC matrix.
     b - Column block index.
     c - Column index.
*/
#define CSC_FROW(a,b,c) (a)->rowtab[(a)->cscftab[b].coltab[c]]
/*
   Macro: CSC_FROW

   Accessor to the first value of the column $c$ in
   the column block $b$.

   Parameters:
     a - Pointer to the CSC matrix.
     b - Column block index.
     c - Column index.
*/
#define CSC_FVAL(a,b,c) (a)->valtab[(a)->cscftab[b].coltab[c]]
/*
  Macro: CSC_VALNBR

  Compute the Number of element on the matrix.

  Parameters:
    a - Pointer to the CSC matrix.
*/
#define CSC_VALNBR(a)   (a)->cscftab[(a)->cscfnbr\
				     -1].coltab[(a)->cscftab[(a)->cscfnbr \
							     -1].colnbr]

/* Section: Structures */
/*
   Structure: CscFormat_

   Internal block column structure.

   Contains:
     colnbr - Number of columns in the block column.
     coltab - Array of indexes of the start of each column in
	      the row and value arrays.
*/
struct z_CscFormat_ {
  pastix_int_t   colnbr;
  pastix_int_t * coltab;
};

/*
   Type: CscFormat

   See <CscFormat_> structure.
*/
typedef struct z_CscFormat_ z_CscFormat;

/*
  Structure: z_CscMatrix_

  Internal column block distributed CSC matrix.

  Contains:
    cscfnbr - Number of column block.
    cscftab - Array of Block column structures. (<CscFormat>)
    rowtab  - Array of rows in the matrix.
    valtab  - Array of values of the matrix.
    type    - 'S' for symmetric, 'H' for hermitian, U for unsymmetric.
*/
struct z_CscMatrix_ {
  pastix_int_t         cscfnbr;
  z_CscFormat          * cscftab;
  pastix_int_t       * rowtab;
  pastix_complex64_t     * valtab;
  char                 type;
};
/*
  Type: z_CscMatrix

  See <z_CscMatrix_> structure.
*/
typedef struct z_CscMatrix_ z_CscMatrix;

#endif /* CSC_H */
