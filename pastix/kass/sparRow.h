/************************************************************/
/**                                                        **/
/**   NAME       : sparRow.h                               **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 17 May 2005     **/
/**                                                        **/
/**                                                        **/
/************************************************************/

/*
**  The function prototypes.
*/

#ifndef SPARSE_ROW_
#define SPARSE_ROW_

typedef struct SparRow *csptr;
typedef struct SparRow {
  /*--------------------------------------------- 
    | C-style CSR format - used internally
    | for all matrices in CSR format 
    |---------------------------------------------*/

  PASTIX_INT      n;
  PASTIX_INT     *nnzrow; /* length of each row                               */
  double **ma;     /* pointer-to-pointer to store nonzero entries      */
  PASTIX_INT    **ja;     /* pointer-to-pointer to store column indices       */
  PASTIX_INT      inarow; /* This flag means the matrix has been allocated as 
		      a single block of memory; it must be desallocated 
		      in consequence                                   */
  PASTIX_INT     *jatab;  /* Used if inarow == 1 to store the the matrix in 
		      two contigue block of memory                     */
  double  *matab;
} SparMat;


PASTIX_INT initCS (csptr amat, PASTIX_INT len);
PASTIX_INT cleanCS(csptr amat);
PASTIX_INT CSnnz  (csptr mat);
PASTIX_INT CS_Perm(csptr mat, PASTIX_INT *perm);

#endif
