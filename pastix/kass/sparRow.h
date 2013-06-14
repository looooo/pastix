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

  pastix_int_t      n;
  pastix_int_t     *nnzrow; /* length of each row                               */
  double **ma;     /* pointer-to-pointer to store nonzero entries      */
  pastix_int_t    **ja;     /* pointer-to-pointer to store column indices       */
  pastix_int_t      inarow; /* This flag means the matrix has been allocated as 
		      a single block of memory; it must be desallocated 
		      in consequence                                   */
  pastix_int_t     *jatab;  /* Used if inarow == 1 to store the the matrix in 
		      two contigue block of memory                     */
  double  *matab;
} SparMat;


pastix_int_t initCS (csptr amat, pastix_int_t len);
pastix_int_t cleanCS(csptr amat);
pastix_int_t CSnnz  (csptr mat);
pastix_int_t CS_Perm(csptr mat, pastix_int_t *perm);

#endif
