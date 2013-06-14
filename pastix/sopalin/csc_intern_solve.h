/*
  File: csc_intern_solve.h

  Functions to copy internal CSCd data 
  onto solver matrix coeftab.
*/

#ifndef CSC_INTERN_SOLVE_H
#define CSC_INTERN_SOLVE_H

/*
  Function: Csc2solv_cblk

  Copy the part of the internal CSCd corresponding to
  the column bloc itercblk into the SolverMatrix structure 
  coeftab which will be used to compute the decomposition.

  Used in NUMA mode.
 
  Parameters:
    cscmtx   - The internal CSCd matrix.
    datacode - The SolverMatrix structure used during decomposition.
    trandcsc - The internal CSCd transpose used in LU decomposition.
    itercblk - Column bloc number in which we had the internal CSCd. 
  
  
*/
void Csc2solv_cblk(const CscMatrix *cscmtx, 
		   SolverMatrix    *solvmtx, 
		   PASTIX_FLOAT           *trandcsc, 
		   PASTIX_INT              itercblk);

#endif /* CSC_INTERN_SOLVE_H */
