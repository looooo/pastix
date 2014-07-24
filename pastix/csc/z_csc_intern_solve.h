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
  File: z_csc_intern_solve.h

  Functions to copy internal CSCd data 
  onto solver matrix coeftab.
*/

#ifndef Z_CSC_INTERN_SOLVE_H
#define Z_CSC_INTERN_SOLVE_H

#include "z_solver.h"
#include "z_csc.h"

/*
  Function: z_Csc2solv_cblk

  Copy the part of the internal CSCd corresponding to
  the column bloc itercblk into the z_SolverMatrix structure 
  coeftab which will be used to compute the decomposition.

  Used in NUMA mode.
 
  Parameters:
    cscmtx   - The internal CSCd matrix.
    datacode - The z_SolverMatrix structure used during decomposition.
    trandcsc - The internal CSCd transpose used in LU decomposition.
    itercblk - Column bloc number in which we had the internal CSCd. 
  
  
*/
void z_Csc2solv_cblk(const z_CscMatrix  *cscmtx,
                     z_SolverMatrix     *solvmtx,
                     pastix_complex64_t *trandcsc,
                     pastix_int_t        itercblk);

#endif /* Z_CSC_INTERN_SOLVE_H */
