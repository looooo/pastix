/**
 *  PaStiX distributed CSC management routines.
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
#include "z_cscd_utils.h"
#ifndef CSCD_UTILS_INTERN_H
#define CSCD_UTILS_INTERN_H
int z_cscd_addlocal_int(pastix_int_t   n   , const pastix_int_t *  ia   , const pastix_int_t *  ja   , const pastix_complex64_t *  a   , const pastix_int_t * l2g,
                        pastix_int_t   addn, const pastix_int_t *  addia,       pastix_int_t *  addja, const pastix_complex64_t *  adda, const pastix_int_t * addl2g,
                        pastix_int_t * newn, pastix_int_t ** newia, pastix_int_t ** newja, pastix_complex64_t ** newa,
                        CSCD_OPERATIONS_t OP, int dof, int malloc_flag);

/*
 * Function: cscd_redispatch_int
 *
 * Redistribute the first cscd into a new one using *dl2g*.
 *
 * - gather all new loc2globs on all processors.
 * - allocate *dia*, *dja* and *da*.
 * - Create new CSC for each processor and send it.
 * - Merge all new CSC to the new local CSC with <cscd_addlocal_int>.
 *
 * If communicator size is one, check that n = dn and
 * l2g = dl2g and simply create a copy of the first cscd.
 *
 * Parameters:
 *   n           - Number of local columns
 *   ia          - First cscd starting index of each column in *ja* and *a*
 *   ja          - Row of each element in first CSCD
 *   a           - value of each cscd in first CSCD (can be NULL)
 *   rhs         - right-hand-side member corresponding to the first CSCD (can be NULL)
 *   nrhs        - number of right-hand-side.
 *   l2g         - local 2 global column numbers for first cscd
 *   dn          - Number of local columns
 *   dia         - New cscd starting index of each column in *ja* and *a*
 *   dja         - Row of each element in new CSCD
 *   da          - value of each cscd in new CSCD
 *   rhs         - right-hand-side member corresponding to the new CSCD
 *   dl2g        - local 2 global column numbers for new cscd
 *   malloc_flag - Internal (API_YES) or external (API_NO) malloc use.
 *   comm        - MPI communicator
 *
 * Returns:
 *   EXIT_SUCCESS - If all goes well
 *   EXIT_FAILURE - If commsize = 1 and *n* != *dn* or *l2g* != *dl2g*.
 */
int z_cscd_redispatch_int(pastix_int_t   n, pastix_int_t *   ia, pastix_int_t *   ja, pastix_complex64_t *   a, pastix_complex64_t *  rhs,  pastix_int_t nrhs, pastix_int_t *   l2g,
                        pastix_int_t  dn, pastix_int_t ** dia, pastix_int_t ** dja, pastix_complex64_t ** da, pastix_complex64_t ** drhs, pastix_int_t *  dl2g,
                        int  malloc_flag, MPI_Comm comm, pastix_int_t dof);

int z_cscd_symgraph_int(pastix_int_t      n, const pastix_int_t *ia, const pastix_int_t *ja, const pastix_complex64_t *a,
                        pastix_int_t * newn, pastix_int_t **  newia, pastix_int_t **  newja, pastix_complex64_t ** newa,
                        pastix_int_t *  l2g, MPI_Comm comm, int malloc_flag);


/*
  Function: cscd_build_g2l

  Construct global to local tabular containing local number of global columns
  if one column is local, and -owner if column is not local.

  For i in 0, gN
     g2l[i] = i local number if i is local
     g2l[i] = -p if p is the owner of the column i

  Parameters:
    n        - Number of local columns
    colptr   - Starting index of each columns in *ja*
    rows     - Row of each element.
    values   - Value of each element.
    l2g      - global number of each local column.
    correct  - Flag indicating if we can correct the symmetry.
    dof      - Number of degrees of freedom.
    comm     - MPI communicator
 */
int z_cscd_build_g2l( pastix_int_t   ncol,
                      pastix_int_t  *loc2glob,
                      MPI_Comm       comm,
                      pastix_int_t  *gN,
                      pastix_int_t **g2l);
/*
   Function: cscd_noDiag

   Removes diagonal elements from a CSCD.
   *ja* and *a* can be reallocated to
   ia[n]-1 elements after this call.

   Parameters:
     n           - Number of local columns
     ia          - First cscd starting index of each column in *ja* and *a*
     ja          - Row of each element in first CSCD
     a           - value of each cscd in first CSCD (can be NULL)
     l2g         - local 2 global column numbers for first cscd

*/
int z_cscd_noDiag(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_complex64_t * a, const pastix_int_t * l2g);


/*
  Function: cscd_checksym

  Check if the CSCD graph is symetric.

  Parameters:
    n   - Number of local columns
    ia  - Starting index of each columns in *ja*
    ja  - Row of each element.
    l2g - global number of each local column.
    correct  - Flag indicating if we can correct the symmetry.
    alloc    - indicate if allocation on CSC uses internal malloc.
    dof      - Number of degrees of freedom.
    comm     - MPI communicator
*/
int z_cscd_checksym(pastix_int_t      n,
                  pastix_int_t     *colptr,
                  pastix_int_t    **rows,
                  pastix_complex64_t  **values,
                  pastix_int_t     *l2g,
                  int      correct,
                  int      alloc,
                  int      dof,
                  MPI_Comm comm);


/**
 *   Function: cscd2csc_int
 *
 *   Transform a cscd to a csc.
 *   colptr2, row2, avals2, rhs2, perm2, invp2 are allocated here.
 *
 *   External function, allocation are not of the internal type.
 *
 *   Parameters:
 *      lN          - number of local column.
 *      lcolptr     - starting index of each local column in row and avals.
 *      lrow        _ row number of each local element.
 *      lavals      - values of each local element.
 *      lrhs        - local part of the right hand side.
 *      lperm       - local part of the permutation tabular.
 *      linvp       - Means nothing, to suppress.
 *      gN          - global number of columns (output).
 *      gcolptr     - starting index of each column in row2 and avals2 (output).
 *      grow        - row number of each element (output).
 *      gavals      - values of each element (output).
 *      grhs        - global right hand side (output).
 *      gperm       - global permutation tabular (output).
 *      ginvp       - global reverse permutation tabular (output).
 *      loc2glob    - global number of each local column.
 *      pastix_comm - PaStiX MPI communicator.
 *      ndof        - Number of degree of freedom per node.
 *      intern_flag - Decide if malloc will use internal or external macros.
 */
void  z_cscd2csc_int(pastix_int_t  lN, pastix_int_t *  lcolptr, pastix_int_t * lrow, pastix_complex64_t * lavals,
                   pastix_complex64_t * lrhs, pastix_int_t * lperm, pastix_int_t * linvp,
                   pastix_int_t *gN, pastix_int_t ** gcolptr, pastix_int_t **grow, pastix_complex64_t **gavals,
                   pastix_complex64_t **grhs, pastix_int_t **gperm, pastix_int_t **ginvp,
                   pastix_int_t *loc2glob, MPI_Comm pastix_comm, pastix_int_t ndof, int intern_flag);

/*
  Function: cscd_redispatch_scotch

  Renumber the columns to have first columns on first proc for Scotch

  Parameters:
    n           - Number of local columns
    ia          - First cscd starting index of each column in *ja* and *a*
    ja          - Row of each element in first CSCD
    a           - value of each cscd in first CSCD (can be NULL)
    l2g         - local 2 global column numbers for first cscd
    dn          - Number of local columns
    dia         - First cscd starting index of each column in *ja* and *a*
    dja         - Row of each element in first CSCD
    da          - value of each cscd in first CSCD (can be NULL)
    l2g         - local 2 global column numbers for first cscd
    comm        - MPI communicator

  Returns:
    EXIT_SUCCESS if already well distributed, 2 if redistributed

 */
int z_cscd_redispatch_scotch(pastix_int_t   n, pastix_int_t *   ia, pastix_int_t *   ja, pastix_complex64_t *   a, pastix_int_t *   l2g,
                           pastix_int_t *dn, pastix_int_t ** dia, pastix_int_t ** dja, pastix_complex64_t ** da, pastix_int_t ** dl2g,
                           MPI_Comm comm);
#endif /* CSCD_UTILS_INTERN_H */
