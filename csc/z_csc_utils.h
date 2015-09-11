/**
 *  PaStiX CSC management routines.
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
#ifndef CSC_UTILS_H
#define CSC_UTILS_H


/**
 *
 * @ingroup csc_utils
 *
 * Modify the CSC to a symetric graph one.
 * Don't use it on a lower symetric CSC
 * it would give you all the CSC upper + lower.
 *
 * @param[in]  n     Number of columns/vertices
 * @param[in]  ia    Starting index of each column in *ja* and *a*
 * @param[in]  ja    Row index of each element
 * @param[in]  a     Value of each element,can be NULL
 * @param[out] newn  New number of column
 * @param[out] newia Starting index of each column in *ja* and *a*
 * @param[out] newja Row index of each element
 * @param[out] newa  Value of each element,can be NULL
 **/
int z_csc_symgraph ( pastix_int_t               n,
                     const pastix_int_t        *ia,
                     const pastix_int_t        *ja,
                     const pastix_complex64_t  *a,
                     pastix_int_t              *newn,
                     pastix_int_t             **newia,
                     pastix_int_t             **newja,
                     pastix_complex64_t       **newa);




/**
 *
 * @ingroup csc_utils
 *
 * Modify the CSC to a symetric graph one.
 * Don't use it on a lower symetric CSC
 * it would give you all the CSC upper + lower.
 *
 * Internal function.
 *
 * @param[in]  n           Number of columns/vertices
 * @param[in]  ia          Starting index of each column in *ja* and *a*
 * @param[in]  ja          Row index of each element
 * @param[in]  a           Value of each element,can be NULL
 * @param[out] newn        New number of column
 * @param[out] newia       Starting index of each column in *ja* and *a*
 * @param[out] newja       Row index of each element
 * @param[out] newa        Value of each element,can be NULL
 * @param[in]  malloc_flag flag to indicate if function call is intern to pastix
 *                         or extern.
 **/
int z_csc_symgraph_int ( pastix_int_t               n,
                         const pastix_int_t        *ia,
                         const pastix_int_t        *ja,
                         const pastix_complex64_t  *a,
                         pastix_int_t              *newn,
                         pastix_int_t             **newia,
                         pastix_int_t             **newja,
                         pastix_complex64_t       **newa,
                         int                        malloc_flag);



/**
 *
 * @ingroup csc_utils
 *
 * Supress diagonal term.
 * After this call, *ja* can be reallocated to *ia[n] -1*.
 *
 * @param[in]      baseval     Initial numbering value (0 or 1).
 * @param[in]      n           Number of columns/vertices
 * @param[in,out]  ia          Starting index of each column in *ja* and *a*
 * @param[in,out]  ja          Row index of each element
 * @param[in,out]  a           Value of each element,can be NULL
 **/
void z_csc_noDiag( pastix_int_t        baseval,
                   pastix_int_t        n,
                   pastix_int_t       *ia,
                   pastix_int_t       *ja,
                   pastix_complex64_t *a);

/**
 *
 * @ingroup csc_utils
 *
 * Check if the csc contains doubles and if correct if asked
 *
 * Assumes that the CSC is sorted.
 *
 * Assumes that the CSC is Fortran numeroted (base 1)
 *
 * @param[in]     n          Size of the matrix.
 * @param[in,out] colptr     Index in *rows* and *values* of the first element
 *                           of each column
 * @param[in,out] rows       row of each element
 * @param[in,out] values     value of each element
 * @param[in]     dof        Number of degrees of freedom
 * @param[in]     flag       Indicate if user wants correction (<API_BOOLEAN>)
 * @param[in]     flagalloc  indicate if allocation on CSC uses internal malloc.
 *
 *
 * @Returns:
 *   API_YES - If the matrix contained no double or was successfully corrected.
 *   API_NO  - Otherwise.
*/

int z_csc_check_doubles(pastix_int_t         n,
                        pastix_int_t        *colptr,
                        pastix_int_t       **rows,
                        pastix_complex64_t **values,
                        int                  dof,
                        int                  flag,
                        int                  flagalloc);


/**
 *
 * @ingroup csc_utils
 *
 * Check if the CSC graph is symetric.
 *
 *   For all local column C,
 *
 *   For all row R in the column C,
 *
 *   We look in column R if we have the row number C.
 *
 *   If we can correct we had missing non zeros.
 *
 *   Assumes that the CSC is Fortran numbered (1 based).
 *
 *   Assumes that the matrix is sorted.
 *
 * @param[in]     n         Number of local columns
 * @param[in,out] colptr    Starting index of each columns in *ja*
 * @param[in,out] rows      Row of each element.
 * @param[in,out] values    Value of each element.
 * @param[in]     correct   Flag indicating if we can correct the symmetry.
 * @param[in]     alloc     indicate if allocation on CSC uses internal malloc.
 * @param[in]     dof       Number of degrees of freedom.
*/
int z_csc_checksym(pastix_int_t         n,
                   pastix_int_t        *colptr,
                   pastix_int_t       **rows,
                   pastix_complex64_t **values,
                   int                  correct,
                   int                  alloc,
                   int                  dof);

void z_csc_colPerm(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_complex64_t *a, pastix_int_t *cperm);
void z_csc_colScale(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_complex64_t *a, pastix_complex64_t *dcol);
void z_csc_rowScale(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_complex64_t *a, pastix_complex64_t *drow);

void z_csc_sort(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_complex64_t *a, pastix_int_t ndof);
void z_csc_Fnum2Cnum(pastix_int_t *ja, pastix_int_t *ia, pastix_int_t n);
void z_csc_Cnum2Fnum(pastix_int_t *ja, pastix_int_t *ia, pastix_int_t n);

/*
  Function: CSC_buildZerosAndNonZerosGraphs

  Separate a graph in two graphs, following
  wether the diagonal term of a column is null or not.

  Parameters:
    n, colptr, rows, values  - The initial CSC
    n_nz, colptr_nz, rows_nz - The graph of the non-null diagonal part.
    n_z, colptr_z, rows_z    - The graph of the null diagonal part.
    perm                     - Permutation to go from the first graph to
                               the one composed of the two graph concatenated.
    revperm                  - Reverse permutation tabular.
    criteria                 - Value beside which a number is said null.
*/
int z_csc_buildZerosAndNonZerosGraphs(pastix_int_t     n,
                                    pastix_int_t    *colptr,
                                    pastix_int_t    *rows,
                                    pastix_complex64_t  *values,
                                    pastix_int_t    *n_nz,
                                    pastix_int_t   **colptr_nz,
                                    pastix_int_t   **rows_nz,
                                    pastix_int_t    *n_z,
                                    pastix_int_t   **colptr_z,
                                    pastix_int_t   **rows_z,
                                    pastix_int_t    *perm,
                                    pastix_int_t    *revperm,
                                    double  criteria);

/*
  Function: CSC_isolate

  Isolate a list of unknowns at the end of the CSC.

  Parameters:
    n            - Number of columns.
    colptr       - Index of first element of each column in *ia*.
    rows         - Rows of each non zeros.
    n_isolate    - Number of unknow to isolate.
    isolate_list - List of unknown to isolate.
*/
int z_csc_isolate(pastix_int_t     n,
                pastix_int_t    *colptr,
                pastix_int_t    *rows,
                pastix_int_t     n_isolate,
                pastix_int_t    *isolate_list,
                pastix_int_t    *perm,
                pastix_int_t    *revperm);


/*
  Function: csc_save

  Save a csc on disk.

  Parameters:
    n       - number of columns
    colptr  - First cscd starting index of each column in *ja* and *a*
    rows    - Row of each element in first CSCD
    values  - value of each cscd in first CSCD (can be NULL)
    dof     - Number of degrees of freedom
    outfile - Output stream.

  Return:
    NO_ERR

*/
int z_csc_save(pastix_int_t      n,
             pastix_int_t    * colptr,
             pastix_int_t    * rows,
             pastix_complex64_t  * values,
             int      dof,
             FILE   * outfile);
/*
  Function: csc_load

  Load a csc from disk.

  Fill *n*, *colptr*, *rows*, *values* and *dof* from *infile*.

  Parameters:
    n       - number of columns
    colptr  - First cscd starting index of each column in *ja* and *a*
    rows    - Row of each element in first CSCD
    values  - value of each cscd in first CSCD (can be NULL)
    dof     - Number of degrees of freedom
    outfile - Output stream.

  Return:
    NO_ERR

*/
int z_csc_load(pastix_int_t    *  n,
             pastix_int_t    ** colptr,
             pastix_int_t    ** rows,
             pastix_complex64_t  ** values,
             int    *  dof,
             FILE   *  infile);

#endif /* CSC_UTILS_H */
