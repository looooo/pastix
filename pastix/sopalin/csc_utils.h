#ifndef CSC_UTILS_H
#define CSC_UTILS_H

/*
  Function: csc_symgraph

  
  Modify the CSC to a symetric graph one.
  Don't use it on a lower symetric CSC 
  it would give you all the CSC upper + lower.
  
  External function

  Parameters: 
    n     - Number of columns/vertices
    ia	  - Starting index of each column in *ja* and *a*
    ja	  - Row index of each element
    a 	  - Value of each element,can be NULL    
    newn  - New number of column
    newia - Starting index of each column in *ja* and *a* 
    newja - Row index of each element
    newa  - Value of each element,can be NULL    

 */
int csc_symgraph(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a, 
		 pastix_int_t *newn, pastix_int_t **newia, pastix_int_t **newja, pastix_float_t **newa);



/*
  Function: csc_symgraph_int

  
  Modify the CSC to a symetric graph one.
  Don't use it on a lower symetric CSC 
  it would give you all the CSC upper + lower.
  
  Parameters: 
    n           - Number of columns/vertices
    ia	        - Starting index of each column in *ja* and *a*
    ja	        - Row index of each element
    a 	        - Value of each element,can be NULL    
    newn        - New number of column
    newia       - Starting index of each column in *ja* and *a* 
    newja       - Row index of each element
    newa        - Value of each element,can be NULL    
    malloc_flag - flag to indicate if function call is intern to pastix or extern.
 */
int csc_symgraph_int (pastix_int_t n,     pastix_int_t * ia,    pastix_int_t * ja,    pastix_float_t * a, 
		      pastix_int_t *newn, pastix_int_t **newia, pastix_int_t **newja, pastix_float_t **newa, 
		      int malloc_flag);



/** 
    Function: csc_noDiag
    
    Supress diagonal term.              
    After this call, *ja* can be reallocated to *ia[n] -1*.
    
    Parameters:
      n  - size of the matrix.
      ia - Index in *ja* and *a* of the first element of each column
      ja - row of each element
      a  - value of each element, can be set to NULL

    Returns:
      ia and ja tabulars modified.
*/
void csc_noDiag(pastix_int_t baseval, pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a);

/*
  Function: csc_check_doubles
  
  Check if the csc contains doubles and if correct if asked

  Assumes that the CSC is sorted.

  Assumes that the CSC is Fortran numeroted (base 1)

  Parameters:
    n      - Size of the matrix.
    colptr - Index in *rows* and *values* of the first element of each column
    rows   - row of each element
    values - value of each element
    dof    - Number of degrees of freedom
    flag   - Indicate if user wants correction (<API_BOOLEAN>)
    flagalloc - indicate if allocation on CSC uses internal malloc. 

    
  Returns:
    API_YES - If the matrix contained no double or was successfully corrected.
    API_NO  - Otherwise.
*/
int csc_check_doubles(pastix_int_t      n,
		      pastix_int_t   *  colptr,
		      pastix_int_t   ** rows,
		      pastix_float_t ** values, 
		      int      dof,
		      int      flag,
		      int      flagalloc);

/*
  Function: csc_checksym

    Check if the CSC graph is symetric.
    
    For all local column C, 
    
    For all row R in the column C,
    
    We look in column R if we have the row number C.
       
    If we can correct we had missing non zeros.
    
    Assumes that the CSC is Fortran numbered (1 based).
    
    Assumes that the matrix is sorted.

  Parameters:
    n        - Number of local columns
    colptr   - Starting index of each columns in *ja*
    rows     - Row of each element.
    values   - Value of each element.
    correct  - Flag indicating if we can correct the symmetry.
    alloc    - indicate if allocation on CSC uses internal malloc. 
    dof      - Number of degrees of freedom.
*/
int csc_checksym(pastix_int_t      n, 
		 pastix_int_t     *colptr, 
		 pastix_int_t    **rows, 
		 pastix_float_t  **values, 
		 int      correct,
		 int      alloc,
		 int      dof);

void CSC_colPerm(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a, pastix_int_t *cperm);
void CSC_colScale(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a, pastix_float_t *dcol);
void CSC_rowScale(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a, pastix_float_t *drow);

void CSC_sort(pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a);
void CSC_Fnum2Cnum(pastix_int_t *ja, pastix_int_t *ia, pastix_int_t n);
void CSC_Cnum2Fnum(pastix_int_t *ja, pastix_int_t *ia, pastix_int_t n);

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
int CSC_buildZerosAndNonZerosGraphs(pastix_int_t     n,
				    pastix_int_t    *colptr,
				    pastix_int_t    *rows,
				    pastix_float_t  *values,
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
int CSC_isolate(pastix_int_t     n,
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
int csc_save(pastix_int_t      n,
	     pastix_int_t    * colptr,
	     pastix_int_t    * rows,
	     pastix_float_t  * values,
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
int csc_load(pastix_int_t    *  n,
	     pastix_int_t    ** colptr,
	     pastix_int_t    ** rows,
	     pastix_float_t  ** values,
	     int    *  dof,
	     FILE   *  infile);

#endif /* CSC_UTILS_H */
