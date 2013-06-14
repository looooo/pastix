/*
  File: laplacian.h

  Generate a laplacian 
  
  Example :
  >  2 -1  0  0
  > -1  2 -1  0 
  >  0  1  2 -1
  >  0  0 -1  2

*/

/*
  Function: genlaplacian

  Generate a laplacian of size *n*

  Parameters:
    n  	    - Size of the wanted matrix
    nnzeros - Number of non zeros in the produced matrice
    ia      - Index of first element of each column in *row* and *val*
    ja 	    - Row of eah element				       
    avals   - Value of each element				       
    rhs     - Right-hand-side member
    type    - Type of the matrix				       	 
    rhstype - Type of the right hand side.			       	 
 */
int genlaplacian(pastix_int_t     n, 
		 pastix_int_t    *nnzeros,
		 pastix_int_t   **ia, 
		 pastix_int_t   **ja,
		 pastix_float_t **avals,
		 pastix_float_t **rhs,
		 char           **type,
		 char           **rhstype);
