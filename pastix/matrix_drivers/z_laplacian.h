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
  File: z_laplacian.h

  Generate a laplacian 
  
  Example :
  >  2 -1  0  0
  > -1  2 -1  0 
  >  0  1  2 -1
  >  0  0 -1  2

*/

/*
  Function: z_genlaplacian

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
int z_genlaplacian(pastix_int_t     n, 
		 pastix_int_t    *nnzeros,
		 pastix_int_t   **ia, 
		 pastix_int_t   **ja,
		 pastix_complex64_t **avals,
		 pastix_complex64_t **rhs,
		 char           **type,
		 char           **rhstype);
