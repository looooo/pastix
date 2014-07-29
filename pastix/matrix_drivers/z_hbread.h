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
  File: z_hbread.h

  Interface to the Harwell-Boeing driver in C (iohb.c)
 */

/*
  Function: z_HBRead
   
  Interface to the Harwell-Boeing driver in C (iohb.c)

  Parameters
    filename - Path to the file to read from
    Nrow     - Number of rows
    Ncol     - Number of columns
    Nnzero   - Number of non zeros
    col      - Index of first element of each column in *row* and *val*
    row      - Row of eah element				       
    val      - Value of each element				       
    Type     - Type of the matrix				       
    RhsType  - Type of the right hand side.			       
 */
void z_HBRead(char const      *filename, 
	    pastix_int_t    *Nrow, 
	    pastix_int_t    *Ncol, 
	    pastix_int_t    *Nnzero, 
	    pastix_int_t   **col, 
	    pastix_int_t   **row, 
	    pastix_complex64_t **val, 
	    char           **Type, 
	    char           **RhsType);
