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
  Function: z_driverFdupros

  Interface to the fortran driver from Fabrice Dupros.
  Binary format.

*/
void z_driverFdupros(char const      *filename, 
		   pastix_int_t    *Nrow, 
		   pastix_int_t    *Ncol, 
		   pastix_int_t    *Nnzero, 
		   pastix_int_t   **col, 
		   pastix_int_t   **row, 
		   pastix_complex64_t **val,
		   pastix_complex64_t **rhs, 
		   char           **Type, 
		   char           **RhsType);

/*
  Function: z_driverFdupros_dist

  Interface to the fortran driver from Fabrice Dupros.
  Binary format, distributed.

*/
void z_driverFdupros_dist(char const      *filename, 
			pastix_int_t    *Nrow, 
			pastix_int_t    *Ncol, 
			pastix_int_t    *Nnzero, 
			pastix_int_t   **col, 
			pastix_int_t   **row,
			pastix_int_t   **loc2glob,
			pastix_complex64_t **val,
			pastix_complex64_t **rhs, 
			char           **Type, 
			char           **RhsType,
			MPI_Comm         pastix_comm);

