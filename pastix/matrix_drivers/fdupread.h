
/*
  Function: driverFdupros

  Interface to the fortran driver from Fabrice Dupros.
  Binary format.

*/
void driverFdupros(char const      *filename, 
		   pastix_int_t    *Nrow, 
		   pastix_int_t    *Ncol, 
		   pastix_int_t    *Nnzero, 
		   pastix_int_t   **col, 
		   pastix_int_t   **row, 
		   pastix_float_t **val,
		   pastix_float_t **rhs, 
		   char           **Type, 
		   char           **RhsType);

/*
  Function: driverFdupros_dist

  Interface to the fortran driver from Fabrice Dupros.
  Binary format, distributed.

*/
void driverFdupros_dist(char const      *filename, 
			pastix_int_t    *Nrow, 
			pastix_int_t    *Ncol, 
			pastix_int_t    *Nnzero, 
			pastix_int_t   **col, 
			pastix_int_t   **row,
			pastix_int_t   **loc2glob,
			pastix_float_t **val,
			pastix_float_t **rhs, 
			char           **Type, 
			char           **RhsType,
			MPI_Comm         pastix_comm);

