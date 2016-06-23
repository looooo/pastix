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
  File: z_coefinit.h

  Allocation and initialisation of the coeficient of the z_solver matrix.  
*/
#ifndef Z_COEFINIT_H
#define Z_COEFINIT_H

/* Section: Functions declarations*/

/*
  Function: z_CoefMatrix_Allocate
  
  Allocate matrix coefficients in coeftab and ucoeftab.

  Should be first called with me = -1 to allocated coeftab.
  Then, should be called with me set to thread ID 
  to allocate column blocks coefficients arrays.
  
  Parameters
 
     datacode  - solverMatrix 
     factotype - factorization type (LU, LLT ou LDLT)
     me        - thread number. (-1 for first call, 
                 from main thread. >=0 to allocate column blocks 
		 assigned to each thread.)
 
*/
void z_CoefMatrix_Allocate (z_SopalinParam    *sopar,
			  z_SolverMatrix    *datacode,
			  pthread_mutex_t *mutex,
			  pastix_int_t              factotype, 
			  pastix_int_t              me);

/*
  Function: z_CoefMatrix_Init

  Init coeftab and ucoeftab coefficients.

  Parameters:
     datacode     - solverMatrix 
     barrier      - Barrier used for thread synchronisation.
     me           - Thread ID 
     iparm        - Integer parameters array.
     transcsc     - vecteur transcsc
     sopalin_data - <z_Sopalin_Data_t> structure for NUMA version.
*/
void z_CoefMatrix_Init     (z_SolverMatrix         *datacode, 
			  sopthread_barrier_t  *barrier, 
			  pastix_int_t                   me,
			  pastix_int_t                  *iparm, 
			  pastix_complex64_t               **transcsc, 
			  z_Sopalin_Data_t       *sopalin_data);

/*
  Function: z_CoefMatrix_Free
  
  Free the z_solver matrix coefficient tabular : coeftab and ucoeftab.
  
  WARNING: Call it with one unnique thread. 

  Parameters:
    datacode   - solverMatrix 
    factotype  - factorisation type (<API_FACT>)
    
*/  
void z_CoefMatrix_Free     (z_SopalinParam *sopar,
			  z_SolverMatrix *datacode, 
			  pastix_int_t           factotype);


#endif /* Z_COEFINIT_H */
