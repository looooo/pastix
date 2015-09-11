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
  File: z_csc_intern_updown.h

  Build UpDownVector from user right-hand-side and CSCd.
  Retrieve soltion from UpDownVector.
  Construct UpDownVector such as X[i] = 1, or X[i] = i.

*/
#ifndef Z_CSC_INTERN_UPDOWN_H
#define Z_CSC_INTERN_UPDOWN_H

/*
  Function: z_CscdUpdownRhs

  Fill-in UpDownVector structure from user right-hand-side member.

  Parameters:
    updovct - UpDownVector structure to fill-in.
    symbmtx - Solver matrix.
    rhs     - Right-hand-side member.
    perm    - reverse permutation tabular.
    dof      - Number of degree of freedom.
 */
void z_CscUpdownRhs(UpDownVector       *updovct,
		  const SolverMatrix *symbmtx, 
		  const pastix_complex64_t        *rhs, 
		  const pastix_int_t          *perm,
		  int                 dof);

/*
  Function: z_CscdUpdownRhs

  Fill-in UpDownVector structure from user distributed right-hand-side member.

  Parameters:
    updovct - UpDownVector structure to fill-in.
    symbmtx - Solver matrix.
    rhs     - Right-hand-side member.
    invp    - reverse permutation tabular.
    g2l     - local numbers of global nodes, if not local contains -owner
    dof      - Number of degree of freedom.
 */
void z_CscdUpdownRhs(UpDownVector       *updovct,
		   const SolverMatrix *symbmtx, 
		   const pastix_complex64_t        *rhs, 
		   const pastix_int_t          *invp,
		   const pastix_int_t          *g2l,
		   const pastix_int_t           ln,
		   int                 dof);

/*
  Function:z_CscdRhsUpdown

  Builds solution from UpDownVector structure

  Parameters:
    updovct  - UpDownVector structure containing the solution.
    symbmtx  - Solver matrix structure.
    rhs      - Solution to fill.
    ncol     - Number of columns in local matrix.
    dof      - Number of degree of freedom.
    comm     - MPI communicator.
  
 */
void z_CscRhsUpdown(const UpDownVector *updovct, 
		  const SolverMatrix *symbmtx, 
		  pastix_complex64_t              *rhs, 
		  const pastix_int_t           ncol,
		  const pastix_int_t          *invp,
		  const int           dof, 
		  const int           rhsmaking, 
		  MPI_Comm            comm);

/*
  Function:z_CscdRhsUpdown

  Builds distributed solution from
  UpDownVector structure

  Parameters:
    updovct  - UpDownVector structure containing the solution.
    symbmtx  - Solver matrix structure.
    x        - Solution to fill.
    ncol     - Number of columns in local matrix.
    g2l      - local numbers of global nodes, if not local contains -owner
    ord      - ordering
    dof      - Number of degree of freedom.
    comm     - MPI communicator.
  
 */
void z_CscdRhsUpdown(const UpDownVector *updovct, 
		   const SolverMatrix *symbmtx, 
		   pastix_complex64_t              *x,
		   const pastix_int_t           ncol, 
		   const pastix_int_t          *g2l,
		   const pastix_int_t          *invp,
		   int                 dof,  
		   MPI_Comm            comm);

/*
  Function: z_Csc2updown

  Fill-in UpDownVector structure such as the solution of
  the system Ax=b is x_i=1 (API_RHS_1) or x_i=i (API_RHS_I). 

  Parameters:
    cscmtx   - internal CSCd matrix.
    updovct  - UpDownVector structure to fill-in.
    symbmtx  - Solver matrix.
    mode     - wanted solution API_RHS_1 or API_RHS_I.
    comm     - MPI communicator.
*/
void z_Csc2updown(const z_CscMatrix    *cscmtx, 
		UpDownVector       *updovct,
		const SolverMatrix *symbmtx, 
		int                 mode,
		MPI_Comm            comm);


/*
  Function: z_Csc2updown_X0

  Fill-in initial X0 for refinement if we don't want to use
  Solve step.

  (iparm[IPARM_ONLY_RAFF] == API_YES)

  Parameters:
    updovct - UpDownVector structure were to copy B as the first X0 used for raffinement.
    symbmtx - Solver matrix.
    mode    - Rule to construct X0 (API_RHS_1 : X0[i] = 1, API_RHS_I : X0[i] = i).
    comm    - MPI_Communicator.
*/
void z_Csc2updown_X0(UpDownVector *updovct, 
		   /*const*/ SolverMatrix *symbmtx, 
		   int mode, 
		   MPI_Comm comm);


/*
 Function: z_buildUpDoVect

 Build z_UpDownVector from user B vector or
 computes its to obtain $$ \forall i X[i] = 1 $$ or $$ \forall i X[i] = i $$
 as the solution. (depending on iparm)

 Parameters:
 pastix_data - PaStiX global data structure.
 loc2glob2   - Global  column number of local columns.
 b           - User right-hand-side member.
 pastix_comm - MPI communicator.
 */
int z_buildUpdoVect(z_pastix_data_t    *pastix_data,
                    pastix_int_t       *loc2glob,
                    pastix_complex64_t *b,
                    MPI_Comm            pastix_comm);
#endif /* Z_CSC_INTERN_UPDOWN_H */
