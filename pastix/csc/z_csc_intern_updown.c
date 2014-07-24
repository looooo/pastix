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
  File: csc_intern_updown.c

  Build z_UpDownVector from user right-hand-side and CSCd.
  Retrieve soltion from z_UpDownVector.
  Construct z_UpDownVector such as X[i] = 1, or X[i] = i.

*/
#include "common.h"
#include <pthread.h>
#include "tools.h"
#include "order.h"
#include "z_csc.h"
#include "z_updown.h"
#include "z_ftgt.h"
#include "z_updown.h"
#include "queue.h"
#include "bulles.h"
#include "z_solver.h"
#include "z_csc_intern_updown.h"

#ifdef DEBUG_RAFF
#define CSC_LOG
#endif

#ifdef CPLX
#define SMX_SOL (1.0+0.0*I)
#else
#define SMX_SOL 1.0
#endif

/*
  Function: z_CscdUpdownRhs

  Fill-in z_UpDownVector structure from user right-hand-side member.

  Parameters:
    updovct - z_UpDownVector structure to fill-in.
    solvmtx - Solver matrix.
    rhs     - Right-hand-side member.
    invp    - reverse permutation tabular.
    dof     - Number of degree of freedom.
 */
void z_CscUpdownRhs(z_UpDownVector       *updovct,
                  const z_SolverMatrix *solvmtx,
                  const pastix_complex64_t        *rhs,
                  const pastix_int_t          *invp,
                  int                 dof)
{
  pastix_int_t itercblk;
  pastix_int_t itercol;
  pastix_int_t itersm2x;
  pastix_int_t indice;
  pastix_int_t i;

  print_debug(DBG_CSC_LOG, "-> z_CscUpdownRhs \n");

  for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
    {
      itersm2x = updovct->cblktab[itercblk].sm2xind;
      for (itercol=solvmtx->cblktab[itercblk].fcolnum;
           itercol<solvmtx->cblktab[itercblk].lcolnum+1;
           itercol++)
        {
          indice = invp[(itercol-itercol%dof)/dof] *dof + itercol%dof;
          for (i=0; i<updovct->sm2xnbr; i++)
            {
              updovct->sm2xtab[itersm2x + i*updovct->sm2xsze] =
                rhs[indice + i*updovct->gnodenbr];
            }
          itersm2x++;
        }
    }

  print_debug(DBG_CSC_LOG, "<- z_CscUpdownRhs \n");
}

/*
  Function: z_CscdUpdownRhs

  Fill-in z_UpDownVector structure from user distributed right-hand-side member.

  Parameters:
    updovct - z_UpDownVector structure to fill-in.
    solvmtx - Solver matrix.
    rhs     - Right-hand-side member.
    invp    - reverse permutation tabular.
    g2l     - local numbers of global nodes, if not local contains -owner
    dof     - Number of degree of freedom.
 */
void z_CscdUpdownRhs(z_UpDownVector       *updovct,
                   const z_SolverMatrix *solvmtx,
                   const pastix_complex64_t        *rhs,
                   const pastix_int_t          *invp,
                   const pastix_int_t          *g2l,
                   const pastix_int_t           ln,
                   int                 dof)
{
  pastix_int_t itercblk;
  pastix_int_t itercol;
  pastix_int_t itersm2x;
  pastix_int_t indice;
  pastix_int_t i;

  print_debug(DBG_CSC_LOG, "-> z_CscdUpdownRhs \n");

  for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
    {
      itersm2x = updovct->cblktab[itercblk].sm2xind;
      for (itercol=solvmtx->cblktab[itercblk].fcolnum;
           itercol<solvmtx->cblktab[itercblk].lcolnum+1;
           itercol++)
        {
          /* No problem, Column is local */
          indice = (g2l[invp[(itercol - itercol%dof)/dof]]-1)*dof + itercol%dof;
          for (i=0; i<updovct->sm2xnbr; i++)
            {
              updovct->sm2xtab[itersm2x + i*updovct->sm2xsze] =
                rhs[indice + i*ln];
            }
          itersm2x++;
        }
    }

  print_debug(DBG_CSC_LOG, "<- z_CscdUpdownRhs \n");
}

/*
  Function:z_CscRhsUpdown

  Builds solution from z_UpDownVector structure

  Parameters:
    updovct  - z_UpDownVector structure containing the solution.
    solvmtx  - Solver matrix structure.
    rhs      - Solution to fill.
    ncol     - Number of columns in local matrix.
    dof      - Number of degree of freedom.
    comm     - MPI communicator.

 */
void z_CscRhsUpdown(const z_UpDownVector *updovct,
                  const z_SolverMatrix *solvmtx,
                  pastix_complex64_t              *rhs,
                  const pastix_int_t           ncol,
                  const pastix_int_t          *invp,
                  const int           dof,
                  const int           rhsmaking,
                  MPI_Comm            comm)
{
  pastix_int_t    iter;
  pastix_int_t    itercblk;
  pastix_int_t    itersm2x;
  pastix_int_t    itercol;
  pastix_int_t    indice, i;
  pastix_int_t    size = updovct->sm2xnbr*ncol*dof;
  pastix_complex64_t *rhs2 = NULL;
  (void)comm;

#ifdef INOUT_ALLREDUCE
  rhs2 = rhs;
#else
  MALLOC_INTERN(rhs2, size, pastix_complex64_t);
#endif

  print_debug(DBG_CSC_LOG, "-> z_CscRhsUpdown \n");

  for (iter=0; iter<size; iter++)
    {
      rhs2[iter] = 0.0;
#ifndef INOUT_ALLREDUCE
      rhs[iter]  = 0.0;
#endif
    }

  if (rhsmaking == API_RHS_B)
    {
      for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
        {
          itersm2x = updovct->cblktab[itercblk].sm2xind;
          for (itercol=solvmtx->cblktab[itercblk].fcolnum;
               itercol<solvmtx->cblktab[itercblk].lcolnum+1;
               itercol++)
            {
              indice = invp[(itercol-itercol%dof)/dof] *dof + itercol%dof;
              for (i=0; i<updovct->sm2xnbr; i++)
                {
                  rhs2[indice + i*updovct->gnodenbr] =
                    updovct->sm2xtab[itersm2x + i*updovct->sm2xsze];
                }
              itersm2x++;
            }
        }
    }
  else /* Rhs making */
    {
      for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
        {
          itersm2x = updovct->cblktab[itercblk].sm2xind;
          for (itercol=solvmtx->cblktab[itercblk].fcolnum;
               itercol<solvmtx->cblktab[itercblk].lcolnum+1;
               itercol++)
            {
              for (i=0; i<updovct->sm2xnbr; i++)
                {
                  rhs2[itercol + i*updovct->gnodenbr] =
                    updovct->sm2xtab[itersm2x + i*updovct->sm2xsze];
                }
              itersm2x++;
            }
        }
    }
  MPI_Allreduce((void *) rhs2, (void *) rhs, size,
                COMM_FLOAT, COMM_SUM, comm);

#ifndef INOUT_ALLREDUCE
  memFree_null(rhs2);
#endif

  print_debug(DBG_CSC_LOG, "<- z_CscRhsUpdown \n");
}

/*
  Function:z_CscdRhsUpdown

  Builds distributed solution from
  z_UpDownVector structure

  Parameters:
    updovct  - z_UpDownVector structure containing the solution.
    solvmtx  - Solver matrix structure.
    x        - Solution to fill.
    ncol     - Number of columns in local matrix.
    g2l      - local numbers of global nodes, if not local contains -owner
    ord      - ordering
    dof      - Number of degree of freedom.
    comm     - MPI communicator.

 */
void z_CscdRhsUpdown(const z_UpDownVector *updovct,
                   const z_SolverMatrix *solvmtx,
                   pastix_complex64_t              *x,
                   const pastix_int_t           ncol,
                   const pastix_int_t          *g2l,
                   const pastix_int_t          *invp,
                   int                 dof,
                   MPI_Comm            comm)
{
  pastix_int_t iter;
  pastix_int_t itercblk;
  pastix_int_t itersm2x;
  pastix_int_t itercol;
  pastix_int_t size = updovct->sm2xnbr*ncol*dof;
  (void)comm;

  print_debug(DBG_CSC_LOG, "-> z_CscdRhsUpdown \n");

  for (iter=0; iter<size; iter++)
    {
      x[iter]  = 0.0;
    }

  for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
    {
      itersm2x = updovct->cblktab[itercblk].sm2xind;
      for (itercol=solvmtx->cblktab[itercblk].fcolnum;
           itercol<solvmtx->cblktab[itercblk].lcolnum+1;
           itercol++)
        {
          pastix_int_t i;
          for (i=0; i<updovct->sm2xnbr; i++)
            {
              x[(g2l[invp[(itercol-itercol%dof)/dof]] - 1)*dof +
                itercol%dof + i*ncol] = updovct->sm2xtab[itersm2x+i*
                                                         updovct->sm2xsze];
            }
          itersm2x++;
        }
    }

  print_debug(DBG_CSC_LOG, "<- z_CscdRhsUpdown \n");
}

/*
  Function: z_Csc2updown

  Fill-in z_UpDownVector structure such as the solution of
  the system Ax=b is x_i=1 (API_RHS_1) or x_i=i (API_RHS_I).

  Parameters:
    cscmtx   - internal CSCd matrix.
    updovct  - z_UpDownVector structure to fill-in.
    solvmtx  - Solver matrix.
    mode     - wanted solution API_RHS_1 or API_RHS_I.
    comm     - MPI communicator.
*/
void z_Csc2updown(const z_CscMatrix    *cscmtx,
                z_UpDownVector       *updovct,
                const z_SolverMatrix *solvmtx,
                int                 mode,
                MPI_Comm            comm)
{
  pastix_int_t    itercblk;
  pastix_int_t    itercol;
  pastix_int_t    itertempy;
  pastix_int_t    iterval;
  pastix_int_t    itersmx;
  pastix_int_t    cblknbr;
  pastix_complex64_t *smb   = NULL;
  pastix_complex64_t *tempy = NULL;
  (void)comm;

  print_debug(DBG_CSC_LOG, "-> z_Csc2updown \n");

  cblknbr = solvmtx->cblknbr;

  MALLOC_INTERN(tempy, updovct->gnodenbr, pastix_complex64_t);
#ifdef INOUT_ALLREDUCE
  smb = tempy;
#else
  MALLOC_INTERN(smb,   updovct->gnodenbr, pastix_complex64_t);
#endif

  print_debug(DBG_CSC_LOG, "nodenbr=%ld\n",(long)updovct->gnodenbr);

  for (itersmx=0; itersmx<updovct->sm2xnbr; itersmx++)
    {
      for (itertempy=0; itertempy<updovct->gnodenbr; itertempy++)
        {
          tempy[itertempy] = 0.0;
#ifndef INOUT_ALLREDUCE
          smb[itertempy]   = 0.0;
#endif
        }

      for (itercblk=0; itercblk < cblknbr; itercblk++)
        {
          pastix_int_t colnbr   = solvmtx->cblktab[itercblk].lcolnum - solvmtx->cblktab[itercblk].fcolnum +1;

          for (itercol=0; itercol < colnbr; itercol++)
            {
              pastix_int_t colvalidx  = CSC_COL(cscmtx,itercblk,itercol);
              pastix_int_t ncolvalidx = CSC_COL(cscmtx,itercblk,itercol+1);

              for (iterval=colvalidx; iterval<ncolvalidx; iterval++)
                {
                  switch (mode)
                    {
                    case API_RHS_1:
                      tempy[CSC_ROW(cscmtx,iterval)] += (pastix_complex64_t)(itersmx+SMX_SOL)*CSC_VAL(cscmtx,iterval);
                      break;
                    case API_RHS_I:
                      tempy[CSC_ROW(cscmtx,iterval)] +=
                        ((pastix_complex64_t)(solvmtx->cblktab[itercblk].fcolnum+itercol))*
                        CSC_VAL(cscmtx,iterval);
                      break;
                    }
                }
            }
        }

      MPI_Allreduce((void *) tempy, (void *) smb, updovct->gnodenbr,
                    COMM_FLOAT, COMM_SUM, comm);

      for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
        {
          pastix_int_t iterdval = updovct->cblktab[itercblk].sm2xind+itersmx*updovct->sm2xsze;

          for (iterval=0; iterval<CSC_COLNBR(cscmtx,itercblk); iterval++)
            {
              updovct->sm2xtab[iterdval+iterval] = smb[solvmtx->cblktab[itercblk].fcolnum+iterval];
            }
        }
    }

  memFree_null(tempy);
#ifndef INOUT_ALLREDUCE
  memFree_null(smb);
#endif

  print_debug(DBG_CSC_LOG, "<- z_Csc2updown \n");
}


/*
  Function: z_Csc2updown_X0

  Fill-in initial X0 for reffinement if we don't want to use
  Solve step.

  (iparm[IPARM_ONLY_RAFF] == API_YES)

  Parameters:
    updovct - z_UpDownVector structure were to copy B as the first X0 used for raffinement.
    solvmtx - Solver matrix.
    mode    - Rule to construct X0 (API_RHS_0 : X0[i] = 0, API_RHS_1 : X0[i] = 1, API_RHS_I : X0[i] = i).
    comm    - MPI_Communicator.
*/
void z_Csc2updown_X0(z_UpDownVector *updovct,
                   /*const*/ z_SolverMatrix *solvmtx,
                   int mode,
                   MPI_Comm comm)
{
  pastix_int_t  itercblk;
  pastix_int_t  iterval;
  pastix_int_t  itersmx;
  pastix_int_t  cblknbr = solvmtx->cblknbr;
  (void)comm;

  print_debug(DBG_CSC_LOG, "-> z_Csc2updown_X0 \n");
  print_debug(DBG_CSC_LOG, "nodenbr=%ld\n",(long)updovct->gnodenbr);

  for (itersmx=0; itersmx<updovct->sm2xnbr; itersmx++)
    {
      for (itercblk=0; itercblk < cblknbr; itercblk++)
        {
          pastix_int_t colnbr   = solvmtx->cblktab[itercblk].lcolnum - solvmtx->cblktab[itercblk].fcolnum +1;
          pastix_int_t iterdval = updovct->cblktab[itercblk].sm2xind+itersmx*updovct->sm2xsze;

          for (iterval=0; iterval < colnbr; iterval++)
            {
              switch (mode)
                {
                case API_RHS_0:
                  updovct->sm2xtab[iterdval+iterval] = (pastix_complex64_t)0.0;
                  break;
                case API_RHS_1:
                  updovct->sm2xtab[iterdval+iterval] = (pastix_complex64_t)SMX_SOL;
                  break;
                case API_RHS_I:
                  updovct->sm2xtab[iterdval+iterval] = ((pastix_complex64_t)(solvmtx->cblktab[itercblk].fcolnum+iterval));
                  break;

                }
            }
        }
    }

  print_debug(DBG_CSC_LOG, "<- z_Csc2updown_X0 \n");
}
