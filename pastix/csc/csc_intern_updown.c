/*
  File: csc_intern_updown.c

  Build d_UpDownVector from user right-hand-side and CSCd.
  Retrieve soltion from d_UpDownVector.
  Construct d_UpDownVector such as X[i] = 1, or X[i] = i.

*/
#include "common.h"
#include <pthread.h>
#include "tools.h"
#include "order.h"
#include "csc.h"
#include "d_updown.h"
#include "d_ftgt.h"
#include "d_updown.h"
#include "queue.h"
#include "bulles.h"
#include "d_solver.h"
#include "csc_intern_updown.h"

#ifdef DEBUG_RAFF
#define CSC_LOG
#endif

#ifdef CPLX
#define SMX_SOL (1.0+0.0*I)
#else
#define SMX_SOL 1.0
#endif

/*
  Function: CscdUpdownRhs

  Fill-in d_UpDownVector structure from user right-hand-side member.

  Parameters:
    updovct - d_UpDownVector structure to fill-in.
    solvmtx - Solver matrix.
    rhs     - Right-hand-side member.
    invp    - reverse permutation tabular.
    dof     - Number of degree of freedom.
 */
void CscUpdownRhs(d_UpDownVector       *updovct,
                  const d_SolverMatrix *solvmtx,
                  const pastix_float_t        *rhs,
                  const pastix_int_t          *invp,
                  int                 dof)
{
  pastix_int_t itercblk;
  pastix_int_t itercol;
  pastix_int_t itersm2x;
  pastix_int_t indice;
  pastix_int_t i;

  print_debug(DBG_CSC_LOG, "-> CscUpdownRhs \n");

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

  print_debug(DBG_CSC_LOG, "<- CscUpdownRhs \n");
}

/*
  Function: CscdUpdownRhs

  Fill-in d_UpDownVector structure from user distributed right-hand-side member.

  Parameters:
    updovct - d_UpDownVector structure to fill-in.
    solvmtx - Solver matrix.
    rhs     - Right-hand-side member.
    invp    - reverse permutation tabular.
    g2l     - local numbers of global nodes, if not local contains -owner
    dof     - Number of degree of freedom.
 */
void CscdUpdownRhs(d_UpDownVector       *updovct,
                   const d_SolverMatrix *solvmtx,
                   const pastix_float_t        *rhs,
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

  print_debug(DBG_CSC_LOG, "-> CscdUpdownRhs \n");

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

  print_debug(DBG_CSC_LOG, "<- CscdUpdownRhs \n");
}

/*
  Function:CscRhsUpdown

  Builds solution from d_UpDownVector structure

  Parameters:
    updovct  - d_UpDownVector structure containing the solution.
    solvmtx  - Solver matrix structure.
    rhs      - Solution to fill.
    ncol     - Number of columns in local matrix.
    dof      - Number of degree of freedom.
    comm     - MPI communicator.

 */
void CscRhsUpdown(const d_UpDownVector *updovct,
                  const d_SolverMatrix *solvmtx,
                  pastix_float_t              *rhs,
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
  pastix_float_t *rhs2 = NULL;
  (void)comm;

#ifdef INOUT_ALLREDUCE
  rhs2 = rhs;
#else
  MALLOC_INTERN(rhs2, size, pastix_float_t);
#endif

  print_debug(DBG_CSC_LOG, "-> CscRhsUpdown \n");

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

  print_debug(DBG_CSC_LOG, "<- CscRhsUpdown \n");
}

/*
  Function:CscdRhsUpdown

  Builds distributed solution from
  d_UpDownVector structure

  Parameters:
    updovct  - d_UpDownVector structure containing the solution.
    solvmtx  - Solver matrix structure.
    x        - Solution to fill.
    ncol     - Number of columns in local matrix.
    g2l      - local numbers of global nodes, if not local contains -owner
    ord      - ordering
    dof      - Number of degree of freedom.
    comm     - MPI communicator.

 */
void CscdRhsUpdown(const d_UpDownVector *updovct,
                   const d_SolverMatrix *solvmtx,
                   pastix_float_t              *x,
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

  print_debug(DBG_CSC_LOG, "-> CscdRhsUpdown \n");

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

  print_debug(DBG_CSC_LOG, "<- CscdRhsUpdown \n");
}

/*
  Function: Csc2updown

  Fill-in d_UpDownVector structure such as the solution of
  the system Ax=b is x_i=1 (API_RHS_1) or x_i=i (API_RHS_I).

  Parameters:
    cscmtx   - internal CSCd matrix.
    updovct  - d_UpDownVector structure to fill-in.
    solvmtx  - Solver matrix.
    mode     - wanted solution API_RHS_1 or API_RHS_I.
    comm     - MPI communicator.
*/
void Csc2updown(const CscMatrix    *cscmtx,
                d_UpDownVector       *updovct,
                const d_SolverMatrix *solvmtx,
                int                 mode,
                MPI_Comm            comm)
{
  pastix_int_t    itercblk;
  pastix_int_t    itercol;
  pastix_int_t    itertempy;
  pastix_int_t    iterval;
  pastix_int_t    itersmx;
  pastix_int_t    cblknbr;
  pastix_float_t *smb   = NULL;
  pastix_float_t *tempy = NULL;
  (void)comm;

  print_debug(DBG_CSC_LOG, "-> Csc2updown \n");

  cblknbr = solvmtx->cblknbr;

  MALLOC_INTERN(tempy, updovct->gnodenbr, pastix_float_t);
#ifdef INOUT_ALLREDUCE
  smb = tempy;
#else
  MALLOC_INTERN(smb,   updovct->gnodenbr, pastix_float_t);
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
                      tempy[CSC_ROW(cscmtx,iterval)] += (pastix_float_t)(itersmx+SMX_SOL)*CSC_VAL(cscmtx,iterval);
                      break;
                    case API_RHS_I:
                      tempy[CSC_ROW(cscmtx,iterval)] +=
                        ((pastix_float_t)(solvmtx->cblktab[itercblk].fcolnum+itercol))*
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

  print_debug(DBG_CSC_LOG, "<- Csc2updown \n");
}


/*
  Function: Csc2updown_X0

  Fill-in initial X0 for reffinement if we don't want to use
  Solve step.

  (iparm[IPARM_ONLY_RAFF] == API_YES)

  Parameters:
    updovct - d_UpDownVector structure were to copy B as the first X0 used for raffinement.
    solvmtx - Solver matrix.
    mode    - Rule to construct X0 (API_RHS_0 : X0[i] = 0, API_RHS_1 : X0[i] = 1, API_RHS_I : X0[i] = i).
    comm    - MPI_Communicator.
*/
void Csc2updown_X0(d_UpDownVector *updovct,
                   /*const*/ d_SolverMatrix *solvmtx,
                   int mode,
                   MPI_Comm comm)
{
  pastix_int_t  itercblk;
  pastix_int_t  iterval;
  pastix_int_t  itersmx;
  pastix_int_t  cblknbr = solvmtx->cblknbr;
  (void)comm;

  print_debug(DBG_CSC_LOG, "-> Csc2updown_X0 \n");
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
                  updovct->sm2xtab[iterdval+iterval] = (pastix_float_t)0.0;
                  break;
                case API_RHS_1:
                  updovct->sm2xtab[iterdval+iterval] = (pastix_float_t)SMX_SOL;
                  break;
                case API_RHS_I:
                  updovct->sm2xtab[iterdval+iterval] = ((pastix_float_t)(solvmtx->cblktab[itercblk].fcolnum+iterval));
                  break;

                }
            }
        }
    }

  print_debug(DBG_CSC_LOG, "<- Csc2updown_X0 \n");
}
