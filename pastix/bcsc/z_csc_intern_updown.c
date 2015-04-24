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

  Build UpDownVector from user right-hand-side and CSCd.
  Retrieve soltion from UpDownVector.
  Construct UpDownVector such as X[i] = 1, or X[i] = i.

*/
#include "common.h"
#include <pthread.h>
#include "z_tools.h"
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

#ifdef TYPE_COMPLEX
#define SMX_SOL (1.0+0.0*I)
#else
#define SMX_SOL 1.0
#endif

/*
  Function: z_CscdUpdownRhs

  Fill-in UpDownVector structure from user right-hand-side member.

  Parameters:
    updovct - UpDownVector structure to fill-in.
    solvmtx - Solver matrix.
    rhs     - Right-hand-side member.
    invp    - reverse permutation tabular.
    dof     - Number of degree of freedom.
 */
void z_CscUpdownRhs(UpDownVector           *updovct,
                    const SolverMatrix     *solvmtx,
                    const void               *rhs,
                    const pastix_int_t       *invp,
                    int                       dof)
{
  pastix_complex64_t *rhsptr = rhs;
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
                rhsptr[indice + i*updovct->gnodenbr];
            }
          itersm2x++;
        }
    }

  print_debug(DBG_CSC_LOG, "<- z_CscUpdownRhs \n");
}

/*
  Function: z_CscdUpdownRhs

  Fill-in UpDownVector structure from user distributed right-hand-side member.

  Parameters:
    updovct - UpDownVector structure to fill-in.
    solvmtx - Solver matrix.
    rhs     - Right-hand-side member.
    invp    - reverse permutation tabular.
    g2l     - local numbers of global nodes, if not local contains -owner
    dof     - Number of degree of freedom.
 */
void z_CscdUpdownRhs(UpDownVector           *updovct,
                     const SolverMatrix     *solvmtx,
                     const void               *rhs,
                     const pastix_int_t       *invp,
                     const pastix_int_t       *g2l,
                     const pastix_int_t        ln,
                     int                       dof)
{
  pastix_complex64_t *rhsptr = rhs;
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
                rhsptr[indice + i*ln];
            }
          itersm2x++;
        }
    }

  print_debug(DBG_CSC_LOG, "<- z_CscdUpdownRhs \n");
}

/*
  Function:z_CscRhsUpdown

  Builds solution from UpDownVector structure

  Parameters:
    updovct  - UpDownVector structure containing the solution.
    solvmtx  - Solver matrix structure.
    rhs      - Solution to fill.
    ncol     - Number of columns in local matrix.
    dof      - Number of degree of freedom.
    comm     - MPI communicator.

 */
void z_CscRhsUpdown(const UpDownVector *updovct,
                    const SolverMatrix *solvmtx,
                    void                 *rhs,
                    const pastix_int_t    ncol,
                    const pastix_int_t   *invp,
                    const int             dof,
                    const int             rhsmaking,
                    MPI_Comm              comm)
{
  pastix_complex64_t *rhsptr = rhs;
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
      rhsptr[iter]  = 0.0;
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
  UpDownVector structure

  Parameters:
    updovct  - UpDownVector structure containing the solution.
    solvmtx  - Solver matrix structure.
    x        - Solution to fill.
    ncol     - Number of columns in local matrix.
    g2l      - local numbers of global nodes, if not local contains -owner
    ord      - ordering
    dof      - Number of degree of freedom.
    comm     - MPI communicator.

 */
void z_CscdRhsUpdown(const UpDownVector *updovct,
                     const SolverMatrix *solvmtx,
                     void                 *x,
                     const pastix_int_t    ncol,
                     const pastix_int_t   *g2l,
                     const pastix_int_t   *invp,
                     int                   dof,
                     MPI_Comm              comm)
{
  pastix_complex64_t *xptr = x;
  pastix_int_t iter;
  pastix_int_t itercblk;
  pastix_int_t itersm2x;
  pastix_int_t itercol;
  pastix_int_t size = updovct->sm2xnbr*ncol*dof;
  (void)comm;

  print_debug(DBG_CSC_LOG, "-> z_CscdRhsUpdown \n");

  for (iter=0; iter<size; iter++)
    {
      xptr[iter]  = 0.0;
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
              xptr[(g2l[invp[(itercol-itercol%dof)/dof]] - 1)*dof +
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

  Fill-in UpDownVector structure such as the solution of
  the system Ax=b is x_i=1 (API_RHS_1) or x_i=i (API_RHS_I).

  Parameters:
    cscmtx   - internal CSCd matrix.
    updovct  - UpDownVector structure to fill-in.
    solvmtx  - Solver matrix.
    mode     - wanted solution API_RHS_1 or API_RHS_I.
    comm     - MPI communicator.
*/
void z_Csc2updown(const z_CscMatrix    *cscmtx,
                  UpDownVector       *updovct,
                  const SolverMatrix *solvmtx,
                  int                   mode,
                  MPI_Comm              comm)
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
    updovct - UpDownVector structure were to copy B as the first X0 used for raffinement.
    solvmtx - Solver matrix.
    mode    - Rule to construct X0 (API_RHS_0 : X0[i] = 0, API_RHS_1 : X0[i] = 1, API_RHS_I : X0[i] = i).
    comm    - MPI_Communicator.
*/
void z_Csc2updown_X0(UpDownVector           *updovct,
                     /*const*/ SolverMatrix *solvmtx,
                     int                       mode,
                     MPI_Comm                  comm)
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
                    MPI_Comm            pastix_comm)
{
    pastix_int_t           *iparm    = pastix_data->iparm;
    z_SolverMatrix         *solvmatr = &(pastix_data->solvmatr);
    Order                  *ordemesh = pastix_data->ordemesh;
    pastix_int_t            procnum  = pastix_data->procnum;
    pastix_int_t           *invp     = ordemesh->peritab;
    (void)loc2glob;

    /* Rhs taking from b */
    if (iparm[IPARM_RHS_MAKING]==API_RHS_B)
    {
        /* Using b */
        if (solvmatr->updovct.sm2xsze > 0 &&  b==NULL)
        {
            errorPrint("b must be allocated.");
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }
        else
        {
            /* Permuter b avec la permutation inverse */
            if (iparm[IPARM_GRAPHDIST] == API_NO )
            {
                z_CscUpdownRhs(&(solvmatr->updovct),
                             solvmatr,
                             b,
                             invp,
                             (int)iparm[IPARM_DOF_NBR]);
            }
#ifdef PASTIX_DISTRIBUTED
            else
            {
                z_CscdUpdownRhs(&(solvmatr->updovct),
                              solvmatr,
                              b,
                              invp,
                              pastix_data->glob2loc,
                              pastix_data->n2,
                              (int)iparm[IPARM_DOF_NBR]);
            }
#endif /* PASTIX_DISTRIBUTED */
        }
    }
    /* Generate rhs */
    else
    {
        if (procnum == 0)
        {
            if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
            {
                if (iparm[IPARM_RHS_MAKING]==API_RHS_1)
                    fprintf(stdout,GEN_RHS_1);
                else
                    fprintf(stdout,GEN_RHS_I);
            }
        }

        /* In updo step, if we only want to use reffinement.
         Set first solution.
         */
        if ((iparm[IPARM_ONLY_RAFF] == API_YES) && (iparm[IPARM_START_TASK] < API_TASK_REFINE))
        {
            fprintf(stdout,GEN_SOL_0);
            z_Csc2updown_X0(&(solvmatr->updovct),
                          solvmatr,
                          API_RHS_0,
                          pastix_comm);
        }
        else
        {
            z_Csc2updown(&(pastix_data->cscmtx),
                         &(solvmatr->updovct),
                         solvmatr,
                         iparm[IPARM_RHS_MAKING],
                         pastix_comm);
        }
    }
    return PASTIX_SUCCESS;
}
