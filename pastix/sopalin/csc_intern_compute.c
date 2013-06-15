/*
  File: csc_intern_compute.c

  Functions computing operations on the CSC.

*/
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#ifdef FORCE_NOMPI
#else
#include <mpi.h>
#endif

#include "common.h"
#include "tools.h"
#include "order.h"
#include "csc.h"
#include "ftgt.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"

#include "sopalin_define.h"
#include "sopalin_time.h"
#include "sopalin_thread.h"
#include "stack.h"
#include "sopalin3d.h"
#include "sopalin_acces.h"
#include "sopalin_compute.h"
#include "csc_intern_compute.h"


/*
  Variables:

  Integer: iun

  Integer constant equal to one.

  Float: fun

  Floating point constant equal to one.
*/
#ifdef SMP_RAFF
static pastix_int_t   iun   = 1;
static pastix_float_t fun   = 1.0;
#endif

#ifdef DEBUG_RAFF
#define CSC_LOG
#endif
#define SYNCHRO_THREAD SYNCHRO_X_THREAD(SOLV_THRDNBR, sopalin_data->barrier)

/* gradient complex sans CONJ */
/* Section: Macros */
/*
  Macro: CONJ_JJP

  Computes the conjugate of a coefficient.
*/

#ifndef HERMITIAN
/* JJP except in hermitian */
#  define JJP
#endif
#ifdef JJP
#define CONJ_JJP(x) (x)
#else
#define CONJ_JJP(x) CONJ_FLOAT(x)
#endif
/* Section: Functions */
/*
  Function: CscNorm1

  Computes the norm 1 of the internal CSCd

  Norm 1 is equal to the maximum value of the sum of the
  Absolute values of the elements of each columns.

  Parameters:
    cscmtx - Internal CSCd to compute the norm of.
    comm   - MPI communicator.

  Returns:
    The norm 1 of the internal CSCd.
*/
double CscNorm1(const CscMatrix *cscmtx,
                MPI_Comm         comm)
{
  pastix_int_t itercblk;
  pastix_int_t itercol;
  pastix_int_t iterval;
  double themax=0;
#ifndef INOUT_ALLREDUCE
  double themax2;
#endif
  (void)comm;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscNorm1 \n");
#endif

  for (itercblk=0;
       itercblk<CSC_FNBR(cscmtx);
       itercblk++)
    {
      for (itercol=0;
           itercol<CSC_COLNBR(cscmtx, itercblk);
           itercol++)
        {
          pastix_int_t colvalidx = CSC_COL(cscmtx, itercblk, itercol);
          pastix_int_t ncolvalidx = CSC_COL(cscmtx, itercblk, itercol+1);
          double sum41=0;

          for (iterval=colvalidx;
               iterval<ncolvalidx;
               iterval++)
            {
              sum41 += ABS_FLOAT(CSC_VAL(cscmtx, iterval));
            }
          themax = (themax>sum41 ? themax : sum41);
        }
    }

#ifdef INOUT_ALLREDUCE
  MPI_Allreduce(&themax, &themax, 1, MPI_DOUBLE, MPI_MAX, comm);
#else
  MPI_Allreduce(&themax, &themax2, 1, MPI_DOUBLE, MPI_MAX, comm);
#endif

#ifdef CSC_LOG
  fprintf(stdout, "<- CscNorm1 \n");
#endif

#ifdef INOUT_ALLREDUCE
  return themax;
#else
  return themax2;
#endif
}

/*
 * Function: CscGather_X
 *
 * Inlined function to gather a MPI distributed vector.
 *
 *
 */
#define CscGather_X PASTIX_PREFIX_F(CscGather_X)
static inline
int CscGather_X(Sopalin_Data_t        *sopalin_data,
                int                    me,
                const CscMatrix       *cscmtx,
                const UpDownVector    *updovct,
                const SolverMatrix    *datacode,
                MPI_Comm               comm,
                const volatile pastix_float_t * x,
                pastix_float_t                * tmpx,
                pastix_float_t                **tmpx2)
{
  pastix_int_t itercblk;
  pastix_int_t itercol;
  pastix_int_t indcblk, indcol;
#ifdef SMP_RAFF
  pastix_int_t iterttsk;
  pastix_int_t iter;
#endif
  (void)comm;

  /* in transpose mode we need to gather full RHS */
#ifdef SMP_RAFF
  for (iterttsk=0; iterttsk<SOLV_TTSKNBR; iterttsk++)
    {
      itercblk = TASK_CBLKNUM(SOLV_TTSKTAB(iterttsk));
#else /* SMP_RAFF */
      for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
        {
#endif /* SMP_RAFF */
          indcblk = UPDOWN_SM2XIND(itercblk);
          for (itercol=0; itercol<CSC_COLNBR(cscmtx,itercblk); itercol++)
            {
              pastix_int_t indcolglob = SYMB_FCOLNUM(itercblk)+itercol;
              indcol   = indcblk+itercol;
              tmpx[indcolglob] = x[indcol];
            }
#ifdef SMP_RAFF
        }
#else
    }
#endif

#ifdef SMP_RAFF
  /* Reduction sur les threads */
  sopalin_data->ptr_csc[me] = (void *)tmpx;

  SYNCHRO_THREAD;

  MONOTHREAD_BEGIN;

  for (iter = 1; iter < SOLV_THRDNBR; iter ++)
    SOPALIN_AXPY(updovct->gnodenbr, fun,
                 sopalin_data->ptr_csc[iter], iun,
                 tmpx, iun);
#endif /* SMP_RAFF */

  /* Reduction sur les processus MPI */
  MyMPI_Allreduce((void *)tmpx, (void *) (*tmpx2), updovct->gnodenbr,
                  COMM_FLOAT, COMM_SUM, comm);
#ifdef SMP_RAFF
  /* Recopie dans le vecteur tempy2 */
  sopalin_data->ptr_csc[0] = (void *)(*tmpx2);
  MONOTHREAD_END;
  SYNCHRO_THREAD;
  (*tmpx2) = (pastix_float_t *)sopalin_data->ptr_csc[0];
#endif /* SMP_RAFF */
  return PASTIX_SUCCESS;
}
/*
 * Function: CscAtx_norm_thread
 *
 * Inlined function computin ||A^t|| \times ||X||.
 *
 * At the end we have atx containing the
 * "by thread" par of the result, need reduction.
 *
 * Parameters:
 *     sopalin_data - Sopalin_Data_t structure, common to all threads.
 *     me           - thread number.
 *     x            - (global) vector to multiply.
 *     atx          - Result vector.
 *     cscmtx       - Compress Sparse Column matrix *A*
 *     updovct      - x, multiplied vector
 *     datacode     - solver matrix to know tho local structure of the matrix
 *     comm         - MPI communicator
 */
#define CscAtx_norm_thread PASTIX_PREFIX_F(CscAtx_norm_thread)
static inline
int CscAtx_norm_thread(Sopalin_Data_t       *sopalin_data,
                       int                   me,
                       const volatile pastix_float_t *x,
                       pastix_float_t                *atx,
                       const CscMatrix      *cscmtx,
                       const UpDownVector   *updovct,
                       const SolverMatrix   *datacode,
                       MPI_Comm              comm)
{
  pastix_int_t itercblk;
  pastix_int_t itercol;
  pastix_int_t iterval;
  pastix_int_t colstart, colend;
#ifdef SMP_RAFF
  pastix_int_t iterttsk;
#endif

  pastix_float_t * tmpx;
  pastix_float_t * tmpx2 = NULL;
  MALLOC_INTERN(tmpx, updovct->gnodenbr, pastix_float_t);
  memset(tmpx, 0, updovct->gnodenbr*sizeof(pastix_float_t));
  MONOTHREAD_BEGIN;
#ifdef INOUT_ALLREDUCE
  tmpx2 = tmpx;
#else
  MALLOC_INTERN(tmpx2, updovct->gnodenbr, pastix_float_t);
#endif
  MONOTHREAD_END;

  CscGather_X(sopalin_data, me, cscmtx, updovct, datacode,
              comm, x,
              tmpx, &tmpx2);

  /* Compute A x tmpx2 */
#ifdef SMP_RAFF
  for (iterttsk=0; iterttsk<SOLV_TTSKNBR; iterttsk++)
    {
      itercblk = TASK_CBLKNUM(SOLV_TTSKTAB(iterttsk));
#else /* SMP_RAFF */
      for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
        {
#endif /* SMP_RAFF */
          for (itercol=0; itercol<CSC_COLNBR(cscmtx,itercblk); itercol++)
            {
              colstart = CSC_COL(cscmtx,itercblk,itercol);
              colend   = CSC_COL(cscmtx,itercblk,itercol+1);


              for (iterval=colstart; iterval<colend; iterval++)
                {
                  pastix_int_t indcolglob = SYMB_FCOLNUM(itercblk)+itercol;
                  atx[indcolglob] +=
                    ABS_FLOAT(CSC_VAL(cscmtx,iterval)) *
                    ABS_FLOAT(tmpx2[CSC_ROW(cscmtx,iterval)]);
                }
            }
#ifdef SMP_RAFF
        }
#else
    }
#endif
  SYNCHRO_THREAD;
#ifdef INOUT_ALLREDUCE
  MONOTHREAD_BEGIN;
  memFree_null(tmpx2);
  MONOTHREAD_END;
#endif
  memFree_null(tmpx);
  return PASTIX_SUCCESS;
}

/*
 * Function: CscAtx_thread
 *
 * Inlined function computin A^t \times X.
 *
 * At the end we have atx containing the
 * "by thread" par of the result, need reduction.
 *
 * Parameters:
 *     sopalin_data - Sopalin_Data_t structure, common to all threads.
 *     me           - thread number.
 *     x            - (global) vector to multiply.
 *     atx          - Result vector.
 *     cscmtx       - Compress Sparse Column matrix *A*
 *     updovct      - x, multiplied vector
 *     datacode     - solver matrix to know tho local structure of the matrix
 *     comm         - MPI communicator
 */
#define CscAtx_thread PASTIX_PREFIX_F(CscAtx_thread)
static inline
int CscAtx_thread(Sopalin_Data_t       *sopalin_data,
                  int                   me,
                  const volatile pastix_float_t *x,
                  pastix_float_t                *atx,
                  const CscMatrix      *cscmtx,
                  const UpDownVector   *updovct,
                  const SolverMatrix   *datacode,
                  MPI_Comm              comm)
{
  pastix_int_t itercblk;
  pastix_int_t itercol;
  pastix_int_t iterval;
  pastix_int_t colstart, colend;
#ifdef SMP_RAFF
  pastix_int_t iterttsk;
#endif
  pastix_float_t * tmpx;
  pastix_float_t * tmpx2 = NULL;
  MALLOC_INTERN(tmpx, updovct->gnodenbr, pastix_float_t);
  memset(tmpx, 0, updovct->gnodenbr*sizeof(pastix_float_t));
  MONOTHREAD_BEGIN;
#ifdef INOUT_ALLREDUCE
  tmpx2 = tmpx;
#else
  MALLOC_INTERN(tmpx2, updovct->gnodenbr, pastix_float_t);
#endif
  MONOTHREAD_END;

  CscGather_X(sopalin_data, me, cscmtx, updovct, datacode,
              comm, x,
              tmpx, &tmpx2);


  /* Compute A x tmpx2 */
#ifdef SMP_RAFF
  for (iterttsk=0; iterttsk<SOLV_TTSKNBR; iterttsk++)
    {
      itercblk = TASK_CBLKNUM(SOLV_TTSKTAB(iterttsk));
#else /* SMP_RAFF */
      for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
        {
#endif /* SMP_RAFF */
          for (itercol=0; itercol<CSC_COLNBR(cscmtx,itercblk); itercol++)
            {
              colstart = CSC_COL(cscmtx,itercblk,itercol);
              colend   = CSC_COL(cscmtx,itercblk,itercol+1);

              for (iterval=colstart; iterval<colend; iterval++)
                {
                  pastix_int_t indcolglob = SYMB_FCOLNUM(itercblk)+itercol;
                  atx[indcolglob] +=
                    CSC_VAL(cscmtx,iterval) *
                    tmpx2[CSC_ROW(cscmtx,iterval)];
                }
            }
#ifdef SMP_RAFF
        }
#else
    }
#endif
  SYNCHRO_THREAD;
#ifndef INOUT_ALLREDUCE
  MONOTHREAD_BEGIN;
  memFree_null(tmpx2);
  MONOTHREAD_END;
#endif
  memFree_null(tmpx);

  return PASTIX_SUCCESS;
}
/*
 * Function: CscbMAx
 *
 * Computes $$r = b-Ax$$.
 *
 * Parameters:
 *   sopalin_data - Internal structure containing global datas used by PaStiX.
 *   me           - thread ID.
 *   r            - will contains $$b-Ax$$
 *   b            - Vector $$b$$.
 *   cscmtx       - Internal CSCd matrix containing $$A$$.
 *   updovct      - UpDownVector structure containing $$x$$.
 *   solvmtx      - Solver matrix.
 *   comm         - MPI communicator.
 *   transpose    - Indicate if we want to transpose A.
 */
void CscbMAx(Sopalin_Data_t       *sopalin_data,
             int                   me,
             volatile pastix_float_t       *r,
             const volatile pastix_float_t *b,
             const CscMatrix      *cscmtx,
             const UpDownVector   *updovct,
             const SolverMatrix   *solvmtx,
             MPI_Comm              comm,
             pastix_int_t                   transpose)
{
  SolverMatrix * datacode;
  pastix_float_t * tempy  = NULL;
  pastix_float_t * tempy2 = NULL;


  pastix_int_t itertempy = 0;
  pastix_int_t itercblk  = 0;
  pastix_int_t itercol   = 0;
  pastix_int_t iterval   = 0;
  pastix_int_t colstart, colend;
  pastix_int_t itersmx;
  pastix_int_t indcblk, indcol;
#ifdef SMP_RAFF
  pastix_int_t iter;
  pastix_int_t iterttsk;
#endif

#ifdef CSC_LOG
  fprintf(stdout, "-> CscbMAx \n");
#endif
  datacode = sopalin_data->datacode;

  /* Vecteur local a chaque thread */
  MALLOC_INTERN(tempy, updovct->gnodenbr, pastix_float_t);

  MONOTHREAD_BEGIN;
#ifdef INOUT_ALLREDUCE
  tempy2 = tempy;
#else
  MALLOC_INTERN(tempy2, updovct->gnodenbr, pastix_float_t);
#endif
  MONOTHREAD_END;

  for (itersmx=0; itersmx<updovct->sm2xnbr; itersmx++)
    {
      for (itertempy=0; itertempy<updovct->gnodenbr; itertempy++)
        {
          tempy[itertempy] = 0.0;
        }
#ifndef INOUT_ALLREDUCE
      MONOTHREAD_BEGIN;
      for (itertempy=0; itertempy<updovct->gnodenbr; itertempy++)
        {
          tempy2[itertempy] = 0.0;
        }
      MONOTHREAD_END;
#endif

      if (transpose == API_YES)
        {
          CscAtx_thread(sopalin_data, me,
                        &(updovct->sm2xtab[itersmx*updovct->sm2xsze]), tempy,
                        cscmtx, updovct, solvmtx,
                        comm);
        }
      else
        {
          /* Boucle de calcul de Ax local */
#ifdef SMP_RAFF
          for (iterttsk=0; iterttsk<SOLV_TTSKNBR; iterttsk++)
            {
              itercblk = TASK_CBLKNUM(SOLV_TTSKTAB(iterttsk));
#else /* SMP_RAFF */
              for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
                {
#endif /* SMP_RAFF */

                  indcblk = updovct->cblktab[itercblk].sm2xind +
                    itersmx*updovct->sm2xsze;

                  for (itercol=0; itercol<CSC_COLNBR(cscmtx,itercblk);
                       itercol++)
                    {
                      colstart = CSC_COL(cscmtx,itercblk,itercol);
                      colend   = CSC_COL(cscmtx,itercblk,itercol+1);
                      indcol   = indcblk+itercol;

                      for (iterval=colstart; iterval<colend; iterval++)
                        {
                          tempy[CSC_ROW(cscmtx,iterval)] +=
                            CSC_VAL(cscmtx,iterval) * updovct->sm2xtab[indcol];
                        }
                    }
                }
#ifdef _UNUSED_
            }
#endif
        }

      /* Reduction sur les threads */
#ifdef SMP_RAFF
      sopalin_data->ptr_csc[me] = (void *)tempy;
      SYNCHRO_THREAD;
      MONOTHREAD_BEGIN;
      for (iter = 1; iter < SOLV_THRDNBR; iter ++)
        SOPALIN_AXPY(updovct->gnodenbr, fun,
                     sopalin_data->ptr_csc[iter], iun, tempy, iun);
#endif /* SMP_RAFF */

      /* Reduction sur les processus MPI */
      MyMPI_Allreduce((void *)tempy, (void *) tempy2, updovct->gnodenbr,
                      COMM_FLOAT, COMM_SUM, comm);

#ifdef SMP_RAFF
      /* Recopie dans le vecteur tempy2 */
      sopalin_data->ptr_csc[0] = (void *)tempy2;
      MONOTHREAD_END;
      SYNCHRO_THREAD;
      tempy2 = (pastix_float_t *)sopalin_data->ptr_csc[0];

      for (iterttsk=0; iterttsk<SOLV_TTSKNBR; iterttsk++)
        {
          itercblk = TASK_CBLKNUM(SOLV_TTSKTAB(iterttsk));
#else /* SMP_RAFF */
          for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
            {
#endif /* SMP_RAFF */

              pastix_int_t iterdval = updovct->cblktab[itercblk].sm2xind +
                itersmx*updovct->sm2xsze;

              for (iterval=0; iterval<CSC_COLNBR(cscmtx,itercblk); iterval++)
                {
                  /* je ne protege pas les ecritures dans r car elles
                   * se font a des endroits differents pour chaque thread
                   * Par contre les vecteurs sont globaux a tous les threads
                   */
                  r[iterdval+iterval] = b[iterdval+iterval]
                    - tempy2[solvmtx->cblktab[itercblk].fcolnum+iterval];
                }
            }
#ifdef _UNUSED_
        }
#endif

#ifdef SMP_RAFF
      SYNCHRO_THREAD;
#endif
    }/* fin boucle multi-membre */

  memFree_null(tempy);
  MONOTHREAD_BEGIN
#ifndef INOUT_ALLREDUCE
    memFree_null(tempy2);
#endif
  MONOTHREAD_END;

#ifdef CSC_LOG
  fprintf(stdout, "<- CscbMAx \n");
#endif
}


/*
 * function: CscAxPb
 *
 *
 *  Compute the operation $$r=|A||x|+|b|$$
 *
 *  Parameters:
 *     sopalin_data - Sopalin_Data_t structure, common to all threads.
 *     me           - thread number
 *     r            - solution (vector commont to all threads)
 *     b            - Added vector (vector commont to all threads)
 *     cscmtx       - Compress Sparse Column matrix *A*
 *     updovct      - x, multiplied vector
 *     solvmtx      - solver matrix to know tho local structure of the matrix
 *     comm         - MPI communicator
 *     transpose    - Indicate if we want to transpose A.
 */
void CscAxPb(Sopalin_Data_t     *sopalin_data,
             int                 me,
             pastix_float_t              *r,
             const pastix_float_t        *b,
             const CscMatrix    *cscmtx,
             const UpDownVector *updovct,
             const SolverMatrix *solvmtx,
             MPI_Comm            comm,
             pastix_int_t                 transpose)
{
  SolverMatrix * datacode;
  pastix_float_t * tempy  = NULL;
  pastix_float_t * tempy2 = NULL;

  pastix_int_t itertempy = 0;
  pastix_int_t itercblk  = 0;
  pastix_int_t itercol   = 0;
  pastix_int_t iterval   = 0;
  pastix_int_t colstart, colend;
  pastix_int_t itersmx;
  pastix_int_t indcblk, indcol;
#ifdef SMP_RAFF
  pastix_int_t iterttsk;
  pastix_int_t iter;
#endif

#ifdef CSC_LOG
  fprintf(stdout, "-> CscAxPb \n");
#endif

  datacode = sopalin_data->datacode;

  /* vecteurs temporaires locaux */
  MALLOC_INTERN(tempy, updovct->gnodenbr, pastix_float_t);

  MONOTHREAD_BEGIN;
#ifdef INOUT_ALLREDUCE
  tempy2 = tempy;
#else
  MALLOC_INTERN(tempy2, updovct->gnodenbr, pastix_float_t);
#endif
  MONOTHREAD_END;


  /* Boucle sur les seconds membres */
  for (itersmx=0; itersmx<updovct->sm2xnbr; itersmx++)
    {
      /* Mise a zero des vecteurs de stockages du résultat */
      for (itertempy=0; itertempy<updovct->gnodenbr; itertempy++)
        {
          tempy[itertempy] = 0.0;
        }
#ifndef INOUT_ALLREDUCE
      MONOTHREAD_BEGIN;
      for (itertempy=0; itertempy<updovct->gnodenbr; itertempy++)
        {
          tempy2[itertempy] = 0.0;
        }
      MONOTHREAD_END;
#endif
      if (transpose == API_YES)
        {

          CscAtx_norm_thread(sopalin_data, me,
                             &(updovct->sm2xtab[itersmx*updovct->sm2xsze]),
                             tempy,
                             cscmtx, updovct, solvmtx,
                             comm);

        }
      else
        {
          /* Boucle de calcul de Ap local */
#ifdef SMP_RAFF
          for (iterttsk=0; iterttsk<SOLV_TTSKNBR; iterttsk++)
            {
              itercblk = TASK_CBLKNUM(SOLV_TTSKTAB(iterttsk));
#else /* SMP_RAFF */
              for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
                {
#endif /* SMP_RAFF */
                  indcblk = updovct->cblktab[itercblk].sm2xind +
                    itersmx*updovct->sm2xsze;

                  for (itercol=0; itercol<CSC_COLNBR(cscmtx,itercblk);
                       itercol++)
                    {
                      colstart = CSC_COL(cscmtx,itercblk,itercol);
                      colend   = CSC_COL(cscmtx,itercblk,itercol+1);
                      indcol   = indcblk+itercol;

                      for (iterval=colstart; iterval<colend; iterval++)
                        {
                          tempy[CSC_ROW(cscmtx,iterval)] +=
                            ABS_FLOAT(CSC_VAL(cscmtx,iterval)) *
                            ABS_FLOAT(updovct->sm2xtab[indcol]);
                        }
                    }
                }
#ifdef _UNUSED_
            }
#endif
        }
      /* Reduction sur les threads */
#ifdef SMP_RAFF
      /* Reduction sur les threads */
      sopalin_data->ptr_csc[me] = (void *)tempy;

      SYNCHRO_THREAD;

      MONOTHREAD_BEGIN;

      for (iter = 1; iter < SOLV_THRDNBR; iter ++)
        SOPALIN_AXPY(updovct->gnodenbr, fun,
                     sopalin_data->ptr_csc[iter], iun,
                     tempy, iun);
#endif /* SMP_RAFF */

      /* Reduction sur les processus MPI */
      MyMPI_Allreduce((void *)tempy, (void *) tempy2, updovct->gnodenbr,
                      COMM_FLOAT, COMM_SUM, comm);

      /* Recopie dans le vecteur X */
#ifdef SMP_RAFF
      /* Recopie dans le vecteur tempy2 */
      sopalin_data->ptr_csc[0] = (void *)tempy2;
      MONOTHREAD_END;
      SYNCHRO_THREAD;
      tempy2 = (pastix_float_t *)sopalin_data->ptr_csc[0];

      for (iterttsk=0; iterttsk<SOLV_TTSKNBR; iterttsk++)
        {
          itercblk = TASK_CBLKNUM(SOLV_TTSKTAB(iterttsk));
#else /* SMP_RAFF */
          for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
            {
#endif /* SMP_RAFF */

              pastix_int_t iterdval = updovct->cblktab[itercblk].sm2xind +
                itersmx*updovct->sm2xsze;

              for (iterval=0; iterval<CSC_COLNBR(cscmtx, itercblk); iterval++)
                {

                  r[iterdval+iterval] = ABS_FLOAT(b[iterdval+iterval])+
                    tempy2[solvmtx->cblktab[itercblk].fcolnum+iterval];
                }
            }

#ifdef _UNUSED_
        }
#endif

#ifdef SMP_RAFF
      SYNCHRO_THREAD;
#endif
    }/* fin boucle multi-membre */

  memFree_null(tempy);

  MONOTHREAD_BEGIN;
#ifndef INOUT_ALLREDUCE
  memFree_null(tempy2);
#endif
  MONOTHREAD_END;

#ifdef CSC_LOG
  fprintf(stdout, "<- CscAxPb \n");
#endif
}

/*
   Function: CscBerr

   Compute the operation $$ berr= max_{i}(\\frac{|r_{i}|}{|s_{i}|}) $$.

   Parameters :

   sopalin_data - Sopalin_Data_t structure, common to all threads.
   me           - thread number
   r            - vector(s) (common to all threads)
   s            - vector(s) (common to all threads)
   colnbr       - size of the vectors
   smvnbr       - number of vectors in r and s
   berr         - berr (local variable)
   comm         - MPI communicator
*/
void CscBerr(Sopalin_Data_t *sopalin_data,
             int            me,
             const pastix_float_t   *r,
             const pastix_float_t   *s,
             const pastix_int_t      colnbr,
             const pastix_int_t      smxnbr,
             double        *berr,
             MPI_Comm       comm)
{
  SolverMatrix *  datacode;
  pastix_int_t first,  last;
  pastix_int_t first2, last2;
  pastix_int_t step;
  pastix_int_t iter;
  pastix_int_t itersmx;
  double berrmax = 0.0;
  double berr2;
  (void)comm;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscBerr \n");
#endif

  datacode = sopalin_data->datacode;

#ifdef SMP_RAFF
  MONOTHREAD_BEGIN;
  sopalin_data->common_flt[0] = 0.0;
  MONOTHREAD_END;
  SYNCHRO_THREAD;
  step  = (pastix_int_t)ceil((double)(colnbr)/(double)(SOLV_THRDNBR));
#else
  step  = colnbr;
#endif
  first = me * step;
  last  = MIN(colnbr, (me+1) * step);

  for (itersmx=0; itersmx<smxnbr; itersmx++)
    {
      first2 = first + itersmx*colnbr;
      last2  = last  + itersmx*colnbr;
      berr2  = ABS_FLOAT(r[first]/s[first]);

      for (iter=first2+1; iter<last2; iter++)
        {
          if (ABS_FLOAT(r[iter]/s[iter]) > berr2)
            berr2 = ABS_FLOAT(r[iter]/s[iter]);
        }

      /* calcul du max entre les threads */
#ifdef SMP_RAFF

      MUTEX_LOCK(&sopalin_data->mutex_raff);
      if (berr2 > sopalin_data->common_flt[0])
  sopalin_data->common_flt[0] = berr2;
      MUTEX_UNLOCK(&sopalin_data->mutex_raff);

      /* Le thread 0 attend que ttes les contributions soient là */
      SYNCHRO_THREAD;
      MONOTHREAD_BEGIN;
      berr2 = sopalin_data->common_flt[0];
#endif /* SMP_RAFF */

      MyMPI_Allreduce(&berr2, berr, 1, MPI_DOUBLE, MPI_MAX, comm);

      /* On s'assure que les threads on tous la mÃªme valeur */
#ifdef SMP_RAFF
      sopalin_data->common_flt[0] = *berr;
      MONOTHREAD_END;
      SYNCHRO_THREAD;
      *berr = sopalin_data->common_flt[0];
#endif /* SMP_RAFF */

      if (*berr > berrmax)
        berrmax = *berr;
    }

  *berr = berrmax;

#ifdef CSC_LOG
  fprintf(stdout, "<- CscBerr \n");
#endif
}

/*
  Function: CscNormErr

  Computes the norm 2 of r and the norm 2 of b and return the quotient of these
  two vectors.

  This Function is multithreaded, each thread will compute a part of the norm,
  it will be gathered between threads, then between MPI processors.

  Parameters:
    sopalin_data - global PaStix informations.
    me           - Thread ID.
    r            - first vector from which the norm 2 is computed.
    b            - second vector from which the norm 2 is computed.
    colnbr       - Size of the vectors.
    smxnbr       - Number of vectors (multi-right-hand-side method)
    comm         - PaStiX MPI communicator.
*/
double CscNormErr(Sopalin_Data_t       *sopalin_data,
                  int                   me,
                  const volatile pastix_float_t *r,
                  const volatile pastix_float_t *b,
                  const pastix_int_t             colnbr,
                  const pastix_int_t             smxnbr,
                  MPI_Comm              comm)
{
  SolverMatrix *datacode;
  pastix_int_t first,  last;
  pastix_int_t first2, last2;
  pastix_int_t iter;
  pastix_int_t itersmx;
  pastix_int_t step;
  double rnormax = 0.0;
  double bnormax = 1.0;
  double rnorm = 0.0;
  double bnorm = 0.0;
#ifndef INOUT_ALLREDUCE
  double rnorm2 = 0.0;
  double bnorm2 = 0.0;
#endif
  (void)comm;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscNormErr \n");
#endif

  datacode  = sopalin_data->datacode;

#ifdef SMP_RAFF
  step  = (pastix_int_t)ceil((double)(colnbr)/(double)(SOLV_THRDNBR));
#else
  step  = colnbr;
#endif
  first = me * step;
  last  = MIN(colnbr, (me+1) * step);

  for (itersmx=0; itersmx<smxnbr; itersmx++)
    {
      /* Produit scalaire sur les donnÃ©es locales au thread */
      first2 = first + itersmx*colnbr;
      last2  = last  + itersmx*colnbr;
      for (iter=first2; iter<last2; iter++)
        {
#ifdef CPLX
          rnorm += (double)(r[iter]*conj(r[iter]));
          bnorm += (double)(b[iter]*conj(b[iter]));
#else  /* CPLX */
          rnorm += r[iter]*r[iter];
          bnorm += b[iter]*b[iter];
#endif /* CPLX */
        }

      /* En SMP reduction du resultat sur les threads */
#ifdef SMP_RAFF
      sopalin_data->common_flt[me]              = rnorm;
      sopalin_data->common_flt[me+SOLV_THRDNBR] = bnorm;

      SYNCHRO_THREAD;
      MONOTHREAD_BEGIN;
      for (iter = 1; iter < SOLV_THRDNBR; iter++)
        {
          rnorm += sopalin_data->common_flt[iter];
          bnorm += sopalin_data->common_flt[iter+SOLV_THRDNBR];
        }
#endif /* SMP_RAFF */

#ifdef INOUT_ALLREDUCE
      MyMPI_Allreduce(&rnorm, &rnorm,  1, MPI_DOUBLE, MPI_SUM, comm);
      MyMPI_Allreduce(&bnorm, &bnorm,  1, MPI_DOUBLE, MPI_SUM, comm);
#else
      MyMPI_Allreduce(&rnorm, &rnorm2, 1, MPI_DOUBLE, MPI_SUM, comm);
      MyMPI_Allreduce(&bnorm, &bnorm2, 1, MPI_DOUBLE, MPI_SUM, comm);
      rnorm = rnorm2;
      bnorm = bnorm2;
#endif /* INOUT_ALLREDUCE */

      /* on broadast rnorm et bnorm sur tous les threads */
#ifdef SMP_RAFF
      sopalin_data->common_flt[0] = rnorm;
      sopalin_data->common_flt[1] = bnorm;
      MONOTHREAD_END;
      SYNCHRO_THREAD;
      rnorm = sopalin_data->common_flt[0];
      bnorm = sopalin_data->common_flt[1];
#endif /* SMP_RAFF */

      if ((rnorm/bnorm)>(rnormax/bnormax))
        {
          rnormax = rnorm;
          bnormax = bnorm;
        }
    }

  rnorm = rnormax;
  bnorm = bnormax;

#ifdef CSC_LOG
  fprintf(stdout, "<- CscNormErr \n");
#endif

  return sqrt(rnorm/bnorm);
}

/*
 * Function: CscAx
 *
 * Computes *A* times *p* and store the result in *x*.
 *
 * When compiled with SMP_RAFF, this operation is multi-threaded.
 *
 * Parameters:
 *   sopalin_data - Gloabl PaStiX data.
 *   me           - thread ID.
 *   cscmtx       - Internal CSCd matrix, A.
 *   p            - Vector which will be multiplied by the CSCd.
 *   x            - vector which will contain the computation result.
 *   solvmtx      - Solver matrix.
 *   updovct      - Structure used for updown step, it contains information
 *                  about the vectors.
 *   comm         - MPI Communicator.
 *     transpose    - Indicate if we want to transpose A.
 */
void CscAx(Sopalin_Data_t       *sopalin_data,
           int                   me,
           const CscMatrix      *cscmtx,
           const volatile pastix_float_t *p,
           volatile pastix_float_t       *x,
           const SolverMatrix   *solvmtx,
           const UpDownVector   *updovct,
           MPI_Comm              comm,
           pastix_int_t                   transpose)
{
  SolverMatrix * datacode;
  pastix_float_t * tempy  = NULL;
  pastix_float_t * tempy2 = NULL;

  pastix_int_t itertempy = 0;
  pastix_int_t itercblk  = 0;
  pastix_int_t itercol   = 0;
  pastix_int_t iterval   = 0;
  pastix_int_t colstart, colend;
  pastix_int_t itersmx;
  pastix_int_t indcblk, indcol;
#ifdef SMP_RAFF
  pastix_int_t iterttsk;
  pastix_int_t iter;
#endif

#ifdef CSC_LOG
  fprintf(stdout, "-> CscAx \n");
#endif
  datacode = sopalin_data->datacode;

  /* vecteurs temporaires locaux */
  MALLOC_INTERN(tempy, updovct->gnodenbr, pastix_float_t);
  MONOTHREAD_BEGIN;
#ifdef INOUT_ALLREDUCE
  tempy2 = tempy;
#else
  MALLOC_INTERN(tempy2, updovct->gnodenbr, pastix_float_t);
#endif
  MONOTHREAD_END;


  /* Boucle sur les seconds membres */
  for (itersmx=0; itersmx<updovct->sm2xnbr; itersmx++)
    {
      /* Mise a zero des vecteurs de stockages du résultat */
      for (itertempy=0; itertempy<updovct->gnodenbr; itertempy++)
        {
          tempy[itertempy] = 0.0;
        }
#ifndef INOUT_ALLREDUCE
      MONOTHREAD_BEGIN;
      for (itertempy=0; itertempy<updovct->gnodenbr; itertempy++)
        {
          tempy2[itertempy] = 0.0;
        }
      MONOTHREAD_END;
#endif

      if (transpose == API_YES)
        {

          CscAtx_thread(sopalin_data, me,
                        &(p[itersmx*updovct->sm2xsze]), tempy,
                        cscmtx, updovct, solvmtx, comm);
        }
      else
        {
          /* Boucle de calcul de Ap local */
#ifdef SMP_RAFF
          for (iterttsk=0; iterttsk<SOLV_TTSKNBR; iterttsk++)
            {
              itercblk = TASK_CBLKNUM(SOLV_TTSKTAB(iterttsk));
#else /* SMP_RAFF */
              for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
                {
#endif /* SMP_RAFF */

                  indcblk = updovct->cblktab[itercblk].sm2xind +
                    itersmx*updovct->sm2xsze;

                  for (itercol=0; itercol<CSC_COLNBR(cscmtx,itercblk);
                       itercol++)
                    {
                      colstart = CSC_COL(cscmtx,itercblk,itercol);
                      colend   = CSC_COL(cscmtx,itercblk,itercol+1);
                      indcol   = indcblk+itercol;

                      for (iterval=colstart; iterval<colend; iterval++)
                        {
                          tempy[CSC_ROW(cscmtx,iterval)] +=
                            CSC_VAL(cscmtx,iterval) * p[indcol];
                        }
                    }
                }
#ifdef _UNUSED_
            }
#endif
        }


#ifdef DEBUG_RAFF
      {
        FILE *rafffile;
        char  rafffilename[30];
        static pastix_int_t toto = 0;
        sprintf(rafffilename, "CscAxtempy%ld.%ld",(long) toto,me);
        rafffile = fopen(rafffilename, "w");
        dump7((pastix_float_t *)tempy,rafffile,updovct->gnodenbr);
        fclose(rafffile);

        SYNCHRO_THREAD;
        MONOTHREAD_BEGIN;
        toto++;
        MONOTHREAD_END;
      }
#endif
#ifdef SMP_RAFF
      /* Reduction sur les threads */
      sopalin_data->ptr_csc[me] = (void *)tempy;

      SYNCHRO_THREAD;

      MONOTHREAD_BEGIN;

      for (iter = 1; iter < SOLV_THRDNBR; iter ++)
        SOPALIN_AXPY(updovct->gnodenbr, fun,
                     sopalin_data->ptr_csc[iter], iun, tempy, iun);
#endif
      /* Reduction sur les processus MPI */
      MyMPI_Allreduce((void *)tempy, (void *) tempy2, updovct->gnodenbr,
                      COMM_FLOAT, COMM_SUM, comm);

#ifdef DEBUG_RAFF
      {
        FILE *rafffile;
        char  rafffilename[10];
        static pastix_int_t toto = 0;
        sprintf(rafffilename, "CscAx%ld.%ld",(long) toto,(long) SOLV_PROCNUM);
        rafffile = fopen(rafffilename, "w");
        dump7((pastix_float_t *)tempy2,rafffile,updovct->gnodenbr);
        fclose(rafffile);
        toto++;
      }
#endif
#ifdef SMP_RAFF
      /* Recopie dans le vecteur tempy2 */
      sopalin_data->ptr_csc[0] = (void *)tempy2;
      MONOTHREAD_END;
      SYNCHRO_THREAD;
      tempy2 = (pastix_float_t *)sopalin_data->ptr_csc[0];

      for (iterttsk=0; iterttsk<SOLV_TTSKNBR; iterttsk++)
        {
          itercblk = TASK_CBLKNUM(SOLV_TTSKTAB(iterttsk));
#else /* SMP_RAFF */
          for (itercblk=0; itercblk<CSC_FNBR(cscmtx); itercblk++)
            {
#endif /* SMP_RAFF */
              pastix_int_t iterdval = updovct->cblktab[itercblk].sm2xind +
                itersmx*updovct->sm2xsze;

              for (iterval=0; iterval<CSC_COLNBR(cscmtx, itercblk); iterval++)
                {
                  x[iterdval+iterval] =
                    tempy2[solvmtx->cblktab[itercblk].fcolnum+iterval];
                }
            }
#ifdef _UNUSED_
        }
#endif

#ifdef SMP_RAFF
      SYNCHRO_THREAD;
#endif
    }/* fin boucle multi-membre */
  memFree_null(tempy);
  MONOTHREAD_BEGIN;
#ifndef INOUT_ALLREDUCE
  memFree_null(tempy2);
#endif
  MONOTHREAD_END;

#ifdef CSC_LOG
  fprintf(stdout, "<- CscAx \n");
#endif
}

/*
  Function: CscGradAlpha

  Computes the scalar product of *r* with *z*,
  then the scalar product of *x* with *p*
  and finaly store the quotient in *alpha*.

  Multi-threaded in SMP_RAFF mode.

  Parameters:
    sopalin_data - Gloabal PaStiX data structure.
    me           - Thread ID.
    r            - A vector of size *colnbr* times *smxnbr*.
    z            - A vector of size *colnbr* times *smxnbr*.
    x            - A vector of size *colnbr* times *smxnbr*.
    p            - A vector of size *colnbr* times *smxnbr*.
    colnbr       - Number of unkowns.
    smxnbr       - Number of right-hand-side members.
    alpha        - Float which will store the computation result.
    comm         - MPI communicator.
*/
void CscGradAlpha(Sopalin_Data_t       *sopalin_data,
                  int                   me,
                  const volatile pastix_float_t *r,
                  const volatile pastix_float_t *z,
                  const volatile pastix_float_t *x,
                  const volatile pastix_float_t *p,
                  pastix_int_t                   colnbr,
                  pastix_int_t                   smxnbr,
                  double               *alpha,
                  MPI_Comm              comm)
{
  SolverMatrix *datacode;
  pastix_int_t first,  last;
  pastix_int_t step;
  pastix_int_t itersmx;
  pastix_int_t iter     = 0;
  double up    = 0.0;
  double down  = 0.0;
#ifndef INOUT_ALLREDUCE
  double up2   = 0.0;
  double down2 = 0.0;
#endif
  (void)comm;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscGradAlpha \n");
#endif

  datacode  = sopalin_data->datacode;

#ifdef NOSMP_RAFF
  step  = colnbr;
#else
  step  = (pastix_int_t)ceil((double)(colnbr)/(double)(SOLV_THRDNBR));
#endif
  first = me * step;
  last  = MIN(colnbr, (me+1) * step);

  for (itersmx=0; itersmx<smxnbr; itersmx++)
    {
      /* Produit scalaire sur les donnees locales au thread */
      /* up   = <r,z> */
      /* down = <x,p> */
      for (iter=first; iter<last; iter++)
        {
          up   += (double)(r[iter]*CONJ_JJP(z[iter]));
          down += (double)(x[iter]*CONJ_JJP(p[iter]));
        }

      /* En SMP reduction du resultat sur les threads */
#ifdef SMP_RAFF
      sopalin_data->common_flt[me]              = up;
      sopalin_data->common_flt[SOLV_THRDNBR+me] = down;

      SYNCHRO_THREAD;
      MONOTHREAD_BEGIN;
      for (iter = 1; iter < SOLV_THRDNBR; iter++)
        {
          up   += sopalin_data->common_flt[iter];
          down += sopalin_data->common_flt[iter+SOLV_THRDNBR];
        }
#endif /* SMP_RAFF */

#ifdef INOUT_ALLREDUCE
      MyMPI_Allreduce((void*)&up,  (void*)&up,   1, MPI_DOUBLE, MPI_SUM, comm);
      MyMPI_Allreduce((void*)&down,(void*)&down, 1, MPI_DOUBLE, MPI_SUM, comm);
#else
      MyMPI_Allreduce((void*)&up,  (void*)&up2,  1, MPI_DOUBLE, MPI_SUM, comm);
      MyMPI_Allreduce((void*)&down,(void*)&down2,1, MPI_DOUBLE, MPI_SUM, comm);
      up   = up2;
      down = down2;
#endif

      alpha[itersmx] = (up/down);

#ifdef SMP_RAFF
      MONOTHREAD_END;
      SYNCHRO_THREAD;
#endif
      up   = 0;
      down = 0;
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscGradAlpha \n");
#endif
}

/*
  Function: CscGradBeta

  Computes the scalar product between *r* and *z*
  and store the result in *beta*.

  At the end, beta is only on thread 0.

  Parameters:
    sopalin_data - PaStiX data structure.
    me           - Thread ID.
    r            - first vector of size *colnbr* times *smxnbr*.
    z            - second vector of size *colnbr* times *smxnbr*.a
    colnbr       - Number of unknowns.
    smxnbr       - Number of right-hand-side members.
    beta         - Float which will store the solution.
    comm         - MPI communicator.
*/
void CscGradBeta(Sopalin_Data_t       *sopalin_data,
                 int                   me,
                 const volatile pastix_float_t *r,
                 const volatile pastix_float_t *z,
                 pastix_int_t                   colnbr,
                 pastix_int_t                   smxnbr,
                 pastix_float_t                *beta,
                 MPI_Comm              comm)
{
  SolverMatrix *datacode;
  pastix_int_t first,  last;
  pastix_int_t first2, last2;
  pastix_int_t step;
  pastix_int_t itersmx;
  pastix_int_t iter  = 0;
  pastix_float_t up  = 0.0;
#ifndef INOUT_ALLREDUCE
  pastix_float_t up2 = 0.0;
#endif
  (void)comm;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscGradBeta \n");
#endif

  datacode  = sopalin_data->datacode;

#ifdef SMP_RAFF
  step  = (pastix_int_t)ceil((double)(colnbr)/(double)(SOLV_THRDNBR));
#else
  step  = colnbr;
#endif
  first = me * step;
  last  = MIN(colnbr, (me+1) * step);

  for (itersmx=0; itersmx<smxnbr; itersmx++)
    {
      /* Produit scalaire sur les donnÃ©es locales au thread */
      first2 = first + itersmx*colnbr;
      last2  = last  + itersmx*colnbr;
      for (iter=first2; iter<last2; iter++)
        {
          up += (r[iter]*CONJ_JJP(z[iter]));
        }

      /* En SMP reduction du resultat sur les threads */
#ifdef SMP_RAFF
      sopalin_data->common_flt[me] = up;

      SYNCHRO_THREAD;
      MONOTHREAD_BEGIN;
      for (iter = 1; iter < SOLV_THRDNBR; iter++)
        up += sopalin_data->common_flt[iter];
#endif /* SMP_RAFF */

#ifdef INOUT_ALLREDUCE
      MyMPI_Allreduce((void*)&up, (void*)&up,  1, COMM_FLOAT, MPI_SUM, comm);
#else
      MyMPI_Allreduce((void*)&up, (void*)&up2, 1, COMM_FLOAT, MPI_SUM, comm);
      up = up2;
#endif

      beta[itersmx] = up;

#ifdef SMP_RAFF
      MONOTHREAD_END;
      SYNCHRO_THREAD;
#endif
      up=0;
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscGradBeta \n");
#endif
}

/*
  Function: CscGmresBeta

  Computes the scalar product between *r* and *z*
  and store the result in *beta*.

  Parameters:
    sopalin_data - PaStiX data structure.
    me           - Thread ID.
    r            - first vector of size *colnbr* times *smxnbr*.
    z            - second vector of size *colnbr* times *smxnbr*.a
    colnbr       - Number of unknowns.
    smxnbr       - Number of right-hand-side members.
    beta         - Float which will store the solution.
    comm         - MPI communicator.
*/
void CscGmresBeta(Sopalin_Data_t       *sopalin_data,
                  int                   me,
                  const volatile pastix_float_t *r,
                  const volatile pastix_float_t *z,
                  pastix_int_t                   colnbr,
                  pastix_int_t                   smxnbr,
                  double               *beta,
                  MPI_Comm              comm)
{
  SolverMatrix *  datacode;
  pastix_int_t   first,  last;
  pastix_int_t   first2, last2;
  pastix_int_t   step;
  pastix_int_t   itersmx;
  pastix_int_t   iter = 0;
  double up   = 0.0;
#ifndef INOUT_ALLREDUCE
  double up2  = 0.0;
#endif
  (void)comm;

#ifdef CSC_LOG
  fprintf(stdout, "-> CscGmresBeta \n");
#endif

  datacode  = sopalin_data->datacode;
#ifdef DEBUG_RAFF
  {
    MONOTHREAD_BEGIN;
    FILE *rafffile;
    char  rafffilename[30];
    static pastix_int_t toto = 0;
    sprintf(rafffilename, "CscGmresBetar%ld.%ld",
            (long) toto,(long) SOLV_PROCNUM);
    rafffile = fopen(rafffilename, "w");
    dump7((pastix_float_t *)r,rafffile,colnbr);
    fclose(rafffile);

    sprintf(rafffilename, "CscGmresBetaz%ld.%ld",
            (long) toto,(long) SOLV_PROCNUM);
    rafffile = fopen(rafffilename, "w");
    dump7((pastix_float_t *)z,rafffile,colnbr);
    fclose(rafffile);
    toto++;
    MONOTHREAD_END;
  }
#endif

#ifdef SMP_RAFF
  step  = (pastix_int_t)ceil((double)(colnbr)/(double)(SOLV_THRDNBR));
#else
  step  = colnbr;
#endif
  first = me * step;
  last  = MIN(colnbr, (me+1) * step);

  for (itersmx=0; itersmx<smxnbr; itersmx++)
    {
      /* Produit scalaire sur les donnÃ©es locales au thread */
      first2 = first + itersmx*colnbr;
      last2  = last  + itersmx*colnbr;
      for (iter=first2; iter<last2; iter++)
        {
          up += (double)(r[iter]*CONJ_FLOAT(z[iter]));
        }

      /* En SMP reduction du resultat sur les threads */
#ifdef SMP_RAFF
      sopalin_data->common_flt[me] = up;
      /*      fprintf(stdout,"%d: line %d up %.20g\n",me,__LINE__,up);  */

      SYNCHRO_THREAD;
      MONOTHREAD_BEGIN;
      for (iter = 1; iter < SOLV_THRDNBR; iter++)
        up += (double)(sopalin_data->common_flt[iter]);


#endif /* SMP_RAFF */

      /* Reduction en MPI */
#ifdef INOUT_ALLREDUCE
      MyMPI_Allreduce((void*)&up, (void*)&up,  1, MPI_DOUBLE, MPI_SUM, comm);
#else
      MyMPI_Allreduce((void*)&up, (void*)&up2, 1, MPI_DOUBLE, MPI_SUM, comm);
      up = up2;
#endif
      /* on s'assure que tous les threads aient la bonne valeur de up */
#ifdef SMP_RAFF
      sopalin_data->common_flt[0]  = up;
      MONOTHREAD_END;
      SYNCHRO_THREAD;
      up  = sopalin_data->common_flt[0];
#endif /* SMP_RAFF */

      beta[itersmx] = up;

#ifdef SMP_RAFF
      SYNCHRO_THREAD;
#endif
      up = 0;
    }

#ifdef CSC_LOG
  fprintf(stdout, "<- CscGmresBeta \n");
#endif
}
