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
 * File: z_raff_functions.c
 *
 * Functions computing operations for reffinement methods
 *
 */

#include "common.h"
#include "z_csc.h"
#include "z_bcsc.h"

//#include "z_tools.h"
//#ifdef PASTIX_EZTRACE
//#  include "pastix_eztrace.h"
//#else
//#  include "trace.h"
//#endif
//#include "sopalin_define.h"
//#include "symbol.h"
//#include "z_csc.h"
//#include "z_updown.h"
//#include "queue.h"
//#include "bulles.h"
//#include "z_ftgt.h"
//#include "z_solver.h"
//#include "sopalin_thread.h"
//#include "stack.h"
//#include "z_sopalin3d.h"
//#include "z_sopalin_init.h"
//#include "perf.h"
//#include "out.h"
//#include "z_coefinit.h"
//#include "z_ooc.h"
//#include "order.h"
//#include "z_debug_dump.h"
//#include "sopalin_acces.h"
//#include "z_csc_intern_compute.h"
//#ifdef PASTIX_WITH_STARPU
//#  include "starpu_zsubmit_tasks.h"
//#endif

static pastix_complex64_t fun   = 1.0;
// #define z_up_down_smp API_CALL(z_up_down_smp)
void* z_up_down_smp ( void *arg );
// #define z_sopalin_updo_comm API_CALL(z_sopalin_updo_comm)
void *z_sopalin_updo_comm ( void *arg );


#include "z_raff_functions.h"


/*** ALLOCATIONS ET SYNCHRONISATIONS ***/

/* Synchronise le vecteur x dans la nb-ieme variable de la structure */
pastix_complex64_t *z_Pastix_Synchro_Vect(void *arg, void *x, int nb)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  pastix_int_t        me           = argument->me;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
//   MONOTHREAD_BEGIN;
  sopalin_data->ptr_raff[nb] = x;
//   MONOTHREAD_END;
//   SYNCHRO_THREAD;
  return (pastix_complex64_t*) sopalin_data->ptr_raff[nb];
}

/* Alloue un vecteur de taille size octets */
void *z_Pastix_Malloc(void *arg, size_t size)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  pastix_int_t        me           = argument->me;
  void *x = NULL;
//   MONOTHREAD_BEGIN;
  MALLOC_INTERN(x, size, char);
  memset(x, 0, size);
//   MONOTHREAD_END;
  return x;
}

/* Libere un vecteur */
void z_Pastix_Free(void *arg, void *x)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  pastix_int_t        me           = argument->me;
//   MONOTHREAD_BEGIN;
  memFree_null(x);
//   MONOTHREAD_END;
}


/*** GESTION DE L'INTERFACE ***/

/* Affichage à chaque itération et communication de certaines informations à la structure */
void z_Pastix_Verbose(void *arg, double t0, double t3, double tmp, pastix_int_t nb_iter)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  z_SopalinParam     *sopar        = sopalin_data->sopar;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;
  sopalin_data->count_iter = nb_iter;
  sopalin_data->stop = tmp;
//   MONOTHREAD_BEGIN;
  if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    {
      double rst = 0.0;
      double stt, rtt;
      double err, stop = tmp;

      stt = t3 - t0;
      MyMPI_Reduce(&stop, &err, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
      MyMPI_Reduce(&stt,  &rtt, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);

      if (SOLV_PROCNUM == 0)
        {
          fprintf(stdout, OUT_ITERRAFF_ITER, (int)sopalin_data->count_iter);
          if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
            fprintf(stdout, OUT_ITERRAFF_TTS, rst);
          fprintf(stdout, OUT_ITERRAFF_TTT, stt);
          fprintf(stdout, OUT_ITERRAFF_ERR, err);
        }
    }
//   MONOTHREAD_END;
}

/* Affichage final */
void z_Pastix_End(void* arg, pastix_complex64_t tmp, pastix_int_t nb_iter, double t, pastix_complex64_t *x)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SopalinParam     *sopar        = sopalin_data->sopar;
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;

  sopalin_data->stop = tmp;
//   MULTITHREAD_BEGIN;
  z_CscCopy(sopalin_data, me, x, UPDOWN_SM2XTAB,
          UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   MULTITHREAD_END(0);
//   SYNCHRO_THREAD;

  sopar->rberror = tmp;
  sopar->itermax = nb_iter;

  if (sopar->iparm[IPARM_PRODUCE_STATS] == API_YES) {
    pastix_complex64_t *r, *s;

//     MONOTHREAD_BEGIN;
    MALLOC_INTERN(r, UPDOWN_SM2XSZE, pastix_complex64_t);
    MALLOC_INTERN(s, UPDOWN_SM2XSZE, pastix_complex64_t);
    sopalin_data->ptr_raff[0] = (void *)r;
    sopalin_data->ptr_raff[1] = (void *)s;
//     MONOTHREAD_END;
//     SYNCHRO_THREAD;

    r = (pastix_complex64_t *)sopalin_data->ptr_raff[0];
    s = (pastix_complex64_t *)sopalin_data->ptr_raff[1];
//     MULTITHREAD_BEGIN;
    /* compute r = b - Ax */
//     z_CscbMAx(sopalin_data, me, r, sopar->b, sopar->cscmtx,
//             &(datacode->updovct), datacode, PASTIX_COMM,
//             sopar->iparm[IPARM_TRANSPOSE_SOLVE]);
    z_bcscGemv( -1.0, sopar->cscmtx, &(datacode->updovct), 1.0, sopar->b)
    /* |A||x| + |b| */
//     z_CscAxPb( sopalin_data, me, s, sopar->b, sopar->cscmtx,
//              &(datacode->updovct), datacode, PASTIX_COMM,
//              sopar->iparm[IPARM_TRANSPOSE_SOLVE]);
    s = z_bcscAxpb();
    z_CscBerr(sopalin_data, me, r, s, UPDOWN_SM2XSZE,
            1, &(sopalin_data->sopar->dparm[DPARM_SCALED_RESIDUAL]),
            PASTIX_COMM);
//     MULTITHREAD_END(1);
  }

//   MONOTHREAD_BEGIN;
// 
//   if (THREAD_COMM_ON)
//     {
//       if (sopar->iparm[IPARM_END_TASK] >= API_TASK_REFINE)
//         {
//           MUTEX_LOCK(&(sopalin_data->mutex_comm));
//           sopalin_data->step_comm = COMMSTEP_END;
//           print_debug(DBG_THCOMM, "%s:%d END\n", __FILE__, __LINE__);
//           MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
//           pthread_cond_broadcast(&(sopalin_data->cond_comm));
//         }
//     }
// #ifdef OOC
//   z_ooc_stop_thread(sopalin_data);
// #endif
  sopar->dparm[DPARM_RAFF_TIME] = t;
//   MONOTHREAD_END;
//   SYNCHRO_THREAD;
}

/* Vecteur solution X */
void Pastix_X(void *arg, pastix_complex64_t *x)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  z_SopalinParam     *sopar        = sopalin_data->sopar;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;
  pastix_int_t        i;

  if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
    for (i=0;i<UPDOWN_SM2XSZE*UPDOWN_SM2XNBR;i++)
      UPDOWN_SM2XTAB[i]=0.0;
//   MULTITHREAD_BEGIN;
  z_CscCopy(sopalin_data, me, UPDOWN_SM2XTAB, x,
          UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   MULTITHREAD_END(1);
//   SYNCHRO_THREAD;
}

/* Taille d'un vecteur */
pastix_int_t z_Pastix_n(void *arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  return UPDOWN_SM2XSZE;
}

/* Nombre de second membres */
pastix_int_t z_Pastix_m(void *arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  return UPDOWN_SM2XNBR;
}

/* Second membre */
void z_Pastix_B(void *arg, pastix_complex64_t *b)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SopalinParam     *sopar        = sopalin_data->sopar;
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;
//   MULTITHREAD_BEGIN;
  z_CscCopy(sopalin_data, me, sopar->b, b,
          UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   MULTITHREAD_END(0);
//   SYNCHRO_THREAD;
}

/* Epsilon */
pastix_complex64_t z_Pastix_Eps(void *arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SopalinParam     *sopar        = sopalin_data->sopar;
  return sopar->epsilonraff;
}

/* Itermax */
pastix_int_t z_Pastix_Itermax(void *arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SopalinParam     *sopar        = sopalin_data->sopar;
  return sopar->itermax;
}


/* Itermax */
pastix_int_t z_Pastix_Krylov_Space(void *arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SopalinParam     *sopar        = sopalin_data->sopar;
  return sopar->gmresim;
}

/*** OPERATIONS DE BASE ***/
/* Multiplication pour plusieurs second membres */
void z_Pastix_Mult(void *arg, pastix_complex64_t *alpha, pastix_complex64_t *beta, pastix_complex64_t *zeta, int flag)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  pastix_int_t        me           = argument->me;
//   MONOTHREAD_BEGIN;
#ifdef MULT_SMX_RAFF
  {
    pastix_int_t itersmx;
    for(itersmx=0; itersmx<UPDOWN_SM2XNBR;itersmx++)
      {
        zeta[itersmx]=alpha[itersmx]*beta[itersmx];
      }
  }
#else
  zeta[0]=alpha[0]*beta[0];
#endif
//   MONOTHREAD_END;
//   if (flag)
//     SYNCHRO_THREAD;
}

/* Division pour plusieurs second membres */
void z_Pastix_Div(void *arg, pastix_complex64_t *alpha, pastix_complex64_t *beta, pastix_complex64_t *zeta, int flag)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  pastix_int_t        me           = argument->me;
//   MONOTHREAD_BEGIN;
#ifdef MULT_SMX_RAFF
  {
    pastix_int_t itersmx;
    for(itersmx=0; itersmx<UPDOWN_SM2XNBR;itersmx++)
      {
        zeta[itersmx]=alpha[itersmx]/beta[itersmx];
      }
  }
#else
  zeta[0]=alpha[0]/beta[0];
#endif
//   MONOTHREAD_END;
//   if (flag)
//     SYNCHRO_THREAD;
}

/* Calcul de la norme de frobenius */
pastix_complex64_t z_Pastix_Norm2(void* arg, pastix_complex64_t *x)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;
  double            normx;
//   MULTITHREAD_BEGIN;
  normx = z_CscNormFro(sopalin_data, me, x,
                     UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   MULTITHREAD_END(1);
//   NOSMP_SYNC_COEF(normx);
  return normx;
}

/* Copie d'un vecteur */
void z_Pastix_Copy(void *arg, pastix_complex64_t *s, pastix_complex64_t *d, int flag)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;
//   MULTITHREAD_BEGIN;
  z_CscCopy(sopalin_data, me, s, d,
          UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   MULTITHREAD_END(0);

//   if (flag)
//     SYNCHRO_THREAD;
}

/* Application du préconditionneur */
void z_Pastix_Precond(void *arg, pastix_complex64_t *s, pastix_complex64_t *d, int flag)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  z_SopalinParam     *sopar        = sopalin_data->sopar;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;

//   MULTITHREAD_BEGIN;
  z_CscCopy(sopalin_data, me, s, UPDOWN_SM2XTAB,
          UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   MULTITHREAD_END(1);
  /* M-1 updo -> updo */
#ifdef PRECOND
  if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
    {
//       SYNCHRO_THREAD;
      API_CALL(z_up_down_smp)(arg);
//       SYNCHRO_THREAD;
    }
#endif
  MULTITHREAD_BEGIN;
  z_CscCopy(sopalin_data, me, UPDOWN_SM2XTAB, d,
          UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   MULTITHREAD_END(0);
//   if (flag)
//     SYNCHRO_THREAD;
}

/* Calcul de alpha * x */
void z_Pastix_Scal(void *arg, pastix_complex64_t alpha, pastix_complex64_t *x, int flag)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;
//   MULTITHREAD_BEGIN;
  z_CscScal(sopalin_data, me, alpha, x,
          UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   MULTITHREAD_END(0);
//   if (flag)
//     SYNCHRO_THREAD;
}

/* Calcul du produit scalaire */
void z_Pastix_Dotc(void *arg, pastix_complex64_t *x, pastix_complex64_t *y, pastix_complex64_t *r, int flag)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;
//   MULTITHREAD_BEGIN;
  z_CscGradBeta(sopalin_data, me, x, y,
              UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, r, pastix_comm);
//   MULTITHREAD_END(0);
//   if (flag)
//     SYNCHRO_THREAD;
}

void z_Pastix_Dotc_Gmres(void *arg, pastix_complex64_t *x, pastix_complex64_t *y, pastix_complex64_t *r, int flag)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;
//   MULTITHREAD_BEGIN;
  z_CscGmresBeta(sopalin_data, me, x, y,
               UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, r, pastix_comm);
//   MULTITHREAD_END(0);
//   SYNC_COEF(*r);
}

/* Produit matrice vecteur */
void z_Pastix_Ax(void *arg, pastix_complex64_t *x, pastix_complex64_t *r)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  z_SopalinParam     *sopar        = sopalin_data->sopar;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;
//   MULTITHREAD_BEGIN;
  z_CscAx(sopalin_data, me, sopalin_data->sopar->cscmtx, x, r,
        datacode, &(datacode->updovct), pastix_comm,
        sopar->iparm[IPARM_TRANSPOSE_SOLVE]);
//   MULTITHREAD_END(1);
}


/*** A MODIFIER! ***/
void z_Pastix_bMAx(void *arg, pastix_complex64_t *b, pastix_complex64_t *x, pastix_complex64_t *r)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  z_SopalinParam     *sopar        = sopalin_data->sopar;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;

//   MULTITHREAD_BEGIN;
  z_CscCopy(sopalin_data, me, x, UPDOWN_SM2XTAB,
          UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   MULTITHREAD_END(0);
//   SYNCHRO_THREAD;
//   MULTITHREAD_BEGIN;
  z_CscbMAx(sopalin_data, me, r, b, sopalin_data->sopar->cscmtx,
          &(datacode->updovct), datacode, pastix_comm,
          sopar->iparm[IPARM_TRANSPOSE_SOLVE]);
//   MULTITHREAD_END(1);
}

void z_Pastix_BYPX(void *arg, pastix_complex64_t *beta, pastix_complex64_t *y, pastix_complex64_t *x, int flag)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;

#ifdef MULT_SMX_RAFF
  {
    pastix_int_t itersmx;
    for (itersmx=0; itersmx<UPDOWN_SM2XNBR; itersmx++)
      {
//         MULTITHREAD_BEGIN;
        z_CscScal(sopalin_data, me, beta[itersmx], x+(itersmx*UPDOWN_SM2XSZE),
                UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//         MULTITHREAD_END(0);
//         SYNCHRO_THREAD;
      }
  }
//   MONOTHREAD_BEGIN;
  SOPALIN_GEAM("N","N",UPDOWN_SM2XSZE,UPDOWN_SM2XNBR, fun,
               y, UPDOWN_SM2XSZE, x, UPDOWN_SM2XSZE);
//   MONOTHREAD_END;
#else
//   MULTITHREAD_BEGIN;
  z_CscScal(sopalin_data, me, beta[0], x,
          UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   MULTITHREAD_END(0);
//   SYNCHRO_THREAD;
//   MULTITHREAD_BEGIN;
  z_CscAXPY(sopalin_data, me, fun, y, x,
          UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   MULTITHREAD_END(0);
#endif
//   if (flag)
//     SYNCHRO_THREAD;
}


void z_Pastix_AXPY(void *arg, double coeff, pastix_complex64_t *alpha, pastix_complex64_t *x, pastix_complex64_t *y, int flag)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t        me           = argument->me;
  pastix_complex64_t      tmp_flt;
#ifdef MULT_SMX_RAFF
  {
    pastix_int_t itersmx;
    for(itersmx=0; itersmx<UPDOWN_SM2XNBR; itersmx++)
      {
        tmp_flt = (pastix_complex64_t) alpha[itersmx] * coeff;
//         MULTITHREAD_BEGIN;
        z_CscAXPY(sopalin_data, me, tmp_flt, y+(itersmx*UPDOWN_SM2XSZE), x+(itersmx*UPDOWN_SM2XSZE),
                UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//         MULTITHREAD_END(1);
      }
  }
#else
  tmp_flt = (pastix_complex64_t) alpha[0] * coeff;
//   MULTITHREAD_BEGIN;
  z_CscAXPY(sopalin_data, me, tmp_flt, y, x,
          UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   MULTITHREAD_END(0);
#endif
//   if (flag)
//     SYNCHRO_THREAD;
}


pastix_int_t z_Pastix_me(void *arg)
{
  sopthread_data_t *argument = (sopthread_data_t *)arg;
  pastix_int_t        me       = argument->me;
  return me;
}

void z_Pastix_Solveur(struct z_solver *solveur)
{
  /*** ALLOCATIONS ET SYNCHRONISATIONS ***/
  solveur->Synchro     = &z_Pastix_Synchro_Vect;
  solveur->Malloc      = &z_Pastix_Malloc;
  solveur->Free        = &z_Pastix_Free;

  /*** GESTION DE L'INTERFACE ***/
  solveur->Verbose = &z_Pastix_Verbose;
  solveur->End     = &z_Pastix_End;
  solveur->X       = &Pastix_X;
  solveur->N       = &z_Pastix_n;
  solveur->B       = &z_Pastix_B;
  solveur->Eps     = &z_Pastix_Eps;
  solveur->Itermax = &z_Pastix_Itermax;
  solveur->me      = &z_Pastix_me;
  solveur->Krylov_Space = &z_Pastix_Krylov_Space;

  /*** OPERATIONS DE BASE ***/
  solveur->Mult     = &z_Pastix_Mult;
  solveur->Div      = &z_Pastix_Div;
  solveur->Dotc_Gmres = &z_Pastix_Dotc_Gmres;

  solveur->Norm    = &z_Pastix_Norm2;
  solveur->Copy    = &z_Pastix_Copy;
  solveur->Precond = &z_Pastix_Precond;

  solveur->Scal    = &z_Pastix_Scal;
  solveur->Dotc    = &z_Pastix_Dotc;
  solveur->Ax      = &z_Pastix_Ax;

  solveur->AXPY    = &z_Pastix_AXPY;
  solveur->bMAx    = &z_Pastix_bMAx;
  solveur->BYPX    = &z_Pastix_BYPX;
}

/*
 ** Section: Function creating threads
 */
/*
 Function: method)

 Launch sopaparam->nbthrdcomm threads which will compute
 <method_smp)>.

 Parameters:
 datacode  - PaStiX <z_SolverMatrix> structure.
 sopaparam - <z_SopalinParam> parameters structure.
 */
void z_raff_thread(z_SolverMatrix *datacode, z_SopalinParam *sopaparam, void*(*method)(void *))
{
  z_Sopalin_Data_t *sopalin_data = NULL;
  BackupSolve_t b;

  MALLOC_INTERN(sopalin_data, 1, z_Sopalin_Data_t);

  z_solve_backup(datacode,&b);
  z_sopalin_init(sopalin_data, datacode, sopaparam, 0);

  sopalin_launch_thread(sopalin_data,
                        SOLV_PROCNUM,          SOLV_PROCNBR,                datacode->btree,
                        sopalin_data->sopar->iparm[IPARM_VERBOSE],
                        SOLV_THRDNBR,          method,                      sopalin_data,
                        sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm), sopalin_data,
                        OOC_THREAD_NBR,        z_ooc_thread,                  sopalin_data);

  z_sopalin_clean(sopalin_data, 2);
  z_solve_restore(datacode,&b);
  memFree_null(sopalin_data);
}
