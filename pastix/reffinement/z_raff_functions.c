/**
 *
 * @file: z_raff_functions.c
 *
 *  Functions computing operations for reffinement methods
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
#include "common.h"
#include "z_spm.h"
#include "bcsc.h"
#include "z_bcsc.h"
#include "sopalin_thread.h"
#include "sopalin_data.h"
#include "solver.h"
#include "z_raff_functions.h"

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

// static pastix_complex64_t fun   = 1.0;
// #define z_up_down_smp API_CALL(z_up_down_smp)
// void* z_up_down_smp ( void *arg );
// #define z_sopalin_updo_comm API_CALL(z_sopalin_updo_comm)
// void *z_sopalin_updo_comm ( void *arg );


/*** ALLOCATIONS ET SYNCHRONISATIONS ***/

/* Synchronise le vecteur x dans la nb-ieme variable de la structure */
// pastix_complex64_t *z_Pastix_Synchro_Vect(void *arg, void *x, int nb)
// {
//   sopthread_data_t *argument     = (sopthread_data_t *)arg;
//   pastix_int_t        me           = argument->me;
//   sopalin_data_t   *sopalin_data = (sopalin_data_t *)(argument->data);
//   SolverMatrix     *datacode     = sopalin_data->datacode;
// //   MONOTHREAD_BEGIN;
//   sopalin_data->ptr_raff[nb] = x;
// //   MONOTHREAD_END;
// //   SYNCHRO_THREAD;
//   return (pastix_complex64_t*) sopalin_data->ptr_raff[nb];
// }

/* Alloue un vecteur de taille size octets */
void *z_Pastix_Malloc(size_t size)
{
  void *x = NULL;
  MALLOC_INTERN(x, size, char);
  memset(x, 0, size);
  return x;
}

/* Libere un vecteur */
void z_Pastix_Free( void *x)
{
  memFree_null(x);
}


/*** GESTION DE L'INTERFACE ***/

/* Affichage à chaque itération et communication de certaines informations à la structure */
void z_Pastix_Verbose(double t0, double t3, double tmp, pastix_int_t nb_iter)
{
//   sopthread_data_t *argument     = (sopthread_data_t *)arg;
//   sopalin_data_t   *sopalin_data = (sopalin_data_t *)(argument->data);
//   SolverMatrix     *datacode     = sopalin_data->datacode;
//   SopalinParam     *sopar        = sopalin_data->sopar;
//   MPI_Comm          pastix_comm  = PASTIX_COMM;
//   pastix_int_t        me           = argument->me;
//   sopalin_data->count_iter = nb_iter;
//   sopalin_data->stop = tmp;
//   MONOTHREAD_BEGIN;
//   if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
//     {
      double rst = 0.0;
      double stt, rtt;
      double err, stop = tmp;

      stt = t3 - t0;
      err = stop;
      rtt = stt;
//       MyMPI_Reduce(&stop, &err, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
//       MyMPI_Reduce(&stt,  &rtt, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);

//       if (SOLV_PROCNUM == 0)
//         {
          fprintf(stdout, OUT_ITERRAFF_ITER, (int)nb_iter);
//           if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
            fprintf(stdout, OUT_ITERRAFF_TTS, rst);
          fprintf(stdout, OUT_ITERRAFF_TTT, stt);
          fprintf(stdout, OUT_ITERRAFF_ERR, err);
//         }
//     }
//   MONOTHREAD_END;
}

/* Affichage final */
void z_Pastix_End(SopalinParam *sopar, pastix_complex64_t tmp, pastix_int_t nb_iter, double t, void *x, pastix_complex64_t *gmresx)
{
//   sopthread_data_t *argument     = (sopthread_data_t *)arg;
//   Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
//   SopalinParam     *sopar        = sopalin_data->sopar;
//   SolverMatrix     *datacode     = sopalin_data->datacode;
//   MPI_Comm          pastix_comm  = PASTIX_COMM;
//   PASTIX_INT        me           = argument->me;
    pastix_complex64_t *xptr = (pastix_complex64_t *)x;
    pastix_int_t        n = sopar->gN;
    pastix_int_t i;

//   sopalin_data->stop = tmp;
//   CscCopy(sopalin_data, me, x, UPDOWN_SM2XTAB,
//           UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
//   memcpy(x,gmresx,sopar->gN * sizeof(pastix_complex64_t));
    for (i=0; i<n; i++)
      xptr[i] = gmresx[i];

//   if (pastix_data->iparm[IPARM_ONLY_RAFF] == API_NO)
//   {
//   MULTITHREAD_END(0);
//   SYNCHRO_THREAD;

  sopar->rberror = tmp;
  sopar->itermax = nb_iter;

//   if (sopar->iparm[IPARM_PRODUCE_STATS] == API_YES) {
//     PASTIX_FLOAT *r, *s;
// 
//     MONOTHREAD_BEGIN;
//     MALLOC_INTERN(r, UPDOWN_SM2XSZE, PASTIX_FLOAT);
//     MALLOC_INTERN(s, UPDOWN_SM2XSZE, PASTIX_FLOAT);
//     sopalin_data->ptr_raff[0] = (void *)r;
//     sopalin_data->ptr_raff[1] = (void *)s;
//     MONOTHREAD_END;
//     SYNCHRO_THREAD;
// 
//     r = (PASTIX_FLOAT *)sopalin_data->ptr_raff[0];
//     s = (PASTIX_FLOAT *)sopalin_data->ptr_raff[1];
//     MULTITHREAD_BEGIN;
//     /* compute r = b - Ax */
//     CscbMAx(sopalin_data, me, r, sopar->b, sopar->cscmtx,
//             &(datacode->updovct), datacode, PASTIX_COMM,
//             sopar->iparm[IPARM_TRANSPOSE_SOLVE]);
//     /* |A||x| + |b| */
//     CscAxPb( sopalin_data, me, s, sopar->b, sopar->cscmtx,
//              &(datacode->updovct), datacode, PASTIX_COMM,
//              sopar->iparm[IPARM_TRANSPOSE_SOLVE]);
//     CscBerr(sopalin_data, me, r, s, UPDOWN_SM2XSZE,
//             1, &(sopalin_data->sopar->dparm[DPARM_SCALED_RESIDUAL]),
//             PASTIX_COMM);
//     MULTITHREAD_END(1);
//   }
}

/* Vecteur solution X */
void z_Pastix_X(pastix_data_t *pastix_data, void *x, pastix_complex64_t *gmresx)
{
  pastix_int_t        i;
  pastix_int_t        n = pastix_data->bcsc->gN;
  pastix_complex64_t *xptr = (pastix_complex64_t *)x;

//   if (pastix_data->iparm[IPARM_ONLY_RAFF] == API_NO)
//   {
    for (i=0; i<n; i++, xptr++)
      gmresx[i]= *xptr;
//   }
//   else
//   {
//     for (i=0; i<n; i++, xptr++)
//       gmresx[i]=0.0;
//   }
}

/* Taille d'un vecteur */
pastix_int_t z_Pastix_n(SopalinParam *sopar)
{
  return sopar->gN;
}

/* Nombre de second membres */
// pastix_int_t z_Pastix_m(void *arg)
// {
//   sopthread_data_t *argument     = (sopthread_data_t *)arg;
//   sopalin_data_t   *sopalin_data = (sopalin_data_t *)(argument->data);
//   SolverMatrix     *datacode     = sopalin_data->datacode;
//   return UPDOWN_SM2XNBR;
// }

/* Second membre */
void z_Pastix_B(void *b, pastix_complex64_t *raffb, pastix_int_t n)
{
  pastix_complex64_t *bptr = (pastix_complex64_t *)b;
  pastix_int_t i;

  for (i=0; i<n; i++, bptr++)
  {
      raffb[i]= *bptr;
  }
//   memcpy(raffb, b, n * sizeof( pastix_complex64_t ));
}

/* Epsilon */
pastix_complex64_t z_Pastix_Eps(SopalinParam *sopar)
{
  return sopar->epsilonraff;
}

/* Itermax */
pastix_int_t z_Pastix_Itermax(SopalinParam *sopar)
{
  return sopar->itermax;
}


/* Itermax */
pastix_int_t z_Pastix_Krylov_Space(SopalinParam *sopar)
{
  return sopar->gmresim;
}

/*** OPERATIONS DE BASE ***/
/* Multiplication pour plusieurs second membres */
// void z_Pastix_Mult(void *arg, pastix_complex64_t *alpha, pastix_complex64_t *beta, pastix_complex64_t *zeta, int flag)
// {
//   sopthread_data_t *argument     = (sopthread_data_t *)arg;
//   sopalin_data_t   *sopalin_data = (sopalin_data_t *)(argument->data);
//   SolverMatrix     *datacode     = sopalin_data->datacode;
//   pastix_int_t        me           = argument->me;
// //   MONOTHREAD_BEGIN;
// #ifdef MULT_SMX_RAFF
//   {
//     pastix_int_t itersmx;
//     for(itersmx=0; itersmx<UPDOWN_SM2XNBR;itersmx++)
//       {
//         zeta[itersmx]=alpha[itersmx]*beta[itersmx];
//       }
//   }
// #else
//   zeta[0]=alpha[0]*beta[0];
// #endif
// //   MONOTHREAD_END;
// //   if (flag)
// //     SYNCHRO_THREAD;
// }

/* Division pour plusieurs second membres */
// void z_Pastix_Div(void *arg, pastix_complex64_t *alpha, pastix_complex64_t *beta, pastix_complex64_t *zeta, int flag)
// {
//   sopthread_data_t *argument     = (sopthread_data_t *)arg;
//   sopalin_data_t   *sopalin_data = (sopalin_data_t *)(argument->data);
//   SolverMatrix     *datacode     = sopalin_data->datacode;
//   pastix_int_t        me           = argument->me;
// //   MONOTHREAD_BEGIN;
// #ifdef MULT_SMX_RAFF
//   {
//     pastix_int_t itersmx;
//     for(itersmx=0; itersmx<UPDOWN_SM2XNBR;itersmx++)
//       {
//         zeta[itersmx]=alpha[itersmx]/beta[itersmx];
//       }
//   }
// #else
//   zeta[0]=alpha[0]/beta[0];
// #endif
// //   MONOTHREAD_END;
// //   if (flag)
// //     SYNCHRO_THREAD;
// }

/* Calcul de la norme de frobenius */
pastix_complex64_t z_Pastix_Norm2(pastix_complex64_t *x, pastix_int_t n)
{
  double normx;
  void *xptr = (void*)x;
  normx = z_vectFrobeniusNorm(xptr, n);
  return normx;
}

/* Copie d'un vecteur */
// void z_Pastix_Copy(void *arg, pastix_complex64_t *s, pastix_complex64_t *d, int flag)
// {
//   sopthread_data_t *argument     = (sopthread_data_t *)arg;
//   sopalin_data_t   *sopalin_data = (sopalin_data_t *)(argument->data);
//   SolverMatrix     *datacode     = sopalin_data->datacode;
//   MPI_Comm          pastix_comm  = PASTIX_COMM;
//   pastix_int_t        me           = argument->me;
// //   MULTITHREAD_BEGIN;
//   z_CscCopy(sopalin_data, me, s, d,
//           UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
// //   MULTITHREAD_END(0);
// 
// //   if (flag)
// //     SYNCHRO_THREAD;
// }

/* Application du préconditionneur */
void z_Pastix_Precond(pastix_data_t *pastix_data, pastix_complex64_t *s, pastix_complex64_t *d, int flag)
{
  pastix_int_t n = pastix_data->bcsc->gN;
  pastix_int_t nrhs = 1;
  void* bptr = (void*)d;

  memcpy(d, s, n * sizeof( pastix_complex64_t ));
  if (pastix_data->iparm[IPARM_ONLY_RAFF] == API_NO)
    {
        sopalin_data_t sopalin_data;
        sopalin_data.solvmtx = pastix_data->solvmatr;

        switch ( pastix_data->iparm[IPARM_FACTORIZATION] ){
        case PastixFactLLT:
            sequential_ztrsm( PastixLeft, PastixLower, PastixNoTrans,   PastixNonUnit, &sopalin_data, nrhs, bptr, n );
            sequential_ztrsm( PastixLeft, PastixLower, PastixConjTrans, PastixNonUnit, &sopalin_data, nrhs, bptr, n );
            break;

        case PastixFactLDLT:
            sequential_ztrsm( PastixLeft, PastixLower, PastixNoTrans, PastixUnit, &sopalin_data, nrhs, bptr, n );
            sequential_zdiag( &sopalin_data, nrhs, bptr, n );
            sequential_ztrsm( PastixLeft, PastixLower, PastixTrans,   PastixUnit, &sopalin_data, nrhs, bptr, n );
            break;

        case PastixFactLDLH:
            sequential_ztrsm( PastixLeft, PastixLower, PastixNoTrans,   PastixUnit, &sopalin_data, nrhs, bptr, n );
            sequential_zdiag( &sopalin_data, nrhs, bptr, n );
            sequential_ztrsm( PastixLeft, PastixLower, PastixConjTrans, PastixUnit, &sopalin_data, nrhs, bptr, n );
            break;

        case PastixFactLU:
        default:
            sequential_ztrsm( PastixLeft, PastixLower, PastixNoTrans, PastixUnit,    &sopalin_data, nrhs, bptr, n );
            sequential_ztrsm( PastixLeft, PastixUpper, PastixNoTrans, PastixNonUnit, &sopalin_data, nrhs, bptr, n );
            break;
        }
    }
}

/* Calcul de alpha * x */
void z_Pastix_Scal(pastix_int_t n, pastix_complex64_t alpha, pastix_complex64_t *x, int flag)
{
    z_bcscScal( x, alpha, n, 1);
}

/* Calcul du produit scalaire */
void z_Pastix_Dotc(pastix_int_t n, pastix_complex64_t *x, pastix_complex64_t *y, pastix_complex64_t *r, int flag)
{
  *r = z_bcscDotc(x, y, n);
}

void z_Pastix_Dotc_Gmres(pastix_int_t n, pastix_complex64_t *x, pastix_complex64_t *y, pastix_complex64_t *r, int flag)
{
  (*r) = z_bcscDotc(x, y, n);
}

/* Produit matrice vecteur */
void z_Pastix_Ax(pastix_bcsc_t *bcsc, pastix_complex64_t *x, pastix_complex64_t *r)
{
    pastix_int_t alpha = 1.0;
    pastix_int_t beta = 0.0;
    void* xptr = (void*)x;
    void* yptr = (void*)r;

    switch (bcsc->mtxtype) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PastixHermitian:
        z_bcscHemv(alpha, bcsc, xptr, beta, yptr );
        break;
#endif
    case PastixSymmetric:
        z_bcscSymv(alpha, bcsc, xptr, beta, yptr );
        break;
    case PastixGeneral:
    default:
        z_bcscGemv(PastixNoTrans, alpha, bcsc, xptr, beta, yptr );
    }
}


/*** A MODIFIER! ***/
void z_Pastix_bMAx(pastix_bcsc_t *bcsc, pastix_complex64_t *b, pastix_complex64_t *x, pastix_complex64_t *r)
{
    pastix_int_t alpha = -1.0;
    pastix_int_t beta = 1.0;
    void* xptr = (void*)x;
    void* yptr = (void*)r;

    memcpy(r, b, bcsc->gN * sizeof( pastix_complex64_t ));
    switch (bcsc->mtxtype) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PastixHermitian:
        z_bcscHemv(alpha, bcsc, xptr, beta, yptr );
        break;
#endif
    case PastixSymmetric:
        z_bcscSymv(alpha, bcsc, xptr, beta, yptr );
        break;
    case PastixGeneral:
    default:
        z_bcscGemv(PastixNoTrans, alpha, bcsc, xptr, beta, yptr );
    }
}

// void z_Pastix_BYPX(void *arg, pastix_complex64_t *beta, pastix_complex64_t *y, pastix_complex64_t *x, int flag)
// {
//   sopthread_data_t *argument     = (sopthread_data_t *)arg;
//   sopalin_data_t   *sopalin_data = (sopalin_data_t *)(argument->data);
//   SolverMatrix     *datacode     = sopalin_data->datacode;
//   MPI_Comm          pastix_comm  = PASTIX_COMM;
//   pastix_int_t        me           = argument->me;
// 
// #ifdef MULT_SMX_RAFF
//   {
//     pastix_int_t itersmx;
//     for (itersmx=0; itersmx<UPDOWN_SM2XNBR; itersmx++)
//       {
// //         MULTITHREAD_BEGIN;
//         z_CscScal(sopalin_data, me, beta[itersmx], x+(itersmx*UPDOWN_SM2XSZE),
//                 UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
// //         MULTITHREAD_END(0);
// //         SYNCHRO_THREAD;
//       }
//   }
// //   MONOTHREAD_BEGIN;
//   SOPALIN_GEAM("N","N",UPDOWN_SM2XSZE,UPDOWN_SM2XNBR, fun,
//                y, UPDOWN_SM2XSZE, x, UPDOWN_SM2XSZE);
// //   MONOTHREAD_END;
// #else
// //   MULTITHREAD_BEGIN;
//   z_CscScal(sopalin_data, me, beta[0], x,
//           UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
// //   MULTITHREAD_END(0);
// //   SYNCHRO_THREAD;
// //   MULTITHREAD_BEGIN;
//   z_CscAXPY(sopalin_data, me, fun, y, x,
//           UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);
// //   MULTITHREAD_END(0);
// #endif
// //   if (flag)
// //     SYNCHRO_THREAD;
// }


void z_Pastix_AXPY(pastix_int_t n, double coeff, pastix_complex64_t *alpha, pastix_complex64_t *x, pastix_complex64_t *y, int flag)
{
    void *yptr = (void*)y;
    void *xptr = (void*)x;
    z_bcscAxpy( coeff*(*alpha), yptr, n, xptr, 1 );
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
//   solveur->Synchro     = &z_Pastix_Synchro_Vect;
  solveur->Malloc      = &z_Pastix_Malloc;
  solveur->Free        = &z_Pastix_Free;

  /*** GESTION DE L'INTERFACE ***/
  solveur->Verbose = &z_Pastix_Verbose;
  solveur->End     = &z_Pastix_End;
  solveur->X       = &z_Pastix_X;
  solveur->N       = &z_Pastix_n;
  solveur->B       = &z_Pastix_B;
  solveur->Eps     = &z_Pastix_Eps;
  solveur->Itermax = &z_Pastix_Itermax;
  solveur->me      = &z_Pastix_me;
  solveur->Krylov_Space = &z_Pastix_Krylov_Space;

  /*** OPERATIONS DE BASE ***/
//   solveur->Mult     = &z_Pastix_Mult;
//   solveur->Div      = &z_Pastix_Div;
  solveur->Dotc_Gmres = &z_Pastix_Dotc_Gmres;

  solveur->Norm    = &z_Pastix_Norm2;
//   solveur->Copy    = &z_Pastix_Copy;
  solveur->Precond = &z_Pastix_Precond;

  solveur->Scal    = &z_Pastix_Scal;
  solveur->Dotc    = &z_Pastix_Dotc;
  solveur->Ax      = &z_Pastix_Ax;

  solveur->AXPY    = &z_Pastix_AXPY;
  solveur->bMAx    = &z_Pastix_bMAx;
//   solveur->BYPX    = &z_Pastix_BYPX;
}

/*
 ** Section: Function creating threads
 */
/*
 Function: method)

 Launch sopaparam->nbthrdcomm threads which will compute
 <method_smp)>.

 Parameters:
 datacode  - PaStiX <SolverMatrix> structure.
 sopaparam - <SopalinParam> parameters structure.
 */
// void z_raff_thread(SolverMatrix *datacode, SopalinParam *sopaparam, void*(*method)(void *))
// {
//   sopalin_data_t *sopalin_data = NULL;
//   BackupSolve_t b;
// 
//   MALLOC_INTERN(sopalin_data, 1, sopalin_data_t);
// 
//   z_solve_backup(datacode,&b);
//   z_sopalin_init(sopalin_data, datacode, sopaparam, 0);
// 
//   sopalin_launch_thread(sopalin_data,
//                         SOLV_PROCNUM,          SOLV_PROCNBR,                datacode->btree,
//                         sopalin_data->sopar->iparm[IPARM_VERBOSE],
//                         SOLV_THRDNBR,          method,                      sopalin_data,
//                         sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm), sopalin_data,
//                         OOC_THREAD_NBR,        z_ooc_thread,                  sopalin_data);
// 
//   z_sopalin_clean(sopalin_data, 2);
//   z_solve_restore(datacode,&b);
//   memFree_null(sopalin_data);
// }
