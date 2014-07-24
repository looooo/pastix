#include "csc_intern_compute.h"

#define RAFF_CLOCK_INIT {clockInit(raff_clk);clockStart(raff_clk);}
#define RAFF_CLOCK_STOP {clockStop((raff_clk));}
#define RAFF_CLOCK_GET  clockGet()
/*clockVal((raff_clk))*/

/* #define DEBUG_RAFF */

/*
** Section: Functions declarations
*/

/* Raffinement du second membre */
#define pivotstatique_smp API_CALL(pivotstatique_smp)
#define gmres_smp         API_CALL(gmres_smp)
#define grad_smp          API_CALL(grad_smp)
#define pivot_thread      API_CALL(pivot_thread)
#define gmres_thread      API_CALL(gmres_thread)
#define grad_thread       API_CALL(grad_thread)

void* pivotstatique_smp(void *arg);
void* gmres_smp        (void *arg);
void* grad_smp         (void *arg);

/* Lancement d'une des fonctions seules */
void pivot_thread(d_SolverMatrix *datacode, SopalinParam *sopaparam);
void gmres_thread(d_SolverMatrix *datacode, SopalinParam *sopaparam);
void grad_thread (d_SolverMatrix *datacode, SopalinParam *sopaparam);
/*
** Section: Threads routines
*/

/*
 * Function: API_CALL(pivotstatique_smp)
 *
 * Refine the solution.
 *
 * Computes :
 *
 *   $r   = b-Ax$
 *
 *   $r^{\\prime} = |A||x| + |b|$
 *
 *   $err = max_{i = 0..n} (r_i/r^{\\prime}_i)$
 *
 *   $rberror = ||r|| / ||b||$
 *
 * While the maximum number of iterations is not reached and the solution is
 * not precise enough, iterates :
 *
 *   Copy Up-down vector into $r^{\\prime}$.
 *
 *   Copy $r$ into Up-down vector.
 *
 *   Solves $Ax_1 = r$
 *
 *   Adds $x_1$ to previous $x$ (stored in $r^{\\prime}$)
 *
 * Parameters:
 *   arg - Pointer to a <sopthread_data_t> structure containing
 *         the <Sopalin_Data_t> structure and the thread number ID.
 */
void* pivotstatique_smp ( void *arg )
{
  double            raff_clk;
  double            t0           = 0;
  double            t1           = 0;
  double            t2           = 0;
  double            t3           = 0;
  pastix_float_t * volatile  lub          = NULL;
  pastix_float_t * volatile  lur          = NULL;
  pastix_float_t * volatile  lur2         = NULL;
  double            tmp_berr     = 0.0;
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  d_SolverMatrix     *datacode     = sopalin_data->datacode;
  SopalinParam     *sopar        = sopalin_data->sopar;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t               me           = argument->me;
  int               iter         = 0;

  MONOTHREAD_BEGIN;
  if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    {
      if (SOLV_PROCNUM == 0)
        {
          fprintf(stdout, OUT_ITERRAFF_PIVOT);
        }
    }
  MALLOC_INTERN(lub,  UPDOWN_SM2XNBR*UPDOWN_SM2XSZE, pastix_float_t);
  MALLOC_INTERN(lur,  UPDOWN_SM2XNBR*UPDOWN_SM2XSZE, pastix_float_t);
  MALLOC_INTERN(lur2, UPDOWN_SM2XNBR*UPDOWN_SM2XSZE, pastix_float_t);

  SOPALIN_COPY(UPDOWN_SM2XNBR*UPDOWN_SM2XSZE,sopar->b,iun,lub,iun);

  sopalin_data->ptr_raff[0] = (void *)lur;
  sopalin_data->ptr_raff[1] = (void *)lub;
  sopalin_data->ptr_raff[2] = (void *)lur2;

  MONOTHREAD_END;
  SYNCHRO_THREAD;

  lur  = (pastix_float_t *)sopalin_data->ptr_raff[0];
  lub  = (pastix_float_t *)sopalin_data->ptr_raff[1];
  lur2 = (pastix_float_t *)sopalin_data->ptr_raff[2];

  RAFF_CLOCK_INIT;

  while(sopalin_data->flag_gmres)
    {
      iter++;
      RAFF_CLOCK_STOP;
      t0 = RAFF_CLOCK_GET;
#ifndef SMP_RAFF
      MONOTHREAD_BEGIN;
#endif /* SMP_RAFF */
      /* r=b-ax */
      CscbMAx(sopalin_data, me, lur, lub, sopalin_data->sopar->cscmtx,
              &(datacode->updovct), datacode, pastix_comm,
              sopar->iparm[IPARM_TRANSPOSE_SOLVE]);


      /* r'=|A||x|+|b| */
      CscAxPb( sopalin_data, me, lur2, lub, sopalin_data->sopar->cscmtx,
               &(datacode->updovct), datacode, pastix_comm,
               sopar->iparm[IPARM_TRANSPOSE_SOLVE]);



      /* tmp_berr =  max_i(|lur_i|/|lur2_i|)*/
      CscBerr(sopalin_data, me, lur, lur2, UPDOWN_SM2XSZE,
              1, &tmp_berr , pastix_comm);

      MONOTHREAD_BEGIN;
      sopalin_data->berr = tmp_berr;
      if (sopalin_data->lberr == 0)
        /* force le premier raffinement */
        sopalin_data->lberr = 3*sopalin_data->berr;

      if (SOLV_PROCNUM == 0)
        {
          print_debug(DBG_RAFF_PIVOT, "RAFF : berr lberr %6g %6g\n",
                      sopalin_data->berr, sopalin_data->lberr);
        }
      MONOTHREAD_END;

      /* Calcul de ||r|| et ||r||/||b|| */
      tmp_berr = CscNormErr(sopalin_data, me, lur,lub,
                            UPDOWN_SM2XSZE,UPDOWN_SM2XNBR, pastix_comm);
      MONOTHREAD_BEGIN;
      sopar->rberror = tmp_berr;
      print_debug(DBG_RAFF_PIVOT, "RAFF : rberror %6g\n", sopar->rberror);
      MONOTHREAD_END;

#ifndef SMP_RAFF
      MONOTHREAD_END;
#endif
      SYNCHRO_THREAD;

      if ((sopalin_data->raffnbr < sopar->itermax)
          && (sopalin_data->berr > sopar->epsilonraff)
          && (sopalin_data->berr <= (sopalin_data->lberr/2)))
        {

          MONOTHREAD_BEGIN;
          /* LU dx = r */
          /* lur2 <= updo_vect (ie X_i)
           * updo_vect <= lur (ie B-AX_i)
           */
          SOPALIN_COPY(UPDOWN_SM2XSZE*UPDOWN_SM2XNBR,UPDOWN_SM2XTAB,
                       iun,lur2,iun);
          SOPALIN_COPY(UPDOWN_SM2XSZE*UPDOWN_SM2XNBR,lur,iun,
                       UPDOWN_SM2XTAB,iun);
          MONOTHREAD_END;

#ifdef PRECOND
          if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
            {
              SYNCHRO_THREAD;

              RAFF_CLOCK_STOP;
              t1 = RAFF_CLOCK_GET;

              API_CALL(up_down_smp)(arg);

              SYNCHRO_THREAD;

              RAFF_CLOCK_STOP;
              t2 = RAFF_CLOCK_GET;
            }
#endif

          MONOTHREAD_BEGIN;

          /* updo_vect <= updo_vect (ie PRECOND(B-AX_i)) + lur2 (ie X_i) */
#ifdef MULT_SMX_RAFF
          SOPALIN_GEAM("N","N",UPDOWN_SM2XSZE,UPDOWN_SM2XNBR,fun,lur2,
                       UPDOWN_SM2XSZE,UPDOWN_SM2XTAB,UPDOWN_SM2XSZE);
#else
          SOPALIN_AXPY(UPDOWN_SM2XSZE,fun,lur2,iun,UPDOWN_SM2XTAB,iun);
#endif


          /* lastberr = berr */
          sopalin_data->lberr = sopalin_data->berr;
          sopalin_data->raffnbr++;

          MONOTHREAD_END;
        }
      else
        {
          MONOTHREAD_BEGIN;

          sopalin_data->flag_gmres = 0;

          MONOTHREAD_END;
        }

      SYNCHRO_THREAD;

      RAFF_CLOCK_STOP;
      t3 = RAFF_CLOCK_GET;

      MONOTHREAD_BEGIN;

      if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        {
          double sst, rst = 0.0;
          double stt, rtt;
          double err, berr = sopalin_data->berr;

          stt = t3 - t0;
          if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
            {
              sst = t2-t1;
              MyMPI_Reduce(&sst, &rst, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
            }

          MyMPI_Reduce(&berr, &err, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
          MyMPI_Reduce(&stt,  &rtt, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
          if (SOLV_PROCNUM == 0)
            {
              fprintf(stdout, OUT_ITERRAFF_ITER, (int)sopalin_data->raffnbr);
              if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
                fprintf(stdout, OUT_ITERRAFF_TTS, rst);
              fprintf(stdout, OUT_ITERRAFF_TTT, rtt);
              fprintf(stdout, OUT_ITERRAFF_ERR, err);
            }
        }
      MONOTHREAD_END;

      t0 = t3;
    }

  MONOTHREAD_BEGIN;
  memFree_null(lub);
  memFree_null(lur);
  memFree_null(lur2);
  sopar->itermax = sopalin_data->raffnbr;

  if (THREAD_COMM_ON)
    {
      if (sopar->iparm[IPARM_END_TASK] >= API_TASK_REFINE)
        {
          MUTEX_LOCK(&(sopalin_data->mutex_comm));
          sopalin_data->step_comm = COMMSTEP_END;
          print_debug(DBG_THCOMM, "%s:%d END\n", __FILE__, __LINE__);
          MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
          pthread_cond_broadcast(&(sopalin_data->cond_comm));
        }
    }
#ifdef OOC
  ooc_stop_thread(sopalin_data);
#endif

  RAFF_CLOCK_STOP;
  print_debug(DBG_SOPALIN_RAFF, "%d : refinement time %lf\n", (int)me, clockVal((raff_clk)));
  set_dparm(sopar->dparm, DPARM_RAFF_TIME, clockVal((raff_clk)));

  MONOTHREAD_END;
  SYNCHRO_THREAD;

  return 0;
}

/*
  Function: API_CALL(gmres_smp)

  Function computing GMRES iterative reffinement.

  Parameters:
  arg - Pointer to a <sopthread_data_t> structure containing
  the <Sopalin_Data_t> structure and the thread number ID.
*/
void* API_CALL(gmres_smp)(void *arg)
{
  double            raff_clk;
  double            t0 = 0;
  double            t1 = 0;
  double            t2 = 0;
  double            t3 = 0;
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  d_SolverMatrix     *datacode     = sopalin_data->datacode;
  SopalinParam     *sopar        = sopalin_data->sopar;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t               me           = argument->me;
  pastix_float_t          *  gmrestemp    = NULL;
  volatile pastix_int_t      gmresim      = 0;
  volatile pastix_int_t      gmresmaxits  = 0;
  pastix_float_t          *  gmresb       = NULL;
  pastix_float_t          ** gmresvv      = NULL;
  pastix_float_t          ** gmreshh      = NULL;
  pastix_float_t          *  gmresc       = NULL;
  pastix_float_t          *  gmress       = NULL;
  pastix_float_t          *  gmresrs      = NULL;
  pastix_float_t          ** gmresw       = NULL;
  pastix_float_t             gmresalpha;
  volatile pastix_int_t      gmresincx    = 0;
  volatile pastix_int_t      gmresiters   = 0;
  pastix_float_t          *  gmreswk1;
  pastix_float_t          *  gmreswk2     = NULL;
  pastix_float_t             gmrest;
  volatile double   gmreseps     = 0;
  volatile double   gmresnormb;
  double            gmresnormb2;
  volatile pastix_int_t      gmresi1      = 0;
  volatile pastix_int_t i=0;
  pastix_int_t j,ii,k;
  double beta;

  gmresim     = sopar->gmresim;     /* ? */
  if (gmresim < 1)
    {
      errorPrintW("The Krylov space size must be greater than 1,"
                  " forced to 1");
      gmresim = 1;
    }
  gmresmaxits = sopar->itermax;     /* ? */
  if (gmresmaxits < 1)
    {
      errorPrintW("Iterations number must be greater than 1,"
                  " forced to 1");
      gmresmaxits = 1;
    }
  gmreseps    = sopar->epsilonraff; /* ? */
  gmresim     = MIN(gmresim, gmresmaxits);

  MONOTHREAD_BEGIN;
  if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    {
      if (SOLV_PROCNUM == 0)
        {
          fprintf(stdout, OUT_ITERRAFF_GMRES);
        }
    }
  MALLOC_INTERN(gmrestemp,  UPDOWN_SM2XSZE, pastix_float_t);
  MALLOC_INTERN(gmresb,     UPDOWN_SM2XSZE, pastix_float_t);
  MALLOC_INTERN(gmresc,     gmresim,        pastix_float_t);
  MALLOC_INTERN(gmress,     gmresim,        pastix_float_t);
  MALLOC_INTERN(gmresrs,    gmresim+1,      pastix_float_t);
  MALLOC_INTERN(gmresvv,    gmresim+1,      pastix_float_t*);
  MALLOC_INTERN(gmreshh,    gmresim,        pastix_float_t*);
  MALLOC_INTERN(gmresw,     gmresim,        pastix_float_t*);

  SOPALIN_COPY(UPDOWN_SM2XSZE, sopar->b, iun, gmresb, iun);

  /* Allocation des tableaux */
  for(i=0; i<(gmresim+1); i++)
    {
      MALLOC_INTERN(gmresvv[i], UPDOWN_SM2XSZE, pastix_float_t);
    }

  for (i=0; i<gmresim; i++)
    {
      MALLOC_INTERN(gmreshh[i], gmresim+1,      pastix_float_t);
      MALLOC_INTERN(gmresw[i],  UPDOWN_SM2XSZE, pastix_float_t);
    }

  sopalin_data->gmresout_flag = 1;

  /* mettre SM2XTAB a zero (pour initialiser le gmres avec eps=norm(b)) !!! */
  if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
    for (i=0;i<UPDOWN_SM2XSZE*UPDOWN_SM2XNBR;i++)
      UPDOWN_SM2XTAB[i]=0.0;

  gmresnormb=0.0;
  for (i=0;i<UPDOWN_SM2XSZE;i++)
    gmresnormb+=ABS_FLOAT(gmresb[i])*ABS_FLOAT(gmresb[i]) ;

  MyMPI_Allreduce((void *) &gmresnormb, (void *) &gmresnormb2, 1,
                  MPI_DOUBLE, MPI_SUM, pastix_comm);
  gmresnormb=sqrt(gmresnormb2);

  print_debug(DBG_RAFF_GMRES, "gmresnormb = %e\n", gmresnormb);
  /* On s'assure que tous les threads possèdent l'adresse des
     tableaux alloues */
  sopalin_data->ptr_raff[0] = (void *)gmrestemp;
  sopalin_data->ptr_raff[1] = (void *)gmresb;
  sopalin_data->ptr_raff[2] = (void *)gmresvv;
  sopalin_data->ptr_raff[3] = (void *)gmreshh;
  sopalin_data->ptr_raff[4] = (void *)gmresc;
  sopalin_data->ptr_raff[5] = (void *)gmress;
  sopalin_data->ptr_raff[6] = (void *)gmresrs;
  sopalin_data->ptr_raff[7] = (void *)gmresw;
  sopalin_data->ptr_raff[8] = (void *)&gmresnormb;

  MONOTHREAD_END;
  SYNCHRO_THREAD;
  gmrestemp  = (pastix_float_t *) sopalin_data->ptr_raff[0];
  gmresb     = (pastix_float_t *) sopalin_data->ptr_raff[1];
  gmresvv    = (pastix_float_t **)sopalin_data->ptr_raff[2];
  gmreshh    = (pastix_float_t **)sopalin_data->ptr_raff[3];
  gmresc     = (pastix_float_t *) sopalin_data->ptr_raff[4];
  gmress     = (pastix_float_t *) sopalin_data->ptr_raff[5];
  gmresrs    = (pastix_float_t *) sopalin_data->ptr_raff[6];
  gmresw     = (pastix_float_t **)sopalin_data->ptr_raff[7];
  /* On s'assure que tous les threads aient bien gmresnormb correct */
  gmresnormb = (double)(*((double*)sopalin_data->ptr_raff[8]));

  gmresalpha = -fun;
  gmresincx  = 1;
  gmresiters = 0;

  RAFF_CLOCK_INIT;
  while (sopalin_data->gmresout_flag)
    {
      RAFF_CLOCK_STOP;
      t0 = RAFF_CLOCK_GET;

#ifndef SMP_RAFF
      MONOTHREAD_BEGIN;
#endif

      gmreswk2 = gmresvv[0];

      /* wk2 = A*x */
      CscAx(sopalin_data, me, sopalin_data->sopar->cscmtx,
            UPDOWN_SM2XTAB, gmreswk2,
            datacode, &(datacode->updovct), pastix_comm,
            sopar->iparm[IPARM_TRANSPOSE_SOLVE]);

      MONOTHREAD_BEGIN;
      /* wk2 = b - A*x */
      for (j=0; j<UPDOWN_SM2XSZE; j++)
        {
          gmresvv[0][j] = gmresb[j] - gmresvv[0][j];
        }
      MONOTHREAD_END;

#ifdef SMP_RAFF
      SYNCHRO_THREAD;
#endif /* SMP_RAFF */

      /* ro = vv[0].vv[0] */
      CscGmresBeta(sopalin_data, me, gmresvv[0], gmresvv[0],
                   UPDOWN_SM2XSZE, 1, &beta, pastix_comm);

      MONOTHREAD_BEGIN;
      print_debug(DBG_RAFF_GMRES, "%ld: line %d beta %.20g\n",
                  (long)me, __LINE__, beta);

      sopalin_data->gmresro = sqrt(beta);

      print_debug(DBG_RAFF_GMRES, "gmresro = %e\n",
                  (double)(fabs(sopalin_data->gmresro)));

      MONOTHREAD_END;
#ifndef SMP_RAFF
      MONOTHREAD_END;
#endif

      SYNCHRO_THREAD;
      if ((double)ABS_FLOAT((pastix_float_t)sopalin_data->gmresro) <=
          sopar->epsilonraff)
        {
          sopalin_data->gmresout_flag = 0;
          break;
        }

      gmrest = (pastix_float_t)(1.0/sopalin_data->gmresro);

      MONOTHREAD_BEGIN;

      SOPALIN_SCAL(UPDOWN_SM2XSZE, gmrest, gmresvv[0], gmresincx);

      gmresrs[0] = (pastix_float_t)sopalin_data->gmresro;
      sopalin_data->gmresin_flag = 1;

      MONOTHREAD_END;

      i=-1;
      SYNCHRO_THREAD;

      print_debug(DBG_RAFF_GMRES, "%ld: gmresin_flag %d\n",
                  (long)me, (int)sopalin_data->gmresin_flag);

      while(sopalin_data->gmresin_flag)
        {
          i++;
          gmresi1 = i+1;

          gmreswk1 = gmresvv[i];
          gmreswk2 = gmresw[i];

          MONOTHREAD_BEGIN;

          /* wk2 = M-1wk1 */
          SOPALIN_COPY(UPDOWN_SM2XSZE,UPDOWN_SM2XTAB,iun,gmrestemp,iun);
          SOPALIN_COPY(UPDOWN_SM2XSZE,gmreswk1,iun,UPDOWN_SM2XTAB,iun);

          MONOTHREAD_END;

#ifdef PRECOND
          if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
            {
              RAFF_CLOCK_STOP;
              t1 = RAFF_CLOCK_GET;
              SYNCHRO_THREAD;
              API_CALL(up_down_smp)(arg);
              SYNCHRO_THREAD;
              RAFF_CLOCK_STOP;
              t2 = RAFF_CLOCK_GET;
            }
#endif /* PRECOND */

          MONOTHREAD_BEGIN;

          SOPALIN_COPY(UPDOWN_SM2XSZE,UPDOWN_SM2XTAB,iun,gmreswk2,iun);
          SOPALIN_COPY(UPDOWN_SM2XSZE,gmrestemp,iun,UPDOWN_SM2XTAB,iun);

#ifdef SMP_RAFF
          MONOTHREAD_END;
          SYNCHRO_THREAD;
#endif /* SMP_RAFF */

          gmreswk1 = gmresvv[gmresi1];

          /* vv[i1] = A*wk2 */
          CscAx(sopalin_data, me, sopalin_data->sopar->cscmtx,
                gmreswk2, gmreswk1,
                datacode, &(datacode->updovct), pastix_comm,
                sopar->iparm[IPARM_TRANSPOSE_SOLVE]);


          /* classical gram - schmidt */
          for (j=0; j<=i; j++)
            {
              /* vv[j]*vv[i1] */
              CscGmresBeta(sopalin_data, me, gmresvv[gmresi1], gmresvv[j],
                           UPDOWN_SM2XSZE, 1, &beta, pastix_comm);

              MONOTHREAD_BEGIN;
              print_debug(DBG_RAFF_GMRES, "%ld: line %d beta %.20g\n",
                          (long)me, __LINE__, beta);
              MONOTHREAD_END;

              gmreshh[i][j] = (pastix_float_t)beta;
            }

#ifdef SMP_RAFF
          SYNCHRO_THREAD;
          MONOTHREAD_BEGIN;
#endif /* SMP_RAFF */

          for (j=0;j<=i;j++)
            {
              gmresalpha = -gmreshh[i][j];
#ifdef TYPE_COMPLEX
              print_debug(DBG_RAFF_GMRES,
                          "   %ld: gmresalpha %.20g + %20g * I\n",
                          (long)me, creal(gmresalpha), cimag(gmresalpha));
#else
              print_debug(DBG_RAFF_GMRES,
                          "   %ld: gmresalpha %.20g\n", (long)me, gmresalpha);
#endif
              SOPALIN_AXPY(UPDOWN_SM2XSZE,gmresalpha,gmresvv[j],
                           gmresincx, gmresvv[gmresi1], gmresincx);
            }

#ifdef SMP_RAFF
          MONOTHREAD_END;
          SYNCHRO_THREAD;
#endif /* SMP_RAFF */

          CscGmresBeta(sopalin_data, me, gmresvv[gmresi1],gmresvv[gmresi1],
                       UPDOWN_SM2XSZE, 1, &beta, pastix_comm);

          MONOTHREAD_BEGIN;
          print_debug(DBG_RAFF_GMRES,"   %ld: line %d beta %.20g\n", (long)me,
                      __LINE__, beta);
          MONOTHREAD_END;

#ifdef SMP_RAFF
          MONOTHREAD_BEGIN;
#endif

          gmrest = (pastix_float_t)sqrt(beta);

          gmreshh[i][gmresi1] = gmrest;

          if (ABS_FLOAT(gmrest) > 10e-50)
            {
              gmrest = fun / gmrest;
              SOPALIN_SCAL(UPDOWN_SM2XSZE,gmrest,
                           gmresvv[gmresi1],gmresincx);

            }

          if (i != 0)
            {
              for (j=1; j<=i;j++)
                {
                  gmrest = gmreshh[i][j-1];
#ifdef CPLX
                  gmreshh[i][j-1] = (pastix_float_t)conj(gmresc[j-1])*gmrest +
                    (pastix_float_t)conj(gmress[j-1])*gmreshh[i][j];
#else /* CPLX */
                  gmreshh[i][j-1] =  gmresc[j-1]*gmrest +
                    gmress[j-1]*gmreshh[i][j];
#endif /* CPLX */
                  gmreshh[i][j]   = -gmress[j-1]*gmrest +
                    gmresc[j-1]*gmreshh[i][j];
                }
            }
#ifdef CPLX
          gmrest = (pastix_float_t)csqrt(ABS_FLOAT(gmreshh[i][i]*gmreshh[i][i])+
                                gmreshh[i][gmresi1]*gmreshh[i][gmresi1]);
#else
          gmrest = (pastix_float_t)sqrt(gmreshh[i][i]*gmreshh[i][i]+
                               gmreshh[i][gmresi1]*gmreshh[i][gmresi1]);
#endif
          if (ABS_FLOAT(gmrest) <= sopar->epsilonraff)
            gmrest = (pastix_float_t)sopar->epsilonraff;

          gmresc[i] = gmreshh[i][i]/gmrest;
          gmress[i] = gmreshh[i][gmresi1]/gmrest;
          gmresrs[gmresi1] = -gmress[i]*gmresrs[i];

#ifdef CPLX
          gmresrs[i] = (pastix_float_t)conj(gmresc[i])*gmresrs[i];
          gmreshh[i][i] = (pastix_float_t)conj(gmresc[i])*gmreshh[i][i] +
            gmress[i]*gmreshh[i][gmresi1];
#else
          gmresrs[i] = gmresc[i]*gmresrs[i];
          gmreshh[i][i] = gmresc[i]*gmreshh[i][i] +
            gmress[i]*gmreshh[i][gmresi1];
#endif
          sopalin_data->gmresro = ABS_FLOAT(gmresrs[gmresi1]);

          if ((i+1 >= gmresim) || (sopalin_data->gmresro/gmresnormb <= gmreseps) || (gmresiters >= gmresmaxits))
            {
              sopalin_data->gmresin_flag = 0;
            }
          MONOTHREAD_END;

          gmresiters++;


          RAFF_CLOCK_STOP;
          t3 = RAFF_CLOCK_GET;

          MONOTHREAD_BEGIN;
          if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
            {
              double serr,rerr;
              double sst, rst = 0.0;
              double stt, rtt;

              serr = sopalin_data->gmresro / gmresnormb;
              stt  = t3 - t0;
              if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
                {
                  sst = t2 - t1;
                  MyMPI_Reduce(&sst, &rst, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
                }

              MyMPI_Reduce(&serr, &rerr, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
              MyMPI_Reduce(&stt,  &rtt,  1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);

              if (SOLV_PROCNUM == 0)
                {
                  fprintf(stdout, OUT_ITERRAFF_ITER, (int)gmresiters);
                  if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
                    fprintf(stdout, OUT_ITERRAFF_TTS, rst);
                  fprintf(stdout, OUT_ITERRAFF_TTT, rtt);
                  fprintf(stdout, OUT_ITERRAFF_ERR, rerr);

                  if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
                    {
                      fprintf(stdout, OUT_ITERRAFF_NORMR, sopalin_data->gmresro);
                      fprintf(stdout, OUT_ITERRAFF_NORMB, gmresnormb);
                      fprintf(stdout, OUT_ITERRAFF_BDIVR, sopalin_data->gmresro/gmresnormb);
                    }

                }
            }
          MONOTHREAD_END;
          t0 = t3;
          SYNCHRO_THREAD;
        }

      MONOTHREAD_BEGIN;

      gmresrs[i] = gmresrs[i]/gmreshh[i][i];
      for (ii=2; ii<=i+1; ii++)
        {
          k = i-ii+1;
          gmrest = gmresrs[k];
          for (j=k+1; j<=i; j++)
            {
              gmrest = gmrest - gmreshh[j][k]*gmresrs[j];
            }
          gmresrs[k] = gmrest/gmreshh[k][k];
        }

      for (j=0; j<=i;j++)
        {
          gmrest = gmresrs[j];
          SOPALIN_AXPY(UPDOWN_SM2XSZE, gmrest, gmresw[j], gmresincx, UPDOWN_SM2XTAB, gmresincx);
        }

      if ((sopalin_data->gmresro/gmresnormb<= gmreseps) || (gmresiters >= gmresmaxits))
        {
          sopalin_data->gmresout_flag = 0;
        }

      MONOTHREAD_END;
      SYNCHRO_THREAD;
    }

  MONOTHREAD_BEGIN;
  sopar->rberror = sopalin_data->gmresro/gmresnormb;
  sopar->itermax = gmresiters;

  for (i=0; i<gmresim+1; i++)
    memFree_null(gmresvv[i]);
  for (i=0; i<gmresim; i++)
    {
      memFree_null(gmreshh[i]);
      memFree_null(gmresw[i]);
    }
  memFree_null(gmrestemp);
  memFree_null(gmresb);
  memFree_null(gmresvv);
  memFree_null(gmreshh);
  memFree_null(gmresc);
  memFree_null(gmress);
  memFree_null(gmresrs);
  memFree_null(gmresw);

  if (THREAD_COMM_ON)
    {
      if (sopar->iparm[IPARM_END_TASK] >= API_TASK_REFINE)
        {
          MUTEX_LOCK(&(sopalin_data->mutex_comm));
          sopalin_data->step_comm = COMMSTEP_END;
          print_debug(DBG_THCOMM, "%s:%d END\n", __FILE__, __LINE__);
          MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
          pthread_cond_broadcast(&(sopalin_data->cond_comm));
        }
    }
#ifdef OOC
  ooc_stop_thread(sopalin_data);
#endif

  RAFF_CLOCK_STOP;
  print_debug(DBG_RAFF_GMRES, "%d : refinement time %lf, %d iters, norm %lg =%lg/%lg \n",
              (int)me, clockVal((raff_clk)), (int)gmresiters,  (double)(sopalin_data->gmresro/gmresnormb),
              (double)(sopalin_data->gmresro), (double)gmresnormb);
  set_dparm(sopar->dparm, DPARM_RAFF_TIME, clockVal((raff_clk)));

  MONOTHREAD_END;

  return 0;
}

/*
  Function: API_CALL(grad_smp)

  Refine the solution using conjugate gradian method.

  Parameters:
  arg - Pointer to a <sopthread_data_t> structure containing
  the <Sopalin_Data_t> structure and the thread number ID.

*/
void* API_CALL(grad_smp)(void *arg)
{
  double            raff_clk;
  double            t0 = 0;
  double            t1 = 0;
  double            t2 = 0;
  double            t3 = 0;
  double            tmp = 0.0;
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  d_SolverMatrix     *datacode     = sopalin_data->datacode;
  SopalinParam     *sopar        = sopalin_data->sopar;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  pastix_int_t               me           = argument->me;
  pastix_int_t               i;

  pastix_float_t * gradb = NULL;
  pastix_float_t * gradr = NULL;
  pastix_float_t * gradp = NULL;
  pastix_float_t * gradz = NULL;
  pastix_float_t * grad2 = NULL;

  /* Alpha et Beta ne sont utilisés que par le thread 0 */
  pastix_float_t *alpha = NULL;
  pastix_float_t *beta  = NULL;
  pastix_float_t   tmp_flt;
  MONOTHREAD_BEGIN;
  if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    {
      if (SOLV_PROCNUM == 0)
        {
          fprintf(stdout, OUT_ITERRAFF_GRAD);
        }
    }
  MALLOC_INTERN(gradb, UPDOWN_SM2XNBR*UPDOWN_SM2XSZE, pastix_float_t);
  MALLOC_INTERN(gradr, UPDOWN_SM2XNBR*UPDOWN_SM2XSZE, pastix_float_t);
  MALLOC_INTERN(gradp, UPDOWN_SM2XNBR*UPDOWN_SM2XSZE, pastix_float_t);
  MALLOC_INTERN(gradz, UPDOWN_SM2XNBR*UPDOWN_SM2XSZE, pastix_float_t);
  MALLOC_INTERN(grad2, UPDOWN_SM2XNBR*UPDOWN_SM2XSZE, pastix_float_t);
  MALLOC_INTERN(alpha, UPDOWN_SM2XNBR               , pastix_float_t);
  MALLOC_INTERN(beta,  UPDOWN_SM2XNBR               , pastix_float_t);

  SOPALIN_COPY(UPDOWN_SM2XSZE*UPDOWN_SM2XNBR,sopar->b,iun,gradb,iun);

  /* mettre SM2XTAB a zero !!! */
  if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
    for (i=0;i<UPDOWN_SM2XSZE*UPDOWN_SM2XNBR;i++)
      UPDOWN_SM2XTAB[i]=0.0;

  /* On s'assure que tous les threads connaissent tous les pointeurs */
  sopalin_data->ptr_raff[0] = (void *)gradb;
  sopalin_data->ptr_raff[1] = (void *)gradp;
  sopalin_data->ptr_raff[2] = (void *)gradr;
  sopalin_data->ptr_raff[3] = (void *)gradz;
  sopalin_data->ptr_raff[4] = (void *)grad2;
  MONOTHREAD_END;
  SYNCHRO_THREAD;
  gradb = (pastix_float_t *)sopalin_data->ptr_raff[0];
  gradp = (pastix_float_t *)sopalin_data->ptr_raff[1];
  gradr = (pastix_float_t *)sopalin_data->ptr_raff[2];
  gradz = (pastix_float_t *)sopalin_data->ptr_raff[3];
  grad2 = (pastix_float_t *)sopalin_data->ptr_raff[4];

  RAFF_CLOCK_INIT;

#ifndef SMP_RAFF
  MONOTHREAD_BEGIN;
#endif

  /* r=b-ax */
  CscbMAx(sopalin_data, me, gradr, gradb, sopalin_data->sopar->cscmtx,
          &(datacode->updovct), datacode, pastix_comm,
          sopar->iparm[IPARM_TRANSPOSE_SOLVE]);

  tmp = CscNormErr( sopalin_data, me, gradr, gradb, UPDOWN_SM2XSZE,
                    UPDOWN_SM2XNBR, pastix_comm);

#ifdef SMP_RAFF
  MONOTHREAD_BEGIN;
#endif

  sopalin_data->stop = tmp;
  print_debug(DBG_SOPALIN_RAFF, "raff.c:%d: %ld stop %g\n",__LINE__, (long)me, sopalin_data->stop);

  /* LUz=r */
  /* updo->grad2, p=x */
  /* gradr->updo, x=r */
  SOPALIN_COPY(UPDOWN_SM2XNBR*UPDOWN_SM2XSZE,UPDOWN_SM2XTAB,iun,gradp,iun);
  SOPALIN_COPY(UPDOWN_SM2XNBR*UPDOWN_SM2XSZE,gradr,iun,UPDOWN_SM2XTAB,iun);

  MONOTHREAD_END;

  /* M-1 updo -> updo */
#ifdef PRECOND
  if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
    {
      RAFF_CLOCK_STOP;
      t1 = RAFF_CLOCK_GET;
      SYNCHRO_THREAD;
      API_CALL(up_down_smp)(arg);
      SYNCHRO_THREAD;
      RAFF_CLOCK_STOP;
      t2 = RAFF_CLOCK_GET;
      MONOTHREAD_BEGIN;
      print_debug(DBG_RAFF_GRAD, "up down time : %.4g\n", t2-t1);
      MONOTHREAD_END;
    }
#endif

  MONOTHREAD_BEGIN;

  /* updo->gradz, z = M-1 r */
  /* gradp->updo, x = x */
  /* gradz->gradp, p = z */
  SOPALIN_COPY(UPDOWN_SM2XNBR*UPDOWN_SM2XSZE,UPDOWN_SM2XTAB,iun,gradz,iun);
  SOPALIN_COPY(UPDOWN_SM2XNBR*UPDOWN_SM2XSZE,gradp,iun,UPDOWN_SM2XTAB,iun);
  SOPALIN_COPY(UPDOWN_SM2XNBR*UPDOWN_SM2XSZE,gradz,iun,gradp,iun);

  MONOTHREAD_END;
  SYNCHRO_THREAD;

  while ((sopalin_data->stop > sopar->epsilonraff) && (sopalin_data->count_iter < sopar->itermax))
    {
      RAFF_CLOCK_STOP;
      t0 = RAFF_CLOCK_GET;

      MONOTHREAD_BEGIN;

      sopalin_data->count_iter++;

      if (SOLV_PROCNUM==0)
        {
          print_debug(DBG_SOPALIN_RAFF,"CG : iter %ld stop %e\n",
                      (long) sopalin_data->count_iter, sopalin_data->stop);
        }

#ifdef SMP_RAFF
      MONOTHREAD_END;
#endif
      /* grad2 = Ap */
      CscAx( sopalin_data, me, sopalin_data->sopar->cscmtx, gradp, grad2,
             datacode, &(datacode->updovct), pastix_comm,
             sopar->iparm[IPARM_TRANSPOSE_SOLVE]);
      /* alpha = <r,z>/<Ap,p> */

      CscGradBeta(sopalin_data, me, gradr, gradz,
                  UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, beta, pastix_comm);
      CscGradBeta(sopalin_data, me, grad2, gradp,
                  UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, alpha, pastix_comm);
      /*
        CscGradAlpha(sopalin_data, me, gradr, gradz, grad2, gradp,
        UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, alpha, pastix_comm);
      */

#ifdef SMP_RAFF
      MONOTHREAD_BEGIN;
#endif

#ifdef MULT_SMX_RAFF
      {
        pastix_int_t itersmx;
        for(itersmx=0; itersmx<UPDOWN_SM2XNBR;itersmx++)
          {
            alpha[itersmx]=beta[itersmx]/alpha[itersmx];
          }
      }
#else
      alpha[0]=beta[0]/alpha[0];
#endif

#ifdef TYPE_COMPLEX
      print_debug(DBG_RAFF_GRAD, "%ld alpha %g %g\n", (long)me,
                  creal(alpha[0]), cimag(alpha[0]));
#else
      print_debug(DBG_RAFF_GRAD, "%ld alpha %g\n", (long)me, alpha[0]);
#endif
      /* x = x+alphap */
#ifdef MULT_SMX_RAFF
      {
        pastix_int_t itersmx;
        for(itersmx=0; itersmx<UPDOWN_SM2XNBR;itersmx++)
          {
            tmp_flt = (pastix_float_t) alpha[itersmx];
            SOPALIN_AXPY(UPDOWN_SM2XSZE, tmp_flt,
                         gradp+(itersmx*UPDOWN_SM2XSZE),iun,
                         UPDOWN_SM2XTAB+(itersmx*UPDOWN_SM2XSZE),iun);
          }
      }
#else
      tmp_flt = (pastix_float_t) alpha[0];
      SOPALIN_AXPY(UPDOWN_SM2XSZE, tmp_flt,gradp,iun,UPDOWN_SM2XTAB,iun);
#endif

      /* r = r-alphaAp */
#ifdef MULT_SMX_RAFF
      {
        pastix_int_t itersmx;

        for(itersmx=0; itersmx<UPDOWN_SM2XNBR;itersmx++)
          {
            tmp_flt = (pastix_float_t)-alpha[itersmx];
            SOPALIN_AXPY(UPDOWN_SM2XSZE, tmp_flt,
                         grad2+(itersmx*UPDOWN_SM2XSZE),iun,
                         gradr+(itersmx*UPDOWN_SM2XSZE),iun);
          }
      }
#else
      tmp_flt = (pastix_float_t)-alpha[0];
      SOPALIN_AXPY(UPDOWN_SM2XSZE, tmp_flt,grad2,iun,gradr,iun);
#endif

      /* z = M-1 r */
      /* updo->grad2, r-> updo */
      SOPALIN_COPY(UPDOWN_SM2XNBR*UPDOWN_SM2XSZE,UPDOWN_SM2XTAB,iun,grad2,iun);
      SOPALIN_COPY(UPDOWN_SM2XNBR*UPDOWN_SM2XSZE,gradr,iun,UPDOWN_SM2XTAB,iun);

      MONOTHREAD_END;

#ifdef PRECOND
      if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
        {

          RAFF_CLOCK_STOP;
          t1 = RAFF_CLOCK_GET;
          SYNCHRO_THREAD;
          API_CALL(up_down_smp)(arg);
          SYNCHRO_THREAD;
          RAFF_CLOCK_STOP;
          t2 = RAFF_CLOCK_GET;
          MONOTHREAD_BEGIN;
          print_debug(DBG_RAFF_GRAD, "up down time : %.4g\n", t2-t1);
          MONOTHREAD_END;
        }
#endif

      MONOTHREAD_BEGIN;

      /* grad2->updo */
      SOPALIN_COPY(UPDOWN_SM2XNBR*UPDOWN_SM2XSZE,UPDOWN_SM2XTAB,iun,gradz,iun);
      SOPALIN_COPY(UPDOWN_SM2XNBR*UPDOWN_SM2XSZE,grad2,iun,UPDOWN_SM2XTAB,iun);

#ifdef SMP_RAFF
      MONOTHREAD_END;
      SYNCHRO_THREAD;
#endif

      CscGradBeta(sopalin_data, me, gradr, gradz,
                  UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, alpha, pastix_comm);

#ifdef SMP_RAFF
      MONOTHREAD_BEGIN;
#endif

#ifdef MULT_SMX_RAFF
      {
        pastix_int_t itersmx;
        for(itersmx=0; itersmx<UPDOWN_SM2XNBR;itersmx++)
          {
            beta[itersmx]=alpha[itersmx]/beta[itersmx];
          }
      }
#else
      beta[0]=alpha[0]/beta[0];
#endif

      /* p = z+beta p */
#ifdef MULT_SMX_RAFF
      {
        pastix_int_t itersmx;
        for (itersmx=0; itersmx<UPDOWN_SM2XNBR; itersmx++)
          {
            tmp_flt = beta[itersmx];
            SOPALIN_SCAL(UPDOWN_SM2XSZE,tmp,
                         gradp+(itersmx*UPDOWN_SM2XSZE),iun);
          }
      }
      SOPALIN_GEAM("N","N",UPDOWN_SM2XSZE,UPDOWN_SM2XNBR, fun,
                   gradz, UPDOWN_SM2XSZE, gradp, UPDOWN_SM2XSZE);
#else
      tmp_flt = (pastix_float_t)beta[0];
      SOPALIN_SCAL(UPDOWN_SM2XSZE,beta[0],gradp,iun);
      SOPALIN_AXPY(UPDOWN_SM2XSZE,fun,gradz,iun,gradp,iun);
#endif

#ifdef SMP_RAFF
      MONOTHREAD_END;
      SYNCHRO_THREAD;
#endif

      tmp = CscNormErr(sopalin_data, me, gradr, gradb, UPDOWN_SM2XSZE, UPDOWN_SM2XNBR, pastix_comm);

#ifndef SMP_RAFF
      MONOTHREAD_END;
#endif

      RAFF_CLOCK_STOP;
      t3 = RAFF_CLOCK_GET;

      MONOTHREAD_BEGIN;
      sopalin_data->stop = tmp;
      print_debug(DBG_RAFF_GRAD, "raff.c:%d: %ld stop %g\n",__LINE__, (long)me, sopalin_data->stop);

      if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        {
          double sst, rst = 0.0;
          double stt, rtt;
          double err, stop = sopalin_data->stop;

          if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
            {
              sst = t2-t1;
              MyMPI_Reduce(&sst, &rst, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
            }

          stt = t3-t0;
          MyMPI_Reduce(&stop, &err, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
          MyMPI_Reduce(&stt,  &rtt, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);

          if (SOLV_PROCNUM == 0)
            {
              fprintf(stdout, OUT_ITERRAFF_ITER, (int)sopalin_data->count_iter);
              if (sopar->iparm[IPARM_ONLY_RAFF] == API_NO)
                fprintf(stdout, OUT_ITERRAFF_TTS, rst);
              fprintf(stdout, OUT_ITERRAFF_TTT, rtt);
              fprintf(stdout, OUT_ITERRAFF_ERR, err);
            }
        }
      MONOTHREAD_END;
      t0 = t3;
      SYNCHRO_THREAD;
    }

  MONOTHREAD_BEGIN;
  sopar->rberror = sopalin_data->stop;
  sopar->itermax = sopalin_data->count_iter;

  memFree_null(gradb);
  memFree_null(gradr);
  memFree_null(gradp);
  memFree_null(gradz);
  memFree_null(grad2);
  memFree_null(alpha);
  memFree_null(beta);

  if (THREAD_COMM_ON)
    {
      if (sopar->iparm[IPARM_END_TASK] >= API_TASK_REFINE)
        {
          MUTEX_LOCK(&(sopalin_data->mutex_comm));
          sopalin_data->step_comm = COMMSTEP_END;
          print_debug(DBG_THCOMM, "%s:%d END\n", __FILE__, __LINE__);
          MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
          pthread_cond_broadcast(&(sopalin_data->cond_comm));
        }
    }
#ifdef OOC
  ooc_stop_thread(sopalin_data);
#endif

  RAFF_CLOCK_STOP;
  set_dparm(sopar->dparm, DPARM_RAFF_TIME, clockVal((raff_clk)));

  MONOTHREAD_END;

  return 0;
}
/*
** Section: Function creating threads
*/

/*
  Function: API_CALL(pivot_thread)

  Launch sopaparam->nbthrdcomm threads which will compute
  <API_CALL(pivotstatique_smp)>.

  Parameters:
  datacode  - PaStiX <d_SolverMatrix> structure.
  sopaparam - <SopalinParam> parameters structure.
*/
void API_CALL(pivot_thread)(d_SolverMatrix *datacode,
                            SopalinParam *sopaparam)
{
  Sopalin_Data_t *sopalin_data = NULL;
  BackupSolve_t b;

  MALLOC_INTERN(sopalin_data, 1, Sopalin_Data_t);

  solve_backup(datacode,&b);
  sopalin_init(sopalin_data, datacode, sopaparam, 0);

  sopalin_launch_thread(sopalin_data,
                        SOLV_PROCNUM,          SOLV_PROCNBR,                datacode->btree,
                        sopalin_data->sopar->iparm[IPARM_VERBOSE],
                        SOLV_THRDNBR,          API_CALL(pivotstatique_smp), sopalin_data,
                        sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm), sopalin_data,
                        OOC_THREAD_NBR,        ooc_thread,                  sopalin_data);

  sopalin_clean(sopalin_data, 2);
  solve_restore(datacode,&b);
  memFree_null(sopalin_data);
}

/*
  Function: API_CALL(gmres_thread)

  Launch sopaparam->nbthrdcomm threads which will compute
  <API_CALL(gmres_smp)>.

  Parameters:
  datacode  - PaStiX <d_SolverMatrix> structure.
  sopaparam - <SopalinParam> parameters structure.
*/
void API_CALL(gmres_thread)(d_SolverMatrix *datacode,
                            SopalinParam *sopaparam)
{
  Sopalin_Data_t *sopalin_data = NULL;
  BackupSolve_t b;

  MALLOC_INTERN(sopalin_data, 1, Sopalin_Data_t);

  solve_backup(datacode,&b);
  sopalin_init(sopalin_data, datacode, sopaparam, 0);

  sopalin_launch_thread(sopalin_data,
                        SOLV_PROCNUM,          SOLV_PROCNBR,                datacode->btree,
                        sopalin_data->sopar->iparm[IPARM_VERBOSE],
                        SOLV_THRDNBR,          API_CALL(gmres_smp),         sopalin_data,
                        sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm), sopalin_data,
                        OOC_THREAD_NBR,        ooc_thread,                  sopalin_data);

  sopalin_clean(sopalin_data, 2);
  solve_restore(datacode,&b);
  memFree_null(sopalin_data);
}

/*
  Function: API_CALL(grad_thread)

  Launch sopaparam->nbthrdcomm threads which will compute
  <API_CALL(grad_smp)>.

  Parameters:
  datacode  - PaStiX <d_SolverMatrix> structure.
  sopaparam - <SopalinParam> parameters structure.
*/
void API_CALL(grad_thread)(d_SolverMatrix *datacode,
                           SopalinParam *sopaparam)
{
  Sopalin_Data_t *sopalin_data = NULL;
  BackupSolve_t b;

  MALLOC_INTERN(sopalin_data, 1, Sopalin_Data_t);

  solve_backup(datacode,&b);
  sopalin_init(sopalin_data, datacode, sopaparam, 0);

  sopalin_launch_thread(sopalin_data,
                        SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree,
                        sopalin_data->sopar->iparm[IPARM_VERBOSE],
                        SOLV_THRDNBR,          API_CALL(grad_smp),          sopalin_data,
                        sopaparam->nbthrdcomm, API_CALL(sopalin_updo_comm), sopalin_data,
                        OOC_THREAD_NBR,        ooc_thread,                  sopalin_data);

  sopalin_clean(sopalin_data, 2);
  solve_restore(datacode,&b);
  memFree_null(sopalin_data);
}
