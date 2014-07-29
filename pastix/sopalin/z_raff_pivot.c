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
** Section: Functions declarations
*/

/* Raffinement du second membre */
#define z_pivotstatique_smp API_CALL(z_pivotstatique_smp)
#define z_pivot_thread      API_CALL(z_pivot_thread)

void* z_pivotstatique_smp(void *arg);

/* Lancement d'une des fonctions seules */
void z_pivot_thread(z_SolverMatrix *datacode, z_SopalinParam *sopaparam);

/*
** Section: Threads routines
*/

/*
 * Function: API_CALL(z_pivotstatique_smp)
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
 *         the <z_Sopalin_Data_t> structure and the thread number ID.
 */
void* z_pivotstatique_smp ( void *arg )
{
  Clock             raff_clk;
  double            t0           = 0;
  double            t1           = 0;
  double            t2           = 0;
  double            t3           = 0;
  pastix_complex64_t * volatile  lub          = NULL;
  pastix_complex64_t * volatile  lur          = NULL;
  pastix_complex64_t * volatile  lur2         = NULL;
  double            tmp_berr     = 0.0;
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  z_SolverMatrix     *datacode     = sopalin_data->datacode;
  z_SopalinParam     *sopar        = sopalin_data->sopar;
  MPI_Comm          pastix_comm  = PASTIX_COMM;
  PASTIX_INT               me           = argument->me;
  int               iter         = 0;

  MONOTHREAD_BEGIN;
  if (sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    {
      if (SOLV_PROCNUM == 0)
        {
          fprintf(stdout, OUT_ITERRAFF_PIVOT);
        }
    }
  MALLOC_INTERN(lub,  UPDOWN_SM2XNBR*UPDOWN_SM2XSZE, pastix_complex64_t);
  MALLOC_INTERN(lur,  UPDOWN_SM2XNBR*UPDOWN_SM2XSZE, pastix_complex64_t);
  MALLOC_INTERN(lur2, UPDOWN_SM2XNBR*UPDOWN_SM2XSZE, pastix_complex64_t);

  SOPALIN_COPY(UPDOWN_SM2XNBR*UPDOWN_SM2XSZE,sopar->b,iun,lub,iun);

  sopalin_data->ptr_raff[0] = (void *)lur;
  sopalin_data->ptr_raff[1] = (void *)lub;
  sopalin_data->ptr_raff[2] = (void *)lur2;

  MONOTHREAD_END;
  SYNCHRO_THREAD;

  lur  = (pastix_complex64_t *)sopalin_data->ptr_raff[0];
  lub  = (pastix_complex64_t *)sopalin_data->ptr_raff[1];
  lur2 = (pastix_complex64_t *)sopalin_data->ptr_raff[2];

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
      z_CscbMAx(sopalin_data, me, lur, lub, sopalin_data->sopar->cscmtx,
              &(datacode->updovct), datacode, pastix_comm,
              sopar->iparm[IPARM_TRANSPOSE_SOLVE]);


      /* r'=|A||x|+|b| */
      z_CscAxPb( sopalin_data, me, lur2, lub, sopalin_data->sopar->cscmtx,
               &(datacode->updovct), datacode, pastix_comm,
               sopar->iparm[IPARM_TRANSPOSE_SOLVE]);



      /* tmp_berr =  max_i(|lur_i|/|lur2_i|)*/
      z_CscBerr(sopalin_data, me, lur, lur2, UPDOWN_SM2XSZE,
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
      tmp_berr = z_CscNormErr(sopalin_data, me, lur,lub,
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
          SYNCHRO_THREAD;

          RAFF_CLOCK_STOP;
          t1 = RAFF_CLOCK_GET;

          API_CALL(z_up_down_smp)(arg);

          SYNCHRO_THREAD;

          RAFF_CLOCK_STOP;
          t2 = RAFF_CLOCK_GET;
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
          sst = t2-t1;
          MyMPI_Reduce(&sst, &rst, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);

          MyMPI_Reduce(&berr, &err, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
          MyMPI_Reduce(&stt,  &rtt, 1, MPI_DOUBLE, MPI_MAX, 0, pastix_comm);
          if (SOLV_PROCNUM == 0)
            {
              fprintf(stdout, OUT_ITERRAFF_ITER, (int)sopalin_data->raffnbr);
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
  z_ooc_stop_thread(sopalin_data);
#endif

  RAFF_CLOCK_STOP;
  print_debug(DBG_SOPALIN_RAFF, "%d : refinement time %lf\n", (int)me, RAFF_CLOCK_GET);
  sopar->dparm[DPARM_RAFF_TIME] = RAFF_CLOCK_GET;

  MONOTHREAD_END;
  SYNCHRO_THREAD;

  return 0;
}


/*
** Section: Function creating threads
*/

/*
  Function: API_CALL(z_pivot_thread)

  Launch sopaparam->nbthrdcomm threads which will compute
  <API_CALL(z_pivotstatique_smp)>.

  Parameters:
  datacode  - PaStiX <z_SolverMatrix> structure.
  sopaparam - <z_SopalinParam> parameters structure.
*/
void API_CALL(z_pivot_thread)(z_SolverMatrix *datacode,
                            z_SopalinParam *sopaparam)
{
  z_Sopalin_Data_t *sopalin_data = NULL;
  BackupSolve_t b;

  MALLOC_INTERN(sopalin_data, 1, z_Sopalin_Data_t);

  z_solve_backup(datacode,&b);
  z_sopalin_init(sopalin_data, datacode, sopaparam, 0);

  sopalin_launch_thread(sopalin_data,
                        SOLV_PROCNUM,          SOLV_PROCNBR,                datacode->btree,
                        sopalin_data->sopar->iparm[IPARM_VERBOSE],
                        SOLV_THRDNBR,          API_CALL(z_pivotstatique_smp), sopalin_data,
                        sopaparam->nbthrdcomm, API_CALL(z_sopalin_updo_comm), sopalin_data,
                        OOC_THREAD_NBR,        z_ooc_thread,                  sopalin_data);

  z_sopalin_clean(sopalin_data, 2);
  z_solve_restore(datacode,&b);
  memFree_null(sopalin_data);
}
