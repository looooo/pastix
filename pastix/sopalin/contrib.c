/* Attentes des contributions locale et MPI */

#ifdef RECV_FANIN_OR_BLOCK
#define RECV_ONE_FANIN API_CALL(recv_waitone_fob)(sopalin_data, me)
#define RECV_ONE_BLOCK API_CALL(recv_waitone_fob)(sopalin_data, me)
#else
#define RECV_ONE_FANIN API_CALL(recv_waitone_fanin)(sopalin_data, me, TASK_PRIONUM(i))
#define RECV_ONE_BLOCK API_CALL(recv_waitone_block)(sopalin_data, me, TASK_PRIONUM(i))
#endif

#ifndef X_ARCHi686_pc_linux
#define inline
#endif

static inline void API_CALL(wait_contrib_comp_1d)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT i){

  SolverMatrix  *datacode    = sopalin_data->datacode;
#ifdef TRACE_SOPALIN
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif

#if (defined FORCE_CONSO)
  if (THREAD_FUNNELED_OFF)
    {
      /* Attente en Multiple / force_conso */
      while(TASK_CTRBCNT(i))
        {
          API_CALL(rcsd_testall_fab)(sopalin_data, me);
        }
    }
  else
#endif
    {
      if (THREAD_COMM_OFF)
        {
          /* Attente en multiple sans force conso */
          while(TASK_FTGTCNT(i))
            {
              RECV_ONE_FANIN;
            }
        }
    }
  trace_begin_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                   STATE_WAITLOC, i);

  MUTEX_LOCK(&(sopalin_data->mutex_task[i]));
#if (DBG_PASTIX_DYNSCHED > 0)
  ASSERTDBG(sopalin_data->taskmark[i] == 0, MOD_SOPALIN);
  sopalin_data->taskmark[i]++;
#endif
  while (TASK_CTRBCNT(i))
  {
    COND_WAIT(&(sopalin_data->cond_task[i]), &(sopalin_data->mutex_task[i]));
  }
  MUTEX_UNLOCK(&(sopalin_data->mutex_task[i]));

}


static inline
void API_CALL(wait_contrib_comp_2d)(Sopalin_Data_t *sopalin_data,
                                    PASTIX_INT me, PASTIX_INT i){

  SolverMatrix  *datacode    = sopalin_data->datacode;
#ifdef TRACE_SOPALIN
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif
#ifdef SMP_SOPALIN
  PASTIX_INT            firsttask   = TASK_MASTER(i);
#endif

  /* Attente contribution locale et MPI */
#if (defined FORCE_CONSO)
  if (THREAD_FUNNELED_OFF)
    {
      while ((!(TASK_BTAGPTR(i)))
             || (!(RTASK_COEFTAB(i))))
        {
          API_CALL(rcsd_testall_fab)(sopalin_data, me);
        }
    }
  else
#endif
    {
      if (THREAD_COMM_OFF)
        {
          while ((TASK_BTAGPTR(i) == NULL) && (sopalin_data->taskmark[i] > 0))
            {
              ASSERTDBG(i == firsttask, MOD_SOPALIN);
              RECV_ONE_BLOCK;
            }
        }
    }

  trace_begin_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                   STATE_WAITLOC, i);

#if (DBG_PASTIX_DYNSCHED > 0)
  MUTEX_LOCK(&(sopalin_data->mutex_task[i]));
  ASSERTDBG(((TASK_TASKID(i) == E2) && (sopalin_data->taskmark[i] == 0))
            || ((TASK_TASKID(i) == E1) &&
                (sopalin_data->taskmark[i] == 1)), MOD_SOPALIN);
  ASSERTDBG(TASK_BTAGPTR(i)  != NULL, MOD_SOPALIN);
  ASSERTDBG(RTASK_COEFTAB(i) != NULL, MOD_SOPALIN);
  sopalin_data->taskmark[i]++;
  MUTEX_UNLOCK(&(sopalin_data->mutex_task[i]));
#endif

  MUTEX_LOCK(&(sopalin_data->mutex_task[firsttask]));
  while ((!(TASK_BTAGPTR(i)))
         || (!(RTASK_COEFTAB(i))))
  COND_WAIT(&(sopalin_data->cond_task[firsttask]),
            &(sopalin_data->mutex_task[firsttask]));
  MUTEX_UNLOCK(&(sopalin_data->mutex_task[firsttask]));
}
