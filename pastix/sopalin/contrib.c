/* Attentes des contributions locale et MPI */

#define RECV_ONE_FANIN API_CALL(recv_waitone_fanin)(sopalin_data, me, TASK_PRIONUM(i))

static inline void API_CALL(wait_contrib_comp_1d)(Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t i){

  d_SolverMatrix  *datacode    = sopalin_data->datacode;
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
