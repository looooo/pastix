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
/*************************************/
/*          RECEIVE FUNCTIONS        */
/*************************************/
/* Handle received data */
#define recv_handle_fanin API_CALL(recv_handle_fanin)
void  recv_handle_fanin (z_Sopalin_Data_t *, pastix_int_t, void *buffer,
                         MPI_Status status, int elected);

/* Wait for one reception */
#define z_recv_waitone_fanin API_CALL(z_recv_waitone_fanin)
#define recv_waitone_fob   API_CALL(recv_waitone_fob)
void  z_recv_waitone_fanin(z_Sopalin_Data_t *, pastix_int_t, pastix_int_t tag);
void  recv_waitone_fob  (z_Sopalin_Data_t *, pastix_int_t);

/* Test reception */
#define recv_testone_fob API_CALL(recv_testone_fob)
#define recv_testall_fab API_CALL(recv_testall_fab)
void  recv_testone_fob  (z_Sopalin_Data_t *, pastix_int_t);
void  recv_testall_fab  (z_Sopalin_Data_t *, pastix_int_t);

/*************************************/
/*          SEND FUNCTIONS           */
/*************************************/
/* Send one communication */
#define send_one_fanin API_CALL(send_one_fanin)
int   send_one_fanin    (z_Sopalin_Data_t *, pastix_int_t, pastix_int_t t);

/* Send all available communications */
#define send_all_fanin API_CALL(send_all_fanin)
void  send_all_fanin    (z_Sopalin_Data_t *, pastix_int_t, pastix_int_t dest);

/* Free data structure associate with send request */
#define send_free_fanin API_CALL(send_free_fanin)
void  send_free_fanin   (z_Sopalin_Data_t *, pastix_int_t, pastix_int_t s_index);

/* Test and wait send requests */
#define send_testall       API_CALL(send_testall)
#define send_testall_fanin API_CALL(send_testall_fanin)
#define send_testall_fab   API_CALL(send_testall_fab)
#define send_waitone       API_CALL(send_waitone)
#define send_waitone_fanin API_CALL(send_waitone_fanin)
#define z_send_waitall_fab   API_CALL(z_send_waitall_fab)
void  send_testall      (z_Sopalin_Data_t *,
                         pastix_int_t, void (*funcfree)(z_Sopalin_Data_t*,
                                               pastix_int_t, pastix_int_t));
void  send_testall_fanin(z_Sopalin_Data_t *, pastix_int_t);
void  send_testall_fab  (z_Sopalin_Data_t *, pastix_int_t);
int   send_waitone      (z_Sopalin_Data_t *, pastix_int_t,
                         void (*funcfree)(z_Sopalin_Data_t*, pastix_int_t, pastix_int_t));
int   send_waitone_fanin(z_Sopalin_Data_t *, pastix_int_t);
void  z_send_waitall_fab  (z_Sopalin_Data_t *, pastix_int_t);

/* Test all receive and send requests */
#define z_rcsd_testall_fab API_CALL(z_rcsd_testall_fab)
void  z_rcsd_testall_fab  (z_Sopalin_Data_t *, pastix_int_t);

/* Fonction pour thread de comm */
#define sendrecv_smp API_CALL(sendrecv_smp)
void* sendrecv_smp(void *arg);

#ifdef FORCE_NOMPI
#define z_recv_waitone_fanin API_CALL(z_recv_waitone_fanin)
#define recv_waitone_fob   API_CALL(recv_waitone_fob)
#define recv_testone_fob   API_CALL(recv_testone_fob)
#define recv_testall_fab   API_CALL(recv_testall_fab)
#define send_all_fanin     API_CALL(send_all_fanin)
#define send_free_fanin    API_CALL(send_free_fanin)
#define z_rcsd_testall_fab   API_CALL(z_rcsd_testall_fab)
#define sendrecv_smp       API_CALL(sendrecv_smp)
void  z_recv_waitone_fanin(z_Sopalin_Data_t *sopalin_data,
                         pastix_int_t me, pastix_int_t tag)
{
  (void)sopalin_data; (void)me; (void)tag;
}
void  recv_waitone_fob  (z_Sopalin_Data_t *sopalin_data, pastix_int_t me)
{
  (void)sopalin_data; (void)me;
}
void  recv_testone_fob  (z_Sopalin_Data_t *sopalin_data, pastix_int_t me)
{
  (void)sopalin_data; (void)me;
}
void  recv_testall_fab  (z_Sopalin_Data_t *sopalin_data, pastix_int_t me)
{
  (void)sopalin_data; (void)me;
}
void  send_all_fanin    (z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t dest)
{
  (void)sopalin_data; (void)me; (void)dest;
}
void  send_free_fanin   (z_Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t s_index)
{
  (void)sopalin_data; (void)me; (void)s_index;
}
void  z_rcsd_testall_fab  (z_Sopalin_Data_t *sopalin_data, pastix_int_t me)
{
  (void)sopalin_data; (void)me;
}
void* sendrecv_smp(void *arg){ (void)arg; return NULL;}
#else

/****************************************************************************/
/* RECEIVE ROUTINES                                                         */
/****************************************************************************/

/*
 * Function: recv_handle_fanin
 *
 * Add fanin contribution received in recv_buffer.
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *    recv_buffer  - Received data
 *    status       - MPI communication status
 *    elected      - Index of communication used
 *
 */
void
recv_handle_fanin(z_Sopalin_Data_t *sopalin_data,
                  pastix_int_t             me,
                  void           *recv_buffer,
                  MPI_Status      status,
                  int             elected)
{
  z_SolverMatrix  *datacode    = sopalin_data->datacode;
#  if (defined TRACE_SOPALIN) || (defined TEST_IRECV)
  z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#  endif
  void  *buffer      = recv_buffer;
  pastix_complex64_t *ga,  *gb;
#  ifdef SOPALIN_LU
  pastix_complex64_t *ga2, *gb2;
#  endif
  pastix_int_t    packnbr, pack, ind;
  pastix_int_t    ctrbnbr, taskdst, cblkdst, blokdst, stride, str, dimi, dimj;
  pastix_int_t    fcolnum, lcolnum, frownum, lrownum, ofrownum, olrownum;
#ifdef TEST_IRECV
  pastix_int_t    size;
#endif
  int    flag = 0;

  print_debug(DBG_SOPALIN_RECV,"%ld: receive fanin target\n",(long)me);

#ifdef TEST_IRECV
  size = PACKMAX*(sizeof(pastix_int_t)*MAXINFO)+PACKAREA*sizeof(pastix_complex64_t);
#endif
  taskdst = ((pastix_int_t*)buffer)[FTGT_TASKDST];
  packnbr = ((pastix_int_t*)buffer)[FTGT_PRIONUM];

  trace_begin_task(thread_data->tracefile,
                   SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                   STATE_L2_RECVF, taskdst);
  trace_recv(thread_data->tracefile,
             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, status.MPI_SOURCE,
             COMM_FANIN, taskdst, 0, ((pastix_int_t*)buffer)[FTGT_IDTRACE]);

  /* Loop on all ftgt */
  for (pack=0;pack<packnbr;pack++)
    {
      /* add fanintarget */
      ctrbnbr = ((pastix_int_t*)buffer)[FTGT_CTRBNBR];
      taskdst = ((pastix_int_t*)buffer)[FTGT_TASKDST];
      blokdst = ((pastix_int_t*)buffer)[FTGT_BLOKDST];
      fcolnum = ((pastix_int_t*)buffer)[FTGT_FCOLNUM];
      lcolnum = ((pastix_int_t*)buffer)[FTGT_LCOLNUM];
      frownum = ((pastix_int_t*)buffer)[FTGT_FROWNUM];
      lrownum = ((pastix_int_t*)buffer)[FTGT_LROWNUM];
      dimj    = lcolnum - fcolnum + 1;
      dimi    = lrownum - frownum + 1;

#  ifdef DRUNK_SOPALIN
      if (taskdst == -DRUNK)
        {
          pastix_int_t c;

          /* find cblkdst without taskdst */
          taskdst = SOLV_TASKNBR-1;
          for (c=TASK_CBLKNUM(taskdst); c<SYMB_CBLKNBR; c++)
            if ((blokdst>=SYMB_BLOKNUM(c)) && (blokdst<SYMB_BLOKNUM(c+1)))
              break;

          ASSERTDBG(c!=SYMB_CBLKNBR,MOD_SOPALIN);

          cblkdst = c;

          print_debug(DBG_SOPALIN_DRUNK, "add on DRUNK ctrncnt=%ld\n",
                      (long)TASK_CTRBCNT(SOLV_TASKNBR-1));
          print_debug(DBG_SOPALIN_DRUNK, "cblkdst blokdst DRUNK = %ld %ld \n",
                      (long)c, (long)blokdst);
          print_debug(DBG_SOPALIN_DRUNK, "stride dimi dimj"
                      " DRUNK = %ld %ld %ld \n",
                      (long)stride, (long)dimi, (long)dimj);
        }
      else
#  endif /* DRUNK_SOPALIN */
        {
          cblkdst = TASK_CBLKNUM(taskdst);
        }

      stride  = SOLV_STRIDE(cblkdst);

      /* Load unpredicted cblk */
      z_ooc_hack_load(sopalin_data, cblkdst, me);
      z_ooc_wait_for_cblk(sopalin_data, cblkdst, me);

      ind = SOLV_COEFIND(blokdst)
        + (fcolnum - SYMB_FCOLNUM(cblkdst)) * stride
        +  frownum - SYMB_FROWNUM(blokdst);

      ga  = &(SOLV_COEFTAB(cblkdst)[ind]);
#  ifdef SOPALIN_LU
      ga2 = &(SOLV_UCOEFTAB(cblkdst)[ind]);
#  endif

      gb  =(pastix_complex64_t *)(((char *) buffer) + MAXINFO*sizeof(pastix_int_t));

      print_debug(DBG_SOPALIN_RECV,
                  "%ld: Recv fanintarget\n"
                  "%ld: ctrbnbr ctrbcnt procdst taskdst blokdst prionum"
                  " fcolnum lcolnum frownum lrownum (cblkdst)\n"
                  "%ld: %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld (%ld)\n",
                  (long)me, (long)me, (long)me, (long)ctrbnbr,
                  (long)((pastix_int_t*)buffer)[FTGT_CTRBCNT],
                  (long)((pastix_int_t*)buffer)[FTGT_PROCDST],
                  (long)taskdst, (long)blokdst,
                  (long)((pastix_int_t*)buffer)[FTGT_PRIONUM],
                  (long)fcolnum, (long)lcolnum,
                  (long)frownum, (long)lrownum, (long)cblkdst);

      /* pas de test possible sur prionum (=packnbr) */
#  ifdef EXACT_THREAD
      if (THREAD_COMM_OFF)
        ASSERTDBG(me==((pastix_int_t*)buffer)[FTGT_PROCDST]%SOLV_THRDNBR,MOD_SOPALIN);
#  endif
      ASSERTDBG((SYMB_FCOLNUM(cblkdst)<=fcolnum) &&
                (SYMB_LCOLNUM(cblkdst)>=lcolnum),MOD_SOPALIN);

      str = dimi;

#  ifdef NAPA_SOPALIN
      ofrownum = frownum;
      olrownum = lrownum;
      blokdst--;
      flag = 1;
      do {

        pastix_int_t trace = 0;

        /* il peut y avoir plusieurs cibles partielles */
        if (!flag)
          print_debug(DBG_SOPALIN_NAPA, "ILU: plusieurs cibles distantes\n");

        frownum = ofrownum;
        lrownum = olrownum;
        blokdst++;

        if ((!flag)
            || (SYMB_FROWNUM(blokdst) > frownum)
            || (SYMB_LROWNUM(blokdst) < lrownum))
          {
            trace = 1;
            if (flag)
              print_debug(DBG_SOPALIN_NAPA,
                          "\nILU: debug fanin"
                          " SF=%ld F=%ld SL=%ld L=%ld gb=%ld ga=%ld"
                          " (stride=%ld SFC=%ld FC=%ld)\n",
                          (long)SYMB_FROWNUM(blokdst), (long)frownum,
                          (long)SYMB_LROWNUM(blokdst), (long)lrownum, (long)0,
                          (long)(SOLV_COEFIND(blokdst)
                                 + (fcolnum-SYMB_FCOLNUM(cblkdst))*stride
                                 + frownum-SYMB_FROWNUM(blokdst)),
                          (long)stride, (long)SYMB_FCOLNUM(cblkdst),
                          (long)fcolnum);
          }

        if (SYMB_FROWNUM(blokdst)>frownum)
          {
            frownum = SYMB_FROWNUM(blokdst);
            print_debug(DBG_SOPALIN_NAPA, "ILU: tronque frownum\n");
          }
        if (SYMB_LROWNUM(blokdst)<lrownum)
          {
            lrownum=SYMB_LROWNUM(blokdst);
            print_debug(DBG_SOPALIN_NAPA, "ILU: tronque lrownum\n");
          }
        dimi = lrownum - frownum + 1;

        ind = SOLV_COEFIND(blokdst)
          + (fcolnum - SYMB_FCOLNUM(cblkdst)) * stride
          +  frownum - SYMB_FROWNUM(blokdst);

        ga  = &(SOLV_COEFTAB(cblkdst)[ind]);
#    ifdef SOPALIN_LU
        ga2 = &(SOLV_UCOEFTAB(cblkdst)[ind]);
#    endif
        gb  = (pastix_complex64_t *) ( ((char *) buffer)
                          + MAXINFO*sizeof(pastix_int_t)
                          +(frownum-ofrownum)*sizeof(pastix_complex64_t));

        if (trace)
          {
            print_debug(DBG_SOPALIN_NAPA,
                        "ILU: debug fanin"
                        " SF=%ld F=%ld SL=%ld L=%ld gb=%ld ga=%ld"
                        " (stride=%ld SFC=%ld FC=%ld)\n",
                        (long)SYMB_FROWNUM(blokdst), (long)frownum,
                        (long)SYMB_LROWNUM(blokdst), (long)lrownum,
                        (long)(frownum - ofrownum),
                        (long)(SOLV_COEFIND(blokdst)
                               +(fcolnum-SYMB_FCOLNUM(cblkdst))*stride
                               +frownum-SYMB_FROWNUM(blokdst)),
                        (long)stride, (long)SYMB_FCOLNUM(cblkdst),
                        (long)fcolnum);
          }
#  endif /* NAPA_SOPALIN */

        ASSERTDBG((SYMB_FROWNUM(blokdst) <= frownum) &&
                  (SYMB_LROWNUM(blokdst) >= lrownum),MOD_SOPALIN);

        MUTEX_LOCK(&(sopalin_data->mutex_blok[blokdst]));
        SOPALIN_GEAM("N","N",dimi,dimj,fun,gb,str,ga,stride);
#  ifdef SOPALIN_LU
        gb2 = gb + str*dimj;
        SOPALIN_GEAM("N","N",dimi,dimj,fun,gb2,str,ga2,stride);
#  endif

        z_ooc_save_coef(sopalin_data, -1, cblkdst, me);
        MUTEX_UNLOCK(&(sopalin_data->mutex_blok[blokdst]));

        /* updatecontrib cnt */
        ASSERTDBG(ctrbnbr > 0, MOD_SOPALIN);

#  ifdef NAPA_SOPALIN
        if (flag)
          {
#  endif

            MUTEX_LOCK(&(sopalin_data->mutex_task[taskdst]));
            TASK_CTRBCNT(taskdst)-=ctrbnbr;
            TASK_FTGTCNT(taskdst)--; /*-=ctrbnbr*/

            /* Unlock taskdst if counter is null */
#  ifdef PASTIX_DYNSCHED
            if ( (!TASK_CTRBCNT(taskdst)) &&
                   (sopalin_data->taskmark[taskdst] == -1))
              {
                pastix_int_t i;

                ASSERTDBG(TASK_TASKID(taskdst) != E2, MOD_SOPALIN);

#    if (DBG_PASTIX_DYNSCHED > 0)
                ASSERTDBG(sopalin_data->taskmark[taskdst] == -1, MOD_SOPALIN);
#    endif
                sopalin_data->taskmark[taskdst]++;
                MUTEX_UNLOCK(&(sopalin_data->mutex_task[taskdst]));

                i = TASK_THREADID(taskdst);
                MUTEX_LOCK(&(sopalin_data->tasktab_mutex[i]));
                queueAdd(&(sopalin_data->taskqueue[i]),
                         taskdst, (double)TASK_PRIONUM(taskdst));
                MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[i]));
                pthread_cond_broadcast(&(sopalin_data->tasktab_cond[i]));
              }
            else
              MUTEX_UNLOCK(&(sopalin_data->mutex_task[taskdst]));
#  else
            MUTEX_UNLOCK(&(sopalin_data->mutex_task[taskdst]));
            pthread_cond_broadcast(&(sopalin_data->cond_task[taskdst]));
#  endif

            if (THREAD_FUNNELED_ON)
              {
                /* MUTEX_LOCK(&(sopalin_data->mutex_comm));   */
                SOLV_FTGTCNT--;
                /* MUTEX_UNLOCK(&(sopalin_data->mutex_comm)); */
              }

#ifdef NAPA_SOPALIN
          }
        flag = 0;
      } while ((blokdst+1<SYMB_BLOKNUM(cblkdst+1)) &&
               (
                 ((SYMB_FROWNUM(blokdst+1)<=ofrownum) &&
                  (SYMB_LROWNUM(blokdst+1)>=ofrownum)) ||
                 ((SYMB_LROWNUM(blokdst+1)>=olrownum) &&
                  (SYMB_FROWNUM(blokdst+1)<=olrownum)) ||
                 ((SYMB_FROWNUM(blokdst+1)>=ofrownum) &&
                  (SYMB_LROWNUM(blokdst+1)<=olrownum)) ||
                 ((SYMB_FROWNUM(blokdst+1)<=ofrownum) &&
                  (SYMB_LROWNUM(blokdst+1)>=olrownum))
                ));
#endif

      /* Next contribution */
#ifdef SOPALIN_LU
      buffer = ((char *) buffer)+MAXINFO*sizeof(pastix_int_t)+str*dimj*sizeof(pastix_complex64_t)*2;
#else
      buffer = ((char *) buffer)+MAXINFO*sizeof(pastix_int_t)+str*dimj*sizeof(pastix_complex64_t);
#endif

      print_debug(DBG_SOPALIN_RECV, "%ld: fin ajout fanin target\n",(long)me);
    }


  /* If we use Irecv, we launch a new communication on the elected buffer */
#ifdef TEST_IRECV
  {
    CALL_MPI MPI_Irecv(thread_data->recv_fanin_buffer[elected], size, MPI_BYTE,
                       MPI_ANY_SOURCE, me, PASTIX_COMM,
                       &(thread_data->recv_fanin_request[elected]));
    TEST_MPI("MPI_Irecv");
  }
#endif /* TEST_IRECV */

  trace_end_task(thread_data->tracefile,
                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                 STATE_L2_RECVF, taskdst);
}

/*
 * Function: z_recv_waitone_fanin
 *
 * Wait one fanin communication and call recv_handle_fanin
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *    tag          - tag used for communication
 *
 */
void
z_recv_waitone_fanin(z_Sopalin_Data_t *sopalin_data,
                   pastix_int_t             me,
                   pastix_int_t             tag)
{
  z_SolverMatrix  *datacode    = sopalin_data->datacode;
  z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  MPI_Status     status;
  int            elected     = 0;

  print_debug(DBG_SOPALIN_RECV, "%ld: z_recv_waitone_fanin\n", (long)me);

#ifdef TEST_IRECV

  /* Si on est en irecv, on attend la reception */
  CALL_MPI MPI_Waitany(MAX_R_REQUESTS, thread_data->recv_fanin_request,
                       &elected, &status);
  TEST_MPI("MPI_Waitany");
  thread_data->recv_buffer = thread_data->recv_fanin_buffer[elected];

#else /* TEST_IRECV */

  {
    pastix_int_t size;

    /* sinon on prepare le mpi_recv */
    tag = TAG_FANIN;
    if (THREAD_COMM_OFF)
      {
#  if (defined EXACT_THREAD)
        tag = me;
#  elif (defined EXACT_TAG)
        tag = tag;
#  endif
      }

    size = PACKMAX*(sizeof(pastix_int_t)*MAXINFO)+PACKAREA*sizeof(pastix_complex64_t);
    OOC_RECEIVING;
    CALL_MPI MPI_Recv(thread_data->recv_buffer, size, MPI_BYTE,
                      MPI_ANY_SOURCE, tag, PASTIX_COMM, &status);
    OOC_RECEIVED;
    TEST_MPI("MPI_Recv");
  }
#endif /* TEST_IRECV */

  recv_handle_fanin(sopalin_data, me, thread_data->recv_buffer,
                    status, elected);
}

/*
 * Function: recv_waitone_fob
 *
 * Wait one fanin or one block communication and call
 * associate recv_handle_(fanin|block)
 * Works only without EXACT_TAG or EXACT_THREAD
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *
 */
void
recv_waitone_fob(z_Sopalin_Data_t *sopalin_data,
                 pastix_int_t             me)
{
  print_debug(DBG_SOPALIN_RECV, "%ld: recv_waitone_fob\n",(long)me);

#ifdef TEST_IRECV
  {
    z_SolverMatrix  *datacode    = sopalin_data->datacode;
    z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
    MPI_Status    *statuses    = thread_data->srteststatus;
    int           *indices     = thread_data->srtestindices;
    int            i, elected, outcount1;

    if (thread_data->srteststatus == NULL)
      {
        MALLOC_INTERN(thread_data->srteststatus,  MAX(MAX_R_REQUESTS,
                                                      MAX_S_REQUESTS),
                      MPI_Status);
        MALLOC_INTERN(thread_data->srtestindices, MAX(MAX_R_REQUESTS,
                                                      MAX_S_REQUESTS), int);
        statuses = thread_data->srteststatus;
        indices  = thread_data->srtestindices;
      }

    outcount1 = 0;
    while (!outcount1)
      {
        CALL_MPI MPI_Testsome(MAX_R_REQUESTS, thread_data->recv_fanin_request,
                              &outcount1, indices, statuses);
        TEST_MPI("MPI_Testsome");

        for (i=0;i<outcount1;i++)
          {
            elected = indices[i];
            recv_handle_fanin(sopalin_data, me,
                              thread_data->recv_fanin_buffer[elected],
                              statuses[elected], elected);
          }
      }
  }
#else

#  if (defined EXACT_TAG) || (defined EXACT_THREAD)
  /* Can add some deadlock if we wait on any_tag in smp */
#    ifdef SMP_SOPALIN
  errorPrintW("Tag EXACT or THREAD are incompatible with the function"
              " recv_waitone_fob");
#    else
  {
    z_SolverMatrix  *datacode    = sopalin_data->datacode;
    z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
    MPI_Status     status;
    pastix_int_t            size;
    size = PACKMAX*(sizeof(pastix_int_t)*MAXINFO)+
        PACKAREA*sizeof(pastix_complex64_t);

    /* Test one fanin */
    CALL_MPI MPI_Recv(thread_data->recv_buffer, size, MPI_BYTE,
                      MPI_ANY_SOURCE, MPI_ANY_TAG, PASTIX_COMM, &status);
    TEST_MPI("MPI_Recv");

    assert(status.MPI_TAG < SEPFB);
    recv_handle_fanin(sopalin_data, me,
                      thread_data->recv_buffer,
                      status, 0);
  }
#    endif
#  else
  {
    z_SolverMatrix  *datacode    = sopalin_data->datacode;
    z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
    MPI_Status     status;
    int            elected     = 0;
    pastix_int_t            size;
    size = PACKMAX*(sizeof(pastix_int_t)*MAXINFO)+
        PACKAREA*sizeof(pastix_complex64_t));


    CALL_MPI MPI_Recv(thread_data->recv_buffer, size, MPI_BYTE,
                      MPI_ANY_SOURCE, MPI_ANY_TAG, PASTIX_COMM, &status);
    TEST_MPI("MPI_Recv");
    switch(status.MPI_TAG)
      {
      case TAG_FANIN:
        recv_handle_fanin(sopalin_data, me,
                          thread_data->recv_buffer,
                          status, elected);
        break;
      default:
        print_debug(DBG_SOPALIN_COMM, "tag unknown\n");
        EXIT(MOD_SOPALIN,INTERNAL_ERR);
      }
  }
#  endif
#endif
}

/*
 * Function: recv_testone_fob
 *
 * Test one fanin or one block communication and call associate
 * recv_handle_(fanin|block)
 * Works only without EXACT_TAG or EXACT_THREAD
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *
 */
void recv_testone_fob(z_Sopalin_Data_t *sopalin_data, pastix_int_t me)
{
  print_debug(DBG_SOPALIN_RECV, "%ld: recv_testone_fob\n",(long)me);

#ifdef TEST_IRECV
  {
    z_SolverMatrix  *datacode    = sopalin_data->datacode;
    z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
    MPI_Status     status;
    int            flag;
    int            elected;

    CALL_MPI MPI_Testany(MAX_R_REQUESTS, thread_data->recv_fanin_request,
                         &elected, &flag, &status);
    TEST_MPI("MPI_Testany");
    if (flag)
      {
        recv_handle_fanin(sopalin_data, me,
                          thread_data->recv_fanin_buffer[elected],
                          status, elected);
      }
  }

#else /* Test_Irecv */

#  if (defined EXACT_TAG)
  /* Can add some deadlock if we wait on any_tag in smp */
#    ifdef SMP_SOPALIN
  errorPrintW("Tag EXACT is incompatible with the function recv_testall_fab");
#    else
  {
    MPI_Status    status;
    int           flag;

    flag = 0;
    while ((!flag))
      {
        /* Test one fanin */
        CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, PASTIX_COMM,
                            &flag, &status);
        TEST_MPI("MPI_Iprobe");
        if (flag)
          {
              assert(flag < SEPFB);
              z_recv_waitone_fanin(sopalin_data, me, status.MPI_TAG);
          }
      }
  }
#    endif
#  elif (defined EXACT_THREAD)
  {
#    ifdef SMP_SOPALIN
    z_SolverMatrix *datacode    = sopalin_data->datacode;
#    endif
    MPI_Status    status;
    int           flag1;

    flag1 = 0;
    while (!flag1)
      {
        /* Test one fanin */
        CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, me, PASTIX_COMM,
                            &flag1, &status);
        TEST_MPI("MPI_Iprobe");
        if (flag1)
          z_recv_waitone_fanin(sopalin_data, me, status.MPI_TAG);
      }
  }
#  else
  {
    MPI_Status status;
    int        flag;

    CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, PASTIX_COMM,
                        &flag, &status);
    TEST_MPI("MPI_Iprobe");
    if (flag)
      {
        switch (status.MPI_TAG)
          {
          case TAG_FANIN:
            z_recv_waitone_fanin(sopalin_data, me, status.MPI_TAG);
            break;
          default:
            print_debug(DBG_SOPALIN_COMM, "tag unknown\n");
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
          }
      }
  }
#  endif
#endif /* TEST_IRECV */
}

/*
 * Function: recv_testall_fab
 *
 * Test all active receive communication and call associate
 * recv_handle_(fanin|block)
 * Works only without EXACT_TAG or EXACT_THREAD
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *
 */
void
recv_testall_fab(z_Sopalin_Data_t *sopalin_data,
                 pastix_int_t             me)
{

  print_debug(DBG_SOPALIN_RECV, "%ld: recv_testall_fab\n",(long)me);

  /* We don't care about the tag since there is launched requests */
#ifdef TEST_IRECV
  {
    z_SolverMatrix  *datacode    = sopalin_data->datacode;
    z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
    MPI_Status    *statuses;
    int           *indices;
    int            i, elected, outcount;

    if (thread_data->srteststatus == NULL)
      {
        MALLOC_INTERN(thread_data->srteststatus,  MAX(MAX_R_REQUESTS,
                                                      MAX_S_REQUESTS),
                      MPI_Status);
        MALLOC_INTERN(thread_data->srtestindices, MAX(MAX_R_REQUESTS,
                                                      MAX_S_REQUESTS), int);
      }
    statuses = thread_data->srteststatus;
    indices  = thread_data->srtestindices;

    CALL_MPI MPI_Testsome(MAX_R_REQUESTS, thread_data->recv_fanin_request,
                          &outcount, indices, statuses);
    TEST_MPI("MPI_Testsome");

    for (i=0;i<outcount;i++)
      {
        elected = indices[i];
        recv_handle_fanin(sopalin_data, me,
                          thread_data->recv_fanin_buffer[elected],
                          statuses[elected], elected);
      }
  }
#else
#  if (defined EXACT_TAG)
  /* Can add some deadlock if we wait on any_tag */
  errorPrintW("Tag EXACT is incompatible with the function recv_testall_fab");
#  elif (defined EXACT_THREAD)
  {
#    ifdef SMP_SOPALIN
    z_SolverMatrix *datacode    = sopalin_data->datacode;
#    endif
    MPI_Status    status;
    int           flag1;

    flag1 = 1;
    while (flag1)
      {
        /* Test one fanin */
        CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, me, PASTIX_COMM,
                            &flag1, &status);
        TEST_MPI("MPI_Iprobe");
        if (flag1)
          z_recv_waitone_fanin(sopalin_data, me, status.MPI_TAG);
      }
  }
#  else
  /* Normally don't use in this section exept if I modifiy some code
   in thread_comm */
  {
    MPI_Status status;
    int        flag;

    flag = 1;
    while (flag)
      {
        CALL_MPI MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, PASTIX_COMM,
                            &flag, &status);
        TEST_MPI("MPI_Iprobe");
        if (flag)
          {
            switch (status.MPI_TAG)
              {
              case TAG_FANIN:
                z_recv_waitone_fanin(sopalin_data, me, status.MPI_TAG);
                break;
              default:
                print_debug(DBG_SOPALIN_COMM, "tag unknown\n");
                EXIT(MOD_SOPALIN,INTERNAL_ERR);
              }
          }
      }
  }
#  endif /* tag */
#endif /* TEST_IRECV */
}


/****************************************************************************/
/* SEND ROUTINES MPI                                                        */
/****************************************************************************/

/*
 * Function: send_one_fanin
 *
 * Send all contribution for a same task on a same destination.
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *    t            - First fanin number
 *
 * Returns:
 *    Number of fanin sent
 */
int
send_one_fanin ( z_Sopalin_Data_t *sopalin_data,
                 pastix_int_t             me,
                 pastix_int_t             t)
{
  z_SolverMatrix  *datacode    = sopalin_data->datacode;
  z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  Queue         *sendqueue;
  pastix_int_t            id_req  = 0;
  pastix_int_t            packnbr = 0;
  pastix_int_t            tdeb, tag;
  pastix_int_t           *extra; /* contient la liste supplementaire
                         * des fanintgt a liberer */
  pastix_int_t            sum_pack;
#ifdef NO_MPI_TYPE
  pastix_int_t            iter;
  pastix_int_t            copied;
#else  /* NO_MPI_TYPE */
  MPI_Datatype   newtype;
#endif /* NO_MPI_TYPE  */

  print_debug(DBG_SOPALIN_BLEND,
              "%ld: Send fanintarget\n %ld:"
              " ctrnbr ctrbcnt procdst taskdst blokdst"
              " prionum fcolnum lcolnum frownum lrownum ->"
              " tag\n %ld: %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n",
              (long)me,(long)me,(long)me,(long)FANIN_CTRBNBR(t),
              (long)FANIN_CTRBCNT(t),(long)FANIN_PROCDST(t),
              (long)FANIN_TASKDST(t),(long)FANIN_BLOKDST(t),
              (long)FANIN_PRIONUM(t),(long)FANIN_FCOLNUM(t),
              (long)FANIN_LCOLNUM(t),(long)FANIN_FROWNUM(t),
              (long)FANIN_LROWNUM(t));

  tdeb  = t;
  extra = NULL;

  /*************************/
  /* add first contibution */
  /*************************/

  /* Number of elements for each type */
  thread_data->gtabsize[packnbr]   = MAXINFO;
  thread_data->gtabsize[packnbr+1] = (FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1) *
    (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1);
#ifdef SOPALIN_LU
  thread_data->gtabsize[packnbr+1]*= 2;
#endif

  sum_pack = thread_data->gtabsize[packnbr+1];

  /* Type of each vector */
#ifdef NO_MPI_TYPE
  thread_data->gtabtype[packnbr]   = sizeof(pastix_int_t);
  thread_data->gtabtype[packnbr+1] = sizeof(pastix_complex64_t);
#else /* NO_MPI_TYPE */
  thread_data->gtabtype[packnbr]   = PASTIX_MPI_INT;
  thread_data->gtabtype[packnbr+1] = COMM_FLOAT;
#endif /* NO_MPI_TYPE */

#ifdef OOC_FTGT
  print_debug(DBG_OOC_FTGT, "WAIT %4d %4d\n", (int)t, (int) -1);
  MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
  z_ooc_wait_for_ftgt(sopalin_data, t, me);
  ASSERTDBG(((unsigned long)(*(((double*)FANIN_COEFTAB(t))-1))) ==
            sizeof(pastix_complex64_t)*((FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1) *
                           (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1) *
                           ((sopalin_data->sopar->factotype ==
                             API_FACT_LU)?2:1))
            , MOD_SOPALIN);
  MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
#endif

  print_debug(DBG_SOPALIN_SEND, "%ld: FANIN_INFOTAB : %x\n", (long)me,
              (unsigned int)(intptr_t)FANIN_INFOTAB(t));
  print_debug(DBG_SOPALIN_SEND, "%ld: FANIN_COEFTAB : %x\n", (long)me,
              (unsigned int)(intptr_t)FANIN_COEFTAB(t));
  print_debug(DBG_SOPALIN_SEND, "%ld: FANIN_CTRBCNT : %x\n", (long)me,
              (unsigned int)(intptr_t)FANIN_CTRBCNT(t));

  /* Adress of each vector */
#ifdef NO_MPI_TYPE
  thread_data->gtaboffs[packnbr]   = FANIN_INFOTAB(t);
  thread_data->gtaboffs[packnbr+1] = FANIN_COEFTAB(t);
#else /* NO_MPI_TYPE */
  CALL_MPI MPI_Get_address(FANIN_INFOTAB(t),&(thread_data->gtaboffs[packnbr]));
  TEST_MPI("MPI_Get_address");
  CALL_MPI MPI_Get_address(FANIN_COEFTAB(t),&(thread_data->gtaboffs[packnbr+1]));
  TEST_MPI("MPI_Get_address");
#endif /* NO_MPI_TYPE */

  /* Add other contribution for the same task */
  if (THREAD_FUNNELED_ON)
    {
      sendqueue = sopalin_data->sendqueue;
    }
  else
    {
      sendqueue = &(sopalin_data->fanintgtsendqueue[FANIN_PROCDST(tdeb)]);
    }

  if (queueSize(sendqueue))
    {
      t = queueRead(sendqueue);

      print_debug(DBG_SOPALIN_COMM,
                  "%ld-%ld C: ftgt %ld / dest %ld / key %ld / task %ld\n",
                  (long)SOLV_PROCNUM, (long)me, (long)t, (long)FANIN_PROCDST(t),
                  (long)FANIN_PRIONUM(t), (long)FANIN_TASKDST(t));

      /* allocation du tableau listant les contributions qui seront envoyees */
      /* Attention : FANIN_CTRBCNT est testé en 2 fois et doit donc
       *             etre protégé */
      MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
      if ((FANIN_PROCDST(tdeb) == FANIN_PROCDST(t)) &&
          (FANIN_TASKDST(tdeb) == FANIN_TASKDST(t)) &&
          (!(FANIN_CTRBCNT(t))) && ((packnbr/2)<(PACKMAX-1)))
        MALLOC_INTERN(extra, PACKMAX, pastix_int_t);

      while ((FANIN_PROCDST(tdeb) == FANIN_PROCDST(t)) &&
             (FANIN_TASKDST(tdeb) == FANIN_TASKDST(t)) &&
             (!(FANIN_CTRBCNT(t))) && ((packnbr/2)<(PACKMAX-1)))
        {
          pastix_int_t next;

          MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
          next = (FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1)*
            (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1);

#ifdef SOPALIN_LU
          next*=2;
#endif

          /* Si cela fait trop de donnees a envoyer, on sort de la boucle */
          if (sum_pack+next>PACKAREA) break;

          t = queueGet(sendqueue);

          print_debug(DBG_SOPALIN_SEND, "%ld: Extra fanintarget\n"
                      "%ld: ctrnbr ctrbcnt procdst taskdst blokdst prionum"
                      " fcolnum lcolnum frownum lrownum\n"
                      "%ld: %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n",
                      (long)me,(long)me,(long)me,(long)FANIN_CTRBNBR(t),
                      (long)FANIN_CTRBCNT(t),(long)FANIN_PROCDST(t),
                      (long)FANIN_TASKDST(t),
                      (long)FANIN_BLOKDST(t),(long)FANIN_PRIONUM(t),
                      (long)FANIN_FCOLNUM(t),
                      (long)FANIN_LCOLNUM(t),(long)FANIN_FROWNUM(t),
                      (long)FANIN_LROWNUM(t));

          extra[packnbr/2] = t;
          packnbr         += 2;

          /* create MPI type */
          thread_data->gtabsize[packnbr]   = MAXINFO;
          thread_data->gtabsize[packnbr+1] =
            (FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1)*
            (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1);
#ifdef SOPALIN_LU
          thread_data->gtabsize[packnbr+1]*= 2;
#endif

          sum_pack += thread_data->gtabsize[packnbr+1];

#ifdef NO_MPI_TYPE
          thread_data->gtabtype[packnbr]   = sizeof(pastix_int_t);
          thread_data->gtabtype[packnbr+1] = sizeof(pastix_complex64_t);
#else /* NO_MPI_TYPE */
          thread_data->gtabtype[packnbr]   = PASTIX_MPI_INT;
          thread_data->gtabtype[packnbr+1] = COMM_FLOAT;
#endif /* NO_MPI_TYPE */

#ifdef OOC_FTGT
          print_debug(DBG_OOC_FTGT, "WAIT %4d %4d\n", (int)t, (int) -1);
          MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
          z_ooc_wait_for_ftgt(sopalin_data, t, me);
          ASSERTDBG(((unsigned long)(*(((double*)FANIN_COEFTAB(t))-1))) ==
                    sizeof(pastix_complex64_t)*((FANIN_LROWNUM(t)-FANIN_FROWNUM(t)+1) *
                                   (FANIN_LCOLNUM(t)-FANIN_FCOLNUM(t)+1)   *
                                   ((sopalin_data->sopar->factotype == API_FACT_LU)?2:1))
                    , MOD_SOPALIN);
          MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
#endif

          print_debug(DBG_SOPALIN_SEND, "%ld: FANIN_INFOTAB : %x\n",
                      (long)me, (unsigned int)(intptr_t)FANIN_INFOTAB(t));
          print_debug(DBG_SOPALIN_SEND, "%ld: FANIN_COEFTAB : %x\n",
                      (long)me, (unsigned int)(intptr_t)FANIN_COEFTAB(t));
          print_debug(DBG_SOPALIN_SEND, "%ld: FANIN_CTRBCNT : %x\n",
                      (long)me, (unsigned int)(intptr_t)FANIN_CTRBCNT(t));

#ifdef NO_MPI_TYPE
          thread_data->gtaboffs[packnbr]   = FANIN_INFOTAB(t);
          thread_data->gtaboffs[packnbr+1] = FANIN_COEFTAB(t);
#else /* NO_MPI_TYPE */
          CALL_MPI MPI_Get_address(FANIN_INFOTAB(t),
                               &(thread_data->gtaboffs[packnbr]));
          TEST_MPI("MPI_Get_address");
          CALL_MPI MPI_Get_address(FANIN_COEFTAB(t),
                               &(thread_data->gtaboffs[packnbr+1]));
          TEST_MPI("MPI_Get_address");
#endif /* NO_MPI_TYPE */

          if (queueSize(sendqueue))
            {
              t = queueRead(sendqueue);
              MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
            }
          else
            {
              MUTEX_LOCK(&(sopalin_data->mutex_fanin[t]));
              break;
            }
        }

      MUTEX_UNLOCK(&(sopalin_data->mutex_fanin[t]));
      if (extra)
        extra[packnbr/2] = -1;
    }

  if (packnbr/2 > PACKMAX)
    print_debug(DBG_SOPALIN_SEND, "packnbr/2=%ld<?(%ld)\n",
                (long)(packnbr/2),(long)PACKMAX);
  if (sum_pack > PACKAREA)
    print_debug(DBG_SOPALIN_SEND, "sum_pack=%ld<?(%ld)\n",
                (long)(sum_pack),(long)PACKAREA);
  ASSERTDBG(packnbr/2 <= PACKMAX, MOD_SOPALIN);
  ASSERTDBG(sum_pack  <= PACKAREA,MOD_SOPALIN);

  /*   fprintf(thread_data->tracefile,"\n"); */
  print_debug(DBG_SOPALIN_SEND, "%ld: packnbr/2=%ld sum_pack=%ld\n",
              (long)me,(long)(packnbr/2),(long)sum_pack);

  /* Choix du tag suivant la version */
  if (sopalin_data->sopar->type_comm != 3)
    tag = TAG_FANIN;
  else
    tag = FANIN_INFOTAB(tdeb)[FTGT_PROCDST]%SOLV_THRDNBR;
  if (THREAD_COMM_OFF)
    {
#if (defined EXACT_THREAD)
      tag = FANIN_INFOTAB(tdeb)[FTGT_PROCDST]%SOLV_THRDNBR;
#elif (defined EXACT_TAG)
      tag = FANIN_PRIONUM(tdeb);
#endif
    }

  /* le nombre de paquets est code dans le champ prionum de la 1ere fanin */
  FANIN_PRIONUM(tdeb) = packnbr/2+1;

  /* On recherche la premiere requete disponible */
#ifdef TEST_ISEND
  id_req = send_waitone_fanin(sopalin_data, me);
  ASSERTDBG(id_req<MAX_S_REQUESTS,MOD_SOPALIN);
#endif /* TEST_ISEND */

  /* Envoi des donnees */
  trace_send(thread_data->tracefile,
             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, FANIN_PROCDST(tdeb),
             COMM_FANIN, FANIN_TASKDST(tdeb), thread_data->gtabsize[1],
             &(FANIN_IDTRACE(tdeb)));

#ifdef NO_MPI_TYPE
  thread_data->send_fanin_buffer_size[id_req] = 0;
  for (iter = 0; iter < 2*(packnbr/2+1); iter++) {
    thread_data->send_fanin_buffer_size[id_req] +=
      thread_data->gtabsize[iter]*thread_data->gtabtype[iter];
  }
  MALLOC_INTERN(thread_data->send_fanin_buffer[id_req],
                thread_data->send_fanin_buffer_size[id_req],
                char);
  copied = 0;
  for (iter = 0; iter < 2*(packnbr/2+1); iter++) {
    memcpy(thread_data->send_fanin_buffer[id_req]+copied,
           thread_data->gtaboffs[iter],
           thread_data->gtabtype[iter]*thread_data->gtabsize[iter] );
    copied += thread_data->gtabsize[iter]*thread_data->gtabtype[iter];
  }

#  ifdef TEST_ISEND
  CALL_MPI MPI_Isend(thread_data->send_fanin_buffer[id_req],
                     thread_data->send_fanin_buffer_size[id_req],
                     MPI_CHAR,FANIN_PROCDST(tdeb),tag,
                     PASTIX_COMM,&((thread_data->send_fanin_requests[id_req])));
  TEST_MPI("MPI_Isend");
#  else
  CALL_MPI MPI_Rsend(thread_data->send_fanin_buffer[id_req],
                     thread_data->send_fanin_buffer_size[id_req],
                     MPI_CHAR,FANIN_PROCDST(tdeb),tag,
                     PASTIX_COMM);
  TEST_MPI("MPI_Rsend");
#  endif
#else /* NO_MPI_TYPE */
  CALL_MPI MPI_Type_create_struct(2*(packnbr/2+1), thread_data->gtabsize,
                           thread_data->gtaboffs,
                           thread_data->gtabtype, &newtype);
  TEST_MPI("MPI_Type_create_struct");
  CALL_MPI MPI_Type_commit(&newtype);
  TEST_MPI("MPI_Type_commit");
#  ifdef TEST_ISEND
  thread_data->send_fanin_mpitypes[id_req] = newtype;
  CALL_MPI MPI_Isend(MPI_BOTTOM,1,thread_data->send_fanin_mpitypes[id_req],
                     FANIN_PROCDST(tdeb),tag,
                     PASTIX_COMM,&((thread_data->send_fanin_requests[id_req])));
  TEST_MPI("MPI_Isend");
#  else
  CALL_MPI MPI_Rsend(MPI_BOTTOM,1,newtype,FANIN_PROCDST(tdeb),tag,
                     PASTIX_COMM);
  TEST_MPI("MPI_Rsend");
  CALL_MPI MPI_Type_free(&newtype);
  TEST_MPI("MPI_Type_free");
#  endif
#endif /* NO_MPI_TYPE */

  thread_data->send_fanin_target[id_req]       = tdeb;
  thread_data->send_fanin_target_extra[id_req] = extra;

  /* Try to free some buffer */
#ifdef TEST_ISEND
  send_testall_fanin(sopalin_data, me);
#else
  send_free_fanin(sopalin_data, me, id_req);
#endif

  return (packnbr/2+1);
}

/*
 * Function: send_all_fanin
 *
 * Send all contribution for a different task on a same destination.
 * Can't be called with API_THREAD_FUNNELED
 *
 * Parameters:
 *    sopalin_data - Solver structure
 *    me           - Thread number
 *    dest         - First fanin number
 *
 */
void
send_all_fanin(z_Sopalin_Data_t *sopalin_data,
               pastix_int_t             me,
               pastix_int_t             dest)
{
  if (THREAD_FUNNELED_ON)
    {
      errorPrintW("API_THREAD_FUNNELED is incompatible with the function"
                  " send_all_fanin");
    }
  else
    {
      z_SolverMatrix  *datacode    = sopalin_data->datacode;
#ifdef TRACE_SOPALIN
      z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
#endif
      int            flag;

      print_debug(DBG_SOPALIN_SEND, "%ld : --> send_all_fanin\n", (long)me);

      trace_begin_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                       STATE_L2_SENDF, 0);

      flag = 1;
      if (dest != SOLV_PROCNUM)
        {
          MUTEX_LOCK(&(sopalin_data->mutex_queue_fanin[dest]));
          while ( (flag) &&
                  (queueSize(&(sopalin_data->fanintgtsendqueue[dest]))))
            {
              pastix_int_t t = queueRead(&(sopalin_data->fanintgtsendqueue[dest]));

              print_debug(DBG_SOPALIN_SEND,
                          "send dest %ld fanintarget %ld\n",(long)dest,(long)t);

              /* If the target is ready, we send it */
              /* With -DCOMM_REORDER : Always true  */
              if (!(FANIN_CTRBCNT(t)))
                {
                  double key;

                  t = queueGet2(&(sopalin_data->fanintgtsendqueue[dest]),
                                &key, NULL);

                  ASSERTDBG(FANIN_PROCDST(t) != SOLV_PROCNUM, MOD_SOPALIN);

                  /* send fanin target t */
                  send_one_fanin(sopalin_data, me, t);

                }
              else
                flag = 0;
            }
          MUTEX_UNLOCK(&(sopalin_data->mutex_queue_fanin[dest]));
        }

      trace_end_task(thread_data->tracefile,
                     SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 2,
                     STATE_L2_SENDF, 0);

      print_debug(DBG_SOPALIN_SEND, "%ld : <-- send_all_fanin\n", (long)me);
    }
  return;
}

/*
 * Function: send_free_fanin
 *
 * Free associated structure to fanin sent.
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
void
send_free_fanin ( z_Sopalin_Data_t *sopalin_data,
                  pastix_int_t             me,
                  pastix_int_t             s_index)
{
#if ((!defined OOC_FTGT) || defined PASTIX_DEBUG )
  z_SolverMatrix  *datacode    = sopalin_data->datacode;
#endif
  z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  pastix_int_t f = 0;
  pastix_int_t i;

  i = thread_data->send_fanin_target[s_index];
  print_debug(DBG_SOPALIN_SEND, "->send_free_fanin %ld \n", (long)i);
#ifdef OOC_FTGT
  z_ooc_reset_ftgt(sopalin_data, i, me);
#else
  memFree_null(FANIN_COEFTAB(i));
#endif
#ifndef NO_MPI_TYPE
  /* free MPI type */
  CALL_MPI MPI_Type_free(&(thread_data->send_fanin_mpitypes[s_index]));
  TEST_MPI("MPI_Type_free");
#endif /* NO_MPI_TYPE */

  print_debug(DBG_SOPALIN_ALLOC, "free fanin coeff %x\n",
              (unsigned int)(intptr_t)FANIN_COEFTAB(i));

  STATS_SUB((FANIN_LROWNUM(i)-FANIN_FROWNUM(i)+1)*
            (FANIN_LCOLNUM(i)-FANIN_FCOLNUM(i)+1));

  if (thread_data->send_fanin_target_extra[s_index])
    {
      while((i = thread_data->send_fanin_target_extra[s_index][f]) != -1)
        {
#ifdef OOC_FTGT
          z_ooc_reset_ftgt(sopalin_data,i,me);
#else
          memFree_null(FANIN_COEFTAB(i));
#endif

          print_debug(DBG_SOPALIN_ALLOC, "free fanin coeff %x\n",
                      (unsigned int)(intptr_t)FANIN_COEFTAB(i));

          STATS_SUB((FANIN_LROWNUM(i)-FANIN_FROWNUM(i)+1)*
                    (FANIN_LCOLNUM(i)-FANIN_FCOLNUM(i)+1));

          f++;
        }

      memFree_null(thread_data->send_fanin_target_extra[s_index]);
      thread_data->send_fanin_target_extra[s_index] = NULL;
    }

#ifdef NO_MPI_TYPE
  memFree_null(thread_data->send_fanin_buffer[s_index]);
  thread_data->send_fanin_buffer[s_index] = NULL;
#endif /* NO_MPI_TYPE */
  print_debug(DBG_SOPALIN_SEND, "<-send_free_fanin %ld \n", (long)i);
}

/*********************************/
/*
 * Function: send_testall_fanin
 *
 * Test all fanin sent to progress communications
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
/*********************************/
void send_testall ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me,
                    void (*funcfree)(z_Sopalin_Data_t*, pastix_int_t, pastix_int_t))
{
  z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  int            i;

  print_debug(DBG_SOPALIN_UPDO, "%ld: test_all_downsend\n", (long)me);

#ifndef PASTIX_TERA
  {
#  if (defined TEST_IRECV) || ((defined TEST_ISEND) && (defined SMP_SOPALIN))
    z_SolverMatrix  *datacode = sopalin_data->datacode;
#  endif
    MPI_Status    *statuses;
    int           *indices;
    int            outcount = 0;

    if (thread_data->srteststatus == NULL)
      {
        MALLOC_INTERN(thread_data->srteststatus,
                      MAX(MAX_R_REQUESTS, MAX_S_REQUESTS), MPI_Status);
        MALLOC_INTERN(thread_data->srtestindices,
                      MAX(MAX_R_REQUESTS, MAX_S_REQUESTS), int);
      }
    statuses = thread_data->srteststatus;
    indices  = thread_data->srtestindices;

    if (thread_data->maxsrequest_fanin > 0)
      {
        CALL_MPI MPI_Testsome(thread_data->maxsrequest_fanin,
                              thread_data->send_fanin_requests,
                              &outcount, indices, statuses);
        TEST_MPI("MPI_Testsome");

        for(i=0; i<outcount; i++)
          {
            funcfree(sopalin_data, me, indices[i]);
          }
      }
  }
#else /* PASTIX_TERA */
  /* Can be removed in next release, if the first version is ok on TERA10 */
  {
    MPI_Status s_status;
    int        s_flag;

    for(i=0; i<thread_data->maxsrequest_fanin; i++)
      {
        if (!MPI_Request_is_equal(thread_data->send_fanin_requests[i],
                                  MPI_REQUEST_NULL))
          {
            CALL_MPI MPI_Test(&(thread_data->send_fanin_requests[i]),
                              &s_flag, &s_status);
            TEST_MPI("MPI_Test");
            if (s_flag)
              {
                funcfree(sopalin_data, me, i);
              }
          }
      }
  }
#endif /* PASTIX_TERA */
}

void send_testall_fanin(z_Sopalin_Data_t *sopalin_data, pastix_int_t me)
{
  send_testall(sopalin_data, me, send_free_fanin);
  return;
}

/*********************************/
/*
 * Function: send_testall_fab
 *
 * Test all block sent to progress communications
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
/*********************************/
void send_testall_fab(z_Sopalin_Data_t *sopalin_data, pastix_int_t me)
{
  send_testall_fanin(sopalin_data, me);
}

/*********************************/
/*
 * Function: send_waitone_fanin
 *
 * Test fanin sent to return an id or wait until one fanin finished.
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
/*********************************/
int send_waitone ( z_Sopalin_Data_t *sopalin_data, pastix_int_t me,
                   void (*funcfree)(z_Sopalin_Data_t*, pastix_int_t, pastix_int_t) )
{
#if (defined TEST_ISEND) && (defined SMP_SOPALIN)
  z_SolverMatrix  *datacode    = sopalin_data->datacode;
#endif
  z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  MPI_Status     status;
  int            flag = 0;
  int            req  = 0;

  print_debug(DBG_SOPALIN_RECV, "%ld: send_waitone_fanin\n",(long)me);

  while((!flag) && (thread_data->maxsrequest_fanin > 0))
    {
      CALL_MPI MPI_Testany(thread_data->maxsrequest_fanin,
                           thread_data->send_fanin_requests,
                           &req, &flag, &status);
      TEST_MPI("MPI_Testany");

      if ( flag )
        {
          if (req != MPI_UNDEFINED)
            {
              funcfree(sopalin_data, me, req);
              return req;
            }
          else
            /* Case where all requests are finished */
            return 0;
        }
      else
        {
          if (thread_data->maxsrequest_fanin < MAX_S_REQUESTS)
            {
              req = thread_data->maxsrequest_fanin;
              thread_data->maxsrequest_fanin++;
              return req;
            }
        }
    }

  /* On n'est pas rentré dans la boucle */
  thread_data->maxsrequest_fanin++;
  return req;
}

int send_waitone_fanin(z_Sopalin_Data_t *sopalin_data, pastix_int_t me)
{
  return send_waitone(sopalin_data, me, send_free_fanin);
}

/*********************************/
/*
 * Function: z_send_waitall_fab
 *
 * Wait for all pending communications (fanin and block).
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
/*********************************/
void z_send_waitall_fab(z_Sopalin_Data_t *sopalin_data, pastix_int_t me)
{
  z_Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  MPI_Status     s_status;
  int            i;

#ifdef FORCE_CONSO
  pastix_int_t        nb_envois_fanin = 0;
  int        s_flag          = 0;

  while(nb_envois_fanin != thread_data->maxsrequest_fanin)
    {
      nb_envois_fanin = 0;
      for(i=0; i<thread_data->maxsrequest_fanin; i++)
        {
          if (!MPI_Request_is_equal(thread_data->send_fanin_requests[i], MPI_REQUEST_NULL))
            {
              CALL_MPI MPI_Test(&(thread_data->send_fanin_requests[i]),
                                &s_flag, &s_status);
              TEST_MPI("MPI_Test");
              if (s_flag)
                {
                  send_free_fanin(sopalin_data, me, i);
                  nb_envois_fanin++;
                }
            }
          else
            nb_envois_fanin++;
        }
    }
#endif

  for (i=0;i<thread_data->maxsrequest_fanin;i++)
    if (!MPI_Request_is_equal(thread_data->send_fanin_requests[i],
                              MPI_REQUEST_NULL))
      {
        CALL_MPI MPI_Wait(&thread_data->send_fanin_requests[i], &s_status);
        TEST_MPI("MPI_Wait");

        send_free_fanin(sopalin_data, me, i);
      }
}

/*********************************/
/*
 * Function: z_rcsd_testall_fab
 *
 * Launch threads for solving step.
 *
 * Parameters:
 *
 * Returns:
 *   void
 */
/*********************************/
void z_rcsd_testall_fab(z_Sopalin_Data_t *sopalin_data, pastix_int_t me)
{
#ifdef TEST_ISEND
  /* Test des envois */
  send_testall_fab(sopalin_data, me);
#endif

  /* Test des receptions */
  recv_testall_fab(sopalin_data, me);
}


/*
 * Function: sendrecv_smp
 *
 * fonction de réception des comms dans la version threads séparés
 *
 */
void* sendrecv_smp ( void *arg )
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  z_Sopalin_Data_t   *sopalin_data = (z_Sopalin_Data_t *)(argument->data);
  if (THREAD_COMM_ON)
    {
      z_SolverMatrix     *datacode     = sopalin_data->datacode;
      MPI_Comm          pastix_comm  = PASTIX_COMM;
      int               type_thcomm  = sopalin_data->sopar->type_comm;
      int               me           = argument->me;
#ifdef TRACE_SOPALIN
      z_Thread_Data_t    *thread_data;
#endif
      void            **receive_buffer;
      MPI_Request      *request;
      MPI_Status        status;
      int               type_comm;
      int               tag_fanin;
      int               nbrequest   = 1;
      int               nbrequesttot= 1;
      pastix_int_t      size;
      int               i;
      int               init;
      int               flag, wait;
      int               nbsend, nbsend_fanin;
      int               nbrecv;
      pastix_int_t      save_ftgtcnt = -1;
      pastix_int_t      ftgt, key;
      double            dest;
      int               nb_proc_end = 1;


      print_debug(DBG_SOPALIN_THREADCOMM, "%d - %d : --> SendRecv\n",
                  (int)SOLV_PROCNUM, (int)me);

      if (SOLV_PROCNBR == 1) goto end;

      /* Allocation de la structure de données spécifique a ce thread */
      init = 0;
      if (THREAD_FUNNELED_ON)
        {
          init = init | INIT_SEND;
        }

      z_sopalin_init_smp(sopalin_data, me, 1, init);
#ifdef TRACE_SOPALIN
      thread_data = sopalin_data->thread_data[me];
#endif
      /***********************************/
      /*         Reception               */
      /***********************************/

      /* Taille des buffers de réception */
      size = PACKMAX*(sizeof(pastix_int_t)*MAXINFO)+PACKAREA*sizeof(pastix_complex64_t);

#ifdef THREAD_COMM_MULTIPLE
      nbrequest = MAX(SOLV_PROCNBR-1, 1);
#endif
      nbrequesttot = nbrequest;

      /* Allocation des buffers de reception et des requetes */
      MALLOC_INTERN(receive_buffer, nbrequesttot, void*);
      MALLOC_INTERN(request       , nbrequesttot, MPI_Request);

      for(i=0; i< nbrequest; i++)
        {
          MALLOC_INTERN(receive_buffer[i], size, char);
        }

      /* Initialisation des requêtes */
      for (i=0; i<nbrequesttot; i++)
        request[i] = MPI_REQUEST_NULL;

      /* Exact thread ou Tag classique
       ( cas 1 thread de comm par thread de calcul )*/
      if (type_thcomm == 3)
        {
          tag_fanin = me-SOLV_THRDNBR;
        }
      else
        {
          tag_fanin = TAG_FANIN;
        }

      /* Initialisation des comms persistantes en réception */
      /* Pas performant, pourrait etre vire */
#ifdef THREAD_COMM_MULTIPLE
      /* Proc de rang inferieur au proc local */
      for(i=0; i < SOLV_PROCNUM; i++){
        CALL_MPI MPI_Recv_init(receive_buffer[i], size, MPI_BYTE,
                               i,tag_fanin,pastix_comm,&request[i]);
        TEST_MPI("MPI_Recv_init");
      }
      /* Proc de rang supérieur au proc local */
      for(i=SOLV_PROCNUM+1; i<SOLV_PROCNBR; i++){
        CALL_MPI MPI_Recv_init(receive_buffer[i-1],size, MPI_BYTE,
                               i,tag_fanin,pastix_comm,&request[i-1]);
        TEST_MPI("MPI_Recv_init");
      }
#else
      CALL_MPI MPI_Recv_init(receive_buffer[0], size, MPI_BYTE,
                             MPI_ANY_SOURCE, tag_fanin, pastix_comm,
                             &request[0]);
      TEST_MPI("MPI_Recv_init");
#endif

      /***********************************/
      /*       Boucle Principale         */
      /***********************************/

      print_debug(DBG_SOPALIN_THREADCOMM,
                  "%d - %d : Boucle Emission-Reception\n",
                  (int)SOLV_PROCNUM, (int)me);

      /* Lancement des communications */
      CALL_MPI MPI_Startall(nbrequesttot, request);
      TEST_MPI("MPI_Startall");

      if (THREAD_FUNNELED_OFF)
        {
          /*
           * Pas de test sur le nombre de comm restant a recevoir
           * car ca pourrait poser pb en multi-thread de comm
           */
          while(1){

            trace_begin_task(thread_data->tracefile,
                             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                             STATE_WAITREM, 1);

            /* On attend une comm */
            CALL_MPI MPI_Waitany(nbrequesttot, request, &type_comm, &status);
            TEST_MPI("MPI_Waitany");

            print_debug(DBG_SOPALIN_THREADCOMM,
                        "%d - %d : 1 reception / nb_proc_end %d\n",
                        (int)SOLV_PROCNUM, (int)me, nb_proc_end);

            /* Recuperation des comms de fin */

            trace_begin_task(thread_data->tracefile,
                             SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                             STATE_COMPUTE, 1);

            if ((status.MPI_TAG == tag_fanin)
                && ( ((int *)(receive_buffer[type_comm]))[0] == -1 ))
              {
                nb_proc_end++;
                /* fprintf(stderr, "%d - %d : reçu msg end de %d (%d, %d)\n",
                 SOLV_PROCNUM, me, status.MPI_SOURCE,
                 status.MPI_TAG, type_comm);*/
                if (nb_proc_end < SOLV_PROCNBR)
                  {
                    CALL_MPI MPI_Start(&request[type_comm]);
                    TEST_MPI("MPI_Start");
                    continue;
                  }
                else
                  break;
              }

            /* On fait le calcul associé */
            assert(type_comm%2 == 0);
            recv_handle_fanin(sopalin_data, me,
                              receive_buffer[type_comm],
                              status, 0);

            /* On relance l'attente sur une comm */
            CALL_MPI MPI_Start(&request[type_comm]);
            TEST_MPI("MPI_Start");
          }
        }
      else
        {
          /* THREAD_FUNNELED_ON */
          wait   = 0;
          nbsend = 2;
          nbsend_fanin = SOLV_FTGTNBR;
          save_ftgtcnt = SOLV_FTGTCNT;
          nbrecv = SOLV_FTGTCNT;

          if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
            {
              fprintf(stdout, OUT2_FUN_STATS,
                      (long)SOLV_PROCNUM, (long)nbrecv,
                      (long)(nbsend_fanin));
            }

          trace_begin_task(thread_data->tracefile,
                           SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                           STATE_WAITREM, 1);

          while( nbsend || nbrecv ) {

            /*
             * Attente si rien à faire
             */
            if (wait)
              {
                MUTEX_LOCK(&(sopalin_data->mutex_comm));
                COND_TIMEWAIT(&(sopalin_data->cond_comm),
                              &(sopalin_data->mutex_comm));
                MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
              }
            wait = 1;

            /*
             * Réception des données
             */
            CALL_MPI MPI_Testany(nbrequesttot, request, &type_comm,
                                 &flag, &status);
            TEST_MPI("MPI_Testany");

            if(flag && (type_comm != MPI_UNDEFINED))
              {
                wait = 0;

                trace_begin_task(thread_data->tracefile,
                                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                                 STATE_COMPUTE, 1);

                /* On fait le calcul associé */
                assert(type_comm == 0);
                recv_handle_fanin(sopalin_data, me,
                                  receive_buffer[type_comm],
                                  status, 0);

                nbrecv = SOLV_FTGTCNT;
                fprintf(stderr, "nbrecv = %d\n", nbrecv );

                /* On relance l'attente sur une comm */
                CALL_MPI MPI_Start(&request[type_comm]);
                TEST_MPI("MPI_Start");

                print_debug(DBG_SOPALIN_THREADCOMM,
                            "%d : %d Reception solv_ftgtcnt %d\n",
                            (int)SOLV_PROCNUM, (int)me, (int)SOLV_FTGTCNT);

                trace_begin_task(thread_data->tracefile,
                                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                                 STATE_WAITREM, 1);
              } /* Fin Boucle réception */

            /*
             * Progression des envois
             */
            send_testall_fanin(sopalin_data, me);

            /*
             * Envoi des données
             */
            if (nbsend > 0)
              {
                trace_begin_task(thread_data->tracefile,
                                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                                 STATE_COMPUTE, 2);

                MUTEX_LOCK(&(sopalin_data->mutex_comm));
                if(queueSize(sopalin_data->sendqueue))
                  {
                    ftgt = queueGet2(sopalin_data->sendqueue, &dest, &key);

                    print_debug(DBG_FUNNELED,
                                "%d-%d C : ftgt %ld / dest %ld / key %ld\n",
                                (int)SOLV_PROCNUM, (int)me,
                                (long)ftgt, (long)dest,
                                (long)key);
                    print_debug(DBG_FUNNELED,
                                "%d-%d C :"
                                " ftgt %ld / dest %ld / key %ld / task %ld\n",
                                (int)SOLV_PROCNUM, (int)me, (long)ftgt,
                                (long)FANIN_PROCDST(ftgt),
                                (long)FANIN_PRIONUM(ftgt),
                                (long)FANIN_TASKDST(ftgt));
                    print_debug(DBG_FUNNELED,
                                "%d-%d C : fanin %ld\n",
                                (int)SOLV_PROCNUM, (int)me, (long)nbsend_fanin);
                    assert(dest > 0);
                    nbsend_fanin -= send_one_fanin(sopalin_data,
                                                   me, ftgt);

                    fprintf(stderr, "nbsend = %d\n", nbsend_fanin );
                    wait = 0;
                  }
                MUTEX_UNLOCK(&(sopalin_data->mutex_comm));

                nbsend = nbsend_fanin;

                trace_begin_task(thread_data->tracefile,
                                 SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                                 STATE_WAITREM, 2);
              }
          }
        }

      trace_begin_task(thread_data->tracefile,
                       SOPALIN_CLOCK_TRACE, SOLV_PROCNUM, me, 1,
                       STATE_IDLE, 1);

      print_debug(DBG_SOPALIN_THREADCOMM,
                  "%d - %d : FIN Boucle Emission-Reception\n",
                  (int)SOLV_PROCNUM, (int)me);


      /***********************************/
      /*         Reception               */
      /***********************************/

      /* Annulation des comms inutiles */
      /* et Liberation de la memoire */
      for(i=0; i<nbrequest; i++)
        {
          int fflag;
          if (THREAD_FUNNELED_OFF)
            {
              if (i == type_comm)
                {
                  memFree_null(receive_buffer[i]);
                  continue;
                }
            }
          /* Annulation des comms */
          CALL_MPI MPI_Cancel(&request[i]);
          TEST_MPI("MPI_Cancel");
          /* Test pour rendre inactive la requete */
          CALL_MPI MPI_Test(&request[i], &fflag, &status);
          TEST_MPI("MPI_Test");
          /* Liberation de la requete persistante */
          CALL_MPI MPI_Request_free(&request[i]);
          TEST_MPI("MPI_Request_free");
          /* Liberation du buffer */
          memFree_null(receive_buffer[i]);
        }

      memFree_null(receive_buffer);
      memFree_null(request);

      /***********************************/
      /*         Emission                */
      /***********************************/
      if (THREAD_FUNNELED_ON)
        {
            fprintf(stderr, "nbsend = %d\n", nbsend );
            /* Attente de la fin des communications en envoi */
            if (SOLV_FTGTCNT > 0)
                z_send_waitall_fab(sopalin_data, me);

            fprintf(stderr, "nbsend = %d DONE\n", nbsend );

            CALL_MPI MPI_Barrier(PASTIX_COMM);
            TEST_MPI("MPI_Barrier");

            /* Restoration de SOLV_FTGTCNT */
            SOLV_FTGTCNT = save_ftgtcnt;
        }
      z_sopalin_clean_smp(sopalin_data, me);

      end:
      if (me == SOLV_THRDNBR)
        {
          MUTEX_LOCK(&(sopalin_data->mutex_comm));
          sopalin_data->step_comm = COMMSTEP_FACTOEND;
          print_debug(DBG_THCOMM, "%s:%d FACTOEND\n", __FILE__, __LINE__);
          MUTEX_UNLOCK(&(sopalin_data->mutex_comm));
          pthread_cond_broadcast(&(sopalin_data->cond_comm));
        }
      print_debug(DBG_SOPALIN_THREADCOMM,
                  "%d - %d : <-- SendRecv\n", (int)SOLV_PROCNUM, (int)me);

      return 0;
    }
  else
    {
      errorPrint("sendrecv_smp sould not be called in THREAD_COMM mode.\n");
      return NULL;
    }
}
#endif /* FORCE_NOMPI */
