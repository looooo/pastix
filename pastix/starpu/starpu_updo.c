/*
 * Updown step written using StarPU.
 *
 */
#  ifdef STARPU_USE_DEPRECATED_API
#    undef STARPU_USE_DEPRECATED_API
#  endif
#  include <starpu.h>
#  include "common.h"
#  include "symbol.h"
#  include "d_ftgt.h"
#  include "d_csc.h"
#  include "d_updown.h"
#  include "queue.h"
#  include "bulles.h"
#  include "d_solver.h"
#  include "sopalin_thread.h"
#  include "sopalin_define.h"
#  include "sopalin3d.h"
#  include "sopalin_acces.h"
#  include "starpu_updo.h"
#  include "starpu_updo_kernels.h"
#  define USE_TASK_DEP


#  if (STARPU_MAJOR_VERSION == 1 &&  STARPU_MINOR_VERSION == 0)
#    define STARPU_1_0
#  endif

#  if (STARPU_MAJOR_VERSION == 1 &&  STARPU_MINOR_VERSION == 1)
#    define STARPU_1_0
#    define STARPU_1_1
#  endif

#  if (STARPU_MAJOR_VERSION == 1 &&  STARPU_MINOR_VERSION == 2)
#    define STARPU_1_0
#    define STARPU_1_1
#    define STARPU_1_2
#  endif

#ifdef PASTIX_WITH_MPI
#  ifdef PASTIX_WITH_STARPU
#    define PASTIX_WITH_STARPU_MPI
#  endif
#endif
#ifndef PASTIX_WITH_STARPU_MPI
#  define starpu_mpi_init(a, ...)            starpu_init(a, __VA_ARGS__)
#  define starpu_mpi_insert_task(a, b, ...)  starpu_insert_task(b, __VA_ARGS__)
#  define starpu_mpi_data_register(...)
#  define starpu_data_get_rank(...)          0
#endif
int starpu_register_sm2x(Sopalin_Data_t       * sopalin_data,
                         starpu_data_handle_t * SM2X_handles)
{
  d_SolverMatrix * datacode = sopalin_data->datacode;
  int itercblk;
  for (itercblk=0;itercblk<SYMB_CBLKNBR;itercblk++) {
    starpu_matrix_data_register(&(SM2X_handles[itercblk]), 0,
                                (uintptr_t)&(UPDOWN_SM2XTAB[
                                               UPDOWN_SM2XIND(itercblk)]),
                                (uint32_t)UPDOWN_SM2XSZE,
                                CBLK_COLNBR(itercblk),
                                (uint32_t)UPDOWN_SM2XNBR,
                                sizeof(pastix_float_t));
    starpu_mpi_data_register(SM2X_handles[itercblk],
                             2* SOLV_GCBLKNBR + UPDOWN_LOC2GLOB(itercblk),
                             SOLV_PROCNUM);

  }
  return NO_ERR;
}

#  define updown_TRSM_model API_CALL(updown_TRSM_model)
static struct starpu_perfmodel updo_TRSM_model;

#  define updown_GEMM_model API_CALL(updown_GEMM_model)
static struct starpu_perfmodel updo_GEMM_model;

#  define updown_DIAG_model API_CALL(updown_DIAG_model)
static struct starpu_perfmodel updo_DIAG_model;

#  define updo_trsm_cl API_CALL(updo_trsm_cl)
static struct starpu_codelet updo_trsm_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = updo_trsm_starpu_cpu,
  .model = &updo_TRSM_model,
  .nbuffers = 2,
  .modes = {
    STARPU_R,
    STARPU_RW}
};

#  define updo_up_gemm_cl API_CALL(updo_up_gemm_cl)
static struct starpu_codelet updo_up_gemm_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = updo_up_gemm_starpu_cpu,
  .model = &updo_GEMM_model,
  .nbuffers = 3,
  .modes = {
    STARPU_R,
    STARPU_R,
    STARPU_RW
  }
};

#  define updo_down_gemm_cl API_CALL(updo_down_gemm_cl)
static struct starpu_codelet updo_down_gemm_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = updo_down_gemm_starpu_cpu,
  .model = &updo_GEMM_model,
  .nbuffers = 3,
  .modes = {
    STARPU_R,
    STARPU_R,
    STARPU_RW
  }
};

#  define updo_diag_cl API_CALL(updo_diag_cl)
struct starpu_codelet updo_diag_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = updo_diag_starpu_cpu,
  .model = &updo_DIAG_model,
  .nbuffers = 2,
  .modes = {
    STARPU_R,
    STARPU_RW}
};

#  define DOWN 0
#  define UP   1

static inline
int starpu_submit_up_or_down(Sopalin_Data_t * sopalin_data,
                             starpu_data_handle_t * cblk_handles,
                             starpu_data_handle_t * SM2X_handles,
                             Queue cblreadyqueue,
                             int DOWN_OR_UP,
                             int * sched_ctxs)
{
  pastix_int_t ii;
  char N = 'N', T = 'T', C = 'C', U = 'U';
  d_SolverMatrix * datacode = sopalin_data->datacode;

#ifdef STARPU_1_2
  updo_TRSM_model.type = STARPU_HISTORY_BASED;
  updo_TRSM_model.symbol = "updo_TRSM";
  starpu_perfmodel_init(&updo_TRSM_model);

  updo_GEMM_model.type = STARPU_HISTORY_BASED;
  updo_GEMM_model.symbol = "updo_GEMM";
  starpu_perfmodel_init(&updo_GEMM_model);

  updo_DIAG_model.type = STARPU_HISTORY_BASED;
  updo_DIAG_model.symbol = "updo_DIAG";
  starpu_perfmodel_init(&updo_DIAG_model);
#endif
  for (ii=0;ii<SYMB_CBLKNBR;ii++) {
    int ret;
    pastix_int_t itercblk;
    pastix_int_t iterblok;
    char * transpose, * diag;
    itercblk = queueGet(&cblreadyqueue);
    /* d_Task bloc itercblk */
    if (DOWN_OR_UP == DOWN) {
      transpose = &N;
#  if (defined CHOL_SOPALIN) && (!defined SOPALIN_LU)
      diag = &N;
#  else
      diag = &U;
#  endif
    } else {
#  ifdef CHOL_SOPALIN
      transpose = &T;
      diag = &N;
#  else
      diag = &U;
#    ifdef HERMITIAN
      transpose = &C;
#    else
      transpose = &T;
#    endif
#  endif
    }
    ret =
      starpu_mpi_insert_task(sopalin_data->sopar->pastix_comm,
                             &updo_trsm_cl,
                             STARPU_VALUE, &sopalin_data, sizeof(Sopalin_Data_t*),
                             STARPU_VALUE, &itercblk, sizeof(pastix_int_t),
                             STARPU_VALUE, transpose, sizeof(char),
                             STARPU_VALUE, diag, sizeof(char),
                             STARPU_R, cblk_handles[itercblk],
                             STARPU_RW, SM2X_handles[itercblk],
                             0);

    STARPU_ASSERT(!ret);

    if (DOWN_OR_UP == DOWN) {
      for (iterblok = SYMB_BLOKNUM(itercblk)+1;
           iterblok < SYMB_BLOKNUM(itercblk+1);
           iterblok++) {
        int ret;
        ret =
          starpu_mpi_insert_task(sopalin_data->sopar->pastix_comm,
                                 &updo_down_gemm_cl,
                                 STARPU_VALUE, &sopalin_data, sizeof(Sopalin_Data_t*),
                                 STARPU_VALUE, &itercblk, sizeof(pastix_int_t),
                                 STARPU_VALUE, &iterblok, sizeof(pastix_int_t),
                                 STARPU_VALUE, &N,        sizeof(char),
                                 STARPU_R, cblk_handles[itercblk],
                                 STARPU_R, SM2X_handles[itercblk],
                                 STARPU_RW, SM2X_handles[SYMB_CBLKNUM(iterblok)],
                                 0);
        STARPU_ASSERT(!ret);
      }
    } else {
      pastix_int_t proc;
      for (proc=UPDOWN_BROWPROCNBR(itercblk)-1;proc>=0;proc--) {
        if (UPDOWN_BROWPROCTAB(itercblk)[proc] != SOLV_PROCNUM) {
          assert(0);
        }
        else {
          pastix_int_t count;
          pastix_int_t listptridx = UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB(itercblk));
          for (count=UPDOWN_LISTPTR(listptridx);
               count<UPDOWN_LISTPTR(listptridx+1);
               count++) {
            pastix_int_t cblk;
            int ret;
            cblk      = UPDOWN_LISTCBLK(count);
            iterblok  = UPDOWN_LISTBLOK(count);

            ASSERTDBG((SYMB_FROWNUM(iterblok)>=
                       SYMB_FCOLNUM(itercblk)) &&
                      (SYMB_LROWNUM(iterblok)<=
                       SYMB_LCOLNUM(itercblk)),
                      MOD_SOPALIN);
            ret = starpu_mpi_insert_task(sopalin_data->sopar->pastix_comm,
                                         &updo_up_gemm_cl,
                                         STARPU_VALUE, &sopalin_data, sizeof(Sopalin_Data_t*),
                                         STARPU_VALUE, &cblk, sizeof(pastix_int_t),
                                         STARPU_VALUE, &iterblok, sizeof(pastix_int_t),
#  ifdef HERMITIAN
                                         STARPU_VALUE, &C, sizeof(char),
#  else
                                         STARPU_VALUE, &T, sizeof(char),
#  endif
                                         STARPU_R, cblk_handles[cblk],
                                         STARPU_R, SM2X_handles[cblk],
                                         STARPU_RW, SM2X_handles[
                                           SYMB_CBLKNUM(iterblok)],
                                         0);
            STARPU_ASSERT(!ret);
          }
        }
      }
    }
  }
  return NO_ERR;
}

int starpu_submit_updown(Sopalin_Data_t * sopalin_data,
                         starpu_data_handle_t * L_handles,
                         starpu_data_handle_t * U_handles,
                         starpu_data_handle_t * SM2X_handles,
			 struct starpu_task  ** tasktab,
                         int                  * sched_ctxs){
    return 0;
}

#define starpu_submit_updown_old API_CALL(starpu_submit_updown_old)
int starpu_submit_updown_old(Sopalin_Data_t * sopalin_data,
                             starpu_data_handle_t * L_handles,
                             starpu_data_handle_t * U_handles,
                             starpu_data_handle_t * SM2X_handles,
                             int                  * sched_ctxs)
{
  d_SolverMatrix * datacode = sopalin_data->datacode;
  Queue cblreadyqueue;
  pastix_int_t   itercblk;

  queueInit(&cblreadyqueue,SYMB_CBLKNBR);
  for (itercblk=0;itercblk<SYMB_CBLKNBR;itercblk++) {
    queueAdd(&cblreadyqueue,itercblk,(double)(TASK_PRIONUM(itercblk)));
  }

  /* Down step */
  starpu_submit_up_or_down(sopalin_data,
                           L_handles,
                           SM2X_handles,
                           cblreadyqueue, DOWN, sched_ctxs);

  /* Diag step */
#  ifndef CHOL_SOPALIN
  for (itercblk=0;itercblk<SYMB_CBLKNBR;itercblk++) {
    int ret;
    ret =
      starpu_mpi_insert_task(sopalin_data->sopar->pastix_comm,
                             &updo_diag_cl,
                             STARPU_VALUE, &sopalin_data, sizeof(Sopalin_Data_t*),
                             STARPU_VALUE, &itercblk, sizeof(pastix_int_t),
                             STARPU_R, L_handles[itercblk],
                             STARPU_RW, SM2X_handles[itercblk],
#  ifdef STARPU_GET_TASK_CTX
                             STARPU_SCHED_CTX, sched_ctxs[0],
#  endif

                             0);
    STARPU_ASSERT(!ret);
  }

#  endif /* not CHOL_SOPALIN */
  /* Up step */
  queueInit(&cblreadyqueue,SYMB_CBLKNBR);
  for (itercblk=0;itercblk<SYMB_CBLKNBR;itercblk++) {
    queueAdd(&cblreadyqueue,itercblk,-(double)(TASK_PRIONUM(itercblk)));
  }

  starpu_submit_up_or_down(sopalin_data,
#  ifdef SOPALIN_LU
                           U_handles,
#  else /* not SOPALIN_LU */
                           L_handles,
#  endif /* not SOPALIN_LU */
                           SM2X_handles,
                           cblreadyqueue, UP,
                           sched_ctxs);
  return NO_ERR;
}
