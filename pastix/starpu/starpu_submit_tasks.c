#  ifndef FORCE_NO_CUDA
#    include <cuda.h>
#  endif
#  ifdef STARPU_USE_DEPRECATED_API
#    undef STARPU_USE_DEPRECATED_API
#  endif
#  include <starpu.h>
#  ifdef WITH_STARPU_MPI
#    include <starpu_mpi.h>
#  endif
#  include <starpu_profiling.h>
#  include <pthread.h>
#  include <string.h>
#  ifdef FORCE_NOMPI
#    include "nompi.h"
#  else
#    include <mpi.h>
#  endif
#  include "common.h"
#  include "out.h"
#  include "sopalin_define.h"
#  include "sopalin_acces.h"
#  include "symbol.h"
#  include "ftgt.h"
#  include "csc.h"
#  include "updown.h"
#  include "queue.h"
#  include "bulles.h"
#  include "solver.h"
#  include "sopalin_thread.h"
#  include "sopalin_time.h"
#  include "sopalin3d.h"
#  include "starpu_kernels.h"
#  include "perf.h"
#  include "starpu_submit_tasks.h"
#  include "starpu_updo.h"
#  include "starpu_pastix_sched_policy.h"
#  include "sopalin_init.h"
#  define dump_all API_CALL(dump_all)
void  dump_all                 (SolverMatrix *, CscMatrix * cscmtx, int);

/* #define starpu_mpi_data_register(data, tag, rank) do {          \ */
/*   starpu_data_set_rank(data, rank);                             \ */
/*   starpu_data_set_tag(data, tag);                               \ */
/*   } while(0) */
#  ifdef STARPU_USE_CUDA
#    if ((!defined PREC_DOUBLE) || (!(defined __CUDA_ARCH__) || __CUDA_ARCH__ >= 130))
#      if !(defined PREC_DOUBLE && defined TYPE_COMPLEX && CUDA_SM_VERSION < 20)
#        ifndef FORCE_NO_CUDA
#          define STARPU_USE_CUDA_GEMM_FUNC
#        endif
#      endif
#    endif
#  endif

#ifndef WITH_STARPU_MPI
#  define starpu_mpi_init(a, ...)            starpu_init(a)
#  define starpu_mpi_insert_task(a, b, ...)  starpu_insert_task(b, __VA_ARGS__)
#  define starpu_mpi_data_register(...)
#  define starpu_mpi_shutdown()
#  define starpu_data_get_rank(...)          0
#endif
#  ifdef TYPE_COMPLEX
#    ifdef PREC_DOUBLE
#      define PREFIX  "Z"
#    else
#      define PREFIX  "C"
#    endif
#  else
#    ifdef PREC_DOUBLE
#      define PREFIX  "D"
#    else
#      define PREFIX  "S"
#    endif
#  endif

#  define CUDA_CALL(x) do                                       \
    {                                                           \
      if (cudaSuccess != x)                                     \
        {                                                       \
          errorPrint("%s (%s,%d)\n",#x, __FILE__,__LINE__);     \
          assert(0);                                            \
        }                                                       \
    } while(0)


#define prof_callback API_CALL(prof_callback)
void prof_callback(void *callback_arg)
{
#  if (defined PASTIX_WITH_STARPU_PROFILING)
  int workerid;
  struct starpu_task * task;
  starpu_task_stats_t * tasks_stats = (starpu_task_stats_t*)callback_arg;
  struct starpu_profiling_task_info *info;
  task = starpu_task_get_current();
  info = task->profiling_info;
  /* How much time did it take before the task started ? */
  tasks_stats[info->workerid].delay_sum +=
    starpu_timing_timespec_delay_us(&info->submit_time,
                                    &info->start_time);
  /* How long was the task execution ? */
  tasks_stats[info->workerid].length_sum +=
    starpu_timing_timespec_delay_us(&info->start_time,
                                    &info->end_time);
  tasks_stats[info->workerid].cnt++;

#    ifdef STARPU_1_2
  tasks_stats[info->workerid].ops +=
    task->cl->model->per_arch[STARPU_CPU_WORKER][0][0][0].size_base(task, NULL, 0);
#    else
  tasks_stats[info->workerid].ops +=
    task->cl->model->per_arch[STARPU_CPU_WORKER][0].size_base(task, 0, 0);
#    endif
#  endif /* (defined PASTIX_WITH_STARPU_PROFILING) */
}

#  if (STARPU_MAJOR_VERSION < 1)
#    error "PaStiX requires STARPU >= 1"
#  endif /* STARPU_MAJOR_VERSION */

#  if (STARPU_MAJOR_VERSION == 1 &&  STARPU_MINOR_VERSION == 0)
/* 1.1 : starpu_perf_archtype => starpu_perfmodel_archtype
 *       starpu_profiling_worker_info => starpu_profiling_worker_info
 */
#    define STARPU_1_0
#    define starpu_perfmodel_archtype        starpu_perf_archtype
#    define starpu_profiling_worker_info     starpu_worker_profiling_info
#    define starpu_profiling_task_info       starpu_task_profiling_info
#    define starpu_profiling_worker_get_info starpu_worker_get_profiling_info
#    define starpu_profiling_set_id          starpu_set_profiling_id
#  endif

#  if (STARPU_MAJOR_VERSION == 1 &&  STARPU_MINOR_VERSION == 1)
#    define STARPU_1_0
#    define STARPU_1_1
#  endif

#  if (STARPU_MAJOR_VERSION == 1 &&  STARPU_MINOR_VERSION == 2)
#    define STARPU_1_0
#    define STARPU_1_1
#    define STARPU_1_2
#    define STARPU_COMMUTABLE
#  endif

#ifndef STARPU_COMMUTABLE
#  define STARPU_COMMUTE 0
#endif

#  define ARCH_CPU  0
#  define ARCH_CUDA 1

#define starpu_id API_CALL(starpu_id)
static int starpu_id = 0;

static size_t trf_size(struct starpu_task *task,
#ifdef STARPU_1_2
		       struct starpu_perfmodel_arch *arch,
#else
		       enum starpu_perfmodel_archtype arch,
#endif
		       unsigned nimpl) {
  Sopalin_Data_t    * sopalin_data;
  SolverMatrix      * datacode;
  pastix_int_t          cblknum;
  pastix_int_t          stride;
  pastix_int_t          tasknum;
  size_t              dima;
  starpu_codelet_unpack_args(task->cl_arg, &sopalin_data, &cblknum, &tasknum);
  datacode = sopalin_data->datacode;
  stride   = SOLV_STRIDE(cblknum);
  dima     = CBLK_COLNBR(cblknum);
  return OPS_PPF(dima);

}
static size_t trsm_size(struct starpu_task *task,
#ifdef STARPU_1_2
			struct starpu_perfmodel_arch *arch,
#else
			enum starpu_perfmodel_archtype arch,
#endif
			unsigned nimpl) {
  Sopalin_Data_t    * sopalin_data;
  SolverMatrix      * datacode;
  pastix_int_t          cblknum;
  pastix_int_t          stride;
  pastix_int_t          tasknum;
  size_t              dima;
  starpu_codelet_unpack_args(task->cl_arg, &sopalin_data, &cblknum, &tasknum);
  datacode = sopalin_data->datacode;
  stride   = SOLV_STRIDE(cblknum);
  dima     = CBLK_COLNBR(cblknum);
  return OPS_TRSM(dima, stride);

}
static size_t trf_trsm_size(struct starpu_task *task,
#ifdef STARPU_1_2
			    struct starpu_perfmodel_arch *arch,
#else
			    enum starpu_perfmodel_archtype arch,
#endif
			    unsigned nimpl) {
  return trf_size(task, arch, nimpl) + trsm_size(task, arch, nimpl);
}

static size_t gemm_size(struct starpu_task *task,
#ifdef STARPU_1_2
                        struct starpu_perfmodel_arch *arch,
#else
                        enum starpu_perfmodel_archtype arch,
#endif
                        unsigned nimpl)
{
  Sopalin_Data_t             * sopalin_data;
  pastix_int_t                   bloknum;
  SolverMatrix               * datacode;
  pastix_int_t                   tasknum;
  pastix_int_t                   cblknum;
  pastix_int_t                   fcblknum;
  pastix_int_t                   indblok;
  pastix_int_t                   stride;
  pastix_int_t                   dimi;
  pastix_int_t                   dimj;
  pastix_int_t                   dima;

  starpu_codelet_unpack_args(task->cl_arg, &sopalin_data,
                             &cblknum, &fcblknum, &bloknum, &tasknum);

  datacode = sopalin_data->datacode;
  if (cblknum < 0) {
    /* HALO Gemm*/
    pastix_int_t hcblk = -(cblknum+1);
    indblok   = HBLOCK_COEFIND(bloknum);
    dimj = HBLOCK_ROWNBR(bloknum);
    dima = HCBLK_COLNBR(hcblk);
    stride = HCBLK_STRIDE(hcblk);
  } else {
    indblok   = SOLV_COEFIND(bloknum);
    dimj = BLOK_ROWNBR(bloknum);
    dima = CBLK_COLNBR(cblknum);
    stride = SOLV_STRIDE(cblknum);
  }
  dimi = stride - indblok;

  return OPS_GEMM(dimi,dimj,dima);
}

static struct starpu_perfmodel GEMM_model;

#  define starpu_partition_data                  API_CALL(starpu_partition_data)


static struct starpu_perfmodel XXTRF_TRSM_model;
static struct starpu_perfmodel XXTRF_model;
static struct starpu_perfmodel TRSM_model;


#      define trf_trsm_cl        API_CALL(trf_trsm_cl)
#      define trf_cl             API_CALL(trf_cl)
#      define trsm_cl            API_CALL(trsm_cl)
#      define trsm_cpu_cl        API_CALL(trsm_cpu_cl)
#      define gemm_cl            API_CALL(gemm_cl)
#      define sparse_gemm_cl     API_CALL(sparse_gemm_cl)
#      define sparse_gemm_cpu_cl API_CALL(sparse_gemm_cpu_cl)

struct starpu_codelet trf_trsm_cl =
{
  .where = STARPU_CPU
#ifdef CHOL_SOPALIN
#  ifndef FORCE_NO_CUDA
#    ifdef WITH_MAGMABLAS
  |STARPU_CUDA
#    endif /* WITH_MAGMABLAS */
#  endif /* not FORCE_NO_CUDA */
#endif /* CHOL_SOPALIN */
  ,
  .cpu_funcs[0] = trfsp1d_starpu_cpu,
#ifdef CHOL_SOPALIN
#  ifndef FORCE_NO_CUDA
#    ifdef STARPU_USE_CUDA
#      ifdef WITH_MAGMABLAS
  .cuda_funcs[0] = trfsp1d_starpu_cuda,
#      endif /* WITH_MAGMABLAS */
#    endif /* STARPU_USE_CUDA */
#  endif /* not FORCE_NO_CUDA */
# endif /* CHOL_SOPALIN */
  .model = &XXTRF_TRSM_model,
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
  /* LU */
  .nbuffers = 2,
  .modes = {
    STARPU_RW,
    STARPU_RW
  }
#    else /* not SOPALIN_LU */
  /* LDt */
  .nbuffers = 1,
  .modes = {
    STARPU_RW
  }
#    endif /* w/wo SOPALIN_LUN */
#  else  /* not CHOL_SOPALIN */
  /* LDLt/LDLh */
  .nbuffers = 2,
  .modes = {
    STARPU_RW,
    STARPU_SCRATCH
  }
#  endif /* w/wo CHOL_SOPALIN */
};

struct starpu_codelet trf_cl =
{
  .where = STARPU_CPU
#ifdef CHOL_SOPALIN
#  ifndef FORCE_NO_CUDA
#    ifdef WITH_MAGMABLAS
  |STARPU_CUDA
#    endif /* WITH_MAGMABLAS */
#  endif /* not FORCE_NO_CUDA */
#endif /* CHOL_SOPALIN */
  ,
  .cpu_funcs[0] = xxtrf_starpu_cpu,
#ifdef CHOL_SOPALIN
#  ifndef FORCE_NO_CUDA
#    ifdef STARPU_USE_CUDA
#      ifdef WITH_MAGMABLAS
  .cuda_funcs[0] = xxtrf_starpu_cuda,
#      endif /* WITH_MAGMABLAS */
#    endif /* STARPU_USE_CUDA */
#  endif /* not FORCE_NO_CUDA */
#endif /* CHOL_SOPALIN */
  .model = &XXTRF_model,
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
  /* LU */
  .nbuffers = 2,
  .modes = {
    STARPU_RW,
    STARPU_RW
  }
#    else /* not SOPALIN_LU */
  /* LDt */
  .nbuffers = 1,
  .modes = {
    STARPU_RW
  }
#    endif /* w/wo SOPALIN_LUN */
#  else  /* not CHOL_SOPALIN */
  /* LDLt/LDLh */
  .nbuffers = 2,
  .modes = {
    STARPU_RW,
    STARPU_SCRATCH
  }
#  endif /* w/wo CHOL_SOPALIN */
};

struct starpu_codelet trsm_cl =
{
  .where = STARPU_CPU
#  ifndef FORCE_NO_CUDA
  |STARPU_CUDA
#  endif /* not FORCE_NO_CUDA */
  ,
  .cpu_funcs[0] = trsm_starpu_cpu,
#  ifndef FORCE_NO_CUDA
  .cuda_funcs[0] = trsm_starpu_cuda,
#  endif /* not FORCE_NO_CUDA */
  .model = &TRSM_model,
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
  /* LU */
  .nbuffers = 2,
  .modes = {
    STARPU_RW,
    STARPU_RW
  }
#    else /* not SOPALIN_LU */
  /* LDt */
  .nbuffers = 1,
  .modes = {
    STARPU_RW
  }
#    endif /* w/wo SOPALIN_LUN */
#  else  /* not CHOL_SOPALIN */
  /* LDLt/LDLh */
  .nbuffers = 2,
  .modes = {
    STARPU_RW,
    STARPU_SCRATCH
  }
#  endif /* w/wo CHOL_SOPALIN */
};

struct starpu_codelet trsm_cpu_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = trsm_starpu_cpu,
  .model = &TRSM_model,
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
  /* LU */
  .nbuffers = 2,
  .modes = {
    STARPU_RW,
    STARPU_RW
  }
#    else /* not SOPALIN_LU */
  /* LDt */
  .nbuffers = 1,
  .modes = {
    STARPU_RW
  }
#    endif /* w/wo SOPALIN_LUN */
#  else  /* not CHOL_SOPALIN */
  /* LDLt/LDLh */
  .nbuffers = 2,
  .modes = {
    STARPU_RW,
    STARPU_SCRATCH
  }
#  endif /* w/wo CHOL_SOPALIN */
};

struct starpu_codelet sparse_gemm_cpu_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = trfsp1d_sparse_gemm_starpu_cpu,
  .model = &GEMM_model,
#  if (defined CHOL_SOPALIN && defined SOPALIN_LU)
  .nbuffers = 6,
  .modes = {
    STARPU_R,
    STARPU_RW|STARPU_COMMUTE,
    STARPU_R,
    STARPU_RW|STARPU_COMMUTE,
    STARPU_SCRATCH,
    STARPU_R
  }
#  else /* not LU */
  .nbuffers = 4,
  .modes = {
    STARPU_R,
    STARPU_RW|STARPU_COMMUTE,
    STARPU_SCRATCH,
    STARPU_R
  }
#  endif /* LU / not LU */
};

struct starpu_codelet sparse_gemm_cl =
{
  .where = STARPU_CPU
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  |STARPU_CUDA
#      endif
  ,
  .cpu_funcs[0] = trfsp1d_sparse_gemm_starpu_cpu,
#      ifdef STARPU_USE_CUDA_GEMM_FUNC
  .cuda_funcs[0] = trfsp1d_sparse_gemm_starpu_cuda,
#      endif
  .model = &GEMM_model,
#  if (defined CHOL_SOPALIN && defined SOPALIN_LU)
  .nbuffers = 6,
  .modes = {
    STARPU_R,
    STARPU_RW|STARPU_COMMUTE,
    STARPU_R,
    STARPU_RW|STARPU_COMMUTE,
    STARPU_SCRATCH,
    STARPU_R
  }
#  else /* not LU */
  .nbuffers = 4,
  .modes = {
    STARPU_R,
    STARPU_RW|STARPU_COMMUTE,
    STARPU_SCRATCH,
    STARPU_R
  }
#  endif /* LU / not LU */
};

/*
 Function: starpu_data_partition

 Initialize column blocks handlers.

 Parameters:
 sopalin_data - PaStiX global data structure.
 L_handle     - Handles for L column blocks.
 U_handle     - Handles for U column blocks.
 */
static inline
void starpu_partition_data(Sopalin_Data_t * sopalin_data,
                           starpu_data_handle_t * L_handle,
                           starpu_data_handle_t * U_handle,
                           starpu_data_handle_t * Lhalo_handle,
                           starpu_data_handle_t * Uhalo_handle)
{
  SolverMatrix       * datacode         = sopalin_data->datacode;
  pastix_int_t itercblk;
  for (itercblk=0;itercblk<SYMB_CBLKNBR;itercblk++) {
    pastix_int_t iter;
    starpu_matrix_data_register(&(L_handle[itercblk]), 0,
                                (uintptr_t)SOLV_COEFTAB(itercblk),
                                (uint32_t)SOLV_STRIDE(itercblk),
                                (uint32_t)SOLV_STRIDE(itercblk),
                                CBLK_COLNBR(itercblk),
                                sizeof(pastix_float_t));
    starpu_mpi_data_register(L_handle[itercblk],
                             UPDOWN_LOC2GLOB(itercblk),
                             SOLV_PROCNUM);

#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
    starpu_matrix_data_register(&(U_handle[itercblk]), 0,
                                (uintptr_t)SOLV_UCOEFTAB(itercblk),
                                (uint32_t)SOLV_STRIDE(itercblk),
                                (uint32_t)SOLV_STRIDE(itercblk),
                                CBLK_COLNBR(itercblk),
                                sizeof(pastix_float_t));
    starpu_mpi_data_register(U_handle[itercblk],
                             SOLV_GCBLKNBR + UPDOWN_LOC2GLOB(itercblk),
                             SOLV_PROCNUM);
#    endif
#  endif
  }
  /* register halo */
  for (itercblk=0;itercblk<SOLV_HCBLKNBR;itercblk++) {
    pastix_int_t iter;
    ASSERT(SOLV_GCBLK2HALO(HCBLK_GCBLK(itercblk)) == itercblk,
           MOD_SOPALIN);
    starpu_matrix_data_register(&(Lhalo_handle[itercblk]), -1,
                                (uintptr_t)NULL,
                                (uint32_t)HCBLK_STRIDE(itercblk),
                                (uint32_t)HCBLK_STRIDE(itercblk),
                                HCBLK_COLNBR(itercblk),
                                sizeof(pastix_float_t));
    starpu_mpi_data_register(Lhalo_handle[itercblk],
                             HCBLK_GCBLK(itercblk),
                             HCBLK_OWNER(itercblk));

#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
    starpu_matrix_data_register(&(Uhalo_handle[itercblk]), -1,
                                (uintptr_t)NULL,
                                (uint32_t)HCBLK_STRIDE(itercblk),
                                (uint32_t)HCBLK_STRIDE(itercblk),
                                HCBLK_COLNBR(itercblk),
                                sizeof(pastix_float_t));
    starpu_mpi_data_register(Uhalo_handle[itercblk],
                             SOLV_GCBLKNBR + HCBLK_GCBLK(itercblk),
                             HCBLK_OWNER(itercblk));
#    endif
#  endif
  }
}

/*
 * Function: starpu_init_smp
 *
 * Initialize thread data structure for factorization when using StarPU.
 */
#  define starpu_init_smp API_CALL(starpu_init_smp)
void*
starpu_init_smp (void * arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  int init;
  init = INIT_COMPUTE;
  if (THREAD_FUNNELED_OFF)
    {
      init = init | INIT_SEND;
      if (THREAD_COMM_OFF)
        {
          init = init | INIT_RECV;
        }
    }
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    {
      sopalin_init_smp(sopalin_data, argument->me, API_YES, init);
    }
  else
    {
      sopalin_init_smp(sopalin_data, argument->me, API_NO, init);
    }
  return NULL;
}

#define starpu_init_kernel API_CALL(starpu_init_kernel)
void starpu_init_kernel(void * buffers[], void * _args)
{
  starpu_init_smp(_args);
}
#define sopalin_init_cl API_CALL(sopalin_init_cl)
struct starpu_codelet sopalin_init_cl =
{
  .where = STARPU_CPU,
  .cpu_funcs[0] = starpu_init_kernel,
  .nbuffers = 0,
  .modes = {}
};

/*
 Struct: sopthread_data

 Structure conotaining the thread number and a pointer to data given to thread in parameters.
 */
struct starpu_loop_data_ {
  int                    me;
  starpu_data_handle_t * L_handle;
  starpu_data_handle_t * Lhalo_handle;
#  if (defined CHOL_SOPALIN && defined SOPALIN_LU)
  starpu_data_handle_t * U_handle;
  starpu_data_handle_t * Uhalo_handle;
#  endif
  starpu_data_handle_t   WORK_handle;
  starpu_data_handle_t * blocktab_handles;
  Sopalin_Data_t       * sopalin_data;                 /*+ Data given to thread as argument +*/
  int                    ctx;
  int                    ctx_nbr;
  int                    thread_per_ctx;
  int                  * sched_ctxs;
  int                    first, last;
  int                   * facto_finished;
  pthread_cond_t        * cond_end_facto;
  pthread_mutex_t       * mutex_end_facto;
  int                   * cpu_workerids;
  int                     ncpus;
  int                   * gpu_workerids;
  int                     ngpus;
  int                   * gpu_gemm_count;
};

static inline int
halo_submit(Sopalin_Data_t * sopalin_data) {
  SolverMatrix * datacode = sopalin_data->datacode;
  starpu_loop_data_t  *starpu_loop_data  = sopalin_data->starpu_loop_data;
  pastix_int_t itertask, workerid = -1;
  starpu_data_handle_t *L_handle         = starpu_loop_data->L_handle;
  starpu_data_handle_t *Lhalo_handle     = starpu_loop_data->Lhalo_handle;
#  if (defined CHOL_SOPALIN)
#    ifdef SOPALIN_LU
  starpu_data_handle_t *U_handle         = starpu_loop_data->U_handle;
  starpu_data_handle_t *Uhalo_handle     = starpu_loop_data->Uhalo_handle;
#    endif
#  endif
  struct starpu_codelet * cl;
  pastix_int_t max_deps = 0;
  int ret;

  for (itertask=0;itertask<SOLV_TASKNBR;itertask++) {
    pastix_int_t itercblk = TASK_CBLKNUM(itertask);
    pastix_int_t gcblk2glist = UPDOWN_GCBLK2GLIST(UPDOWN_LOC2GLOB( itercblk ));
    pastix_int_t browk      = ( gcblk2glist != -1 )?
      UPDOWN_GLISTPTR( gcblk2glist    ):-1;
    pastix_int_t browk1     = ( gcblk2glist != -1 )?
      UPDOWN_GLISTPTR( gcblk2glist + 1):-1;
    pastix_int_t ndeps, iter;
    ndeps = browk1-browk;

    max_deps = MAX(max_deps, ndeps);

    /* several deps can be inserted during same loop if involving same CBLKs */
    for (iter = 0; iter < ndeps;) {
      pastix_int_t gcblk= UPDOWN_GLISTCBLK(browk+iter);
      pastix_int_t lcblk = SOLV_GCBLK2LOC(gcblk);
      if (lcblk == 0) {
        errorPrint("This cblk (%d) should not appear on this proc (%d)\n",
                   gcblk, SOLV_PROCNUM);
        ASSERT(lcblk != 0, MOD_SOPALIN);
      } else if (lcblk > 0) {
        /* ignore Local cblk */
        iter++;
      } else {
        /* halo cblk */
        pastix_int_t hcblk = -(lcblk+1);
        pastix_int_t hblock = HCBLK_BLOKNUM(hcblk);
        pastix_int_t gfcblk = UPDOWN_LOC2GLOB(itercblk);
        /* submit GEMM halo task */
#if defined(PASTIX_WITH_CUDA)
        if ( sopalin_data->sopar->iparm[IPARM_CUDA_NBR] > 0 &&
             SOLV_COLOR(itercblk) >= 0) {
          cl = &sparse_gemm_cl;
          workerid = starpu_loop_data->gpu_workerids[SOLV_COLOR(itercblk)];
          starpu_loop_data->gpu_gemm_count[SOLV_COLOR(itercblk)]++;
        }
#endif /* defined(PASTIX_WITH_CUDA) */
        /* search first block facing itercblk */
        while(!HBLOCK_ISFACING(hblock, SYMB_BLOKNUM(itercblk))) {
          hblock++;
          ASSERT(hblock < HCBLK_BLOKNUM(hcblk+1), MOD_SOPALIN);
        }
        /* insert all tasks involving same CBLKs */
        while(hblock < SOLV_HBLOKNBR &&
              HBLOCK_ISFACING(hblock, SYMB_BLOKNUM(itercblk))) {
          ret =
            starpu_mpi_insert_task(sopalin_data->sopar->pastix_comm, &sparse_gemm_cl,
                                   STARPU_VALUE, &sopalin_data, sizeof(Sopalin_Data_t*),
                                   STARPU_VALUE, &lcblk,        sizeof(pastix_int_t),
                                   STARPU_VALUE, &itercblk,     sizeof(pastix_int_t),
                                   STARPU_VALUE, &hblock,       sizeof(pastix_int_t),
                                   STARPU_VALUE, &itertask,     sizeof(pastix_int_t),
                                   STARPU_R,     Lhalo_handle[hcblk],
                                   STARPU_COMMUTE|STARPU_RW,    L_handle[itercblk],
#  if (defined CHOL_SOPALIN)
#    ifdef SOPALIN_LU
                                   STARPU_R,     Uhalo_handle[hcblk],
                                   STARPU_COMMUTE|STARPU_RW,    U_handle[itercblk],
#    endif
#  endif /* CHOL_SOPALIN */
                                   STARPU_SCRATCH, starpu_loop_data->WORK_handle,
                                   STARPU_R, starpu_loop_data->blocktab_handles[SOLV_PROCNUM],
#ifdef STARPU_1_1
                                   STARPU_EXECUTE_ON_WORKER, workerid,
#endif
                                   STARPU_CALLBACK,     prof_callback,
                                   STARPU_CALLBACK_ARG, sopalin_data->hgemm_stats,
#  ifdef STARPU_CONTEXT
                                   STARPU_SCHED_CTX, sched_ctxs[my_ctx],
#  endif
                                   0);
	  STARPU_CHECK_RETURN_VALUE(ret, "starpu_mpi_insert_task");
          TASK_CTRBCNT(itertask)++;
          hblock++;
          iter++;
        }
      }
    }
  }

  return NO_ERR;
}
int
starpu_submit_one_trf(pastix_int_t itertask, Sopalin_Data_t * sopalin_data)
{
  SolverMatrix        *datacode          = sopalin_data->datacode;
  starpu_loop_data_t  *starpu_loop_data  = (starpu_loop_data_t *)sopalin_data->starpu_loop_data;
  int                  me                = starpu_loop_data->me;
  int                 *sched_ctxs        = starpu_loop_data->sched_ctxs;
  pastix_int_t itercblk = TASK_CBLKNUM(itertask);
  int        ret;
  pastix_int_t iterbloc;
  struct starpu_task *task_diag;
  int this_workerid;;
  int workerid = -1;
  char * separate_trsm;
  struct starpu_codelet * cl;
#  ifdef STARPU_CONTEXT
  pastix_int_t threadid = TASK_THREADID(itertask);
  pastix_int_t my_ctx;

  if (threadid > SOLV_THRDNBR)
    my_ctx = 0;
  else
    my_ctx = 1+threadid/starpu_loop_data->thread_per_ctx;
#  endif
  this_workerid = starpu_worker_get_id();
  if (this_workerid == -1) this_workerid = 0;

#  ifdef STARPU_PASTIX_SCHED
  if (TASK_THREADID(itertask) < starpu_loop_data->ncpus) {
    workerid = starpu_loop_data->cpu_workerids[TASK_THREADID(itertask)];
  } else {
    workerid = starpu_loop_data->gpu_workerids[TASK_THREADID(itertask) -
                                               starpu_loop_data->ncpus];
  }
#  endif

  if ((separate_trsm = getenv("PASTIX_STARPU_SEPARATE_TRSM")) && !strcmp(separate_trsm, "1")) {
    ret =
      starpu_mpi_insert_task(sopalin_data->sopar->pastix_comm, &trf_cl,
                             STARPU_VALUE, &sopalin_data, sizeof(Sopalin_Data_t*),
                             STARPU_VALUE, &itercblk, sizeof(pastix_int_t),
                             STARPU_VALUE, &itertask, sizeof(pastix_int_t),
                             STARPU_RW, starpu_loop_data->L_handle[itercblk],
#  if (defined CHOL_SOPALIN && defined SOPALIN_LU)
                             STARPU_RW, starpu_loop_data->U_handle[itercblk],
#  endif
#  ifndef CHOL_SOPALIN
                             STARPU_SCRATCH, starpu_loop_data->WORK_handle,
#  endif
                             STARPU_PRIORITY, TASK_PRIONUM(itertask),
                             STARPU_CALLBACK,     prof_callback,
                             STARPU_CALLBACK_ARG, sopalin_data->xxtrf_stats,
#  ifdef STARPU_CONTEXT
                             STARPU_SCHED_CTX, sched_ctxs[my_ctx],
#  endif
                             0);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_mpi_insert_task");

    ret =
      starpu_mpi_insert_task(sopalin_data->sopar->pastix_comm, &trsm_cl,
                             STARPU_VALUE, &sopalin_data, sizeof(Sopalin_Data_t*),
                             STARPU_VALUE, &itercblk, sizeof(pastix_int_t),
                             STARPU_VALUE, &itertask, sizeof(pastix_int_t),
                             STARPU_RW, starpu_loop_data->L_handle[itercblk],
#  if (defined CHOL_SOPALIN && defined SOPALIN_LU)
                             STARPU_RW, starpu_loop_data->U_handle[itercblk],
#  endif
#  ifndef CHOL_SOPALIN
                             STARPU_SCRATCH, starpu_loop_data->WORK_handle,
#  endif
#ifdef STARPU_1_1
                             STARPU_EXECUTE_ON_WORKER, workerid,
#endif

                             STARPU_PRIORITY, TASK_PRIONUM(itertask),
                             STARPU_CALLBACK,     prof_callback,
                             STARPU_CALLBACK_ARG, sopalin_data->trsm_stats,
#  ifdef STARPU_CONTEXT
                             STARPU_SCHED_CTX, sched_ctxs[my_ctx],
#  endif
                             0);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_mpi_insert_task");
  } else {
    ret =
      starpu_mpi_insert_task(sopalin_data->sopar->pastix_comm, &trf_trsm_cl,
                             STARPU_VALUE, &sopalin_data, sizeof(Sopalin_Data_t*),
                             STARPU_VALUE, &itercblk, sizeof(pastix_int_t),
                             STARPU_VALUE, &itertask, sizeof(pastix_int_t),
                             STARPU_RW, starpu_loop_data->L_handle[itercblk],
#  if (defined CHOL_SOPALIN && defined SOPALIN_LU)
                             STARPU_RW, starpu_loop_data->U_handle[itercblk],
#  endif
#  ifndef CHOL_SOPALIN
                             STARPU_SCRATCH, starpu_loop_data->WORK_handle,
#  endif
                             STARPU_PRIORITY, TASK_PRIONUM(itertask),
                             STARPU_CALLBACK,     prof_callback,
                             STARPU_CALLBACK_ARG, sopalin_data->xxtrf_stats,
#  ifdef STARPU_CONTEXT
                             STARPU_SCHED_CTX, sched_ctxs[my_ctx],
#  endif
                             0);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_mpi_insert_task");
  }
  return NO_ERR;
}


int
starpu_submit_bunch_of_gemm (pastix_int_t itertask, Sopalin_Data_t * sopalin_data)
{
  SolverMatrix        *datacode          = sopalin_data->datacode;
  starpu_loop_data_t  *starpu_loop_data  = (starpu_loop_data_t *)sopalin_data->starpu_loop_data;
  int                  me                = starpu_loop_data->me;
  int                 *sched_ctxs        = starpu_loop_data->sched_ctxs;
  starpu_data_handle_t *L_handle         = starpu_loop_data->L_handle;
  starpu_data_handle_t *Lhalo_handle     = starpu_loop_data->Lhalo_handle;
#  if (defined CHOL_SOPALIN)
#    ifdef SOPALIN_LU
  starpu_data_handle_t *U_handle         = starpu_loop_data->U_handle;
  starpu_data_handle_t *Uhalo_handle     = starpu_loop_data->Uhalo_handle;
#    endif
#  endif

  pastix_int_t itercblk = TASK_CBLKNUM(itertask);
  int        ret;
  pastix_int_t iterbloc;
  pastix_int_t handle_idx;
  int workerid, this_workerid;
  struct starpu_codelet * cl;

#  ifdef STARPU_CONTEXT
  pastix_int_t threadid = TASK_THREADID(itertask);
  pastix_int_t my_ctx;

  if (threadid > SOLV_THRDNBR)
    my_ctx = 0;
  else
    my_ctx = 1+threadid/starpu_loop_data->thread_per_ctx;
#  endif

  for (iterbloc = SYMB_BLOKNUM(itercblk)+1;
       iterbloc < SYMB_BLOKNUM(itercblk+1);
       iterbloc ++)
    {
      /* contribution is local */
      struct starpu_task * task_gemm;
      pastix_int_t blocnbr;
      pastix_int_t fcblknum;
      pastix_int_t n,t;
      int dst_cblk;
      starpu_data_handle_t * L_target_handle;
#  if (defined CHOL_SOPALIN)
#    ifdef SOPALIN_LU
      starpu_data_handle_t * U_target_handle;
#    endif
#  endif
      fcblknum = SYMB_CBLKNUM(iterbloc);
      cl = &sparse_gemm_cpu_cl;
#ifndef STARPU_PASTIX_SCHED
      workerid=-1; /* Let StarPU choose CPU worker */
#endif
      if (fcblknum <  0) {
        /* itercblk udates a remote cblk */
        pastix_int_t gfcblknum = -(fcblknum+1);
        fcblknum = SOLV_GCBLK2HALO(gfcblknum);

        L_target_handle = &(Lhalo_handle[fcblknum]);
#  if (defined CHOL_SOPALIN)
#    ifdef SOPALIN_LU
        U_target_handle = &(Uhalo_handle[fcblknum]);
#    endif
#  endif
      } else {
        /* fcblknum is local */
        L_target_handle = &(L_handle[fcblknum]);
#  if (defined CHOL_SOPALIN)
#    ifdef SOPALIN_LU
        U_target_handle = &(U_handle[fcblknum]);
#    endif
#  endif

#if defined(PASTIX_WITH_CUDA)
        if ( starpu_loop_data->ngpus > 0 &&
             SOLV_COLOR(SYMB_CBLKNUM(iterbloc)) >= 0) {
          cl = &sparse_gemm_cl;
          workerid = starpu_loop_data->gpu_workerids[SOLV_COLOR(SYMB_CBLKNUM(iterbloc))];
          starpu_loop_data->gpu_gemm_count[SOLV_COLOR(SYMB_CBLKNUM(iterbloc))]++;
        }
#endif  /* defined(PASTIX_WITH_CUDA) */
      }
      dst_cblk = starpu_data_get_rank(*L_target_handle);
      ret =
        starpu_mpi_insert_task(sopalin_data->sopar->pastix_comm, cl,
			       STARPU_VALUE, &sopalin_data, sizeof(Sopalin_Data_t*),
			       STARPU_VALUE, &itercblk,     sizeof(pastix_int_t),
			       STARPU_VALUE, &fcblknum,     sizeof(pastix_int_t),
			       STARPU_VALUE, &iterbloc,     sizeof(pastix_int_t),
			       STARPU_VALUE, &itertask,     sizeof(pastix_int_t),
			       STARPU_R,     L_handle[itercblk],
			       STARPU_COMMUTE|STARPU_RW,    *L_target_handle,
#  if (defined CHOL_SOPALIN)
#    ifdef SOPALIN_LU
			       STARPU_R,     U_handle[itercblk],
			       STARPU_COMMUTE|STARPU_RW,    *U_target_handle,
#    endif
#  endif /* CHOL_SOPALIN */
			       STARPU_SCRATCH, starpu_loop_data->WORK_handle,
			       STARPU_R,  starpu_loop_data->blocktab_handles[dst_cblk],
#ifdef STARPU_1_1
			       STARPU_EXECUTE_ON_WORKER, workerid,
#endif
			       STARPU_CALLBACK,     prof_callback,
			       STARPU_CALLBACK_ARG, sopalin_data->gemm_stats,
#  ifdef STARPU_CONTEXT
			       STARPU_SCHED_CTX, sched_ctxs[my_ctx],
#  endif
			       0);
      if (ret != -ENODEV) STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
      STARPU_ASSERT(!ret);
    }

  /* If we are submiting the last cblknum updates we can tell the main thread
   * he can wait for all tasks to complete calling starpu_task_wait_for_all()
   */
  if (itercblk == SYMB_CBLKNBR-1) {
    int rc;
    rc = pthread_mutex_lock(starpu_loop_data->mutex_end_facto);
    if (rc) {
      perror("pthread_mutex_lock");
      exit(1);
    }
    *(starpu_loop_data->facto_finished) = API_YES;
    rc = pthread_cond_signal(starpu_loop_data->cond_end_facto);
    if (rc) {
      pthread_mutex_unlock(starpu_loop_data->mutex_end_facto);
      perror("pthread_cond_signal");
      exit(1);
    }
    rc = pthread_mutex_unlock(starpu_loop_data->mutex_end_facto);
    if (rc) {
      perror("pthread_mutex_unlock");
      exit(1);
    }
  }
  return NO_ERR;
};
/*
 * Function: starpu_submit_loop
 *
 * Submit the tasks.
 */
#  define starpu_submit_loop API_CALL(starpu_submit_loop)
void*
starpu_submit_loop (void * arg)
{
  starpu_loop_data_t  *starpu_loop_data  = (starpu_loop_data_t*)(arg);
  Sopalin_Data_t      *sopalin_data      = (Sopalin_Data_t *)(starpu_loop_data->sopalin_data);
  SolverMatrix        *datacode          = sopalin_data->datacode;
  int                  me                = starpu_loop_data->me;
  int                 *sched_ctxs        = starpu_loop_data->sched_ctxs;
  pastix_int_t itertask;
  pastix_int_t n_cblks = 0, n_tasks = 0;
  char      *prefetch;
  if ((prefetch = getenv("PASTIX_STARPU_PREFETCH_ON_NODE")) && !strcmp(prefetch, "1")) {
    /* Prefetch data on GPUs */
    pastix_int_t   iterworker;
    pastix_int_t * memory_nodes;
    fprintf(stdout, "Prefetching data on GPUs %s\n", prefetch);
    MALLOC_INTERN(memory_nodes, starpu_loop_data->ngpus, pastix_int_t);
    for (iterworker = 0; iterworker < starpu_loop_data->ngpus; iterworker++) {
      memory_nodes[iterworker] = starpu_worker_get_memory_node(starpu_loop_data->gpu_workerids[iterworker]);
    }
    for (itertask=0;itertask<SOLV_TASKNBR;itertask++) {
      pastix_int_t itercblk = TASK_CBLKNUM(itertask);
#if defined(PASTIX_WITH_CUDA)
      if (starpu_loop_data->ngpus > 0 &&
          SOLV_COLOR(itercblk) >= 0) {
        pastix_int_t workerid = starpu_loop_data->gpu_workerids[SOLV_COLOR(itercblk)];
        pastix_int_t node = memory_nodes[workerid];
        starpu_data_prefetch_on_node(starpu_loop_data->L_handle[itercblk],
                                     node, 1);
#  if (defined CHOL_SOPALIN)
#    ifdef SOPALIN_LU
        starpu_data_prefetch_on_node(starpu_loop_data->U_handle[itercblk],
                                     node, 1);
#    endif
#  endif
      }
#endif /* defined(PASTIX_WITH_CUDA) */
    }
    memFree_null(memory_nodes);
  }

  /* For all ready task we submit factorization  */
  for (itertask=0;itertask<SOLV_TASKNBR;itertask++) {
    char * nested;
    if (((nested = getenv("PASTIX_STARPU_NESTED_TASK")) && !strcmp(nested, "1")) &&
	TASK_CTRBCNT((itertask))) continue;
    starpu_submit_one_trf(itertask, sopalin_data);
    if (!((nested = getenv("PASTIX_STARPU_NESTED_TASK")) &&
	  !strcmp(nested, "1"))) {
      starpu_submit_bunch_of_gemm(itertask, sopalin_data);
    }
  }
  return NULL;
}

/*
 * Funciton starpu_clean_smp
 *
 * Clean thread data structures when using starpu.
 */
#  define starpu_clean_smp API_CALL(starpu_clean_smp)
void*
starpu_clean_smp (void * arg)
{
  sopthread_data_t *argument     = (sopthread_data_t *)arg;
  Sopalin_Data_t   *sopalin_data = (Sopalin_Data_t *)(argument->data);
  sopalin_clean_smp ( sopalin_data, argument->me );
  return NULL;
}

/*
 Function: starpu_submit_tasks

 Submit tasks to perform the decomposition of the matrix.

 Parameters:
 sopalin_data - PaStiX global data structure.

 Returns:
 NO_ERR
 */
int
starpu_submit_tasks(Sopalin_Data_t  * sopalin_data) {
  SolverMatrix         * datacode         = sopalin_data->datacode;
  Thread_Data_t        * thread_data;
  starpu_data_handle_t * L_handle;
  starpu_data_handle_t * Lhalo_handle;
  starpu_data_handle_t * SM2X_handles = NULL;
  pastix_int_t             task_number;
  starpu_data_handle_t * U_handle = NULL;
  starpu_data_handle_t * Uhalo_handle = NULL;
  starpu_data_handle_t   WORK_handle;
  pastix_int_t itertask;
  int * blocktab;
  starpu_data_handle_t * blocktab_handles;
  struct starpu_conf     conf;
  int                    ret;
  int                    cuda_nbr = sopalin_data->sopar->iparm[IPARM_CUDA_NBR];
  unsigned int * sched_ctxs = NULL;
  int iter;

#  ifdef STARPU_CONTEXT
  int * devices;
  pastix_int_t thread_per_ctx;
#  endif

#ifndef STARPU_INIT_SMP
  sopalin_launch_thread(sopalin_data,
                        SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree,
                        sopalin_data->sopar->iparm[IPARM_VERBOSE],
                        SOLV_THRDNBR+cuda_nbr, starpu_init_smp, sopalin_data,
                        0, NULL, NULL,
                        0, NULL, NULL);
#endif /* not STARPU_INIT_SMP ¨*/

#ifdef STARPU_PASTIX_SCHED
  {
    int k, bubnbr = datacode->bublnbr;
    int priomin = INT_MAX, priomax = INT_MIN;

    for(k = 0; k<bubnbr; k++)
      {
        priomin = MIN(priomin, datacode->btree->nodetab[k].priomin);
        priomax = MAX(priomax, datacode->btree->nodetab[k].priomax);
      }
    starpu_sched_set_min_priority(priomin);
    starpu_sched_set_max_priority(priomax);
  }
#endif

  starpu_conf_init(&conf);
  /* 1 GB */
#if (STARPU_MAJOR_VERSION > 1 || (STARPU_MAJOR_VERSION == 1 && STARPU_MINOR_VERSION >= 1))
  conf.trace_buffer_size = 1<<30;
#endif

  if (NULL != conf.sched_policy_name)
    {
      if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        fprintf(stdout, OUT_STARPU_TP, conf.sched_policy_name);
    }
  else
    {
#ifdef STARPU_PASTIX_SCHED
      conf.sched_policy_name = NULL;
      conf.sched_policy = &starpu_pastix_sched_policy;
#else
      conf.sched_policy_name = "dmda";
      if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        fprintf(stdout, OUT_STARPU_STP, conf.sched_policy_name);
#endif
    }
  conf.ncpus = SOLV_THRDNBR;
  conf.ncuda = cuda_nbr;
  conf.nopencl = 0;

  //starpu_profiling_set_id(starpu_id++);
  ret = starpu_init(&conf);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_init");
  ret = starpu_mpi_init(NULL, NULL, 0);
  STARPU_CHECK_RETURN_VALUE(ret, "starpu_mpi_init");
#ifdef STARPU_USE_CUDA
  if (cuda_nbr) {
    starpu_cublas_init();
  }
#endif /* STARPU_USE_CUDA */

  GEMM_model.type = STARPU_HISTORY_BASED;
  GEMM_model.symbol = PREFIX "GEMM";
#  ifdef STARPU_1_2
  starpu_perfmodel_init(&GEMM_model);
  GEMM_model.per_arch[STARPU_CPU_WORKER][0][0][0].size_base = gemm_size;
  if(starpu_worker_get_count_by_type(STARPU_CUDA_WORKER) != 0)
    GEMM_model.per_arch[STARPU_CUDA_WORKER][0][0][0].size_base = gemm_size;
#  else
  GEMM_model.per_arch[STARPU_CPU_WORKER][0].size_base = gemm_size;
  GEMM_model.per_arch[STARPU_CUDA_WORKER][0].size_base = gemm_size;
#  endif

  XXTRF_TRSM_model.type = STARPU_REGRESSION_BASED;
  XXTRF_model.type = STARPU_REGRESSION_BASED;
  TRSM_model.type = STARPU_REGRESSION_BASED;
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
/* LU */
  XXTRF_TRSM_model.symbol = PREFIX "GETRF_TRSM";
  XXTRF_model.symbol = PREFIX "GETRF";
#    else /* SOPALIN_LU */
  XXTRF_TRSM_model.symbol = PREFIX "POTRF_TRSM";
  XXTRF_model.symbol = PREFIX "POTRF";
#    endif /* SOPALIN_LU */
#  else  /* CHOL_SOPALIN */
/* LDLT */
  XXTRF_TRSM_model.symbol = PREFIX "HETRF_TRSM";
  XXTRF_model.symbol = PREFIX "HETRF";
#  endif  /* CHOL_SOPALIN */
  TRSM_model.symbol = PREFIX "TRSM";

#  ifdef STARPU_1_2
  starpu_perfmodel_init(&XXTRF_TRSM_model);
  XXTRF_TRSM_model.per_arch[STARPU_CPU_WORKER][0][0][0].size_base = trf_trsm_size;
  if(starpu_worker_get_count_by_type(STARPU_CUDA_WORKER) != 0)
    XXTRF_TRSM_model.per_arch[STARPU_CUDA_WORKER][0][0][0].size_base = trf_trsm_size;
  starpu_perfmodel_init(&XXTRF_model);
  XXTRF_TRSM_model.per_arch[STARPU_CPU_WORKER][0][0][0].size_base = trf_size;
  if(starpu_worker_get_count_by_type(STARPU_CUDA_WORKER) != 0)
    XXTRF_TRSM_model.per_arch[STARPU_CUDA_WORKER][0][0][0].size_base = trf_size;
  starpu_perfmodel_init(&TRSM_model);
  XXTRF_TRSM_model.per_arch[STARPU_CPU_WORKER][0][0][0].size_base = trsm_size;
  if(starpu_worker_get_count_by_type(STARPU_CUDA_WORKER) != 0)
    XXTRF_TRSM_model.per_arch[STARPU_CUDA_WORKER][0][0][0].size_base = trsm_size;
#  else
  XXTRF_TRSM_model.per_arch[STARPU_CPU_WORKER][0].size_base  = trf_trsm_size;
  XXTRF_TRSM_model.per_arch[STARPU_CUDA_WORKER][0].size_base = trf_trsm_size;
  XXTRF_model.per_arch[STARPU_CPU_WORKER][0].size_base  = trf_size;
  XXTRF_model.per_arch[STARPU_CUDA_WORKER][0].size_base = trf_size;
  TRSM_model.per_arch[STARPU_CPU_WORKER][0].size_base  = trsm_size;
  TRSM_model.per_arch[STARPU_CUDA_WORKER][0].size_base = trsm_size;
#  endif
#  ifdef STARPU_CONTEXT
  MALLOC_INTERN(devices, SOLV_THRDNBR+cuda_nbr, int);
  starpu_worker_get_ids_by_type(STARPU_CPU_WORKER, devices, SOLV_THRDNBR);

#    ifdef STARPU_USE_CUDA
  starpu_worker_get_ids_by_type(STARPU_CUDA_WORKER, devices+SOLV_THRDNBR, cuda_nbr);
#    endif
  /*create contexts however you want*/
  fprintf(stdout, "creating %d contexts \n", (int)sopalin_data->sopar->iparm[IPARM_STARPU_CTX_NBR]);
  thread_per_ctx = SOLV_THRDNBR/(sopalin_data->sopar->iparm[IPARM_STARPU_CTX_NBR]-1);
  if (SOLV_THRDNBR%(sopalin_data->sopar->iparm[IPARM_STARPU_CTX_NBR]-1))
    thread_per_ctx++;

  MALLOC_INTERN(sched_ctxs, sopalin_data->sopar->iparm[IPARM_STARPU_CTX_NBR], unsigned);
  sched_ctxs[0] = starpu_sched_ctx_create("dmda", devices+SOLV_THRDNBR, cuda_nbr, "ctx_0");
  for (iter = 1; iter < sopalin_data->sopar->iparm[IPARM_STARPU_CTX_NBR]; iter++)
    {
      char string[128];
      int nthreads = thread_per_ctx;
      if (iter == sopalin_data->sopar->iparm[IPARM_STARPU_CTX_NBR]-1 &&
          SOLV_THRDNBR%(thread_per_ctx) != 0)
        nthreads = SOLV_THRDNBR%(thread_per_ctx);
      sprintf(string, "ctx_%d", iter);
      fprintf(stdout, "creating %s contexts with %d cores %d\n", string, nthreads, thread_per_ctx);

      sched_ctxs[iter] = starpu_sched_ctx_create("dmda",
                                                 devices+(iter-1)*thread_per_ctx,
                                                 nthreads, string);
      starpu_sched_ctx_set_inheritor(sched_ctxs[iter], sched_ctxs[0]);
    }
#endif

#  if (defined PASTIX_WITH_STARPU_PROFILING)
  if ((ret = starpu_profiling_status_set(STARPU_PROFILING_ENABLE) < 0))
    {
      errorPrint("Error %d in starpu_profiling_enable\n", ret);
    }
#  endif /* (defined PASTIX_WITH_STARPU_PROFILING) */


#  ifdef STARPU_INIT_SMP
  {
    int threadid;
    sopthread_data_t * init_arg;
    MALLOC_INTERN(init_arg, SOLV_THRDNBR+cuda_nbr, sopthread_data_t);
    for (threadid = 0; threadid < SOLV_THRDNBR+cuda_nbr; threadid++)
      {
        struct starpu_task * task_init;
#    ifdef STARPU_CONTEXT
        pastix_int_t my_ctx;
        my_ctx = 1+threadid/thread_per_ctx;
#    endif
        task_init = starpu_task_create();
        /* We compute GEMM */
        init_arg[threadid].me   = threadid;
        init_arg[threadid].data = sopalin_data;
        task_init->cl = &sopalin_init_cl;
        task_init->cl_arg = &(init_arg[threadid]);
#    ifdef STARPU_PASTIX_SCHED
        task_init->workerid = threadid;
        task_init->priority = 0; /* No priority needed as the task needs
                                  * to be runned first due to dependancies */
#    endif

#    ifdef STARPU_CONTEXT
        ret = starpu_task_submit_to_ctx(task_init,
                                        sched_ctxs[my_ctx]);
        if (ret != -ENODEV) STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit_to_ctx")

#    else
        ret = starpu_task_submit(task_init);
        if (ret != -ENODEV) STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");

#    endif
        STARPU_ASSERT(!ret);
      }
    /* wait for end of init */
    starpu_task_wait_for_all();
    memFree_null(init_arg);
  }
#  endif
  {
    pastix_int_t itertask;
    for(itertask = 0; itertask < SOLV_TASKNBR; itertask++) {
      TASK_CTRBCNT(itertask) = 0;
    }
    for (itertask=0;itertask<SOLV_TASKNBR;itertask++) {
      pastix_int_t itercblk = TASK_CBLKNUM(itertask);
      pastix_int_t iterbloc;
      /* count updates from local cblk */
      for (iterbloc = SYMB_BLOKNUM(itercblk)+1;
           iterbloc < SYMB_BLOKNUM(itercblk+1);
           iterbloc ++) {
        pastix_int_t fcblknum = SYMB_CBLKNUM(iterbloc);
        if (fcblknum <  0) {
          /* itercblk udates a remote cblk */
          continue;
        } else {
          TASK_CTRBCNT(fcblknum)++;
        }
      }
    }
  }
  MALLOC_INTERN(blocktab_handles, SOLV_PROCNBR, starpu_data_handle_t);
  for (iter = 0; iter < SOLV_PROCNBR; iter++) {
    int size;
    if (SOLV_PROCNUM == iter) {
      /* build blocktab */
      int iterblock;
      MALLOC_INTERN(blocktab, (SYMB_BLOKNBR+SOLV_HBLOKNBR)*2, int);
      if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        fprintf(stdout, "sizeof blocktab : %d integers\n",
                (int)(2*(SYMB_BLOKNBR+SOLV_HBLOKNBR)));
      for (iterblock = 0; iterblock < SYMB_BLOKNBR; iterblock++) {
        blocktab[2*iterblock]   = SYMB_FROWNUM(iterblock);
        blocktab[2*iterblock+1] = SYMB_LROWNUM(iterblock);
      }
      for (iterblock = 0; iterblock < SOLV_HBLOKNBR; iterblock++) {
        blocktab[2*(iterblock+SYMB_BLOKNBR)]   = HBLOCK_FROWNUM(iterblock);
        blocktab[2*(iterblock+SYMB_BLOKNBR)+1] = HBLOCK_LROWNUM(iterblock);
      }
      size = 2*(SYMB_BLOKNBR + SOLV_HBLOKNBR);
      MPI_Bcast(&size, 1, MPI_INT, iter, sopalin_data->sopar->pastix_comm);
      starpu_vector_data_register(&(blocktab_handles[iter]), 0,
                                  (uintptr_t)blocktab, size,
                                  sizeof(int));
    } else {
      MPI_Bcast(&size, 1, MPI_INT, iter, sopalin_data->sopar->pastix_comm);
      starpu_vector_data_register(&(blocktab_handles[iter]), 0,
                                  (uintptr_t)NULL, size,
                                  sizeof(int));
    }
    starpu_mpi_data_register(blocktab_handles[iter], -20-iter, iter);

    starpu_data_set_sequential_consistency_flag(blocktab_handles[iter], 0);
  }

#  ifdef PASTIX_DUMP_FACTO
  dump_all(datacode, sopalin_data->sopar->cscmtx,
           ((datacode->updovct.sm2xtab!=NULL)?
            (DUMP_CSC | DUMP_SOLV | DUMP_SMB):(DUMP_CSC | DUMP_SOLV)));
#  endif

  thread_data = sopalin_data->thread_data[0];
  sopalin_data->sopar->diagchange = 0;
  SOPALIN_CLOCK_INIT;

  {
    int solv_coefmax, solv_coefmaxmax;
    solv_coefmax = SOLV_COEFMAX;
    MPI_Allreduce(&solv_coefmax, &solv_coefmaxmax, 1, MPI_INT, MPI_MAX, sopalin_data->sopar->pastix_comm);
#  ifdef CHOL_SOPALIN
    starpu_vector_data_register(&WORK_handle, -1, (uintptr_t)NULL, solv_coefmaxmax,
                                sizeof(pastix_float_t));
#  else
    fprintf(stdout, "2*solv_coefmaxmax = %d\n", 2*solv_coefmaxmax);
    starpu_vector_data_register(&WORK_handle, -1, (uintptr_t)NULL, 2*solv_coefmaxmax,
                              sizeof(pastix_float_t));
#  endif
  }
  starpu_mpi_data_register(WORK_handle, -3, SOLV_PROCNUM);

  MALLOC_INTERN(L_handle,     SYMB_CBLKNBR, starpu_data_handle_t);
  MALLOC_INTERN(Lhalo_handle, SOLV_HCBLKNBR, starpu_data_handle_t);
  memset(Lhalo_handle, 0, SOLV_HCBLKNBR*sizeof(starpu_data_handle_t));

#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
  MALLOC_INTERN(U_handle,    SYMB_CBLKNBR, starpu_data_handle_t);
  MALLOC_INTERN(Uhalo_handle, SOLV_HCBLKNBR, starpu_data_handle_t);
#    endif
#  endif


  {
    int itercblk;
    int max_cblksize = 0;
    int max_cblkcolnbr = 0;
    for (itercblk = 0; itercblk < SYMB_CBLKNBR; itercblk++)
      {
        max_cblksize   = MAX(max_cblksize,
                             CBLK_COLNBR(itercblk)*SOLV_STRIDE(itercblk));
        max_cblkcolnbr = MAX(max_cblkcolnbr,
                             CBLK_COLNBR(itercblk));
      }

    if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
      fprintf(stdout, "Maximum cblk size %d, maximu cblk colnbr %d\n",
              max_cblksize, max_cblkcolnbr);
  }

  task_number = 0;
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    task_number = SYMB_BLOKNBR;
  if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
    task_number += 2*SYMB_BLOKNBR+SYMB_CBLKNBR;
#  if (defined CHOL_SOPALIN && defined SOPALIN_LU)
  starpu_partition_data(sopalin_data, L_handle, U_handle, Lhalo_handle, Uhalo_handle);
#  else
  starpu_partition_data(sopalin_data, L_handle, NULL, Lhalo_handle, NULL);
#  endif

  if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
    {
      if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_SOLVE)
        if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
          errorPrintW("Raffinement not available with StarPU,"
                      " only performing solve\n");
      MALLOC_INTERN(SM2X_handles, SYMB_CBLKNBR, starpu_data_handle_t);
      starpu_register_sm2x(sopalin_data, SM2X_handles);
    }
  SOPALIN_CLOCK_STOP;
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    fprintf(stdout,"----- Time after data registration %lf s\n",
            SOPALIN_CLOCK_GET);

  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
    {
      starpu_loop_data_t * starpu_loop_data;
      pthread_cond_t  cond_end_facto;
      pthread_mutex_t mutex_end_facto;

      pthread_cond_init(&cond_end_facto, NULL);
      pthread_mutex_init(&mutex_end_facto, NULL);

      MALLOC_INTERN(starpu_loop_data, 1, starpu_loop_data_t);
      sopalin_data->starpu_loop_data = starpu_loop_data;
      starpu_loop_data->me               = 0;
      starpu_loop_data->L_handle         = L_handle;
      starpu_loop_data->Lhalo_handle     = Lhalo_handle;
#    if (defined CHOL_SOPALIN && defined SOPALIN_LU)
      starpu_loop_data->U_handle         = U_handle;
      starpu_loop_data->Uhalo_handle         = Uhalo_handle;
#    endif
      starpu_loop_data->WORK_handle      = WORK_handle;
      starpu_loop_data->blocktab_handles = blocktab_handles;
      starpu_loop_data->sopalin_data     = sopalin_data;
      starpu_loop_data->ctx_nbr          = 1;
      starpu_loop_data->first            = 0;
      starpu_loop_data->last             = SOLV_THRDNBR;
#    ifdef STARPU_CONTEXT
      starpu_loop_data->thread_per_ctx   = thread_per_ctx;
      starpu_loop_data->sched_ctxs       = sched_ctxs;
#    endif
      starpu_loop_data->cond_end_facto  = &cond_end_facto;
      starpu_loop_data->mutex_end_facto = &mutex_end_facto;
      MALLOC_INTERN(starpu_loop_data->cpu_workerids, SOLV_THRDNBR, int);
      starpu_loop_data->ncpus = SOLV_THRDNBR;
      MALLOC_INTERN(starpu_loop_data->gpu_workerids, cuda_nbr, int);
      MALLOC_INTERN(starpu_loop_data->gpu_gemm_count, cuda_nbr, int);
      memset(starpu_loop_data->gpu_gemm_count, 0, cuda_nbr*sizeof(int));
      starpu_loop_data->ngpus = cuda_nbr;
      starpu_worker_get_ids_by_type(STARPU_CPU_WORKER,
                                    starpu_loop_data->cpu_workerids,
                                    SOLV_THRDNBR);
      starpu_worker_get_ids_by_type(STARPU_CUDA_WORKER,
                                    starpu_loop_data->gpu_workerids,
                                    cuda_nbr);
      {
        int j;
        for (j = 0; j < cuda_nbr; j++)
          fprintf(stdout, "cuda_id %d\n", starpu_loop_data->gpu_workerids[j]);
      }
      MALLOC_INTERN(starpu_loop_data->facto_finished, 1, int);
      *(starpu_loop_data->facto_finished) = API_NO;

      MALLOC_INTERN(sopalin_data->gemm_stats,
                    starpu_worker_get_count(),
                    starpu_task_stats_t);
      memset(sopalin_data->gemm_stats, 0,
             starpu_worker_get_count()*sizeof(starpu_task_stats_t));
      MALLOC_INTERN(sopalin_data->hgemm_stats,
                    starpu_worker_get_count(),
                    starpu_task_stats_t);
      memset(sopalin_data->hgemm_stats, 0,
             starpu_worker_get_count()*sizeof(starpu_task_stats_t));
      MALLOC_INTERN(sopalin_data->xxtrf_stats,
                    starpu_worker_get_count(),
                    starpu_task_stats_t);
      memset(sopalin_data->xxtrf_stats, 0,
             starpu_worker_get_count()*sizeof(starpu_task_stats_t));
      MALLOC_INTERN(sopalin_data->trsm_stats,
                    starpu_worker_get_count(),
                    starpu_task_stats_t);
      memset(sopalin_data->trsm_stats, 0,
             starpu_worker_get_count()*sizeof(starpu_task_stats_t));

      halo_submit(sopalin_data);
      starpu_submit_loop (starpu_loop_data);
    }
  else
    sopalin_data->starpu_loop_data = NULL;

  if (sopalin_data->starpu_loop_data != NULL) {
    /* wait for last task to be submitted */
    char * nested;
    if ((nested = getenv("PASTIX_STARPU_NESTED_TASK")) && !strcmp(nested, "1")) {
      while(SYMB_BLOKNBR > 0 &&
	    *(sopalin_data->starpu_loop_data->facto_finished) == API_NO) {
	COND_WAIT(sopalin_data->starpu_loop_data->cond_end_facto,
		  sopalin_data->starpu_loop_data->mutex_end_facto);
      }
    }
    if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
      starpu_submit_updown(sopalin_data, L_handle, U_handle, SM2X_handles,
			   NULL, (int*)sched_ctxs);

    memFree_null(sopalin_data->starpu_loop_data->facto_finished);
    pthread_mutex_unlock(sopalin_data->starpu_loop_data->mutex_end_facto);

    pthread_cond_destroy(sopalin_data->starpu_loop_data->cond_end_facto);
    pthread_mutex_destroy(sopalin_data->starpu_loop_data->mutex_end_facto);
  }

  /* wait for all task to complet */
  starpu_task_wait_for_all();
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT) {
    int i;
    int nworkers = 1;
    nworkers = SOLV_THRDNBR + sopalin_data->sopar->iparm[IPARM_CUDA_NBR];;
    memFree_null(sopalin_data->starpu_loop_data->cpu_workerids);
    memFree_null(sopalin_data->starpu_loop_data->gpu_workerids);
    for (i = 0; i < sopalin_data->sopar->iparm[IPARM_CUDA_NBR]; i++)
      fprintf(stdout, "%d GEMMs forced on GPU %d\n",
              sopalin_data->starpu_loop_data->gpu_gemm_count[i], i);
    memFree_null(sopalin_data->starpu_loop_data->gpu_gemm_count);
    memFree_null(sopalin_data->starpu_loop_data);
  }
  SOPALIN_CLOCK_STOP;
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    fprintf(stdout,"----- submission and wait for all %lf s (%ld tasks)\n",
            SOPALIN_CLOCK_GET, (long)SYMB_BLOKNBR);


  /* Unregister buffers and leave starpu */
  if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT) {
#  if (defined PASTIX_WITH_STARPU_PROFILING)
    if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO) {
      int worker;
      /* Display the occupancy of all workers during the test */
      for (worker = 0; worker < starpu_worker_get_count(); worker++) {
        struct starpu_profiling_worker_info worker_info;
        int ret = starpu_profiling_worker_get_info(worker, &worker_info);
        STARPU_ASSERT(!ret);

        double total_time     = starpu_timing_timespec_to_us(&worker_info.total_time);
        double executing_time = starpu_timing_timespec_to_us(&worker_info.executing_time);
        double sleeping_time  = starpu_timing_timespec_to_us(&worker_info.sleeping_time);
        double overhead_time = total_time - executing_time - sleeping_time;

        float executing_ratio = 100.0*executing_time/total_time;
        float sleeping_ratio  = 100.0*sleeping_time/total_time;
        float overhead_ratio = 100.0 - executing_ratio - sleeping_ratio;

        char workername[128];


        starpu_worker_get_name(worker, workername, 128);
        fprintf(stdout, "Worker %s:\n", workername);
        fprintf(stdout, "\ttotal time : %.2lf ms\n", total_time*1e-3);
        fprintf(stdout, "\texec time  : %.2lf ms (%.2f %%)\n", executing_time*1e-3, executing_ratio);
        fprintf(stdout, "\tblocked time  : %.2lf ms (%.2f %%)\n", sleeping_time*1e-3, sleeping_ratio);
        fprintf(stdout, "\toverhead time: %.2lf ms (%.2f %%)\n", overhead_time*1e-3, overhead_ratio);
        if (sopalin_data->xxtrf_stats[worker].cnt != 0) {
          fprintf(stdout, "\tAvg. delay on XXTRF : %2.2lf us, %d tasks\n",
                  sopalin_data->xxtrf_stats[worker].delay_sum/
                  sopalin_data->xxtrf_stats[worker].cnt,
                  sopalin_data->xxtrf_stats[worker].cnt);
          fprintf(stdout, "\tAvg. length on XXTRF : %2.2lf us, %5g %s\n",
                  sopalin_data->xxtrf_stats[worker].length_sum/
                  sopalin_data->xxtrf_stats[worker].cnt,
                  PRINT_FLOPS(sopalin_data->xxtrf_stats[worker].ops/
                              sopalin_data->xxtrf_stats[worker].length_sum),
                  PRINT_FLOPS_UNIT(sopalin_data->xxtrf_stats[worker].ops/
                                   sopalin_data->xxtrf_stats[worker].length_sum));
        }
        if (sopalin_data->trsm_stats[worker].cnt != 0) {
          fprintf(stdout, "\tAvg. delay on TRSM : %2.2lf us, %d tasks\n",
                  sopalin_data->trsm_stats[worker].delay_sum/
                  sopalin_data->trsm_stats[worker].cnt,
                  sopalin_data->trsm_stats[worker].cnt);
          fprintf(stdout, "\tAvg. length on TRSM : %2.2lf us, %5g %s\n",
                  sopalin_data->trsm_stats[worker].length_sum/
                  sopalin_data->trsm_stats[worker].cnt,
                  PRINT_FLOPS(sopalin_data->trsm_stats[worker].ops/
                              (sopalin_data->trsm_stats[worker].length_sum*10e-6)),
                  PRINT_FLOPS_UNIT(sopalin_data->trsm_stats[worker].ops/
                                   (sopalin_data->trsm_stats[worker].length_sum*
                                    10e-6)));
        }
        if (sopalin_data->gemm_stats[worker].cnt != 0) {
          fprintf(stdout, "\tAvg. delay on GEMM : %2.2lf us, %d tasks\n",
                  sopalin_data->gemm_stats[worker].delay_sum/
                  sopalin_data->gemm_stats[worker].cnt,
                  sopalin_data->gemm_stats[worker].cnt);
          fprintf(stdout, "\tAvg. length on GEMM : %2.2lf us, %5g %s\n",
                  sopalin_data->gemm_stats[worker].length_sum/
                  sopalin_data->gemm_stats[worker].cnt,
                  PRINT_FLOPS(sopalin_data->gemm_stats[worker].ops/
                              (sopalin_data->gemm_stats[worker].length_sum*10e-6)),
                  PRINT_FLOPS_UNIT(sopalin_data->gemm_stats[worker].ops/
                                   (sopalin_data->gemm_stats[worker].length_sum*
                                    10e-6)));
        }
        if (sopalin_data->hgemm_stats[worker].cnt != 0) {
          fprintf(stdout, "\tAvg. delay on Halo GEMM : %2.2lf us, %d tasks\n",
                  sopalin_data->hgemm_stats[worker].delay_sum/
                  sopalin_data->hgemm_stats[worker].cnt,
                  sopalin_data->hgemm_stats[worker].cnt);
          fprintf(stdout, "\tAvg. length on Halo GEMM : %2.2lf us, %5g %s\n",
                  sopalin_data->hgemm_stats[worker].length_sum/
                  sopalin_data->hgemm_stats[worker].cnt,
                  PRINT_FLOPS(sopalin_data->hgemm_stats[worker].ops/
                              (sopalin_data->hgemm_stats[worker].length_sum*10e-6)),
                  PRINT_FLOPS_UNIT(sopalin_data->hgemm_stats[worker].ops/
                                   (sopalin_data->hgemm_stats[worker].length_sum*
                                    10e-6)));
        }

      }
    }
#  endif /* (defined PASTIX_WITH_STARPU_PROFILING) */
    memFree_null(sopalin_data->gemm_stats);
    memFree_null(sopalin_data->hgemm_stats);
    memFree_null(sopalin_data->xxtrf_stats);
    memFree_null(sopalin_data->trsm_stats);

  }
  for (itertask=0;itertask<SOLV_TASKNBR;itertask++) {
    pastix_int_t itercblk = TASK_CBLKNUM(itertask);
    pastix_int_t iterbloc;
    for (iterbloc = SYMB_BLOKNUM(itercblk);
         iterbloc < SYMB_BLOKNUM(itercblk+1);
         iterbloc ++) {
      pastix_int_t first_task = 0;
      if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT) {
	first_task = SYMB_BLOKNBR;
      }

      /* TODO: old UPDO code require updates */
      if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT) {
        first_task += SYMB_BLOKNBR;
#  ifndef CHOL_SOPALIN
        first_task += SYMB_CBLKNBR;
#  endif /* not CHOL_SOPALIN */
      }
    }
#  ifndef CHOL_SOPALIN
    {
      if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT) {
        pastix_int_t first_task = SYMB_BLOKNBR;
        if (sopalin_data->sopar->iparm[IPARM_START_TASK] <= API_TASK_NUMFACT)
          first_task += SYMB_BLOKNBR;
      }
    }
#  endif /* CHOL_SOPALIN */
    starpu_data_unregister(L_handle[itercblk]);
    if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT) {
      starpu_data_unregister(SM2X_handles[itercblk]);
    }
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
    starpu_data_unregister(U_handle[itercblk]);
#    endif
#  endif
  }
  for (itertask=0;itertask<SOLV_HCBLKNBR;itertask++) {
    starpu_data_unregister(Lhalo_handle[itertask]);
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
    starpu_data_unregister(Uhalo_handle[itertask]);
#    endif
#  endif
  }
  starpu_data_unregister(WORK_handle);

  for (iter = 0; iter < SOLV_PROCNBR; iter++)
    starpu_data_unregister(blocktab_handles[iter]);
  memFree_null(blocktab_handles);
  memFree_null(blocktab);

  /* Reduction on pivot number */
  sopalin_data->sopar->diagchange = 0;
  {
    pastix_int_t me;
    for (me = 0; me < SOLV_THRDNBR+cuda_nbr; me++)
      {
        sopalin_data->sopar->diagchange += sopalin_data->thread_data[me]->nbpivot;
      }
  }

  SOPALIN_CLOCK_STOP;
  if (sopalin_data->sopar->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    fprintf(stdout,"----- sopalin time %lf\n",
            SOPALIN_CLOCK_GET);
  sopalin_data->sopar->dparm[DPARM_FACT_TIME] = SOPALIN_CLOCK_GET;
#  ifdef PASTIX_DUMP_FACTO
  dump_all(datacode, sopalin_data->sopar->cscmtx,
           ((datacode->updovct.sm2xtab!=NULL)?
            (DUMP_CSC | DUMP_SOLV | DUMP_SMB):(DUMP_CSC | DUMP_SOLV)));
#  endif
  memFree_null(L_handle);
  memFree_null(Lhalo_handle);
  if (sopalin_data->sopar->iparm[IPARM_END_TASK] > API_TASK_NUMFACT)
    memFree_null(SM2X_handles);
#  ifdef CHOL_SOPALIN
#    ifdef SOPALIN_LU
  memFree_null(U_handle);
  memFree_null(Uhalo_handle);
#    endif /* CHOL_SOPALIN */
#  endif /* SOPALIN_LU   */
#ifdef STARPU_USE_CUDA
  if (cuda_nbr) {
    starpu_cublas_shutdown();
  }
#endif /* STARPU_USE_CUDA */

  starpu_mpi_shutdown();
  starpu_shutdown();
  sopalin_clean(sopalin_data, 1);

  sopalin_launch_thread(sopalin_data,
                        SOLV_PROCNUM, SOLV_PROCNBR, datacode->btree,
                        sopalin_data->sopar->iparm[IPARM_VERBOSE],
                        SOLV_THRDNBR, starpu_clean_smp, sopalin_data,
                        0, NULL, NULL,
                        0, NULL, NULL);
  return NO_ERR;
}
