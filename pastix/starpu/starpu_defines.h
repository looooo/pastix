#ifndef STARPU_DEFINES_H
#define STARPU_DEFINES_H

#include <starpu.h>
#include "common.h"
#ifdef PASTIX_WITH_MPI
#  include <starpu_mpi.h>
#endif
#include <starpu_profiling.h>

#include "sopalin3d.h"

#if (STARPU_MAJOR_VERSION < 1)
#  error "PaStiX requires STARPU >= 1"
#endif /* STARPU_MAJOR_VERSION */

#ifndef PASTIX_WITH_MPI
#  define starpu_mpi_init(a, ...)            starpu_init(a)
#  define starpu_mpi_insert_task(a, b, ...)  starpu_insert_task(b, __VA_ARGS__)
#  define starpu_mpi_data_register(...)
#  define starpu_mpi_shutdown()
#  define starpu_data_get_rank(...)          0
#endif

#if (STARPU_MAJOR_VERSION == 1 &&  STARPU_MINOR_VERSION == 0)
/* 1.1 : starpu_perf_archtype => starpu_perfmodel_archtype
 *       starpu_profiling_worker_info => starpu_profiling_worker_info
 */
#  define STARPU_1_0
#  define starpu_perfmodel_archtype        starpu_perf_archtype
#  define starpu_profiling_worker_info     starpu_worker_profiling_info
#  define starpu_profiling_task_info       starpu_task_profiling_info
#  define starpu_profiling_worker_get_info starpu_worker_get_profiling_info
#  define starpu_profiling_set_id          starpu_set_profiling_id
#endif

#if (STARPU_MAJOR_VERSION == 1 &&  STARPU_MINOR_VERSION == 1)
#  define STARPU_1_0
#  define STARPU_1_1
#endif

#if (STARPU_MAJOR_VERSION == 1 &&  STARPU_MINOR_VERSION == 2)
#  define STARPU_1_0
#  define STARPU_1_1
#  define STARPU_1_2
#  define STARPU_COMMUTABLE
#endif

#ifndef STARPU_COMMUTABLE
#  define STARPU_COMMUTE 0
#endif

#define ARCH_CPU  0
#define ARCH_CUDA 1

#define SUBMIT_TRF_IF_NEEDED                                            \
    do {                                                                \
        char * nested;                                                  \
        SolverMatrix *datacode = sopalin_data->datacode;                \
        if (cblk_islocal(datacode, fcblk)) {                            \
            pastix_int_t fcblknum = cblk_getnum(datacode, fcblk);       \
            TASK_CTRBCNT(fcblknum)--;                                   \
            if (pastix_starpu_with_nested_task() == API_YES) {          \
                if (TASK_CTRBCNT(fcblknum) == 0) {                      \
                    starpu_submit_one_trf(fcblknum, sopalin_data);      \
                }                                                       \
            } else {                                                    \
                if (pastix_starpu_with_fanin() == API_YES) {            \
                    pastix_int_t faninnum = fcblk_getnum(datacode,      \
                                                         fcblk,         \
                                                         SOLV_PROCNUM); \
                    pastix_int_t *fanin_ctrbcnt;                        \
                    fanin_ctrbcnt = sopalin_data->fanin_ctrbcnt;        \
                    fanin_ctrbcnt[faninnum]--;                          \
                    if (fanin_ctrbcnt[faninnum] == 0) {                 \
                        pastix_int_t hcblknum;                          \
                        SolverCblk * hcblk;                             \
                        hcblknum = SOLV_GCBLK2HALO(fcblk->gcblknum);    \
                        hcblk = datacode->hcblktab + hcblknum;          \
                        starpu_submit_outgoing_fanin( sopalin_data,     \
                                                      fcblk,            \
                                                      hcblk);           \
                    }                                                   \
                }                                                       \
            }                                                           \
        }                                                               \
    } while(0)


#define SUBMIT_GEMMS_IF_NEEDED                                          \
        do {                                                            \
            char * nested;                                              \
            SolverMatrix *datacode = sopalin_data->datacode;            \
            pastix_int_t tasknum = cblk_getnum(datacode, cblk);         \
            if (pastix_starpu_with_nested_task() == API_YES) {          \
                starpu_submit_bunch_of_gemm(tasknum, sopalin_data);     \
            }                                                           \
        } while(0)


struct starpu_loop_data_ {
    int                     me;
    starpu_data_handle_t *  L_handle;
    starpu_data_handle_t *  Lhalo_handle;
    starpu_data_handle_t *  U_handle;
    starpu_data_handle_t *  Uhalo_handle;
    starpu_data_handle_t ** Lfanin_handle;
    starpu_data_handle_t ** Ufanin_handle;
    starpu_data_handle_t    WORK_handle;
    starpu_data_handle_t *  blocktab_handles;
    Sopalin_Data_t       *  sopalin_data;
    int                     ctx;
    int                     ctx_nbr;
    int                     thread_per_ctx;
    int                  *  sched_ctxs;
    int                     first, last;
    int                   * facto_finished;
    pthread_cond_t        * cond_end_facto;
    pthread_mutex_t       * mutex_end_facto;
    int                   * cpu_workerids;
    int                     ncpus;
    int                   * gpu_workerids;
    int                     ngpus;
    int                   * gpu_gemm_count;
};


#endif /* STARPU_DEFINES_H */
