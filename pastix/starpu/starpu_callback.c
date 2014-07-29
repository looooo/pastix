/**
 * @file starpu_callback.c
 *
 * @author Xavier Lacoste
 */
#include "common.h"
#include "starpu_zdefines.h"
void starpu_prof_callback(void *callback_arg)
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
