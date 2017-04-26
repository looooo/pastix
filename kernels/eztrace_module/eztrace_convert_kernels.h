/*
 * Copyright (C) CNRS, INRIA, Université Bordeaux 1, Télécom SudParis
 * See COPYING in top-level directory.
 */

#ifndef __EZTRACE_CONVERT_KERNELS_H__
#define __EZTRACE_CONVERT_KERNELS_H__

#include <stdio.h>
#include <strings.h>
#include <GTG.h>
#include "eztrace_convert.h"
#include "kernels_ev_codes.h"

/* thread-specific structure */
struct kernels_thread_info_t {
    struct thread_info_t *p_thread;

    int nb_calls;
    int size;
    double time;

    float time_start;
    /* TO COMPLETE: You can add per-thread counters here */
};

#define INIT_KERNELS_THREAD_INFO(p_thread, var)                         \
    struct kernels_thread_info_t *var = (struct kernels_thread_info_t*) \
        ezt_hook_list_retrieve_data(&p_thread->hooks, (uint8_t)KERNELS_EVENTS_ID); \
    if(!(var)) {                                                        \
        var = _kernels_register_thread_hook(p_thread);                  \
    }

int eztrace_convert_kernels_init();
int handle_kernels_events(eztrace_event_t *ev);
int handle_kernels_stats(eztrace_event_t *ev);
void print_kernels_stats();

void handle_lralloc_start();
void handle_lralloc_stop();

#endif /* __EZTRACE_CONVERT_KERNELS_H__ */
