/**
 *
 * @file eztrace_convert_kernels.h
 *
 * Module to convert eztrace events
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Gregoire Pichon
 * @date 2017-04-26
 *
 * @addtogroup eztrace_dev
 * @{
 *
 **/

#ifndef __EZTRACE_CONVERT_KERNELS_H__
#define __EZTRACE_CONVERT_KERNELS_H__

#include <stdio.h>
#include <strings.h>
#include <GTG.h>
#include "eztrace_convert.h"
#include "kernels_ev_codes.h"

typedef struct kernels_e {
    char *name;
    gtg_color_t color;
} kernels_t;

typedef struct kernels_thread_info_e {
    struct thread_info_t *p_thread;

    /* Use by KERNELS_STOP */
    float             time_start;
    kernels_ev_code_t current_ev;

    /* Counters per event */
    int    *nb;
    double *flops;
    double *run_time;
} kernels_thread_info_t;

#define INIT_KERNELS_THREAD_INFO(p_thread, var, stats)           \
    kernels_thread_info_t *var = (kernels_thread_info_t *) \
        ezt_hook_list_retrieve_data(&p_thread->hooks, KERNELS_EVENTS_ID); \
    if(!(var)) {                                                        \
        var = kernels_register_thread_hook(p_thread, stats);                  \
    }

void define_kernels_properties();
void handle_start(kernels_ev_code_t ev, int stats);
void handle_stop(int stats);

int eztrace_convert_kernels_init();
int handle_kernels_events(eztrace_event_t *ev);
int handle_kernels_stats(eztrace_event_t *ev);
void print_kernels_stats();

#endif /* __EZTRACE_CONVERT_KERNELS_H__ */
