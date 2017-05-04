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
 **/

#ifndef __EZTRACE_CONVERT_KERNELS_H__
#define __EZTRACE_CONVERT_KERNELS_H__

#include <stdio.h>
#include <strings.h>
#include <GTG.h>
#include "eztrace_convert.h"
#include "kernels_ev_codes.h"

#define MAX_EVENTS 20

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

/* Properties (name/color) of each event */
static kernels_t kernels_properties[MAX_EVENTS];

#define INIT_KERNELS_THREAD_INFO(p_thread, var)                         \
    kernels_thread_info_t *var = (kernels_thread_info_t *) \
        ezt_hook_list_retrieve_data(&p_thread->hooks, (uint8_t)KERNELS_EVENTS_ID); \
    if(!(var)) {                                                        \
        var = kernels_register_thread_hook(p_thread);                  \
    }

void define_kernels_properties();
void handle_start(kernels_ev_code_t ev);
void handle_stop();

int eztrace_convert_kernels_init();
int handle_kernels_events(eztrace_event_t *ev);
void print_kernels_stats();

#endif /* __EZTRACE_CONVERT_KERNELS_H__ */
