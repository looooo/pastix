/**
 *
 * @file context.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_CONTEXT_H_
#define _PLASMA_CONTEXT_H_

#include <stdio.h>
#include "quark.h"

struct plasma_context_struct;

/***************************************************************************//**
 *  PLASMA context
 **/
typedef struct plasma_context_struct {
    /* Initialization flag */
    PLASMA_bool initialized;

    /* Thread control */
    int world_size, group_size;
    int thread_bind[CONTEXT_THREADS_MAX];
    int thread_rank[CONTEXT_THREADS_MAX];
    pthread_attr_t thread_attr;
    pthread_t thread_id[CONTEXT_THREADS_MAX];

    /* Master-worker communication */
    pthread_mutex_t action_mutex;
    pthread_cond_t action_condt;
    volatile PLASMA_enum action;
    void (*parallel_func_ptr)(struct plasma_context_struct*);
    unsigned char args_buff[ARGS_BUFF_SIZE];

    /* Boolean flags */
    PLASMA_bool errors_enabled;
    PLASMA_bool warnings_enabled;
    PLASMA_bool autotuning_enabled;
    PLASMA_bool dynamic_section;

    /* Enum flags */
    PLASMA_enum scheduling;     // static or dynamic scheduling
    PLASMA_enum householder;    // "domino" (flat) or tree-based (reduction) Householder
    PLASMA_enum translation;    // In place or Out of place layout conversion

    /* Matrix tile attributes */
    int nb;
    int ib;
    int nbnbsize;   // tile size in elements (possibly padded)
    int ibnbsize;   // T or L tile size in elements (---''---)
    int rhblock;    // block size for tree-based (reduction) Householder

    /* Barrier implementation */
    /* Busy waiting version */
    volatile int barrier_in[CONTEXT_THREADS_MAX];
    volatile int barrier_out[CONTEXT_THREADS_MAX];

    /* Conditional version */
    int volatile    barrier_id;             /*+ ID of the barrier                     +*/
    int volatile    barrier_nblocked_thrds; /*+ Number of threads lock in the barrier +*/
    pthread_mutex_t barrier_synclock;       /*+ mutex for the barrier                 +*/
    pthread_cond_t  barrier_synccond;       /*+ condition for the barrier             +*/

    /* Static scheduler implementation */
    int ss_ld;                  // static scheduler progress table leading dimension
    volatile int ss_abort;      // static scheduler abort flag
    volatile int *ss_progress;  // static scheduler progress table

    /* Dynamic scheduler */
    Quark *quark;
} plasma_context_t;

/***************************************************************************//**
 *  Threads contexts map
 **/
typedef struct plasma_context_map_struct {
    pthread_t thread_id;        // thread id
    plasma_context_t *context;  // pointer to associated context
} plasma_context_map_t;

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Internal routines
 **/
plasma_context_t *plasma_context_create();
int plasma_context_insert(plasma_context_t *context, pthread_t thread_id);
int plasma_context_remove(plasma_context_t *context, pthread_t thread_id);
plasma_context_t *plasma_context_self();

#ifdef __cplusplus
}
#endif

#endif
