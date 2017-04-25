/*
 * Copyright (C) CNRS, INRIA, Université Bordeaux 1, Télécom SudParis
 * See COPYING in top-level directory.
 */

#include <stdio.h>
#include <strings.h>
#include <GTG.h>
#include "eztrace_convert.h"

#include "kernels_ev_codes.h"

int eztrace_convert_kernels_init();
int handle_kernels_events(eztrace_event_t *ev);
int handle_kernels_stats(eztrace_event_t *ev);
void print_kernels_stats();



/* thread-specific structure */
struct kernels_thread_info_t {
  struct thread_info_t *p_thread;

    int nb_calls;
    int size;
    double time;

    float time_start;
  /* TO COMPLETE: You can add per-thread counters here */
};

static struct kernels_thread_info_t *_kernels_register_thread_hook(
    struct thread_info_t *p_thread) {
  struct kernels_thread_info_t *p_info = (struct kernels_thread_info_t*) malloc(
      sizeof(struct kernels_thread_info_t));

  p_info->p_thread = p_thread;

  /* TO COMPLETE: If you added per-thread counters, initialize them here*/
  p_info->nb_calls = 0;
  p_info->size     = 0;

  /* add the hook in the thread info structure */
  ezt_hook_list_add(&p_info->p_thread->hooks, p_info,
                   (uint8_t) KERNELS_EVENTS_ID);
  return p_info;
}

#define  INIT_KERNELS_THREAD_INFO(p_thread, var)                      \
  struct kernels_thread_info_t *var = (struct kernels_thread_info_t*)        \
    ezt_hook_list_retrieve_data(&p_thread->hooks, (uint8_t)KERNELS_EVENTS_ID); \
  if(!(var)) {                                                        \
    var = _kernels_register_thread_hook(p_thread);                            \
  }


/* Constructor of the plugin.
 * This function registers the current module to eztrace_convert
 */
struct eztrace_convert_module kernels_module;
void libinit(void) __attribute__ ((constructor));
void libinit(void)
{
  kernels_module.api_version = EZTRACE_API_VERSION;

  /* Specify the initialization function.
   * This function will be called once all the plugins are loaded
   * and the trace is started.
   * This function usually declared StateTypes, LinkTypes, etc.
   */
  kernels_module.init = eztrace_convert_kernels_init;

  /* Specify the function to call for handling an event
   */
  kernels_module.handle = handle_kernels_events;

  /* Specify the function to call for handling an event when
   * eztrace_stats is called
   */
  kernels_module.handle_stats = handle_kernels_stats;

  /* Specify the function to call for printinf statistics
   */
  kernels_module.print_stats = print_kernels_stats;

  /* Specify the module prefix */
  kernels_module.module_prefix = KERNELS_EVENTS_ID;

  asprintf(&kernels_module.name, "kernels");
  asprintf(&kernels_module.description, "PaStiX kernels");

  kernels_module.token.data = &kernels_module;

  /* Register the module to eztrace_convert */
  eztrace_convert_register_module(&kernels_module);
}

void libfinalize(void) __attribute__ ((destructor));
void libfinalize(void)
{
  printf("unloading module \n");
}



/*
 * This function will be called once all the plugins are loaded
 * and the trace is started.
 * This function usually declared StateTypes, LinkTypes, etc.
 */
int
eztrace_convert_kernels_init()
{
  if (get_mode() == EZTRACE_CONVERT) {
      addEntityValue("STV_lralloc", "ST_Thread", "LRALLOC", GTG_GREEN);
  }
  return 0;
}

void handle_lralloc_start(eztrace_event_t *ev)
{
  DECLARE_THREAD_ID_STR(thread_id, CUR_INDEX, CUR_THREAD_ID);
  DECLARE_CUR_THREAD(p_thread);
  INIT_KERNELS_THREAD_INFO(p_thread, p_info);

  int size;
  GET_PARAM_PACKED_1(CUR_EV, size);

  p_info->nb_calls++;
  p_info->size += size;

  p_info->time_start = CURRENT;

  pushState(CURRENT, "ST_Thread", thread_id, "STV_lralloc");
}

void handle_lralloc_stop(eztrace_event_t *ev)
{
  DECLARE_THREAD_ID_STR(thread_id, CUR_INDEX, CUR_THREAD_ID);
  DECLARE_CUR_THREAD(p_thread);
  INIT_KERNELS_THREAD_INFO(p_thread, p_info);

  p_info->time += (CURRENT - p_info->time_start);

  pushState(CURRENT, "ST_Thread", thread_id, "STV_Working");
}

/* This function is called by eztrace_convert for each event to
 * handle.
 * It shall return 1 if the event was handled successfully or
 * 0 otherwise.
 */
int
handle_kernels_events(eztrace_event_t *ev)
{

  if(! CUR_TRACE->start)
    return 0;

  switch (LITL_READ_GET_CODE(ev)) {

  case KERNELS_LRALLOC_START:
      handle_lralloc_start(ev);
      break;
  case KERNELS_LRALLOC_STOP:
      handle_lralloc_stop(ev);
      break;
    default:
      /* The event was not handled */
      return 0;
    }
  return 1;
}

/* This function is called by eztrace_stats for each event to
 * handle.
 * It shall return 1 if the event was handled successfully or
 * 0 otherwise.
 */
int
handle_kernels_stats(eztrace_event_t *ev)
{
  /* By default, use the same function as for eztrace_convert */
  return handle_kernels_events(ev);
}


void
print_kernels_stats()
{
  printf("\n:\n");
  printf("-------\n");

  int i;

  int sum = 0;
  /* Browse the list of processes */
  for (i = 0; i < NB_TRACES; i++) {
    struct eztrace_container_t *p_process = GET_PROCESS_CONTAINER(i);
    int j;
    /* For each process, browse the list of threads */
    for(j=0; j<p_process->nb_children; j++) {
      struct eztrace_container_t *thread_container = p_process->children[j];
      struct thread_info_t *p_thread = (struct thread_info_t*) thread_container->container_info;
      if(!p_thread)
       continue;
      INIT_KERNELS_THREAD_INFO(p_thread, p_info);
      printf("\tThread %20s NB CALLS %5d SIZE %8d TIME %.3g\n", thread_container->name, p_info->nb_calls, p_info->size, p_info->time);

      sum+=p_info->nb_calls;
      /* TO COMPLETE: you can print per-thread counters here */
    }
  }
  printf("TOTAL NUMBER OF CALLS %d\n", sum);
}


