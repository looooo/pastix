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
struct _thread_info_t {
  struct thread_info_t *p_thread;

  /* TO COMPLETE: You can add per-thread counters here */
};

static struct _thread_info_t *_kernels_register_thread_hook(
    struct thread_info_t *p_thread) {
  struct _thread_info_t *ptr = (struct _thread_info_t*) malloc(
      sizeof(struct _thread_info_t));

  ptr->p_thread = p_thread;

  /* TO COMPLETE: If you added per-thread counters, initialize them here*/

  /* add the hook in the thread info structure */
  ezt_hook_list_add(&ptr->p_thread->hooks, ptr,
                   (uint8_t) KERNELS_EVENTS_ID);
  return ptr;
}

#define  INIT_KERNELS_THREAD_INFO(p_thread, var)                      \
  struct _thread_info_t *var = (struct _thread_info_t*)        \
    ezt_hook_list_retrieve_data(&p_thread->hooks, (uint8_t)KERNELS_EVENTS_ID); \
  if(!(var)) {                                                         \
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

  //printf("module  loaded\n");
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
      printf("START\n\n");
      break;
  case KERNELS_LRALLOC_STOP:
      printf("STOP\n\n");
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
      INIT_KERNELS_THREAD_INFO(p_thread, ptr);
      printf("\tThread %s\n", thread_container->name);

      /* TO COMPLETE: you can print per-thread counters here */
    }
  }

}


