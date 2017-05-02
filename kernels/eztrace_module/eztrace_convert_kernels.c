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

#include "eztrace_convert_kernels.h"

static kernels_thread_info_t *kernels_register_thread_hook(
    struct thread_info_t *p_thread) {
    kernels_thread_info_t *p_info = (kernels_thread_info_t*) malloc(
        sizeof(kernels_thread_info_t));

    p_info->p_thread = p_thread;

    int    *nb       = malloc(MAX_EVENTS * sizeof(int));
    int    *flops    = malloc(MAX_EVENTS * sizeof(int));
    double *run_time = malloc(MAX_EVENTS * sizeof(double));

    int i;
    for (i=0; i<MAX_EVENTS; i++){
        nb[i]        = 0;
        flops[i]     = 0;
        run_time[i]  = 0;
    }

    p_info->nb       = nb;
    p_info->flops    = flops;
    p_info->run_time = run_time;

    ezt_hook_list_add(&p_info->p_thread->hooks, p_info,
                      (uint8_t) KERNELS_EVENTS_ID);
    return p_info;
}

struct eztrace_convert_module kernels_module;
void libinit(void) __attribute__ ((constructor));
void libinit(void)
{
    kernels_module.api_version = EZTRACE_API_VERSION;
    kernels_module.init = eztrace_convert_kernels_init;
    kernels_module.handle = handle_kernels_events;
    kernels_module.handle_stats = handle_kernels_stats;
    kernels_module.print_stats = print_kernels_stats;
    kernels_module.module_prefix = KERNELS_EVENTS_ID;

    asprintf(&kernels_module.name, "kernels");
    asprintf(&kernels_module.description, "PaStiX kernels");

    kernels_module.token.data = &kernels_module;
    eztrace_convert_register_module(&kernels_module);
}

void libfinalize(void) __attribute__ ((destructor));
void libfinalize(void)
{
    printf("unloading module \n");
}


void define_kernels_properties()
{
    /* Low-rank operations */
    kernels_properties[LR_ALLOC] = (kernels_t) {"lr_alloc", GTG_YELLOW};
    kernels_properties[LR_TRSM]  = (kernels_t) {"lr_trsm" , GTG_DARKBLUE};
    kernels_properties[LR_GEMM]  = (kernels_t) {"lr_gemm" , GTG_LIGHTPINK};

    /* Dense operations */
    kernels_properties[GETRF]      = (kernels_t) {"getrf"  ,    GTG_RED};
    kernels_properties[POTRF]      = (kernels_t) {"potrf"  ,    GTG_RED};
    kernels_properties[DENSE_TRSM] = (kernels_t) {"dense_trsm", GTG_SEABLUE};
    kernels_properties[DENSE_GEMM] = (kernels_t) {"dense_gemm", GTG_GREEN};
}

int
eztrace_convert_kernels_init()
{
    if (NB_EVENTS > MAX_EVENTS + 1){
        printf("FATAL ERROR: static table is not large enough\n");
    }

    if (get_mode() == EZTRACE_CONVERT) {

        define_kernels_properties();

        int k;
        for (k=1; k<NB_EVENTS; k++){
            addEntityValue(kernels_properties[k].name, "ST_Thread",
                           kernels_properties[k].name, kernels_properties[k].color);
        }
    }
    return 0;
}

void handle_start(kernels_ev_code_t ev)
{
    int size;

    DECLARE_THREAD_ID_STR(thread_id, CUR_INDEX, CUR_THREAD_ID);
    DECLARE_CUR_THREAD(p_thread);
    INIT_KERNELS_THREAD_INFO(p_thread, p_info);

    GET_PARAM_PACKED_1(CUR_EV, size);

    p_info->nb[ev]++;
    p_info->flops[ev]+=size;

    p_info->current_ev = ev;
    p_info->time_start = CURRENT;

    pushState(CURRENT, "ST_Thread", thread_id, kernels_properties[ev].name);
}

void handle_stop()
{
    DECLARE_THREAD_ID_STR(thread_id, CUR_INDEX, CUR_THREAD_ID);
    DECLARE_CUR_THREAD(p_thread);
    INIT_KERNELS_THREAD_INFO(p_thread, p_info);

    p_info->run_time[p_info->current_ev] += (CURRENT - p_info->time_start);

    popState(CURRENT, "ST_Thread", thread_id);
}

int
handle_kernels_events(eztrace_event_t *ev)
{

    if(! CUR_TRACE->start)
        return 0;

    switch (LITL_READ_GET_CODE(ev)) {

    case KERNELS_CODE(STOP):
        handle_stop();
        break;
    default:
        if (LITL_READ_GET_CODE(ev) > KERNELS_PREFIX)
            handle_start((kernels_ev_code_t) (LITL_READ_GET_CODE(ev) - KERNELS_PREFIX));
        break;
    }
    return 1;
}

int
handle_kernels_stats(eztrace_event_t *ev)
{
    return handle_kernels_events(ev);
}


void
print_kernels_stats()
{
    define_kernels_properties();

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
            INIT_KERNELS_THREAD_INFO(p_thread, p_info);
            printf("\tThread %20s\n", thread_container->name);

            int k;
            for (k=1; k<NB_EVENTS; k++){
                printf("Kernel %20s was called %5d times, flops=%8d, duration=%.3g\n",
                       kernels_properties[k].name, p_info->nb[k], p_info->flops[k], p_info->run_time[k]);
            }
        }
    }
}


