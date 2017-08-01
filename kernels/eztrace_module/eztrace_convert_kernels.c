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
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "eztrace_convert_kernels.h"

static kernels_t kernels_properties[KERNELS_NB_EVENTS];

/**
 *******************************************************************************
 *
 * @brief Initialize data associated with each kernel for a thread
 *
 *******************************************************************************
 *
 * @param[in] p_thread
 *          The reference to the thread managed by eztrace
 *
 * @param[in] stats
 *          - stats == 0, save event
 *          - stats == 1, save statistics
 *
 *******************************************************************************
 *
 * @return  The structure containing the thread and statistics for each kernel
 *
 *******************************************************************************/
static kernels_thread_info_t *kernels_register_thread_hook(struct thread_info_t *p_thread,
                                                           int stats) {

    kernels_thread_info_t *p_info = (kernels_thread_info_t*) malloc(sizeof(kernels_thread_info_t));

    p_info->p_thread = p_thread;

    if (stats == 1){
        int    *nb       = malloc(KERNELS_NB_EVENTS * sizeof(int));
        double *flops    = malloc(KERNELS_NB_EVENTS * sizeof(double));
        double *run_time = malloc(KERNELS_NB_EVENTS * sizeof(double));

        int i;
        for (i=0; i<KERNELS_NB_EVENTS; i++){
            nb[i]        = 0;
            flops[i]     = 0;
            run_time[i]  = 0;
        }

        p_info->nb       = nb;
        p_info->flops    = flops;
        p_info->run_time = run_time;
    }

    ezt_hook_list_add(&p_info->p_thread->hooks, p_info, KERNELS_EVENTS_ID);
    return p_info;
}

/**
 * @brief Start to use eztrace by loading PaStiX module
 */
struct eztrace_convert_module kernels_module;
void libinit(void) __attribute__ ((constructor));
void libinit(void)
{
    kernels_module.api_version   = EZTRACE_API_VERSION;
    kernels_module.init          = eztrace_convert_kernels_init;
    kernels_module.handle        = handle_kernels_events;
    kernels_module.handle_stats  = handle_kernels_stats;
    kernels_module.print_stats   = print_kernels_stats;
    kernels_module.module_prefix = KERNELS_EVENTS_ID;

    asprintf(&kernels_module.name, "kernels");
    asprintf(&kernels_module.description, "PaStiX kernels");

    kernels_module.token.data = &kernels_module;
    eztrace_convert_register_module(&kernels_module);
}

/**
 * @brief Stop to use eztrace
 */
void libfinalize(void) __attribute__ ((destructor));
void libfinalize(void)
{
    printf("unloading module \n");
}


/**
 * @brief Define events properties such as name or color
 */
void define_kernels_properties()
{
    int i;
    for (i=1; i<KERNELS_NB_EVENTS; i++){
        kernels_properties[i] = (kernels_t) {"unknown", GTG_BLACK};
    }

    /* Low-rank operations */
    kernels_properties[LR_INIT]   = (kernels_t) {"lr_init",   GTG_YELLOW};
    kernels_properties[LR_INIT_Q] = (kernels_t) {"lr_init_q", GTG_YELLOW};
    kernels_properties[LR_TRSM]   = (kernels_t) {"lr_trsm",   GTG_SEABLUE};

    kernels_properties[LR_GEMM_PRODUCT]  = (kernels_t) {"lr_gemm_product",  GTG_GREEN};
    kernels_properties[LR_GEMM_ADD_Q]    = (kernels_t) {"lr_gemm_add_q",    GTG_ORANGE};
    kernels_properties[LR_GEMM_ADD_RRQR] = (kernels_t) {"lr_gemm_add_rrqr", GTG_LIGHTBROWN};

    kernels_properties[UNCOMPRESS] = (kernels_t) {"lr_uncompress", GTG_PINK};

    /* Dense operations */
    kernels_properties[GETRF]      = (kernels_t) {"getrf",      GTG_RED};
    kernels_properties[POTRF]      = (kernels_t) {"potrf",      GTG_RED};
    kernels_properties[DENSE_TRSM] = (kernels_t) {"dense_trsm", GTG_SEABLUE};
    kernels_properties[DENSE_GEMM] = (kernels_t) {"dense_gemm", GTG_GREEN};
}

/**
 *******************************************************************************
 *
 * @brief Init events metainformation such as name or color
 *
 *******************************************************************************
 *
 * @retval 0 The events where correcly initialized
 *
 * @retval 1 The events where not correcly initialized
 *
 *******************************************************************************/
int eztrace_convert_kernels_init()
{
    if (get_mode() == EZTRACE_CONVERT) {

        define_kernels_properties();

        int k;
        for (k=1; k<KERNELS_NB_EVENTS; k++){
            addEntityValue(kernels_properties[k].name, "ST_Thread",
                           kernels_properties[k].name, kernels_properties[k].color);
        }
    }
    return 0;
}

/**
 *******************************************************************************
 *
 * @brief Handle the start of an elemental event
 *
 *******************************************************************************
 *
 * @param[in] ev
 *          The reference of the kernel for which statistics will be updated
 *
 * @param[in] stats
 *          - stats == 0, save event
 *          - stats == 1, save statistics
 *
 *******************************************************************************/
void handle_start(kernels_ev_code_t ev, int stats)
{
    DECLARE_THREAD_ID_STR(thread_id, CUR_INDEX, CUR_THREAD_ID);
    DECLARE_CUR_THREAD(p_thread);
    INIT_KERNELS_THREAD_INFO(p_thread, p_info, stats);

    if (stats == 1){

        p_info->nb[ev]++;
        p_info->current_ev = ev;
        p_info->time_start = CURRENT;
    }
    else{
        pushState(CURRENT, "ST_Thread", thread_id, kernels_properties[ev].name);
    }
}

/**
 *******************************************************************************
 *
 * @brief Handle the end of an elemental event
 *
 *******************************************************************************
 *
 * @param[in] stats
 *          - stats == 0, save event
 *          - stats == 1, save statistics
 *
 *******************************************************************************/
void handle_stop(int stats)
{
    double size;

    DECLARE_THREAD_ID_STR(thread_id, CUR_INDEX, CUR_THREAD_ID);
    DECLARE_CUR_THREAD(p_thread);
    INIT_KERNELS_THREAD_INFO(p_thread, p_info, stats);

    if (stats == 1){
        GET_PARAM_PACKED_1(CUR_EV, size);
        p_info->run_time[p_info->current_ev] += (CURRENT - p_info->time_start);
        p_info->flops[p_info->current_ev]    += size;
    }
    else{
        popState(CURRENT, "ST_Thread", thread_id);
    }
}

/**
 *******************************************************************************
 *
 * @brief Fonction called by eztrace_convert to handle a single event
 *
 *******************************************************************************
 *
 * @param[in] ev
 *          The event to be handled
 *
 *******************************************************************************
 *
 * @retval 0 The event was not correclty handled
 *
 * @retval 1 The event was correclty handled
 *
 *******************************************************************************/
int handle_kernels_events(eztrace_event_t *ev)
{
    if(! CUR_TRACE->start)
        return 0;

    switch (LITL_READ_GET_CODE(ev)) {

    case KERNELS_CODE(STOP):
        handle_stop(0);
        break;
    default:
        if (LITL_READ_GET_CODE(ev) > KERNELS_PREFIX)
            handle_start((kernels_ev_code_t) (LITL_READ_GET_CODE(ev) - KERNELS_PREFIX), 0);
        break;
    }
    return 1;
}

/**
 *******************************************************************************
 *
 * @brief Fonction called by eztrace_stats to handle a single event
 *
 *******************************************************************************
 *
 * @param[in] ev
 *          The event to be handled
 *
 *******************************************************************************
 *
 * @retval 0 The event was not correclty handled
 *
 * @retval 1 The event was correclty handled
 *
 *******************************************************************************/
int handle_kernels_stats(eztrace_event_t *ev)
{
    if(! CUR_TRACE->start)
        return 0;

    switch (LITL_READ_GET_CODE(ev)) {

    case KERNELS_CODE(STOP):
        handle_stop(1);
        break;
    default:
        if (LITL_READ_GET_CODE(ev) > KERNELS_PREFIX)
            handle_start((kernels_ev_code_t) (LITL_READ_GET_CODE(ev) - KERNELS_PREFIX), 1);
        break;
    }
    return 1;
}

/**
 * @brief Print the statistics of each kernel for each thread
 */
void print_kernels_stats()
{
    int64_t i, j, k;
    double total_flops = 0;

    define_kernels_properties();

    /* Browse the list of processes */
    for (i = 0; i < NB_TRACES; i++) {
        struct eztrace_container_t *p_process = GET_PROCESS_CONTAINER(i);

        /* For each process, browse the list of threads */
        for(j=0; j<(int64_t)(p_process->nb_children); j++) {

            struct eztrace_container_t *thread_container;
            struct thread_info_t       *p_thread;

            thread_container = p_process->children[j];
            p_thread = (struct thread_info_t*)(thread_container->container_info);

            if(!p_thread)
                continue;

            INIT_KERNELS_THREAD_INFO(p_thread, p_info, 1);
            printf("\tThread %20s\n", thread_container->name);

            for (k=1; k<KERNELS_NB_EVENTS; k++){
                double perf = 1000 * p_info->flops[k] / p_info->run_time[k];
                total_flops += p_info->flops[k];

                if (p_info->nb[k] > 0)
                    printf("Kernel %20s was called %8d times, flops=%8.3g, time=%7.2g s, perf=%7.2lf %cFlop/s\n",
                           kernels_properties[k].name, p_info->nb[k],
                           p_info->flops[k], p_info->run_time[k] / 1000.,
                           printflopsv( perf ), printflopsu( perf ) );
            }
        }
    }
    printf("\n\n\tTotal number of operations: %5.2lf %cFlops\n",
           printflopsv( total_flops ), printflopsu( total_flops ) );
}
