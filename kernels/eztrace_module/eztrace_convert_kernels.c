/**
 *
 * @file eztrace_convert_kernels.c
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
#define _GNU_SOURCE 1
#include "eztrace_convert_kernels.h"

struct eztrace_convert_module kernels_module;

static kernels_t kernels_properties[PastixKernelsNbr];

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
static kernels_thread_info_t *
kernels_register_thread_hook( struct thread_info_t *p_thread,
                              int stats )
{
    kernels_thread_info_t *p_info = (kernels_thread_info_t*) malloc(sizeof(kernels_thread_info_t));

    p_info->p_thread = p_thread;

    if (stats == 1){
        int    *nb       = malloc(PastixKernelsNbr * sizeof(int));
        double *flops    = malloc(PastixKernelsNbr * sizeof(double));
        double *run_time = malloc(PastixKernelsNbr * sizeof(double));

        int i;
        for (i=0; i<PastixKernelsNbr; i++){
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
void libinit(void) __attribute__ ((constructor));
void libinit(void)
{
    kernels_module.api_version   = EZTRACE_API_VERSION;
    kernels_module.init          = eztrace_convert_kernels_init;
    kernels_module.handle        = handle_kernels_events;
    kernels_module.handle_stats  = handle_kernels_stats;
    kernels_module.print_stats   = print_kernels_stats;
    kernels_module.module_prefix = KERNELS_EVENTS_ID;

    asprintf(&kernels_module.name,        "kernels"       );
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
void
define_kernels_properties()
{
    kernels_t *kernels_lvl1, *kernels_lvl2;
    int i;

    for (i=0; i<PastixKernelsNbr; i++){
        kernels_properties[i] = (kernels_t) {"unknown", GTG_BLACK};
    }

    /* Level 1 kernels */
    kernels_lvl1 = kernels_properties + PastixKernelLvl0Nbr;
    kernels_lvl1[PastixKernelGETRF]     = (kernels_t) {"lvl1_getrf",      GTG_RED};
    kernels_lvl1[PastixKernelHETRF]     = (kernels_t) {"lvl1_hetrf",      GTG_RED};
    kernels_lvl1[PastixKernelPOTRF]     = (kernels_t) {"lvl1_potrf",      GTG_RED};
    kernels_lvl1[PastixKernelPXTRF]     = (kernels_t) {"lvl1_pxtrf",      GTG_RED};
    kernels_lvl1[PastixKernelSYTRF]     = (kernels_t) {"lvl1_sytrf",      GTG_RED};
    kernels_lvl1[PastixKernelSCALOCblk] = (kernels_t) {"lvl1_scalo_cblk", GTG_SEABLUE};
    kernels_lvl1[PastixKernelSCALOBlok] = (kernels_t) {"lvl1_scalo_blok", GTG_GREEN};

    kernels_lvl1[PastixKernelTRSMCblk1d  ] = (kernels_t) {"lvl1_trsm_cblk_1d", GTG_BLUE};
    kernels_lvl1[PastixKernelTRSMCblk2d  ] = (kernels_t) {"lvl1_trsm_cblk_2d", GTG_BLUE};
    kernels_lvl1[PastixKernelTRSMCblkLR  ] = (kernels_t) {"lvl1_trsm_cblk_lr", GTG_BLUE};
    kernels_lvl1[PastixKernelTRSMBlok2d  ] = (kernels_t) {"lvl1_trsm_blok_2d", GTG_BLUE};
    kernels_lvl1[PastixKernelTRSMBlokLR  ] = (kernels_t) {"lvl1_trsm_blok_lr", GTG_BLUE};
    kernels_lvl1[PastixKernelGEMMCblk1d1d] = (kernels_t) {"lvl1_gemm_cblk_1d1d", GTG_GREEN};
    kernels_lvl1[PastixKernelGEMMCblk1d2d] = (kernels_t) {"lvl1_gemm_cblk_1d2d", GTG_GREEN};
    kernels_lvl1[PastixKernelGEMMCblk2d2d] = (kernels_t) {"lvl1_gemm_cblk_2d2d", GTG_GREEN};
    kernels_lvl1[PastixKernelGEMMCblkFRLR] = (kernels_t) {"lvl1_gemm_cblk_frlr", GTG_GREEN};
    kernels_lvl1[PastixKernelGEMMCblkLRLR] = (kernels_t) {"lvl1_gemm_cblk_lrlr", GTG_GREEN};
    kernels_lvl1[PastixKernelGEMMBlok2d2d] = (kernels_t) {"lvl1_gemm_blok_2d2d", GTG_GREEN};
    kernels_lvl1[PastixKernelGEMMBlok2dLR] = (kernels_t) {"lvl1_gemm_blok_2dlr", GTG_GREEN};

    /* Level 2 kernels - Factorization kernels */
    kernels_lvl2 = kernels_lvl1 + PastixKernelLvl1Nbr;
    kernels_lvl2[PastixKernelLvl2GETRF] = (kernels_t) {"lvl2_getrf", GTG_RED};
    kernels_lvl2[PastixKernelLvl2HETRF] = (kernels_t) {"lvl2_hetrf", GTG_RED};
    kernels_lvl2[PastixKernelLvl2POTRF] = (kernels_t) {"lvl2_potrf", GTG_RED};
    kernels_lvl2[PastixKernelLvl2PXTRF] = (kernels_t) {"lvl2_pxtrf", GTG_RED};
    kernels_lvl2[PastixKernelLvl2SYTRF] = (kernels_t) {"lvl2_sytrf", GTG_RED};

    /* Level 2 kernels - Full-rank operations */
    kernels_lvl2[PastixKernelLvl2_FR_TRSM] = (kernels_t) {"lvl2_fr_trsm", GTG_SEABLUE};
    kernels_lvl2[PastixKernelLvl2_FR_GEMM] = (kernels_t) {"lvl2_fr_gemm", GTG_GREEN  };

    /* Level 2 kernels - Low-rank operations */
    kernels_lvl2[PastixKernelLvl2_LR_INIT]          = (kernels_t) {"lvl2_lr_init",          GTG_YELLOW    };
    kernels_lvl2[PastixKernelLvl2_LR_INIT_Q]        = (kernels_t) {"lvl2_lr_init_q",        GTG_YELLOW    };
    kernels_lvl2[PastixKernelLvl2_LR_TRSM]          = (kernels_t) {"lvl2_lr_trsm",          GTG_SEABLUE   };
    kernels_lvl2[PastixKernelLvl2_LR_GEMM_PRODUCT]  = (kernels_t) {"lvl2_lr_gemm_product",  GTG_GREEN     };
    kernels_lvl2[PastixKernelLvl2_LR_GEMM_ADD_Q]    = (kernels_t) {"lvl2_lr_gemm_add_q",    GTG_ORANGE    };
    kernels_lvl2[PastixKernelLvl2_LR_GEMM_ADD_RRQR] = (kernels_t) {"lvl2_lr_gemm_add_rrqr", GTG_LIGHTBROWN};
    kernels_lvl2[PastixKernelLvl2_LR_UNCOMPRESS]    = (kernels_t) {"lvl2_lr_uncompress",    GTG_PINK      };
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
int
eztrace_convert_kernels_init()
{
    if (get_mode() == EZTRACE_CONVERT) {
        int k;

        define_kernels_properties();

        for (k=0; k<PastixKernelsNbr; k++) {
            addEntityValue( kernels_properties[k].name, "ST_Thread",
                            kernels_properties[k].name, kernels_properties[k].color );
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
void
handle_start( int event_id, int stats )
{
    DECLARE_THREAD_ID_STR(thread_id, CUR_INDEX, CUR_THREAD_ID);
    DECLARE_CUR_THREAD(p_thread);
    INIT_KERNELS_THREAD_INFO(p_thread, p_info, stats);

    if (stats == 1)
    {
        p_info->nb[event_id]++;
        p_info->current_ev = event_id;
        p_info->time_start = CURRENT;
    }
    else {
        pushState(CURRENT, "ST_Thread", thread_id, kernels_properties[event_id+1].name);
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
void
handle_stop( int stats )
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
int
handle_kernels_events(eztrace_event_t *ev)
{
    int event_id;

    if( !(CUR_TRACE->start) ||
        !IS_A_KERNELS_EV(ev) )
    {
        return 0;
    }

    event_id = KERNELS_GET_CODE( ev );
    switch (event_id) {
    case PastixKernelStop:
        handle_stop( 0 );
        break;
    default:
        handle_start( event_id-1, 0 );
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
    int event_id;

    if( !(CUR_TRACE->start) ||
        !IS_A_KERNELS_EV(ev) )
    {
        return 0;
    }

    event_id = KERNELS_GET_CODE( ev );
    switch (event_id) {
    case PastixKernelStop:
        handle_stop( 1 );
        break;
    default:
        handle_start( event_id-1, 1 );
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

            if( !p_thread ) {
                continue;
            }

            INIT_KERNELS_THREAD_INFO(p_thread, p_info, 1);
            printf("\tThread %20s\n", thread_container->name);

            for (k=0; k<PastixKernelsNbr; k++)
            {
                if (p_info->nb[k] > 0)
                {
                    double flops = p_info->flops[k];
                    double time  = p_info->run_time[k] * 1.e-3;
                    double perf  = flops / time;

                    total_flops += flops;

                    printf( "Kernel %20s was called %8d times, flops=%8.3g, time=%7.2g s, perf=%7.2lf %cFlop/s\n",
                            kernels_properties[k].name, p_info->nb[k],
                            flops, time,
                            printflopsv( perf ), printflopsu( perf ) );
                }
            }
        }
    }
    printf( "\n\n\tTotal number of operations: %5.2lf %cFlops\n",
            printflopsv( total_flops ), printflopsu( total_flops ) );
}
