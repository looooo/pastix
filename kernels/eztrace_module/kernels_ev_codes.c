/**
 *
 * @file kernels_ev_codes.x
 *
 * Wrappers to trace kernels with eztrace
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
#include "eztrace_module/kernels_ev_codes.h"

extern int eztrace_level;

/**
 * @brief Start eztrace module
 */
void start_eztrace_kernels(){
#if defined(PASTIX_WITH_EZTRACE)
    /* Set tracing level */
    char* LEVEL = pastix_getenv ("EZTRACE_LEVEL");
    if (LEVEL != NULL){
        eztrace_level = atoi(LEVEL);
    }
    else{
        eztrace_level = 1;
    }

    eztrace_start ();
#endif /* defined(PASTIX_WITH_EZTRACE) */
}

/**
 * @brief Stop eztrace module
 */
void stop_eztrace_kernels(){
#if defined(PASTIX_WITH_EZTRACE)
    eztrace_stop ();
#endif /* defined(PASTIX_WITH_EZTRACE) */
}

