/**
 *
 * @file eztrace_wrapper.c
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
 **/

#include "common.h"
#include "kernels_ev_codes.h"

#if defined(PASTIX_WITH_EZTRACE)
#include "eztrace.h"
#include "eztrace_sampling.h"

void start_trace_kernel(enum kernels_ev_code_e state, pastix_int_t flops){
    EZTRACE_EVENT_PACKED_1(KERNELS_CODE(state), flops);
}

void stop_trace_kernel(){
    EZTRACE_EVENT_PACKED_0(KERNELS_CODE(STOP));
}

#else
void start_trace_kernel(enum kernels_ev_code_e state, pastix_int_t flops){
    (void) state;
    (void) flops;
}

void stop_trace_kernel(){
}

#endif /* defined(PASTIX_WITH_EZTRACE) */
