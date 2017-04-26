#include "common.h"
#include "kernels_ev_codes.h"

#ifdef PASTIX_WITH_EZTRACE
#include "eztrace.h"
#include "eztrace_sampling.h"

void trace_kernel(enum kernels_ev_code_e state, pastix_int_t flops){
    EZTRACE_EVENT_PACKED_1(KERNELS_CODE(state), flops);
}

#else
void trace_kernel(enum kernels_ev_code_e state, pastix_int_t flops){
    (void) state;
    (void) flops;
}

#endif
