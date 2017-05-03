/**
 *
 * @file kernels_ev_codes.h
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

#ifndef __KERNELS_EV_CODES_H__
#define __KERNELS_EV_CODES_H__

#include "common.h"

#if defined(PASTIX_WITH_EZTRACE)
#include "eztrace.h"
#include "eztrace_sampling.h"
#include "ev_codes.h"

#define KERNELS_EVENTS_ID    USER_MODULE_ID(0x51)
#define KERNELS_PREFIX       (KERNELS_EVENTS_ID << 3)
#define KERNELS_CODE(event)  (KERNELS_PREFIX | event )

#endif /* defined(PASTIX_WITH_EZTRACE) */

enum kernels_ev_code_e {
    STOP,

    /* Low-rank operations */
    LR_ALLOC,
    LR_TRSM,
    LR_GEMM,

    /* General kernels: similar in low-rank and dense */
    GETRF,
    POTRF,

    /* Dense operations */
    DENSE_TRSM,
    DENSE_GEMM,

    NB_EVENTS,
};
typedef enum kernels_ev_code_e kernels_ev_code_t;

/**
 *******************************************************************************
 *
 * @ingroup pastix_eztrace
 *
 * @brief Start to trace a kernel
 *
 *******************************************************************************
 *
 * @param[in] state
 *          The kernel's name
 *
 * @param[in] flops
 *          The number of operations performed
 *
 *******************************************************************************/
static inline void start_trace_kernel(kernels_ev_code_t state, double flops){
#if defined(PASTIX_WITH_EZTRACE)
    EZTRACE_EVENT_PACKED_1(KERNELS_CODE(state), flops);
#else
    (void) state;
    (void) flops;
#endif /* defined(PASTIX_WITH_EZTRACE) */
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_eztrace
 *
 * @brief Stop to trace a kernel
 *
 *******************************************************************************/
static inline void stop_trace_kernel(){
#if defined(PASTIX_WITH_EZTRACE)
    EZTRACE_EVENT_PACKED_0(KERNELS_CODE(STOP));
#endif /* defined(PASTIX_WITH_EZTRACE) */
}

#endif	/* __KERNELS_EV_CODES_H__ */
