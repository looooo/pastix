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
 * @addtogroup eztrace_dev
 * @{
 *
 **/

#ifndef __KERNELS_EV_CODES_H__
#define __KERNELS_EV_CODES_H__

#include "common.h"
#include "flops.h"

#if defined(PASTIX_WITH_EZTRACE)
#include "eztrace.h"
#include "eztrace_sampling.h"
#include "ev_codes.h"

#define KERNELS_EVENTS_ID    USER_MODULE_ID(0x51)
#define KERNELS_PREFIX       (KERNELS_EVENTS_ID << 5)
#define KERNELS_CODE(event)  (KERNELS_PREFIX | event )

#endif /* defined(PASTIX_WITH_EZTRACE) */

typedef enum kernels_ev_code_e {
    STOP,

    /* Low-rank operations */
    LR_INIT,
    LR_INIT_Q,
    LR_TRSM,
    LR_GEMM,

    /* General kernels: similar in low-rank and dense */
    GETRF,
    POTRF,

    /* Dense operations */
    DENSE_TRSM,
    DENSE_GEMM,

    KERNELS_NB_EVENTS,
} kernels_ev_code_t;

/**
 *******************************************************************************
 *
 * @brief Start to trace a kernel
 *
 *******************************************************************************
 *
 * @param[in] state
 *          The kernel's name
 *
 *******************************************************************************/
static inline void start_trace_kernel(kernels_ev_code_t state){
#if defined(PASTIX_WITH_EZTRACE)
    EZTRACE_EVENT_PACKED_0(KERNELS_CODE(state));
#else
    (void) state;
#endif /* defined(PASTIX_WITH_EZTRACE) */
}

/**
 *******************************************************************************
 *
 * @brief Stop to trace a kernel
 *
 *******************************************************************************
 *
 * @param[in] flops
 *          The number of operations performed during the call
 *
 *******************************************************************************/
static inline void stop_trace_kernel(double flops){
#if defined(PASTIX_WITH_EZTRACE)
    EZTRACE_EVENT_PACKED_1(KERNELS_CODE(STOP), flops);
#else
    (void) flops;
#endif /* defined(PASTIX_WITH_EZTRACE) */
}

/**
 * @brief Start eztrace module
 */
static inline void start_eztrace_kernels(){
#if defined(PASTIX_WITH_EZTRACE)
  eztrace_start ();
#endif /* defined(PASTIX_WITH_EZTRACE) */
}

/**
 * @brief Stop eztrace module
 */
static inline void stop_eztrace_kernels(){
#if defined(PASTIX_WITH_EZTRACE)
  eztrace_stop ();
#endif /* defined(PASTIX_WITH_EZTRACE) */
}

#endif	/* __KERNELS_EV_CODES_H__ */
