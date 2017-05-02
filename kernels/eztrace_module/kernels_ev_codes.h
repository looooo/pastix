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
#include "ev_codes.h"

#define KERNELS_EVENTS_ID    USER_MODULE_ID(0x99)
#define KERNELS_PREFIX       (KERNELS_EVENTS_ID << 3)
#define KERNELS_CODE(event)  (KERNELS_PREFIX | event )

#endif /* defined(PASTIX_WITH_EZTRACE) */

enum kernels_ev_code_e {
    STOP,

    LR_ALLOC,
    LR_TRSM,
    LR_GEMM,

    GETRF,
    POTRF,

    DENSE_TRSM,
    DENSE_GEMM,

    NB_EVENTS,
};

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
void start_trace_kernel(enum kernels_ev_code_e state, pastix_int_t flops);

/**
 *******************************************************************************
 *
 * @ingroup pastix_eztrace
 *
 * @brief Stop to trace a kernel
 *
 *******************************************************************************/
void stop_trace_kernel();

#endif	/* __KERNELS_EV_CODES_H__ */
