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

#ifndef _kernels_ev_codes_h_
#define _kernels_ev_codes_h_

#include "common.h"
#include "flops.h"

#if defined(PASTIX_WITH_EZTRACE)
#include "eztrace.h"
#include "eztrace_sampling.h"
#include "ev_codes.h"

#define KERNELS_EVENTS_ID    USER_MODULE_ID(0x51)
#define KERNELS_PREFIX       (KERNELS_EVENTS_ID << 5)
#define KERNELS_CODE(event)  (KERNELS_PREFIX | event )

extern int pastix_eztrace_level;

#endif /* defined(PASTIX_WITH_EZTRACE) */

/**
 * @brief Kernels
 */
typedef enum kernels_ev_code_e {
    STOP,

    /* Low-rank operations */
    LR_INIT,          /**< try to compress a dense block (RRQR)                  */
    LR_INIT_Q,        /**< form Q when compression succeeded                     */
    LR_TRSM,          /**< trsm on a low-rank block                              */
    LR_GEMM_PRODUCT,  /**< formation of a product of low-rank blocks             */
    LR_GEMM_ADD_Q,    /**< getrf/unmqr during recompression                      */
    LR_GEMM_ADD_RRQR, /**< compression of concatenated matrices in recompression */
    UNCOMPRESS,       /**< uncompress a low-rank block into a dense block        */

    /* General kernels: similar in low-rank and dense */
    GETRF,
    HETRF,
    POTRF,
    PXTRF,
    SYTRF,

    /* Dense operations */
    DENSE_TRSM, /**< trsm on a dense block           */
    DENSE_GEMM, /**< gemm between three dense blocks */

    LVL1_GETRF,
    LVL1_HETRF,
    LVL1_POTRF,
    LVL1_PXTRF,
    LVL1_SYTRF,
    LVL1_SCALO,
    LVL1_TRSM_CBLK_1D,
    LVL1_TRSM_CBLK_2D,
    LVL1_TRSM_CBLK_LR,
    LVL1_TRSM_BLOK_2D,
    LVL1_TRSM_BLOK_LR,
    LVL1_GEMM_CBLK_1D1D,
    LVL1_GEMM_CBLK_1D2D,
    LVL1_GEMM_CBLK_2D2D,
    LVL1_GEMM_CBLK_FRLR,
    LVL1_GEMM_CBLK_LRLR,
    LVL1_GEMM_BLOK_2D2D,
    LVL1_GEMM_BLOK_2DLR,
    LVL1_GEMDM,

    KERNELS_NB_EVENTS,
} kernels_ev_code_t;

void start_eztrace_kernels();
void stop_eztrace_kernels();

/**
 *******************************************************************************
 *
 * @brief Start to trace a kernel
 *
 *******************************************************************************
 *
 * @param[in] level
 *          Level of tracing defined by environment variable EZTRACE_LVL
 *
 * @param[in] state
 *          The kernel's name
 *
 *******************************************************************************/
static inline void start_trace_kernel(int level, kernels_ev_code_t state){
#if defined(PASTIX_WITH_EZTRACE)
    if (level == pastix_eztrace_level){
        EZTRACE_EVENT_PACKED_0(KERNELS_CODE(state));
    }
#else
    (void) level;
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
 * @param[in] level
 *          Level of tracing defined by environment variable EZTRACE_LVL
 *
 * @param[in] flops
 *          The number of operations performed during the call
 *
 *******************************************************************************/
static inline void stop_trace_kernel(int level, double flops){
#if defined(PASTIX_WITH_EZTRACE)
    if (level == pastix_eztrace_level){
        EZTRACE_EVENT_PACKED_1(KERNELS_CODE(STOP), flops);
    }
#else
    (void) level;
    (void) flops;
#endif /* defined(PASTIX_WITH_EZTRACE) */
}

#endif	/* _kernels_ev_codes_h_ */
