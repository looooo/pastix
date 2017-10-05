/**
 *
 * @file kernels_trace.h
 *
 * Wrappers to trace kernels.
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

#ifndef _kernels_trace_h_
#define _kernels_trace_h_

#include "common.h"
#include "flops.h"

#define PastixKernelStop 0

typedef enum pastix_ktype0_e {
    PastixKernelLvl0Facto,
    PastixKernelLvl0Solve,
    PastixKernelLvl0Diag,
    PastixKernelLvl0Nbr
} pastix_ktype0_t;

typedef enum pastix_ktype_e {
    PastixKernelGETRF,
    PastixKernelHETRF,
    PastixKernelPOTRF,
    PastixKernelPXTRF,
    PastixKernelSYTRF,
    PastixKernelSCALOCblk,
    PastixKernelSCALOBlok,
    PastixKernelTRSMCblk1d,
    PastixKernelTRSMCblk2d,
    PastixKernelTRSMCblkLR,
    PastixKernelTRSMBlok2d,
    PastixKernelTRSMBlokLR,
    PastixKernelGEMMCblk1d1d,
    PastixKernelGEMMCblk1d2d,
    PastixKernelGEMMCblk2d2d,
    PastixKernelGEMMCblkFRLR,
    PastixKernelGEMMCblkLRLR,
    PastixKernelGEMMBlok2d2d,
    PastixKernelGEMMBlok2dLR,
    PastixKernelLvl1Nbr
} pastix_ktype_t;

typedef enum pastix_ktype2_e {

    /* General kernels: similar in low-rank and dense */
    PastixKernelLvl2GETRF,
    PastixKernelLvl2HETRF,
    PastixKernelLvl2POTRF,
    PastixKernelLvl2PXTRF,
    PastixKernelLvl2SYTRF,

    /* Dense operations */
    PastixKernelLvl2_FR_TRSM, /**< trsm on a dense block           */
    PastixKernelLvl2_FR_GEMM, /**< gemm between three dense blocks */

    /* Low-rank operations */
    PastixKernelLvl2_LR_INIT,          /**< try to compress a dense block (RRQR)                  */
    PastixKernelLvl2_LR_INIT_Q,        /**< form Q when compression succeeded                     */
    PastixKernelLvl2_LR_TRSM,          /**< trsm on a low-rank block                              */
    PastixKernelLvl2_LR_GEMM_PRODUCT,  /**< formation of a product of low-rank blocks             */
    PastixKernelLvl2_LR_GEMM_ADD_Q,    /**< getrf/unmqr during recompression                      */
    PastixKernelLvl2_LR_GEMM_ADD_RRQR, /**< compression of concatenated matrices in recompression */
    PastixKernelLvl2_LR_UNCOMPRESS,    /**< uncompress a low-rank block into a dense block        */

    PastixKernelLvl2Nbr
} pastix_ktype2_t;

extern volatile double kernels_flops[PastixKernelLvl1Nbr];

#define PastixKernelsNbr (PastixKernelLvl0Nbr + PastixKernelLvl1Nbr + PastixKernelLvl2Nbr)

#if defined(PASTIX_WITH_EZTRACE)

#include "eztrace_module/kernels_ev_codes.h"

extern int pastix_eztrace_level;

#else

static inline void kernel_trace_start_lvl0( pastix_ktype0_t ktype ) { (void)ktype; };
static inline void kernel_trace_stop_lvl0 ( double flops ) { (void)flops; };
static inline void kernel_trace_start_lvl2( pastix_ktype2_t ktype ) { (void)ktype; };
static inline void kernel_trace_stop_lvl2 ( double flops ) { (void)flops; };

#endif

#if defined(PASTIX_GENERATE_MODEL)

typedef struct pastix_model_entry_s {
    pastix_ktype_t ktype;
    int m;
    int n;
    int k;
    double time;
} pastix_model_entry_t;

extern pastix_model_entry_t *model_entries;     /*< Array to all entries                 */
extern volatile int32_t      model_entries_nbr; /*< Index of the last entry in the array */
extern int32_t               model_size;        /*< Allocated size of the array          */

#endif

void   kernelsTraceStart( const pastix_data_t *pastix_data );
double kernelsTraceStop(  const pastix_data_t *pastix_data );

/**
 *******************************************************************************
 *
 * @brief Start the trace of a single kernel
 *
 *******************************************************************************
 *
 * @param[in] ktype
 *          Type of the kernel starting that need to be traced.
 *          With EZTrace mode, this call is empty if the environment variable
 *          PASTIX_EZTRACE_LEVEL is different from 1.
 *
 *******************************************************************************
 *
 * @return the starting time if PASTIX_GENERATE_MODEL is enabled, 0. otherwise.
 *
 *******************************************************************************/
static inline double
kernel_trace_start( pastix_ktype_t ktype )
{
    double time = 0.;

#if defined(PASTIX_WITH_EZTRACE)

    if (pastix_eztrace_level == 1) {
        EZTRACE_EVENT_PACKED_0( KERNELS_LVL1_CODE(ktype) );
    }

#endif

#if defined(PASTIX_GENERATE_MODEL)

    time = clockGet();

#endif

    (void)ktype;
    return time;
}

/**
 *******************************************************************************
 *
 * @brief Stop the trace of a single kernel
 *
 *******************************************************************************
 *
 * @param[in] ktype
 *          Type of the kernel starting that need to be traced.
 *          With EZTrace mode, this call is empty if the environment variable
 *          PASTIX_EZTRACE_LEVEL is different from 1.
 *
 * @param[in] m
 *          The m parameter of the kernel (used by xxTRF, TRSM, and GEMM)
 *
 * @param[in] n
 *          The n parameter of the kernel (used by TRSM, and GEMM)
 *
 * @param[in] k
 *          The k parameter of the kernel (used by GEMM)
 *
 * @param[in] flops
 *          The number of flops of the kernel
 *
 * @param[in] starttime
 *          The stating time of the kernel. Used only if PASTIX_GENERATE_MODEL
 *          is enabled.
 *
 *******************************************************************************/
static inline void
kernel_trace_stop( pastix_ktype_t ktype, int m, int n, int k, double flops, double starttime )
{

#if defined(PASTIX_WITH_EZTRACE)

    if (pastix_eztrace_level == 1) {
        EZTRACE_EVENT_PACKED_1( KERNELS_CODE( PastixKernelStop ), flops );
    }

#endif

#if defined(PASTIX_GENERATE_MODEL)

    {
        double  time = clockGet() - starttime;
        int32_t index = pastix_atomic_inc_32b( &model_entries_nbr );

        if ( index < model_size ) {
            model_entries[index].ktype = ktype;
            model_entries[index].m     = m;
            model_entries[index].n     = n;
            model_entries[index].k     = k;
            model_entries[index].time  = time;
        }
        else {
            fprintf(stderr, "WARNING: too many kernels to log %d\n", index);
        }
    }

#endif

    /* { */
    /*     double oldval, newval; */
    /*     do { */
    /*         oldval = (uint64_t)(kernels_flops[ktype]); */
    /*         newval = oldval + flops; */
    /*     } while( !pastix_atomic_cas_64b( (uint64_t*)(kernels_flops + ktype), */
    /*                                      (uint64_t)oldval, (uint64_t)newval ) ); */
    /* } */

    (void)flops;
    (void)m;
    (void)n;
    (void)k;
    (void)starttime;
    return;
}

#endif/* _kernels_trace_h_ */
