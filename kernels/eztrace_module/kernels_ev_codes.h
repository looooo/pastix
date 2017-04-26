/*
 * Copyright (C) CNRS, INRIA, Université Bordeaux 1, Télécom SudParis
 * See COPYING in top-level directory.
 */

#ifndef __KERNELS_EV_CODES_H__
#define __KERNELS_EV_CODES_H__

#include "common.h"

#ifdef PASTIX_WITH_EZTRACE
#include "ev_codes.h"

#define KERNELS_EVENTS_ID    USER_MODULE_ID(0x99)
#define KERNELS_PREFIX       (KERNELS_EVENTS_ID << 3)

#define KERNELS_CODE(event)  (KERNELS_PREFIX | event )
#endif

enum kernels_ev_code_e {
    KERNELS_LRALLOC_START,
    KERNELS_LRALLOC_STOP,
};

void trace_kernel(enum kernels_ev_code_e state, pastix_int_t flops);

#endif	/* __KERNELS_EV_CODES_H__ */
