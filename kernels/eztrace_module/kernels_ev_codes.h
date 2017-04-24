/*
 * Copyright (C) CNRS, INRIA, Université Bordeaux 1, Télécom SudParis
 * See COPYING in top-level directory.
 */

#ifndef __KERNELS_EV_CODES_H__
#define __KERNELS_EV_CODES_H__

/* This file defines the event codes that are used by the example
 * module.
 */
#include "ev_codes.h"

/* 7-bits event codes prefix. This identifies the module and thus should be
 * unique.
 * The 0x0? prefix is reserved for eztrace internal use. Thus you can
 * use any prefix between 0x80 and 0xff.
 */
#define KERNELS_EVENTS_ID    USER_MODULE_ID(0x99)
#define KERNELS_PREFIX       (KERNELS_EVENTS_ID << 3)

#define KERNELS_LRALLOC_START   (KERNELS_PREFIX | 0x0010)
#define KERNELS_LRALLOC_STOP    (KERNELS_PREFIX | 0x0011)


/* Define various event codes used by the example module
 * The 2 most significant bytes should correspond to the module id,
 * as below:
 */


#endif	/* __KERNELS_EV_CODES_H__ */
