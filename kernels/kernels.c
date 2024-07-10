/**
 *
 * @file kernels.c
 *
 * Non precision dependent routines and variables associated to the kernels.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2024-07-05
 *
 **/
#include "common.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
pthread_mutex_t    pastix_comm_lock = PTHREAD_MUTEX_INITIALIZER;
volatile pthread_t pastix_comm_tid  = (pthread_t)-1;
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
