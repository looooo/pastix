/**
 *
 * @file kernels.c
 *
 * Non precision dependent routines and variables associated to the kernels.
 *
 * @copyright 2015-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Tony Delarue
 * @date 2019-04-11
 *
 **/
#include "common.h"

pastix_atomic_lock_t pastix_mpi_lock = PASTIX_ATOMIC_UNLOCKED;


