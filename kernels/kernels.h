/**
 *
 * @file kernels.h
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
#ifndef _kernels_h_
#define _kernels_h_

extern pthread_mutex_t    pastix_comm_lock;
extern volatile pthread_t pastix_comm_tid;

#endif /* _kernels_h_ */
