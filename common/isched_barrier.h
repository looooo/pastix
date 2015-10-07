/**
 *
 * @file isched_barrier.h
 *
 * Copyright (c) 2014-2015 The University of Bordeaux, IPB, LaBRI, Inria -
 *                         Bordeaux-Sud-Ouest.  All rights reserved.
 *
 * Copyright (c) 2009-2010 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 *  PaStiX internal scheduler barrier routines.
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * Contains basic functions to bind threads.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 */
#ifndef _ISCHED_BARRIER_H_
#define _ISCHED_BARRIER_H_

#include <unistd.h>
#include <pthread.h>

/* The Linux includes are completely screwed up right now. Even if they
 * correctly export a _POSIX_BARRIER define the barrier functions are
 * not correctly defined in the pthread.h. So until we figure out
 * how to correctly identify their availability, we will have to
 * disable them.
 */
#if defined(_POSIX_BARRIERS) && (_POSIX_BARRIERS - 20012L) >= 0 && 0

BEGIN_C_DECLS

typedef pthread_barrier_t isched_barrier_t;
#define isched_barrier_init      pthread_barrier_init
#define isched_barrier_wait      pthread_barrier_wait
#define isched_barrier_destroy   pthread_barrier_destroy
#define ISCHED_IMPLEMENT_BARRIERS 0

#else

typedef struct isched_barrier_s {
    int                 count;
    volatile int        curcount;
    volatile int        generation;
    pthread_mutex_t     mutex;
    pthread_cond_t      cond;
} isched_barrier_t;

int isched_barrier_init(isched_barrier_t *barrier, const void *pthread_mutex_attr, unsigned int count);
int isched_barrier_wait(isched_barrier_t*);
int isched_barrier_destroy(isched_barrier_t*);
#define ISCHED_IMPLEMENT_BARRIERS 1

END_C_DECLS

#endif


#endif  /* _ISCHED_BARRIER_H_ */
