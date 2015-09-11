/*
 * Copyright (c) 2009-2012 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Timing system specifics imported from PaRSEC project.
 */
/*
 *  File: timing.h
 *
 *  Part of a parallel direct block solver.
 *
 *  These file defines the macro to measure computational times of the
 *  different code sections.
 *
 *  Authors:
 *    Mathieu  Faverge    - faverge@labri.fr
 *    David    GOUDIN     - .
 *    Pascal   HENON      - henon@labri.fr
 *    Xavier   LACOSTE    - lacoste@labri.fr
 *    Francois PELLEGRINI - .
 *    Pierre   RAMET      - ramet@labri.fr
 *
 */
#ifndef _TIMING_H_
#define _TIMING_H_

#include <stdint.h>

typedef double Clock;

/** TIMING SYSTEM-SPECIFICS **/
#if defined(HAVE_MPI)
#define clockGet() MPI_Wtime()
#else
#if defined(HAVE_CLOCK_GETTIME)
#include <unistd.h>
#include <time.h>
static inline double clockGet(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ((double) ts.tv_sec + (double) ts.tv_nsec * (double)1.0e-9L);
}

#elif defined(HAVE_GETRUSAGE)
#include <sys/time.h>
static inline double clockGet(void)
{
    struct rusage       data;
    getrusage (RUSAGE_SELF, &data);
    return (((double) data.ru_utime.tv_sec  + (double) data.ru_stime.tv_sec) +
            ((double) data.ru_utime.tv_usec + (double) data.ru_stime.tv_usec) *
            1.0e-6L);
}

#elif defined(__IA64)
static inline double clockGet(void)
{
    uint64_t ret;
    __asm__ __volatile__ ("mov %0=ar.itc" : "=r"(ret));
    return (double)ret;
}

#elif defined(__X86)
static inline double clockGet(void)
{
    uint64_t ret;
    unsigned hi, lo;
    __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
    ret = ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
    return (double)ret;
}

#elif defined(__bgp__)

#include <bpcore/ppc450_inlines.h>
static inline double clockGet(void)
{
    return (double)_bgp_GetTimeBase();
}

#elif defined(HAVE_GETRUSAGE)

#include <sys/time.h>
static inline double clockGet(void)
{
    struct rusage       data;
    getrusage (RUSAGE_SELF, &data);
    return (((double) data.ru_utime.tv_sec  + (double) data.ru_stime.tv_sec) +
            ((double) data.ru_utime.tv_usec + (double) data.ru_stime.tv_usec) *
            1.0e-6L);
}

#else

#include <sys/time.h>
static inline double clockGet(void)
{
    struct timeval tv;
    gettimeofday( &tv, NULL );
    return ((double)tv.tv_sec + (double) tv.tv_usec * (double)1.0e-6L);
}

#endif
#endif /* defined(HAVE_MPI) */

#define clockInit(clk)  do { clk = clockGet(); } while(0)
#define clockVal(clk)   (clk)

#define clockStart(clk) do { clk = clockGet(); } while(0)
#define clockStop(clk)  do { clk = clockGet() - (clk); } while(0)

#if defined(HAVE_MPI)
#define clockSyncStart(clk, comm) do {          \
        MPI_Barrier(comm);                      \
        clk = clockGet(); } while(0)
#define clockSyncStop(clk, comm)  do {          \
        MPI_Barrier(comm);                      \
        clk = clockGet() - (clk); } while(0)
#else
#define clockSyncStart(clk) do { clk = clockGet(); } while(0)
#define clockSyncStop(clk)  do { clk = clockGet() - (clk); } while(0)
#endif

#endif /* _TIMING_H_ */
