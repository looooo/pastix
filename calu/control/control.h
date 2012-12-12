/**
 *
 * @file control.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Jakub Kurzak
 * @date 2010-11-15
 *
 **/
#ifndef _PLASMA_CONTROL_H_
#define _PLASMA_CONTROL_H_

#ifndef __cplusplus
extern int pthread_getconcurrency(void);
extern int pthread_setconcurrency(int);
#endif

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************************//**
 *  Internal routines
 **/
void  plasma_barrier_init(plasma_context_t *plasma);
void  plasma_barrier(plasma_context_t *plasma);
void *plasma_parallel_section(void *plasma);
int   plasma_setaffinity(int rank);
int   plasma_unsetaffinity();
int   plasma_yield();
void  plasma_topology_init();
void  plasma_topology_finalize();
int   plasma_get_numthreads();
int   plasma_get_numthreads_numa();
int   plasma_get_affthreads(int *coresbind);

/***************************************************************************//**
 *  User routines
 **/
int PLASMA_Init(int cores);
int PLASMA_Init_Affinity(int cores, int *bindtab);
int PLASMA_Finalize();

#ifdef __cplusplus
}
#endif

#endif
