/**
 *
 * @file isched.h
 *
 * @copyright 2008-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2010-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 *
 * PaStiX thread binding routines header.
 *
 * @version 6.4.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Alycia Lisito
 * @date 2024-07-05
 *
 **/
#ifndef _isched_h_
#define _isched_h_

#include "isched_barrier.h"

BEGIN_C_DECLS

enum isched_action_e {
    ISCHED_ACT_STAND_BY,
    ISCHED_ACT_PARALLEL,
    ISCHED_ACT_FINALIZE
};

struct isched_s;
typedef struct isched_s isched_t;

/**
 * Thread structure of the execution context of one instance of the scheduler
 */
typedef struct isched_thread_s {
    isched_t        *global_ctx;
    int              rank;
    int              bindto;
} isched_thread_t;

/**
 * Global structure of the execution context of one instance of the scheduler
 */
struct isched_s {
    int              world_size;
    int              socketsnbr;

    isched_barrier_t barrier;
    pthread_mutex_t  statuslock;
    pthread_cond_t   statuscond;
    volatile int     status;

    pthread_t       *tids;
    isched_thread_t *master;

    void           (*pfunc)(isched_thread_t*, void*);
    void            *pargs;
};

#if defined(HAVE_HWLOC)
#include "isched_hwloc.h"
#define isched_topo_init               isched_hwloc_init
#define isched_topo_destroy            isched_hwloc_destroy
#define isched_topo_bind_on_core_index isched_hwloc_bind_on_core_index
#define isched_topo_unbind             isched_hwloc_unbind
#define isched_topo_world_size         isched_hwloc_world_size
#define isched_topo_socketsnbr         isched_hwloc_socketsnbr
#else
#define isched_topo_init               isched_nohwloc_init
#define isched_topo_destroy            isched_nohwloc_destroy
#define isched_topo_bind_on_core_index isched_nohwloc_bind_on_core_index
#define isched_topo_unbind             isched_nohwloc_unbind
#define isched_topo_world_size         isched_nohwloc_world_size
#define isched_topo_socketsnbr         isched_nohwloc_socketsnbr
#endif

int  isched_topo_init( void );
int  isched_topo_destroy( void );
int  isched_topo_bind_on_core_index( int );
int  isched_topo_unbind( void );
int  isched_topo_world_size( void );
int  isched_topo_socketsnbr( void );

static inline void
isched_parallel_call( isched_t *isched, void (*func)(isched_thread_t*, void*), void *args )
{
    isched->pfunc = func;
    isched->pargs = args;
    pthread_mutex_lock(&isched->statuslock);
    isched->status = ISCHED_ACT_PARALLEL;
    pthread_mutex_unlock(&isched->statuslock);
    pthread_cond_broadcast(&isched->statuscond);
    isched_barrier_wait( &(isched->barrier) );
    isched->status = ISCHED_ACT_STAND_BY;
    func( isched->master, args );
    isched_barrier_wait( &(isched->barrier) );
}

isched_t *ischedInit(int cores, const int *coresbind);
int ischedFinalize(isched_t *isched);

END_C_DECLS

#endif /* _isched_h_ */
