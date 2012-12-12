/**
 *
 * @file common.h
 *
 *  PLASMA auxiliary routines
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.4.6
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 **/

/***************************************************************************//**
 *  PLASMA facilities of interest to both PLASMA core developer
 *  and also of interest to PLASMA community contributor.
 **/
#ifndef _PLASMA_COMMON_H_
#define _PLASMA_COMMON_H_

#include "global.h"
#include "context.h"
#include "control.h"
#include "core_blas.h"
#include "core_blas_dag.h"
#include "allocate.h"
#include "auxiliary.h"
#include "descriptor.h"
#include "tile.h"
#include "async.h"
#include "bulge.h"

#if defined( _WIN32 ) || defined( _WIN64 )
# include <io.h>
#else  /* NOT WIN32 */
# include <unistd.h>
#endif  /* IF WIN32 */


/** ****************************************************************************
 *  Determine if weak symbol are allowed 
 */
#if defined(linux) || defined(__linux) || defined(__linux__)
#if defined(__GNUC_EXCL__) || defined(__GNUC__) 
#define PLASMA_HAVE_WEAK 1
#endif
#endif

/** ****************************************************************************
 *  Determine FORTRAN names
 */
#if defined(ADD_)
#define PLASMA_FNAME(lcname, UCNAME)        plasma_##lcname##_
#define PLASMA_TILE_FNAME(lcname, UCNAME)   plasma_##lcname##_tile_
#define PLASMA_ASYNC_FNAME(lcname, UCNAME)  plasma_##lcname##_tile_async_
#define PLASMA_WS_FNAME(lcname, UCNAME)     plasma_alloc_workspace_##lcname##_
#define PLASMA_WST_FNAME(lcname, UCNAME)    plasma_alloc_workspace_##lcname##_tile_
#elif defined(NOCHANGE)
#define PLASMA_FNAME(lcname, UCNAME)        plasma_##lcname
#define PLASMA_TILE_FNAME(lcname, UCNAME)   plasma_##lcname##_tile
#define PLASMA_ASYNC_FNAME(lcname, UCNAME)  plasma_##lcname##_tile_async
#define PLASMA_WS_FNAME(lcname, UCNAME)     plasma_alloc_workspace_##lcname
#define PLASMA_WST_FNAME(lcname, UCNAME)    plasma_alloc_workspace_##lcname##_tile
#elif defined(UPCASE)
#define PLASMA_FNAME(lcname, UCNAME)        PLASMA_##UCNAME
#define PLASMA_TILE_FNAME(lcname, UCNAME)   PLASMA_##UCNAME##_TILE
#define PLASMA_ASYNC_FNAME(lcname, UCNAME)  PLASMA_##UCNAME##_TILE_ASYNC
#define PLASMA_WS_FNAME(lcname, UCNAME)     PLASMA_ALLOC_WORKSPACE_##UCNAME
#define PLASMA_WST_FNAME(lcname, UCNAME)    PLASMA_ALLOC_WORKSPACE_##UCNAME##_TILE
#endif

/***************************************************************************//**
 *  Global shortcuts
 **/
#define PLASMA_RANK        plasma_rank(plasma)
#define PLASMA_SIZE        plasma->world_size
#define PLASMA_GRPSIZE     plasma->group_size
#define PLASMA_NB          plasma->nb
#define PLASMA_IB          plasma->ib
#define PLASMA_NBNBSIZE    plasma->nbnbsize
#define PLASMA_IBNBSIZE    plasma->ibnbsize
#define PLASMA_SCHEDULING  plasma->scheduling
#define PLASMA_RHBLK       plasma->rhblock
#define PLASMA_TRANSLATION plasma->translation

/***************************************************************************//**
 *  IPT internal define
 **/
#define PlasmaIPT_NoDep   0
#define PlasmaIPT_Panel   1
#define PlasmaIPT_All     2

/***************************************************************************//**
 *  Global utilities
 **/
#ifndef max
#define max(a, b) ((a) > (b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif
#ifndef roundup
#define roundup(a, b) (b <= 0) ? (a) : (((a) + (b)-1) & ~((b)-1))
#endif

/***************************************************************************//**
 *  Static scheduler
 **/
#define ss_init(m, n, init_val) \
{ \
    if (PLASMA_RANK == 0) { \
        plasma->ss_progress = (volatile int *)plasma_shared_alloc(plasma, (m)*(n), PlasmaInteger); \
        plasma_memset_int((int*)plasma->ss_progress, (m)*(n), (init_val)); \
    } \
    plasma->ss_ld = (m); \
    plasma->ss_abort = 0; \
    plasma_barrier(plasma); \
}

#define ss_abort()   plasma->ss_abort = 1;
#define ss_aborted() plasma->ss_abort

#define ss_cond_set(m, n, val) \
{ \
    plasma->ss_progress[(m)+plasma->ss_ld*(n)] = (val); \
}

#define ss_cond_wait(m, n, val) \
{ \
    while (!plasma->ss_abort && \
            plasma->ss_progress[(m)+plasma->ss_ld*(n)] != (val)) \
        plasma_yield(); \
    if (plasma->ss_abort) \
        break; \
}

#define ss_finalize() \
{ \
    plasma_barrier(plasma); \
    if (PLASMA_RANK == 0) \
        plasma_shared_free(plasma, (void*)plasma->ss_progress); \
}

/***************************************************************************//**
 *  Global array of LAPACK constants
 **/
extern char *plasma_lapack_constants[];

#ifdef __cplusplus
extern "C" {
#endif

#include "compute_s.h"
#include "compute_d.h"
#define COMPLEX
#include "compute_c.h"
#include "compute_z.h"
#undef COMPLEX

void plasma_pdlag2s(plasma_context_t *plasma);
void plasma_pzlag2c(plasma_context_t *plasma);
void plasma_pslag2d(plasma_context_t *plasma);
void plasma_pclag2z(plasma_context_t *plasma);


#define QUARK_Insert_Task QUARK_Execute_Task

#ifdef __cplusplus
}
#endif

#endif
