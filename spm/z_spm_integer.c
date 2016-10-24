/**
 *
 * @file z_spm_integers.c
 *
 *  PaStiX spm routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 * @precisions normal z -> c d s
 *
 **/
#include <ctype.h>
#include <limits.h>
#include <time.h>
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Sort 2 arrays simultaneously, the first array is an array of
 * pastix_int_t and used as key for sorting.  The second array is an array of
 * pastix_complex64_t.
 *
 *******************************************************************************
 */
static size_t intsortsize[2] = { sizeof(pastix_int_t), sizeof(pastix_complex64_t) };
#define INTSORTNAME            z_spmIntSortAsc
#define INTSORTSIZE(x)         (intsortsize[x])
#define INTSORTNTAB            2
#define INTSORTSWAP(p,q)       do {					\
    pastix_int_t     t;								\
    long    disp_p   = (((pastix_int_t*)p)-((pastix_int_t*)base_ptr));			\
    long    disp_q   = (((pastix_int_t*)q)-((pastix_int_t*)base_ptr));			\
    pastix_complex64_t * floatptr = *(pbase+1);					\
    pastix_complex64_t   f;								\
    /* swap integers */							\
    t = *((pastix_int_t *) (p));							\
    *((pastix_int_t *) (p)) = *((pastix_int_t *) (q));					\
    *((pastix_int_t *) (q)) = t;							\
    /* swap corresponding values */					\
    f = floatptr[disp_p];						\
    floatptr[disp_p] = floatptr[disp_q];				\
    floatptr[disp_q] = f;						\
  } while (0)
#define INTSORTCMP(p,q)             (*((pastix_int_t *) (p)) < *((pastix_int_t *) (q)))
#include "integer_sort_mtypes.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB

