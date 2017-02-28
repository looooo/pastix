/**
 *
 * @file spm_integers.c
 *
 *  PaStiX spm routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Francois Pellegrini
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include <ctype.h>
#include <limits.h>
#include <time.h>
#include <stdlib.h>
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm
 *
 * @brief Convert integer array to the pastix_int_t format if it is not already
 * the case.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of elements in the array.
 *
 * @param[in,out] input
 *          The input array. If the type size is not the same, the array is
 *          freed on exit.
 *
 *******************************************************************************
 *
 * @return The pointer to the new allocated array if size has changed, or to
 *         input if sizes are identical.
 *
 *******************************************************************************/
pastix_int_t *
spmIntConvert( pastix_int_t n, int *input )
{
    if (sizeof(pastix_int_t) != sizeof(int)) {
        pastix_int_t *output, *tmpo;
        int *tmpi, i;

        output = malloc( n * sizeof(pastix_int_t) );
        tmpi = input;
        tmpo = output;
        for(i=0; i<n; i++, tmpi++, tmpo++) {
            *tmpo = (pastix_int_t)(*tmpi);
        }
        free(input);
        return output;
    }
    else {
        return (pastix_int_t*)input;
    }
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_dev
 *
 * @brief Sorts in ascending order array of element composed of one single
 * pastix_int_t with a single key value.
 *
 *******************************************************************************
 */
#define INTSORTNAME                 spmIntSort1Asc1
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define INTSORTSIZE                 (sizeof (pastix_int_t))
#define INTSORTSWAP(p,q)            do {			\
        pastix_int_t t;						\
        t = *((pastix_int_t *) (p));				\
        *((pastix_int_t *) (p)) = *((pastix_int_t *) (q));	\
        *((pastix_int_t *) (q)) = t;				\
    } while (0)
#define INTSORTCMP(p,q)             (*((pastix_int_t *) (p)) < *((pastix_int_t *) (q)))
#include "integer_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_dev
 *
 * @brief Sorts in ascending order array of element composed of two
 * pastix_int_t by ascending order. The first value is used as key.
 *
 *******************************************************************************
 */
#define INTSORTNAME                 spmIntSort2Asc1
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define INTSORTSIZE                 (2 * sizeof (pastix_int_t))
#define INTSORTSWAP(p,q)            do {				\
        pastix_int_t t, u;						\
        t = *((pastix_int_t *) (p));					\
        u = *((pastix_int_t *) (p) + 1);				\
        *((pastix_int_t *) (p)) = *((pastix_int_t *) (q));		\
        *((pastix_int_t *) (p) + 1) = *((pastix_int_t *) (q) + 1);	\
        *((pastix_int_t *) (q)) = t;					\
        *((pastix_int_t *) (q) + 1) = u;				\
    } while (0)
#define INTSORTCMP(p,q)             (*((pastix_int_t *) (p)) < *((pastix_int_t *) (q)))
#include "integer_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_dev
 *
 * @brief Sorts in ascending order array of element composed of three
 * pastix_int_t by ascending order. The first value is used as key.
 *
 *******************************************************************************
 */
#define INTSORTNAME                 spmInt_intSort3Asc1
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define INTSORTSIZE                 (3 * sizeof (pastix_int_t))
#define INTSORTSWAP(p,q)            do {				\
        pastix_int_t t, u, v;						\
        t = *((pastix_int_t *) (p));					\
        u = *((pastix_int_t *) (p) + 1);				\
        v = *((pastix_int_t *) (p) + 2);				\
        *((pastix_int_t *) (p)) = *((pastix_int_t *) (q));		\
        *((pastix_int_t *) (p) + 1) = *((pastix_int_t *) (q) + 1);	\
        *((pastix_int_t *) (p) + 2) = *((pastix_int_t *) (q) + 2);	\
        *((pastix_int_t *) (q)) = t;					\
        *((pastix_int_t *) (q) + 1) = u;				\
        *((pastix_int_t *) (q) + 2) = v;				\
    } while (0)
#define INTSORTCMP(p,q)             (*((pastix_int_t *) (p)) < *((pastix_int_t *) (q)))
#include "integer_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_dev
 *
 * @brief Sorts in ascending order array of element composed of two
 * pastix_int_t by ascending order. Both values are used as key.
 *
 *******************************************************************************
 */
#define INTSORTNAME                 spmIntSort2Asc2
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define INTSORTSIZE                 (2 * sizeof (pastix_int_t))
#define INTSORTSWAP(p,q)            do {				\
        pastix_int_t t, u;						\
        t = *((pastix_int_t *) (p));					\
        u = *((pastix_int_t *) (p) + 1);				\
        *((pastix_int_t *) (p)) = *((pastix_int_t *) (q));		\
        *((pastix_int_t *) (p) + 1) = *((pastix_int_t *) (q) + 1);	\
        *((pastix_int_t *) (q)) = t;					\
        *((pastix_int_t *) (q) + 1) = u;				\
    } while (0)
#define INTSORTCMP(p,q)             ((*((pastix_int_t *) (p)) < *((pastix_int_t *) (q))) || ((*((pastix_int_t *) (p)) == *((pastix_int_t *) (q))) && (*((pastix_int_t *) (p) + 1) < *((pastix_int_t *) (q) + 1))))
#include "integer_sort.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_dev
 *
 * @brief Sort 2 arrays simultaneously, the first array is an array of
 * pastix_int_t and used as primary key for sorting.  The second array is an
 * other array of pastix_int_t used as secondary key.
 *
 *******************************************************************************
 */
#define INTSORTNAME            spmIntMSortIntAsc
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define INTSORTSIZE(x)         (sizeof (pastix_int_t))
#define INTSORTNTAB            2
#define INTSORTSWAP(p,q)       do {                                     \
        pastix_int_t     t;                                             \
        long    disp_p   = (((pastix_int_t*)p)-((pastix_int_t*)base_ptr)); \
        long    disp_q   = (((pastix_int_t*)q)-((pastix_int_t*)base_ptr)); \
        pastix_int_t   * int2ptr  = *(pbase+1);                         \
        /* swap integers */                                             \
        t = *((pastix_int_t *) (p));                                    \
        *((pastix_int_t *) (p)) = *((pastix_int_t *) (q));              \
        *((pastix_int_t *) (q)) = t;                                    \
        /* swap on secont integer array */                              \
        t = int2ptr[disp_p];                                            \
        int2ptr[disp_p] = int2ptr[disp_q];                              \
        int2ptr[disp_q] = t;                                            \
    } while (0)
#define INTSORTCMP(p,q)  ((*((pastix_int_t *) (p)) < *((pastix_int_t *) (q))) || \
                          ((*((pastix_int_t *) (p)) == *((pastix_int_t *) (q))) && \
                           ((( pastix_int_t *)(*(pbase+1)))[(((pastix_int_t*)p)-((pastix_int_t*)base_ptr))] < \
                            (( pastix_int_t *)(*(pbase+1)))[(((pastix_int_t*)q)-((pastix_int_t*)base_ptr))])))
#include "integer_sort_mtypes.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @ingroup pastix_spm_dev
 *
 * @brief Sort 2 arrays simultaneously, the first array is an array of
 * pastix_int_t and used as primary key for sorting.  The second array is an
 * other array of pastix_int_t used as secondary key.
 *
 *******************************************************************************
 */
#define INTSORTNAME            spmIntMSortSmallIntAsc
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define INTSORTSIZE(x)         (sizeof (int))
#define INTSORTNTAB            2
#define INTSORTSWAP(p,q)       do {                             \
        int     t;                                              \
        long    disp_p   = (((int*)p)-((int*)base_ptr));        \
        long    disp_q   = (((int*)q)-((int*)base_ptr));        \
        int   * int2ptr  = *(pbase+1);                          \
        /* swap integers */                                     \
        t = *((int *) (p));                                     \
        *((int *) (p)) = *((int *) (q));                        \
        *((int *) (q)) = t;                                     \
        /* swap on secont integer array */                      \
        t = int2ptr[disp_p];                                    \
        int2ptr[disp_p] = int2ptr[disp_q];                      \
        int2ptr[disp_q] = t;                                    \
        } while (0)
#define INTSORTCMP(p,q)  ((*((int *) (p)) < *((int *) (q))) ||		\
                          ((*((int *) (p)) == *((int *) (q))) &&	\
                           ((( int *)(*(pbase+1)))[(((int*)p)-((int*)base_ptr))] < \
                            (( int *)(*(pbase+1)))[(((int*)q)-((int*)base_ptr))])))
#include "integer_sort_mtypes.c"
#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP
#undef INTSORTNTAB
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
