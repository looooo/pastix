/**
 *
 * @file memory.h
 *
 * @copyright 1998-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * PaStiX memory tracking function.
 *
 * @version 6.4.0
 * @author Francois Pellegrini
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2024-07-05
 *
 */
/*
 *  File: memory.h
 *
 *  Part of a parallel direct block solver.
 *
 *  These lines are the common data
 *  declarations for all modules.
 *
 *  Authors:
 *    Mathieu  Faverge    - faverge@labri.fr
 *    Xavier   LACOSTE    - lacoste@labri.fr
 *    Pierre   RAMET      - ramet@labri.fr
 *
 *  Dates:
 *    Version 0.0 - from 08 may 1998
 *                  to   08 jan 2001
 *    Version 1.0 - from 06 jun 2002
 *                  to   06 jun 2002
 */
#ifndef _memory_h_
#define _memory_h_

#if defined(PASTIX_WITH_CUDA)
#include <cuda_runtime.h>
#endif

/*
 * Function: pastix_protected_malloc
 *
 * PowerPC architectures don't support malloc(0). This function
 * prints a warning when it happens to avoid segfault.
 */
static inline void *pastix_malloc_func( size_t size,
                                        char *filename,
                                        int line )
{
    if (size > 0) {
	return malloc(size);
    }
    else {
	fprintf(stderr, "Pb Alloc 0 %s:%d\n", filename, line);
	return (void *)NULL;
    }
}

#if defined(PASTIX_ARCH_PPC)
#  define memAlloc(size) pastix_malloc_func(size, __FILE__, __LINE__)
#else
#  define memAlloc(size) malloc(size)
#endif

#define memFree(ptr) free((void*)(ptr))
#define memFree_null(ptr) do			\
	{					\
	    memFree( ptr );			\
	    (ptr) = NULL;			\
	} while(0)

#define MALLOC_INTERN(ptr, size, type)                          \
    do {                                                        \
        ptr = (type*)memAlloc((size) * sizeof(type)) ;          \
    } while(0)

#define MALLOC_EXTERN(ptr, size, type)		\
    ptr = (type*)malloc((size) * sizeof(type))

#define __MALLOC_ERROR( _str_, _fname, _line )                                           \
    {                                                                   \
        fprintf(stderr, "%s allocation (line=%d,file=%s)\n",(_str_),(_line),(_fname)); \
        exit(-1);                                                       \
    }

#define MALLOC_ERROR( _str_ ) __MALLOC_ERROR( _str_, __FILE__, __LINE__ )

/*
 * Macro: MALLOC_INTOREXTERN
 *
 * Choose between <MALLOC_INTERN> and <MALLOC_EXTERN>
 * following flag_int.
 *
 * Parameters:
 *   ptr      - address where to allocate.
 *   size     - Number of elements to allocate.
 *   types    - Type of the elements to allocate.
 *   flag_int - 1 for internal allocation, 0 for external.
 */
#define MALLOC_INTOREXTERN(ptr, size, type, flag_int) \
  do {                                                \
    if (flag_int == 1)                          \
      {                                               \
        MALLOC_INTERN(ptr, size, type);               \
      }                                               \
    else                                              \
      {                                               \
        MALLOC_EXTERN(ptr, size, type);               \
      }                                               \
  } while (0)

#define FREE_NULL_INTOREXT(ptr, flag_int)         \
  do {                                            \
    if (flag_int == 1)                      \
      {                                           \
        memFree_null(ptr);                        \
      }                                           \
    else                                          \
      {                                           \
        free(ptr);                                \
        ptr = NULL;                               \
      }                                           \
  } while (0)

#define memRealloc realloc

/**
 *******************************************************************************
 *
 * @brief Allocates pinned memory via cuda and normal memory otherwise.
 *
 *******************************************************************************
 *
 * @param[in] size
 *          The number of bytes to allocate.
 *
 * @param[in] fname
 *          The name of the file, in order to the malloc error.
 *
 * @param[in] line
 *          The number of the line, in order to the malloc error.
 *
 *******************************************************************************/
static inline void *
__pastix_malloc_pinned( size_t size, const char* fname, int line )
{
    void* ptr;
#if defined(PASTIX_WITH_CUDA)
    int rc;
    rc = cudaMallocHost( &ptr, size );
    if( rc )
    {
        __MALLOC_ERROR( "error with pinned-memory", fname, line );
    }
#else
    ptr = malloc( size );
#endif
    (void)fname;
    (void)line;
    return ptr;
}

#define pastix_malloc_pinned( _size_ ) __pastix_malloc_pinned( _size_, __FILE__, __LINE__ )

/**
 *******************************************************************************
 *
 * @brief Frees and sets to NULL memory allocated with pastix_malloc_pinned.
 *
 *******************************************************************************
 *
 * @param[inout] ptr
 *          The address to free.
 *
 *******************************************************************************/
static inline void
__pastix_free_pinned( void *ptr )
{
#if defined(PASTIX_WITH_CUDA)
    cudaFree( ptr );
#else
    free( ptr );
#endif
}

#define pastix_free_pinned( _ptr_ ) do	\
    {                                   \
        __pastix_free_pinned( _ptr_ );  \
        (_ptr_) = NULL;                 \
    } while(0)

#endif /* _memory_h_ */
