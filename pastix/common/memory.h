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
#ifndef _MEMORY_H_
#define _MEMORY_H_

/*
 * Macro: MEMORY_WRITE
 *
 * Adapt the integer given in parameter to the suitable multiple of
 * 2**10 for output.
 */
#define MEMORY_WRITE(mem) ( ((mem) < 1<<10) ?                           \
                            ( (double)(mem) ) :                         \
                            ( ( (mem) < 1<<20 ) ?                       \
                              ( (double)(mem)/(double)(1<<10) ) :       \
                              ( ((mem) < 1<<30 ) ?                      \
                                ( (double)(mem)/(double)(1<<20) ) :     \
                                ( (double)(mem)/(double)(1<<30) ))))
/*
 * Macro: MEMORY_UNIT_WRITE
 *
 * Return the unit adapted to the amount of memory represented by mem.
 */
#define MEMORY_UNIT_WRITE(mem) (((mem) < 1<<10) ?                       \
                                "o" :                                   \
                                ( ( (mem) < 1<<20 ) ?                   \
                                  "Ko" :                                \
                                  ( ( (mem) < 1<<30 ) ?                 \
                                    "Mo" :                              \
                                    "Go" )))

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

#if defined(ARCH_PPC)
#  define memAlloc(size) pastix_protected_malloc(size, __FILE__, __LINE__)
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

#define MALLOC_ERROR( _str_ )                                           \
    {                                                                   \
        fprintf(stderr, "%s allocation (line=%d,file=%s)\n",(_str_),__LINE__,__FILE__); \
        exit(-1);                                                       \
    }

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
 *   flag_int - API_YES for internal allocation, API_NO for external.
 */
#define MALLOC_INTOREXTERN(ptr, size, type, flag_int) \
  do {                                                \
    if (flag_int == API_YES)                          \
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
    if (flag_int == API_YES)                      \
      {                                           \
        memFree_null(ptr);                        \
      }                                           \
    else                                          \
      {                                           \
        free(ptr);                                \
        ptr = NULL;                               \
      }                                           \
  } while (0)

//void *memAlloc(size_t size);
//void  memFree(void *ptr);
//void *memRealloc(void *ptr, size_t size);
#define memRealloc realloc
void *memAllocGroup  (void **, ...);
void *memReallocGroup(void *, ...);

#endif /* _MEMORY_H_ */
