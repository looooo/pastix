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
static inline void *pastix_protected_malloc( size_t size,
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

#define memFree(ptr) free(ptr)
#define memFree_null(ptr) do			\
	{					\
	    memFree( ptr );			\
	    (ptr) = NULL;			\
	} while(0)

/* #ifdef MEMORY_USAGE */
/* void *         memAlloc_func       (size_t,char*,int); */
/* void *         memRealloc_func     (void *, size_t, char*, int); */
/* #  define memRealloc(ptr,size)                          \ */
/*   memRealloc_func(ptr, size, __FILE__, __LINE__) */

/* void           memFree             (void *); */
/* unsigned long  memAllocGetCurrent  (void); */
/* unsigned long  memAllocGetMax      (void); */
/* void           memAllocTraceReset  (void); */
/* #else */
/* void *         mymalloc            (size_t,char*,int); */
/* #endif /\* MEMORY_USAGE *\/ */
/* void *         memAllocGroup       (void **, ...); */
/* void *         memReallocGroup     (void *, ...); */
/* void *         memOffset           (void *, ...); */

#define MALLOC_INTERN(ptr, size, type)		\
    ptr = (type*)memAlloc((size) * sizeof(type))


#endif /* _MEMORY_H_ */
