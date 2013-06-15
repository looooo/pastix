#ifndef _PASTIX_CONFIG_H_
#define _PASTIX_CONFIG_H_

#define PASTIX_VERSION_MAJOR 5
#define PASTIX_VERSION_MINOR 1
#define PASTIX_VERSION_MICRO 0

/* #undef PASTIX_WITH_MPI */
/* #undef PASTIX_WITH_CUDA */

/* system */
#define HAVE_PTHREAD
#define HAVE_SCHED_SETAFFINITY
#define HAVE_CLOCK_GETTIME
/* #undef HAVE_ASPRINTF */
/* #undef HAVE_VASPRINTF */
#define HAVE_STDARG_H
#define HAVE_UNISTD_H
/* #undef HAVE_VA_COPY */
/* #undef HAVE_UNDERSCORE_VA_COPY */
#define HAVE_GETOPT_LONG
#define HAVE_GETRUSAGE
#define HAVE_GETOPT_H
#define HAVE_ERRNO_H
#define HAVE_STDDEF_H
#define HAVE_LIMITS_H
#define HAVE_STRING_H
#define HAVE_COMPLEX_H

/* #undef ARCH_X86 */
#define ARCH_X86_64
/* #undef ARCH_PPC */
/* #undef MAC_OS_X */

/* Optional packages */
#define HAVE_HWLOC
#define HAVE_HWLOC_BITMAP
#define HAVE_HWLOC_PARENT_MEMBER
#define HAVE_HWLOC_CACHE_ATTR
#define HAVE_HWLOC_OBJ_PU

/* Ordering options */
#define PASTIX_ORDERING_SCOTCH
/* #undef PASTIX_ORDERING_METIS */

/* Datatypes used */
#define PASTIX_INT64

#if defined(PASTIX_INT64)
#define FORCE_INT64
#define INTSSIZE64
#endif

#define FORCE_DOUBLE
#define PREC_DOUBLE

#if defined(PASTIX_WITH_MPI)
#define HAVE_MPI
#undef  FORCE_NOMPI
#else
#undef  HAVE_MPI
#define FORCE_NOMPI
#endif

#if defined(PASTIX_ORDERING_SCOTCH)
#define HAVE_SCOTCH
#define WITH_SCOTCH
#endif

#endif  /* _PASTIX_CONFIG_H_ */
