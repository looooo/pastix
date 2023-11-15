/**
 *
 * @file common.h
 *
 * @copyright 2004-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.1
 * @author David Goudin
 * @author Gregoire Pichon
 * @author Mathieu Faverge
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Tony Delarue
 * @author Brieuc Nicolas
 * @author Mohamed Aymane Kherraz
 * @date 2023-08-03
 *
 **/
#ifndef _common_h_
#define _common_h_

#include "pastix.h"
#include "pastix/order.h"
#if !defined(PASTIX_WITH_MPI)
#include "nompi.h"
#endif
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <limits.h>
#include <sched.h>

#if defined(HAVE_MM_SETCSR)
#include <immintrin.h>
#endif
#if defined(HAVE_FTZ_MACROS)
#include <xmmintrin.h>
#endif
#if defined(HAVE_DAZ_MACROS)
#include <pmmintrin.h>
#endif

#if defined(PASTIX_DEBUG_VALGRIND)
#define pastix_yield() sched_yield()
#else
#define pastix_yield() do {} while (0)
#endif

#include "sys/atomic.h"
#include "memory.h"
#include "integer.h"
#include "timing.h"
#include "pastixdata.h"
#include "out.h"
#include "parse_options.h"
#include "pastix_papi.h"

#if defined(HAVE_BUILTIN_EXPECT)
#define pastix_likely( _x_ )   __builtin_expect( (_x_), 1 )
#define pastix_unlikely( _x_ ) __builtin_expect( (_x_), 0 )
#else
#define pastix_likely( _x_ )   ( _x_ )
#define pastix_unlikely( _x_ ) ( _x_ )
#endif

#if defined(PASTIX_OS_WINDOWS)
#include <windows.h>
#define COMMON_RANDOM_RAND 1
#endif

/********************************************************************
 * CBLAS value address
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( a_ ) (&(a_))
#endif

/*
 * Get environment variable
 */
#if defined(PASTIX_OS_WINDOWS)

static inline int
pastix_setenv( const char *var, const char *value, int overwrite ) {
    return !(SetEnvironmentVariable( var, value ));
}

static inline char *
pastix_getenv( const char *var ) {
    char *str;
    int len = 512;
    int rc;
    str = (char*)malloc(len * sizeof(char));
    rc = GetEnvironmentVariable(var, str, len);
    if (rc == 0) {
        free(str);
        str = NULL;
    }
    return str;
}

static inline void
pastix_cleanenv( char *str ) {
    if (str != NULL) free(str);
}

#else /* Other OS systems */

static inline int
pastix_setenv( const char *var, const char *value, int overwrite ) {
    return setenv( var, value, overwrite );
}

static inline char *
pastix_getenv( const char *var ) {
    return getenv( var );
}

static inline void
pastix_cleanenv( char *str ) {
    (void)str;
}

#endif


static inline int
pastix_env_is_set_to(char * str, char * value) {
    char * val;
    if ( (val = pastix_getenv(str)) &&
         !strcmp(val, value))
        return 1;
    return 0;
}

static inline int
pastix_env_is_on(char * str) {
    return pastix_env_is_set_to(str, "1");
}

static inline int
pastix_getenv_get_value_int( char *string, int default_value )
{
    long int ret;
    char *str = pastix_getenv(string);
    if (str == NULL) return default_value;

    if ( sscanf( str, "%ld", &ret ) != 1 ) {
        perror("sscanf");
        return default_value;
    }

    return (int)ret;
}

static inline char *
pastix_getenv_get_value_str( char *varname, const char *default_value )
{
    char *str = pastix_getenv( varname );
    if ( str == NULL ) {
        return strdup( default_value );
    }
    else {
        return strdup( str );
    }
}

/* **************************************** */

static inline void set_iparm(pastix_int_t *iparm, pastix_iparm_t offset, pastix_int_t value)
{
    if (iparm != NULL) iparm[offset] = (pastix_int_t)value;
}

static inline void set_dparm(double *dparm, pastix_dparm_t offset, double value)
{
    if (dparm != NULL) dparm[offset] = (double)value;
}

/*
 * Check if the scheduler have been changed between two steps.
 * If it's the case, check if the new scheduler is in the same family
 * than the previous one.
 * If not, rehabilitate the previous scheduler.
 */
static inline void
pastix_check_and_correct_scheduler( pastix_data_t *pastix_data )
{
    if ( pastix_data->inter_node_procnbr == 1 ) {
        return;
    }
    pastix_int_t *isched = pastix_data->iparm + IPARM_SCHEDULER;
    pastix_int_t  sched  = pastix_data->sched;

    if( (isSchedRuntime(*isched) && (pastix_data->solvmatr != pastix_data->solvglob)) ||
        (isSchedPthread(*isched) && (pastix_data->solvmatr != pastix_data->solvloc )) )
    {
        pastix_print_warning( "Scheduler can't be changed to %s, restore %s scheduler\n",
                              pastix_scheduler_getstr( *isched ), pastix_scheduler_getstr( sched ) );
        *isched = pastix_data->sched;
    }
    /* Backup latest scheduler */
    pastix_data->sched = *isched;
}

void api_dumparm(FILE *stream, pastix_int_t *iparm, double *dparm);

#if !defined(HAVE_GETLINE)
ssize_t getdelim(char **buf, size_t *bufsiz, int delimiter, FILE *fp);
ssize_t getline(char **buf, size_t *bufsiz, FILE *fp);
#endif

/*
 * Set the Flush to Zero bit to 1 so that there is no slow down due to denormals
 */
static inline void
set_ftz( void )
{

/**
 * The following x86 asm solution is proposed in :
 * Mawussi Zounon, Nicholas J. Higham, Craig Lucas, and Françoise Tisseur.
 * "Performance Evaluation of Mixed Precision Algorithms for Solving Sparse Linear Systems"
 * https://www.semanticscholar.org/paper/Performance-Evaluation-of-Mixed-Precision-for-Zounon-Higham/44731cd2cd93aca8299a83ee04029deb7641b1a3
 */
#if defined(__x86_64)
#if defined(HAVE_FTZ_MACROS) && defined(HAVE_DAZ_MACROS)
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#elif defined(HAVE_MM_SET_CSR)
    _mm_setcsr(_mm_getcsr() & ~0x8040);
#else
    asm( "stmxcsr -0x4(%rsp)\n\t"          /* store CSR register on stack */
         "orl     $0x8040, -0x4(%rsp)\n\t" /* set bits 15( FTZ ) and 7( DAZ ) */
         "ldmxcsr -0x4(%rsp)" );           /* load CSR register from stack */
#endif
#endif

#if defined(__arm64__) || defined(__aarch64__)
    uint64_t fpcr;
    asm( "mrs %0,   fpcr" : "=r"( fpcr ));
    asm( "msr fpcr, %0"   :: "r"( fpcr | (1 << 24) ));
#endif

}

#endif /* _common_h_ */
