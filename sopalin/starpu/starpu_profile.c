/**
 *
 * @file starpu_profile.c
 *
 * @copyright 2017-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2023-07-21
 *
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define _GNU_SOURCE
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
#include "common.h"
#if !defined(PASTIX_WITH_STARPU)
#error "This file should not be compiled if Starpu is not enabled"
#endif
#include <stdio.h>
#include "pastix_starpu.h"

#define PROFILING_LOG_FILENAME "profile_log_%d.csv"
#define PROFILING_LOG_HEADERNAME "profile_log_header.csv"

#if defined( PASTIX_STARPU_PROFILING )

/**
 * @brief Chained list of measurement tables and their respective names.
 */
starpu_profile_t *profile_list = NULL;

/**
 *******************************************************************************
 *
 * @brief Registers a measurement table contained in a starpu_profile_t cell
 * into the profile_list.
 *
 *******************************************************************************
 *
 * @param[in] codelet
 *          TODO
 *
 ******************************************************************************/
void
profiling_register_cl( starpu_profile_t *codelet )
{
    assert( codelet->next == NULL );

    codelet->next = profile_list;
    profile_list = codelet;
}

/**
 *******************************************************************************
 *
 * @brief Adds the measures of a codelet into its designated measurement tables.
 *
 *******************************************************************************
 *
 * @param[in] callback_arg
 *          The argument given to the codelet containing at least a
 *          profile_data_t as first field.
 *
 ******************************************************************************/
void
cl_profiling_callback( void *callback_arg )
{
    struct starpu_task                *task = starpu_task_get_current();
    struct starpu_profiling_task_info *info = task->profiling_info;

    if ( info == NULL ) {
        return;
    }

    profile_data_t *arg      = (profile_data_t *) callback_arg;
    double          flops    = arg->flops;
    measure_t      *measures = arg->measures;
    double          duration = starpu_timing_timespec_delay_us( &info->start_time, &info->end_time );
    double          speed    = flops / ( 1000.0 * duration );

    measures[info->workerid].sum  += speed;
    measures[info->workerid].sum2 += speed * speed;
    measures[info->workerid].n    += 1;
}

/**
 *******************************************************************************
 *
 * @brief Display all collected profiling data of a given measurement table.
 *
 *******************************************************************************
 *
 * @param[in] codelet
 *          A measurement table to be displayed.
 *
 ******************************************************************************/
void
profiling_display_info( starpu_profile_t *codelet )
{
    measure_t *measures = codelet->measures;
    unsigned   worker;
    int        header = 0;

    for ( worker = 0; worker < starpu_worker_get_count(); worker++ ) {
        if ( measures[worker].n > 0 ) {
            if ( !header ) {
                printf( "Performance for kernel %s: \n", codelet->name );
                printf( "\tWorker  Gflop/s  delta  Nb\n" );
                header = 1;
            }
            char workername[128];
            starpu_worker_get_name( worker, workername, 128 );

            long   n    = measures[worker].n;
            double sum  = measures[worker].sum;
            double sum2 = measures[worker].sum2;

            double avg = sum / n;
            double sd  = sqrt((sum2 - (sum*sum)/n)/n);

            printf("\t%s\t%.2lf\t%.2lf\t%ld\n", workername, avg, sd, n);
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Displays all profiling data collected into all measurements tables of
 * the profile_list.
 *
 ******************************************************************************/
void
profiling_display_allinfo()
{
    starpu_profile_t *codelet = profile_list;
    while ( codelet != NULL ) {
        profiling_display_info( codelet );
        codelet = codelet->next;
    }
}

#endif /* defined( PASTIX_STARPU_PROFILING ) */

#if defined( PASTIX_STARPU_PROFILING_LOG )

FILE *profiling_log_file[STARPU_NMAXWORKERS];

/**
 *******************************************************************************
 *
 * @brief Initializes the profiling_log functions.
 *
 * The function profiling_log_fini() must be called to end the profiling.
 * Opens a series of FILE* into the given directory for concurrent logging of
 * all individual tasks.
 *
 *******************************************************************************
 *
 * @param[in] dirname
 *          The path to the directory in which the profiling_log will be written.
 *
 ******************************************************************************/
void
profiling_log_init( const char *dirname )
{
    char *tmp;
    int   index = 0;
    for ( ; index < STARPU_NMAXWORKERS; index++ ) {
        asprintf( &tmp, PROFILING_LOG_FILENAME, index );
        profiling_log_file[index] = pastix_fopenw( dirname, tmp, "w" );
        assert( profiling_log_file[index] );
        free( tmp );
    }
    FILE* header = pastix_fopenw( dirname, PROFILING_LOG_HEADERNAME, "w" );
    assert( header );
    fprintf( header, "task_name;cl_name;m;n;k;flops;speed\n" );
    fclose( header );
}

/**
 *******************************************************************************
 *
 * @brief Registers a full observation of an individual task.
 * Concurrent call to this function can be achieved on different STARPU_WORKERS.
 *
 *******************************************************************************
 *
 * @param[in] task_name
 *          The unique name of the task.
 *
 * @param[in] cl_name
 *          The name given to the codelet.
 *
 * @param[in] m
 *          The m parameter of the kernel (used by xxTRF, TRSM, and GEMM).
 *
 * @param[in] n
 *          The n parameter of the kernel (used by TRSM, and GEMM).
 *
 * @param[in] k
 *          The k parameter of the kernel (used by GEMM).
 *
 * @param[in] flops
 *          The number of flops of the kernel.
 *
 * @param[in] speed
 *          The speed of computation achieved by the task (in GFLop/s).
 *
 ******************************************************************************/
void
cl_profiling_log_register( const char *task_name,
                           const char *cl_name,
                           int         m,
                           int         n,
                           int         k,
                           double      flops,
                           double      speed )
{
    struct starpu_task                *task = starpu_task_get_current();
    struct starpu_profiling_task_info *info = task->profiling_info;
    fprintf( profiling_log_file[info->workerid], "%s;%s;%d;%d;%d;%lf;%lf\n", task_name, cl_name, m, n, k, flops, speed );
}

/**
 *******************************************************************************
 *
 * @brief Terminates the profiling_log life cycle.
 *
 ******************************************************************************/
void
profiling_log_fini()
{
    int index = 0;
    for ( ; index < STARPU_NMAXWORKERS; index++ ) {
        fclose( profiling_log_file[index] );
    }
}

#endif /* defined( PASTIX_STARPU_PROFILING_LOG ) */

/**
 * @}
 */
