/**
 *
 * @file models.h
 *
 *  PaStiX sopalin routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @date 2017-02-09
 *
 **/
#ifndef _MODELS_H_
#define _MODELS_H_

#include "common/timing.h"

typedef enum pastix_ktype_e {
    PastixKernelCpuGETRF,
    PastixKernelCpuHETRF,
    PastixKernelCpuPOTRF,
    PastixKernelCpuPXTRF,
    PastixKernelCpuSYTRF,
    PastixKernelCpuTRSM,
    PastixKernelCpuGEMM1D,
    PastixKernelCpuGEMM2D,
    PastixKernelGpuGEMM1D,
    PastixKernelGpuGEMM2D
} pastix_ktype_t;

typedef struct pastix_model_entry_s {
    pastix_ktype_t ktype;
    int m;
    int n;
    int k;
    double time;
} pastix_model_entry_t;

extern pastix_model_entry_t *model_entries;     /*< Array to all entries                 */
extern volatile int32_t      model_entries_nbr; /*< Index of the last entry in the array */
extern int32_t               model_size;        /*< Allocated size of the array          */

#if defined(PASTIX_GENERATE_MODEL)
static inline double
modelGetTime()
{
    return clockGet();
}

static inline void
modelAddEntry( pastix_ktype_t ktype, int m, int n, int k, double starttime )
{
    double  time = clockGet() - starttime;
    int32_t index = pastix_atomic_inc_32b( &model_entries_nbr );

    if ( index < model_size ) {
        model_entries[index].ktype = ktype;
        model_entries[index].m     = m;
        model_entries[index].n     = n;
        model_entries[index].k     = k;
        model_entries[index].time  = time;
    }
    else {
        fprintf(stderr, "WARNING: too many kernels to log %d\n", index);
        //pastix_atomic_dec_32b( &model_entries_nbr );
    }
}

void modelInit( const SolverMatrix *solvmtx );
void modelDumpAndExit( char **directory );

#else
static inline double
modelGetTime() { return 0; }

static inline void
modelAddEntry( pastix_ktype_t ktype, int m, int n, int k, double time )
{
    (void)ktype; (void)m; (void)n; (void)k; (void)time;
}

static inline void
modelInit( const SolverMatrix *solvmtx)
{ (void)solvmtx; };

static inline void
modelDumpAndExit( char **directory )
{ (void)directory; }

#endif

#endif /* _MODELS_H_ */
