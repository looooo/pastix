#include "common.h"
#include "dof.h"
#include "cost.h"
#include "ftgt.h"
#include "symbol.h"
#include "queue.h"
#include "bulles.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "queue.h"
#include "sopalin_acces.h"

#ifndef GPU_MAX_FILL
#  define GPU_MAX_FILL 0.7
#endif
#ifndef GPU_MIN_NBPAGES
#  define GPU_MIN_NBPAGES 0
#endif
#ifndef GPU_MIN_UPDATES
#  define GPU_MIN_UPDATES 0
#endif
#ifndef GPU_MIN_FLOP
#  define GPU_MIN_FLOP 0
#endif


/* Alter threadid of tasks to put some of them on GPUs (threadid >= thread_nbr)
 * cblk are taken from the queue until all GPUs are alocated maxMem memory.
 */
#define blend_distributeOnGPU PASTIX_PREFIX(blend_distributeOnGPU)
int blend_distributeOnGPU(SolverMatrix  * solvmtr,
                          double          maxMem,
                          int             pageSize,
                          int             criterium,
                          enum API_GPU_CRITERIUM nGPUs,
                          enum API_FLOAT  floatType,
                          enum API_FACT   factType);
