/**
 *
 * @file sopalin_zgetrf.c
 *
 *  PaStiX factorization routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "isched.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_PARSEC)
#include <dague.h>
#include <dague/data.h>
#include <dague/data_distribution.h>

int dsparse_zgetrf_1dplus_sp( dague_context_t *dague,
                              sparse_matrix_desc_t *A,
                              sopalin_data_t *sopalin_data );
int dsparse_zgetrf_2d_sp( dague_context_t *dague,
                          sparse_matrix_desc_t *A,
                          sopalin_data_t *sopalin_data );
#endif

void
sequential_zgetrf( pastix_data_t  *pastix_data,
                   sopalin_data_t *sopalin_data )
{
    SolverMatrix       *datacode = pastix_data->solvmatr;
    SolverCblk         *cblk;
    double              threshold = sopalin_data->diagthreshold;
    pastix_complex64_t *work;
    double              tol = datacode->lowrank.tolerance;
    pastix_int_t  i;
    (void)pastix_data;

    MALLOC_INTERN( work, datacode->gemmmax, pastix_complex64_t );

    cblk = datacode->cblktab;
    for (i=0; i<datacode->cblknbr; i++, cblk++){
        /* Compute */
        core_zgetrfsp1d( datacode, cblk, threshold, work, tol );
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( datacode, "getrf_L.txt" );
#endif

    memFree_null( work );
}

void
thread_pzgetrf( int rank, void *args )
{
    sopalin_data_t     *sopalin_data = (sopalin_data_t*)args;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    Task               *t;
    pastix_complex64_t *work;
    double              tol = datacode->lowrank.tolerance;
    pastix_int_t  i, ii;
    pastix_int_t  tasknbr, *tasktab;

    MALLOC_INTERN( work, datacode->gemmmax, pastix_complex64_t );

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        /* Wait */
        do {} while( cblk->ctrbcnt );

        /* Compute */
        core_zgetrfsp1d( datacode, cblk, sopalin_data->diagthreshold, work, tol );
    }

#if defined(PASTIX_DEBUG_FACTO) && 0
    isched_barrier_wait( &(((isched_t*)(sopalin_data->sched))->barrier) );
    if (rank == 0) {
        coeftab_zdump( datacode, "getrf_L.txt" );
    }
    isched_barrier_wait( &(((isched_t*)(sopalin_data->sched))->barrier) );
#endif

    memFree_null( work );
}

void
thread_zgetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    isched_parallel_call( pastix_data->isched, thread_pzgetrf, sopalin_data );
}

#if defined(PASTIX_WITH_PARSEC)
void
parsec_zgetrf_1dplus( pastix_data_t  *pastix_data,
                      sopalin_data_t *sopalin_data )
{
    dague_context_t *ctx;

    /* Start PaRSEC */
    if (pastix_data->parsec == NULL) {
        int argc = 0;
        pastix_parsec_init( pastix_data, &argc, NULL );
    }
    ctx = pastix_data->parsec;

    /* Run the facto */
    dsparse_zgetrf_1dplus_sp( ctx, sopalin_data->solvmtx->parsec_desc, sopalin_data );

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( sopalin_data->solvmtx, "getrf.txt" );
#endif

}

void
parsec_zgetrf_2d( pastix_data_t  *pastix_data,
                  sopalin_data_t *sopalin_data )
{
    dague_context_t *ctx;

    /* Start PaRSEC */
    if (pastix_data->parsec == NULL) {
        int argc = 0;
        pastix_parsec_init( pastix_data, &argc, NULL );
    }
    ctx = pastix_data->parsec;

    /* Run the facto */
    dsparse_zgetrf_2d_sp( ctx, sopalin_data->solvmtx->parsec_desc, sopalin_data );

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( sopalin_data->solvmtx, "getrf.txt" );
#endif

}
#endif

static void (*zgetrf_table[4])(pastix_data_t *, sopalin_data_t *) = {
    sequential_zgetrf,
    thread_zgetrf,
#if defined(PASTIX_WITH_PARSEC)
    parsec_zgetrf_2d,
    parsec_zgetrf_1dplus,
#else
    NULL,
    NULL
#endif
};

void
sopalin_zgetrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zgetrf)(pastix_data_t *, sopalin_data_t *) = zgetrf_table[ sched ];

    if (zgetrf == NULL) {
        zgetrf = thread_zgetrf;
    }
    zgetrf( pastix_data, sopalin_data );
}
