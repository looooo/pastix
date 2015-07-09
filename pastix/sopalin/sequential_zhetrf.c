/**
 *
 * @file sopalin_zhetrf.c
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
 * @precisions normal z -> c
 *
 **/
#include "common.h"
#include "isched.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"

void
sequential_zhetrf( sopalin_data_t *sopalin_data )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk   *cblk;
    pastix_int_t  i;

    cblk = datacode->cblktab;
    for (i=0; i<datacode->cblknbr; i++, cblk++){
        /* Compute */
        core_zhetrfsp1d( datacode, cblk, sopalin_data->diagthreshold );
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( datacode, "hetrf_L.txt" );
#endif
}


void
thread_pzhetrf( int rank, void *args )
{
    sopalin_data_t *sopalin_data = (sopalin_data_t*)args;
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk   *cblk;
    Task         *t;
    pastix_int_t  i, ii;
    pastix_int_t  tasknbr, *tasktab;

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        /* Compute */
        core_zhetrfsp1d( datacode, cblk, sopalin_data->diagthreshold );
    }

#if defined(PASTIX_DEBUG_FACTO)
    isched_barrier_wait( &(((isched_t*)(sopalin_data->sched))->barrier) );
    if (rank == 0) {
        coeftab_zdump( datacode, "hetrf_L.txt" );
    }
    isched_barrier_wait( &(((isched_t*)(sopalin_data->sched))->barrier) );
#endif
}


void
thread_zhetrf( sopalin_data_t *sopalin_data )
{
    isched_parallel_call( sopalin_data->sched, thread_pzhetrf, sopalin_data );
}

static void (*zhetrf_table[4])(sopalin_data_t *) = {
    sequential_zhetrf,
    thread_zhetrf,
    NULL,
    NULL
};

void
sopalin_zhetrf( sopalin_data_t *sopalin_data )
{
    int sched = 0;
    void (*zhetrf)(sopalin_data_t *) = zhetrf_table[ sched ];

    if (zhetrf == NULL) {
        zhetrf = thread_zhetrf;
    }
    zhetrf( sopalin_data );
}
