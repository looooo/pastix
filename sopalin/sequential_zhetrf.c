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

#if defined(PASTIX_WITH_PARSEC)
#include <dague.h>
#include <dague/data.h>
#include <dague/data_distribution.h>
#include "parsec/sparse-matrix.h"
#endif

void
sequential_zhetrf( pastix_data_t  *pastix_data,
                   sopalin_data_t *sopalin_data )
{
    SolverMatrix       *datacode = pastix_data->solvmatr;
    SolverCblk         *cblk;
    double              threshold = sopalin_data->diagthreshold;
    pastix_complex64_t *work1, *work2;
    pastix_int_t  i;
    (void)pastix_data;

    MALLOC_INTERN( work1, pastix_imax(datacode->gemmmax, datacode->diagmax),
                   pastix_complex64_t );
    MALLOC_INTERN( work2, datacode->gemmmax, pastix_complex64_t );

    cblk = datacode->cblktab;
    for (i=0; i<datacode->cblknbr; i++, cblk++){
        /* Compute */
        core_zhetrfsp1d( datacode, cblk, threshold,
                         work1, work2 );
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( datacode, "hetrf_L.txt" );
#endif

    memFree_null( work1 );
    memFree_null( work2 );
}

void
thread_pzhetrf( int rank, void *args )
{
    sopalin_data_t     *sopalin_data = (sopalin_data_t*)args;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    Task               *t;
    pastix_complex64_t *work1, *work2;
    pastix_int_t  i, ii;
    pastix_int_t  tasknbr, *tasktab;

    MALLOC_INTERN( work1, pastix_imax(datacode->gemmmax, datacode->diagmax),
                   pastix_complex64_t );
    MALLOC_INTERN( work2, datacode->gemmmax, pastix_complex64_t );

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        /* Compute */
        core_zhetrfsp1d( datacode, cblk, sopalin_data->diagthreshold,
                         work1, work2 );
    }

#if defined(PASTIX_DEBUG_FACTO) && 0
    isched_barrier_wait( &(((isched_t*)(sopalin_data->sched))->barrier) );
    if (rank == 0) {
        coeftab_zdump( datacode, "hetrf_L.txt" );
    }
    isched_barrier_wait( &(((isched_t*)(sopalin_data->sched))->barrier) );
#endif

    memFree_null( work1 );
    memFree_null( work2 );
}


void
thread_zhetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    isched_parallel_call( pastix_data->isched, thread_pzhetrf, sopalin_data );
}

#if defined(PASTIX_WITH_PARSEC) && 0
void
parsec_zhetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    sparse_matrix_desc_t desc;
    dague_context_t *ctx;
    int argc = 0;

    /* Start PaRSEC */
    if (pastix_data->parsec == NULL) {
        pastix_data->parsec = dague_init( -1, &argc, NULL );
    }
    ctx = pastix_data->parsec;

    /* Create the matrix descriptor */
    sparse_matrix_init( &desc, sopalin_data->solvmtx,
                        pastix_size_of( PastixComplex64 ), 1, 0 );

    /* Run the facto */
    dsparse_zhetrf_sp( ctx, &desc, sopalin_data );

    /* Destroy the decriptor */
    sparse_matrix_destroy( &desc );

    dague_fini( &(pastix_data->parsec) );
}
#endif

static void (*zhetrf_table[4])(pastix_data_t *, sopalin_data_t *) = {
    sequential_zhetrf,
    thread_zhetrf,
#if defined(PASTIX_WITH_PARSEC) && 0
    parsec_zhetrf,
#else
    NULL,
#endif
    NULL
};

void
sopalin_zhetrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zhetrf)(pastix_data_t *, sopalin_data_t *) = zhetrf_table[ sched ];

    if (zhetrf == NULL) {
        zhetrf = thread_zhetrf;
    }
    zhetrf( pastix_data, sopalin_data );
}
