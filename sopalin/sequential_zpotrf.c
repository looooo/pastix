/**
 *
 * @file sopalin_zpotrf.c
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
#include "sopalin/coeftab_z.h"
#include "pastix_zcores.h"

#if defined(PASTIX_WITH_PARSEC)
#include "parsec/pastix_zparsec.h"
#endif

#if defined(PASTIX_WITH_STARPU)
#include "starpu/pastix_zstarpu.h"
#endif

void
sequential_zpotrf( pastix_data_t  *pastix_data,
                   sopalin_data_t *sopalin_data )
{
    SolverMatrix       *datacode = pastix_data->solvmatr;
    SolverCblk         *cblk;
    double              threshold = sopalin_data->diagthreshold;
    pastix_complex64_t *work;
    pastix_int_t  i;
    (void)pastix_data;

    MALLOC_INTERN( work, datacode->gemmmax, pastix_complex64_t );

    cblk = datacode->cblktab;
    for (i=0; i<datacode->cblknbr; i++, cblk++){

        if ( cblk->cblktype & CBLK_IN_SCHUR )
            break;

        /* Compute */
        cpucblk_zpotrfsp1d( datacode, cblk, threshold, work );
    }

    memFree_null( work );
}

void
thread_pzpotrf( isched_thread_t *ctx, void *args )
{
    sopalin_data_t     *sopalin_data = (sopalin_data_t*)args;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    Task               *t;
    pastix_complex64_t *work;
    pastix_int_t  i, ii;
    pastix_int_t  tasknbr, *tasktab;
    int rank = ctx->rank;

    MALLOC_INTERN( work, datacode->gemmmax, pastix_complex64_t );

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        if ( cblk->cblktype & CBLK_IN_SCHUR )
            continue;

        /* Wait */
        do {} while( cblk->ctrbcnt );

        /* Compute */
        cpucblk_zpotrfsp1d( datacode, cblk, sopalin_data->diagthreshold, work );
    }

    memFree_null( work );
}

void
thread_zpotrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    isched_parallel_call( pastix_data->isched, thread_pzpotrf, sopalin_data );
}

static void (*zpotrf_table[4])(pastix_data_t *, sopalin_data_t *) = {
    sequential_zpotrf,
    thread_zpotrf,
#if defined(PASTIX_WITH_PARSEC)
    parsec_zpotrf,
#else
    NULL,
#endif
#if defined(PASTIX_WITH_STARPU)
    starpu_zpotrf
#else
    NULL
#endif
};

void
sopalin_zpotrf( pastix_data_t  *pastix_data,
                sopalin_data_t *sopalin_data )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zpotrf)(pastix_data_t *, sopalin_data_t *) = zpotrf_table[ sched ];

    if (zpotrf == NULL) {
        zpotrf = thread_zpotrf;
    }
    zpotrf( pastix_data, sopalin_data );

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( pastix_data, sopalin_data->solvmtx, "potrf.txt" );
#endif
}
