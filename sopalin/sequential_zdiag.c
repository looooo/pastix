/**
 *
 * @file sequential_zdiag.c
 *
 *  PaStiX routines
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
#include "cblas.h"
#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"

void
sequential_zdiag( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data,
                  int nrhs, pastix_complex64_t *b, int ldb )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk   *cblk;
    pastix_int_t  i;
    (void)pastix_data;

    cblk = datacode->cblktab;
    for (i=0; i<datacode->cblknbr; i++, cblk++){

        if ( (cblk->cblktype & CBLK_IN_SCHUR) && (sopalin_data->solvmode != PastixSolvModeSchur) )
            break;

        solve_zdiag( cblk, nrhs,
                     b + cblk->lcolidx, ldb, NULL );
    }
}

struct args_zdiag_t
{
    sopalin_data_t *sopalin_data;
    int nrhs;
    pastix_complex64_t *b;
    int ldb;
};

void
thread_pzdiag( isched_thread_t *ctx, void *args )
{
    struct args_zdiag_t *arg = (struct args_zdiag_t*)args;
    sopalin_data_t     *sopalin_data = arg->sopalin_data;
    SolverMatrix       *datacode = sopalin_data->solvmtx;
    pastix_complex64_t *b = arg->b;
    int nrhs  = arg->nrhs;
    int ldb   = arg->ldb;
    SolverCblk *cblk;
    Task       *t;
    pastix_int_t i,ii;
    pastix_int_t tasknbr, *tasktab;
    int rank = ctx->rank;

    tasknbr = datacode->ttsknbr[rank];
    tasktab = datacode->ttsktab[rank];

    for (ii=0; ii<tasknbr; ii++) {
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        if ( (cblk->cblktype & CBLK_IN_SCHUR) && (sopalin_data->solvmode != PastixSolvModeSchur) )
            continue;

        solve_zdiag( cblk, nrhs,
                     b + cblk->lcolidx, ldb, NULL );
    }
}

void
thread_zdiag( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data,
              int nrhs, pastix_complex64_t *b, int ldb )
{
    struct args_zdiag_t args_zdiag = {sopalin_data, nrhs, b, ldb};
    isched_parallel_call( pastix_data->isched, thread_pzdiag, &args_zdiag );
}

static void (*zdiag_table[4])(pastix_data_t *, sopalin_data_t *,
                              int, pastix_complex64_t *, int) = {
    sequential_zdiag,
    thread_zdiag,
#if defined(PASTIX_WITH_PARSEC)
    NULL, /* parsec_zdiag not yet implemented */
#else
    NULL,
#endif
    NULL
};

void
sopalin_zdiag( pastix_data_t *pastix_data, sopalin_data_t *sopalin_data,
               int nrhs, pastix_complex64_t *b, int ldb )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*zdiag)(pastix_data_t *, sopalin_data_t *, int, pastix_complex64_t *, int) = zdiag_table[ sched ];

    if (zdiag == NULL) {
        zdiag = thread_zdiag;
    }
    zdiag( pastix_data, sopalin_data, nrhs, b, ldb );
}
