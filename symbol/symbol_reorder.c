/**
 *
 * @file symbol_reorder.c
 *
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Vincent Bridonneau
 * @date 2018-06-04
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include "symbol.h"
#include "pastix/order.h"

void
symbol_reorder_cblk( const symbol_matrix_t *symbptr,
                     const symbol_cblk_t   *cblk,
                     pastix_order_t     *order,
                     const pastix_int_t *levels,
                     pastix_int_t       *depthweight,
                     pastix_int_t        depthmax,
                     pastix_int_t        split_level,
                     pastix_int_t        stop_criteria );

/**
 * Version sequentielle
 */
void
sequential_reorder( pastix_data_t         *pastix_data,
                    pastix_int_t           maxdepth,
                    pastix_int_t          *levels,
                    pastix_int_t          *depthweight )
{
    symbol_matrix_t *symbptr = pastix_data->symbmtx;
    symbol_cblk_t   *cblk;
    pastix_int_t  cblknbr;
    pastix_int_t *iparm = pastix_data->iparm;
    // pastix_solv_mode_t mode = iparm[IPARM_SCHUR_SOLV_MODE];
    pastix_order_t *order = pastix_data->ordemesh;

    cblk = symbptr->cblktab;
    // cblknbr = (mode == PastixSolvModeSchur) ? symbptr->cblknbr : symbptr->cblkschur;
    cblknbr = symbptr->cblknbr;
    
    pastix_int_t itercblk;
    for (itercblk=0; itercblk<cblknbr; itercblk++, cblk++) {

        if (cblk->fcolnum >= symbptr->schurfcol )
            continue;

        memset( depthweight, 0, maxdepth * sizeof(pastix_int_t) );

        symbol_reorder_cblk( symbptr, cblk, order,
                             levels,
                             depthweight, maxdepth,
                             iparm[IPARM_REORDERING_SPLIT],
                             iparm[IPARM_REORDERING_STOP] );
    }
}

/**
 * Version parallel
 */
/* Arguments */
struct args_reorder_t
{

};

/* Fonction appelée par chaque thread */
void
thread_preorder( isched_thread_t *ctx, void *args )
{
    // struct args_reorder_t *arg = (struct args_reorder_t*)args;
    // pastix_data_t      *pastix_data  = arg->pastix_data;
    // sopalin_data_t     *sopalin_data = arg->sopalin_data;
    // SolverMatrix       *datacode = sopalin_data->solvmtx;
    // pastix_complex64_t *b = arg->b;
    // int nrhs  = arg->nrhs;
    // int ldb   = arg->ldb;
    // SolverCblk *cblk;
    // Task       *t;
    // pastix_int_t i, ii, cblknbr;
    // pastix_int_t tasknbr, *tasktab;
    // pastix_solv_mode_t mode = pastix_data->iparm[IPARM_SCHUR_SOLV_MODE];
    // int rank = ctx->rank;

    // tasknbr = datacode->ttsknbr[rank];
    // tasktab = datacode->ttsktab[rank];
    // cblknbr = (mode == PastixSolvModeSchur) ? datacode->cblknbr : datacode->cblkschur;

    // for (ii=0; ii<tasknbr; ii++) {
        // i = tasktab[ii];
        // t = datacode->tasktab + i;

        // if ( t->cblknum >= cblknbr ) {
        //     continue;
        // }
        // cblk = datacode->cblktab + t->cblknum;
        // solve_reorder( cblk, nrhs,
        //              b + cblk->lcolidx, ldb, NULL );
    // }
}

/* Fonction Multi-thread (prototype doit etre le meme que la version seq) */
void
thread_reorder( pastix_data_t         *pastix_data,
                pastix_int_t           maxdepth,
                pastix_int_t          *levels,
                pastix_int_t          *depthweight )
{
    // struct args_reorder_t args_reorder = { XXXX };
    // isched_parallel_call( pastix_data->isched, thread_preorder, &args_reorder );
}

static void (*reorder_table[4])(pastix_data_t *, pastix_int_t , pastix_int_t *, pastix_int_t *) = {
    sequential_reorder,
    thread_reorder,
    NULL,
    NULL
};

void
symbol_reorder( pastix_data_t         *pastix_data,
                pastix_int_t           maxdepth,
                pastix_int_t          *levels,
                pastix_int_t          *depthweight )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*reorder)(pastix_data_t *, pastix_int_t , pastix_int_t *, pastix_int_t *) = reorder_table[ sched ];

    if (reorder == NULL) {
        reorder = sequential_reorder;
    }
    reorder( pastix_data, maxdepth, levels, depthweight );
}
