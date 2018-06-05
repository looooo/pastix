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
                    pastix_int_t          *levels )
{
    symbol_matrix_t *symbptr = pastix_data->symbmtx;
    symbol_cblk_t   *cblk;
    pastix_int_t     cblknbr;
    pastix_int_t    *iparm = pastix_data->iparm;
    pastix_order_t  *order = pastix_data->ordemesh;
    pastix_int_t     itercblk;
    pastix_int_t    *depthweight;

    cblk = symbptr->cblktab;
    cblknbr = symbptr->cblknbr;

    /**
     * Solves the Traveler Salesman Problem on each cblk to minimize the number
     * of off-diagonal blocks per row
     */
    MALLOC_INTERN( depthweight, maxdepth, pastix_int_t );
    
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

    memFree_null( depthweight );
}

/**
 * Version parallel
 */
/* Arguments */
struct args_reorder_t
{
    pastix_data_t         *pastix_data;
    pastix_int_t           maxdepth;
    const pastix_int_t    *levels;
};

/* Fonction appelée par chaque thread */
void
thread_preorder( isched_thread_t *ctx, void *args )
{
    struct args_reorder_t *arg = (struct args_reorder_t*)args;
    pastix_data_t      *pastix_data = arg->pastix_data;
    symbol_matrix_t    *symbptr = pastix_data->symbmtx;
    symbol_cblk_t      *cblk;
    pastix_int_t       *iparm = pastix_data->iparm;
    pastix_order_t     *order = pastix_data->ordemesh;
    pastix_int_t        maxdepth = arg->maxdepth;
    pastix_int_t        ii, cblknbr, tasknbr;
    pastix_int_t       *depthweight;
    const pastix_int_t *levels = arg->levels;
    /* Task       *t;
     * pastix_int_t *tasktab;
     */
    pastix_int_t        rank = (pastix_int_t)ctx->rank;
    pastix_int_t        size = (pastix_int_t)ctx->global_ctx->world_size; // Number of thread

    /* tasktab = symbptr->ttsktab[rank];
     */
    cblknbr = symbptr->cblknbr;
    tasknbr = cblknbr / size;

    if (rank < (cblknbr % size)) {
        tasknbr ++;
    }

    /**
     * Solves the Traveler Salesman Problem on each cblk to minimize the number
     * of off-diagonal blocks per row
     */
    MALLOC_INTERN( depthweight, maxdepth, pastix_int_t );

    cblk = symbptr->cblktab + rank;
    for (ii=0; ii<tasknbr; ii++, cblk += size) {

        if (cblk->fcolnum >= symbptr->schurfcol )
            continue;

        memset( depthweight, 0, maxdepth * sizeof(pastix_int_t) );
        
        symbol_reorder_cblk( symbptr, cblk, order,
                             levels,
                             depthweight, maxdepth,
                             iparm[IPARM_REORDERING_SPLIT],
                             iparm[IPARM_REORDERING_STOP] );
    }

    memFree_null( depthweight );
}

/* Fonction Multi-thread (prototype doit etre le meme que la version seq) */
void
thread_reorder( pastix_data_t         *pastix_data,
                pastix_int_t           maxdepth,
                pastix_int_t          *levels )
{
    struct args_reorder_t args_reorder = { pastix_data, maxdepth, levels };
    isched_parallel_call( pastix_data->isched, thread_preorder, &args_reorder );
}

static void (*reorder_table[4])(pastix_data_t *, pastix_int_t , pastix_int_t *) = {
    sequential_reorder,
    thread_reorder,
    NULL,
    NULL
};

void
symbol_reorder( pastix_data_t         *pastix_data,
                pastix_int_t           maxdepth,
                pastix_int_t          *levels )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*reorder)(pastix_data_t *, pastix_int_t , pastix_int_t * ) = reorder_table[ sched ];

    if (reorder == NULL) {
        reorder = sequential_reorder;
    }
    reorder( pastix_data, maxdepth, levels );
}
