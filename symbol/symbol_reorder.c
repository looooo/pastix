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
#include "blend/queue.h"
#include "blend/extendVector.h"

/**
 * @brief Sequential version for reordering
 */
static inline void
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

    /*
     * Solves the Traveler Salesman Problem on each cblk to minimize the number
     * of off-diagonal blocks per row
     */
    MALLOC_INTERN( depthweight, maxdepth, pastix_int_t );
    
    for (itercblk=0; itercblk<cblknbr; itercblk++, cblk++) {

        if ( cblk->fcolnum >= symbptr->schurfcol )
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

struct args_reorder_t
{
    pastix_data_t         *pastix_data;
    pastix_int_t           maxdepth;
    const pastix_int_t    *levels;
    ExtendVectorINT       *tasktab;
};

/**
 * @brief Parallel basic version for reordering
 */
/* Arguments */
static inline void
thread_preorder_basic_stategy( isched_thread_t *ctx, void *args )
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
    pastix_int_t        rank = (pastix_int_t)ctx->rank;
    pastix_int_t        size = (pastix_int_t)ctx->global_ctx->world_size; // Number of thread

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
 * @brief Parallel improved version for reordering
 */
static inline void
thread_preorder_zigzag_stategy( isched_thread_t *ctx, void *args )
{
    struct args_reorder_t *arg = (struct args_reorder_t*)args;
    pastix_data_t      *pastix_data = arg->pastix_data;
    symbol_matrix_t    *symbptr = pastix_data->symbmtx;
    symbol_cblk_t      *cblk;
    pastix_int_t       *iparm = pastix_data->iparm;
    pastix_order_t     *order = pastix_data->ordemesh;
    pastix_int_t        maxdepth = arg->maxdepth;
    pastix_int_t        ii;
    ExtendVectorINT    *tasktab;
    pastix_int_t        tasknbr;
    pastix_int_t       *depthweight;
    const pastix_int_t *levels = arg->levels;   
    pastix_int_t        rank = (pastix_int_t)ctx->rank;

    /**
     * Solves the Traveler Salesman Problem on each cblk to minimize the number
     * of off-diagonal blocks per row
     */
    MALLOC_INTERN( depthweight, maxdepth, pastix_int_t );

    tasktab = arg->tasktab + rank;
    tasknbr = extendint_Size( tasktab );

    for ( ii=0; ii<tasknbr; ii++ ) {
        
        cblk = symbptr->cblktab + extendint_Read(tasktab, ii);

        memset( depthweight, 0, maxdepth * sizeof(pastix_int_t) );
        
        symbol_reorder_cblk( symbptr, cblk, order,
                             levels,
                             depthweight, maxdepth,
                             iparm[IPARM_REORDERING_SPLIT],
                             iparm[IPARM_REORDERING_STOP] );
    }



    memFree_null( depthweight );
}

/* Fonction appelée par chaque thread */
static inline void
thread_preorder( isched_thread_t *ctx, void *args )
{
#ifdef BASIC_REORDERING_STRATEGY
    thread_preorder_basic_stategy ( ctx, args );
#else
    thread_preorder_zigzag_stategy( ctx, args );
#endif
}

static inline double
cost( symbol_cblk_t *cblk )
{
    double n = (double)(cblk->lcolnum - cblk->fcolnum + 1);
    return n*n * ((double)(cblk[1].bloknum - cblk[0].bloknum) / 2.0 + 1.0);
}

static inline void 
order_tasks( isched_t              *ctx,
             struct args_reorder_t *args )
{
    pastix_data_t   *pastix_data = args->pastix_data;
    symbol_matrix_t *symbmtx = pastix_data->symbmtx;
    pastix_queue_t   cblks;
    pastix_queue_t   procs;
    pastix_int_t     cblknbr = symbmtx->cblknbr;
    pastix_int_t     size = ctx->world_size;
    pastix_int_t     itercblk, iterproc;
    pastix_int_t     cblk_id, proc_id;
    double           cblk_cost, proc_cost;
    symbol_cblk_t   *cblk;

    pqueueInit( &cblks, cblknbr );
    pqueueInit( &procs, size );

    /*
     * Sort the cblks decreasing
     */
    cblk = symbmtx->cblktab;
    for ( itercblk=0; itercblk < cblknbr; itercblk++, cblk++ ) {

        if ( cblk->fcolnum >= symbmtx->schurfcol )
            continue;

        pqueuePush1( &cblks, itercblk, -cost(cblk) );
    }

    for ( iterproc=0; iterproc < size; ++iterproc) {
        pqueuePush1( &procs, iterproc, 0. );
    }

    while ( pqueueSize( &cblks ) > 0 ) {
        cblk_id = pqueuePop1( &cblks, &cblk_cost );
        proc_id = pqueuePop1( &procs, &proc_cost );
        proc_cost += -cblk_cost; // Negative because of reverse sort
        pqueuePush1 ( &procs, proc_id, proc_cost );

        extendint_Add( args->tasktab + proc_id, cblk_id );
    }

    pqueueExit( &cblks );
    pqueueExit( &procs );
}

/* Fonction Multi-thread (prototype doit etre le meme que la version seq) */
static inline void
thread_reorder( pastix_data_t *pastix_data,
                pastix_int_t   maxdepth,
                pastix_int_t  *levels )
{
    struct args_reorder_t args_reorder = { pastix_data, maxdepth, levels, NULL };
    pastix_int_t          size = pastix_data->isched->world_size;
    pastix_int_t          iterproc;
    pastix_int_t          cblknbr = pastix_data->symbmtx->cblknbr;
    pastix_int_t          cblkavg = pastix_imax(1, (pastix_int_t)(cblknbr / size));

    MALLOC_INTERN( args_reorder.tasktab, size, ExtendVectorINT );

    for ( iterproc=0; iterproc < size; ++iterproc ) {
        extendint_Init( args_reorder.tasktab + iterproc, cblkavg );
    }
    order_tasks( pastix_data->isched, &args_reorder );

    isched_parallel_call( pastix_data->isched, thread_preorder, &args_reorder );

    for ( iterproc=0; iterproc < size; ++iterproc ) {
        extendint_Exit( args_reorder.tasktab + iterproc );
    }

    memFree_null( args_reorder.tasktab );
}

static void (*reorder_table[4])(pastix_data_t *, pastix_int_t , pastix_int_t *) = {
    sequential_reorder,
    thread_reorder,
    NULL,
    NULL
};

void
symbol_reorder( pastix_data_t *pastix_data,
                pastix_int_t   maxdepth,
                pastix_int_t  *levels )
{
    int sched = pastix_data->iparm[IPARM_SCHEDULER];
    void (*reorder)(pastix_data_t *, pastix_int_t , pastix_int_t * ) = reorder_table[ sched ];

    if (reorder == NULL) {
        reorder = thread_reorder;
    }
    reorder( pastix_data, maxdepth, levels );
}
