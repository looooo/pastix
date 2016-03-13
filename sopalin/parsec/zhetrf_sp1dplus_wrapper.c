/*
 * Copyright (c) 2010      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> c
 *
 */
#include <dague.h>
#include <dague/data_distribution.h>
#include <dague/private_mempool.h>
#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "parsec/zhetrf_sp1dplus.h"

dague_handle_t*
dsparse_zhetrf_sp_New( sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    dague_zhetrf_sp1dplus_handle_t *dague_zhetrf_sp = NULL;

    dague_zhetrf_sp = dague_zhetrf_sp1dplus_new( (dague_ddesc_t*)A, sopalin_data, NULL, NULL );

    dague_zhetrf_sp->p_work1 = (dague_memory_pool_t*)malloc(sizeof(dague_memory_pool_t));
    dague_private_memory_init( dague_zhetrf_sp->p_work1,
                               pastix_imax(sopalin_data->solvmtx->gemmmax,
                                           sopalin_data->solvmtx->diagmax) * sizeof(pastix_complex64_t) );

    dague_zhetrf_sp->p_work2 = (dague_memory_pool_t*)malloc(sizeof(dague_memory_pool_t));
    dague_private_memory_init( dague_zhetrf_sp->p_work2, sopalin_data->solvmtx->gemmmax * sizeof(pastix_complex64_t) );

    /* dague_matrix_add2arena_rect( dague_zhetrf_sp->arenas[DAGUE_zhetrf_sp1dplus_DEFAULT_ARENA], */
    /*                              dague_datatype_double_complex_t, */
    /*                              sopalin_data->solvmtx->gemmmax, 1, 1 ); */

    return (dague_handle_t*)dague_zhetrf_sp;
}

void
dsparse_zhetrf_sp_Destruct( dague_handle_t *o )
{
    (void)o;
    dague_zhetrf_sp1dplus_handle_t *dague_zhetrf_sp = NULL;
    dague_zhetrf_sp = (dague_zhetrf_sp1dplus_handle_t *)o;

    /*dague_matrix_del2arena( dague_zhetrf_sp->arenas[DAGUE_zhetrf_sp1dplus_DEFAULT_ARENA] );*/

    dague_private_memory_fini( dague_zhetrf_sp->p_work1 );
    free( dague_zhetrf_sp->p_work1 );

    dague_private_memory_fini( dague_zhetrf_sp->p_work2 );
    free( dague_zhetrf_sp->p_work2 );

    /* o->destructor(o); */
    /* o = NULL; */
}

int dsparse_zhetrf_sp( dague_context_t *dague,
                       sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    dague_handle_t *dague_zhetrf_sp = NULL;
    int info = 0;

    dague_zhetrf_sp = dsparse_zhetrf_sp_New( A, sopalin_data );

    if ( dague_zhetrf_sp != NULL )
    {
        dague_enqueue( dague, (dague_handle_t*)dague_zhetrf_sp);
        dague_context_wait( dague );
        dsparse_zhetrf_sp_Destruct( dague_zhetrf_sp );
    }
    return info;
}
