/*
 * Copyright (c) 2010      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#include <dague.h>
#include <dague/data_distribution.h>
#include <dague/private_mempool.h>
#include <dague/arena.h>
#include <data_dist/matrix/matrix.h>
#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "parsec/zpotrf_sp2d.h"

dague_handle_t*
dsparse_zpotrf_sp_New( sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    dague_zpotrf_sp2d_handle_t *dague_zpotrf_sp = NULL;

    dague_zpotrf_sp = dague_zpotrf_sp2d_new( (dague_ddesc_t*)A, sopalin_data, NULL );

    dague_zpotrf_sp->p_work = (dague_memory_pool_t*)malloc(sizeof(dague_memory_pool_t));
    dague_private_memory_init( dague_zpotrf_sp->p_work, sopalin_data->solvmtx->gemmmax * sizeof(pastix_complex64_t) );

    dague_matrix_add2arena_rect( dague_zpotrf_sp->arenas[DAGUE_zpotrf_sp2d_DEFAULT_ARENA],
                                 dague_datatype_double_complex_t,
                                 /*sopalin_data->solvmtx->gemmmax*/ 1, 1, 1 );

    return (dague_handle_t*)dague_zpotrf_sp;
}

void
dsparse_zpotrf_sp_Destruct( dague_handle_t *handle )
{
    dague_zpotrf_sp2d_handle_t *dague_zpotrf_sp = NULL;
    dague_zpotrf_sp = (dague_zpotrf_sp2d_handle_t *)handle;

    dague_matrix_del2arena( dague_zpotrf_sp->arenas[DAGUE_zpotrf_sp2d_DEFAULT_ARENA] );

    dague_private_memory_fini( dague_zpotrf_sp->p_work );
    free( dague_zpotrf_sp->p_work );

    dague_handle_free( handle );
}

int dsparse_zpotrf_sp( dague_context_t *dague,
                       sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    dague_handle_t *dague_zpotrf_sp = NULL;
    int info = 0;

    dague_zpotrf_sp = dsparse_zpotrf_sp_New( A, sopalin_data );

    if ( dague_zpotrf_sp != NULL )
    {
        dague_enqueue( dague, (dague_handle_t*)dague_zpotrf_sp);
        dague_context_wait( dague );
        dsparse_zpotrf_sp_Destruct( dague_zpotrf_sp );
    }
    return info;
}
