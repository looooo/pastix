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
#include "common.h"
#include "sopalin_data.h"
#include "parsec/sparse-matrix.h"
#include "parsec/zgetrf_sp1dplus.h"

dague_handle_t*
dsparse_zgetrf_sp_New( sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    dague_zgetrf_sp1dplus_handle_t *dague_zgetrf_sp = NULL;

    dague_zgetrf_sp = dague_zgetrf_sp1dplus_new( (dague_ddesc_t*)A, sopalin_data, NULL );

    dague_zgetrf_sp->p_work = (dague_memory_pool_t*)malloc(sizeof(dague_memory_pool_t));
    dague_private_memory_init( dague_zgetrf_sp->p_work, sopalin_data->solvmtx->gemmmax * sizeof(pastix_complex64_t) );

    /* dague_matrix_add2arena_rect( dague_zgetrf_sp->arenas[DAGUE_zgetrf_sp1dplus_DEFAULT_ARENA], */
    /*                              dague_datatype_double_complex_t, */
    /*                              sopalin_data->solvmtx->gemmmax, 1, 1 ); */

    return (dague_handle_t*)dague_zgetrf_sp;
}

void
dsparse_zgetrf_sp_Destruct( dague_handle_t *o )
{
    (void)o;
    dague_zgetrf_sp1dplus_handle_t *dague_zgetrf_sp = NULL;
    dague_zgetrf_sp = (dague_zgetrf_sp1dplus_handle_t *)o;

    /*dague_matrix_del2arena( dague_zgetrf_sp->arenas[DAGUE_zgetrf_sp1dplus_DEFAULT_ARENA] );*/

    dague_private_memory_fini( dague_zgetrf_sp->p_work );
    free( dague_zgetrf_sp->p_work );

    /* o->destructor(o); */
    /* o = NULL; */
}

int dsparse_zgetrf_sp( dague_context_t *dague,
                       sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    dague_handle_t *dague_zgetrf_sp = NULL;
    int info = 0;

    dague_zgetrf_sp = dsparse_zgetrf_sp_New( A, sopalin_data );

    if ( dague_zgetrf_sp != NULL )
    {
        dague_enqueue( dague, (dague_handle_t*)dague_zgetrf_sp);
        dague_context_wait( dague );
        dsparse_zgetrf_sp_Destruct( dague_zgetrf_sp );
    }
    return info;
}
