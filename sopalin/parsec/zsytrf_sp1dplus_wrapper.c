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
#include "parsec/zsytrf_sp1dplus.h"

dague_handle_t*
dsparse_zsytrf_sp_New( sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    dague_zsytrf_sp1dplus_handle_t *dague_zsytrf_sp = NULL;

    dague_zsytrf_sp = dague_zsytrf_sp1dplus_new( (dague_ddesc_t*)A, sopalin_data, NULL, NULL );

    dague_zsytrf_sp->p_work1 = (dague_memory_pool_t*)malloc(sizeof(dague_memory_pool_t));
    dague_private_memory_init( dague_zsytrf_sp->p_work1,
                               pastix_imax(sopalin_data->solvmtx->gemmmax,
                                           sopalin_data->solvmtx->diagmax) * sizeof(pastix_complex64_t) );

    dague_zsytrf_sp->p_work2 = (dague_memory_pool_t*)malloc(sizeof(dague_memory_pool_t));
    dague_private_memory_init( dague_zsytrf_sp->p_work2, sopalin_data->solvmtx->gemmmax * sizeof(pastix_complex64_t) );

    /* dague_matrix_add2arena_rect( dague_zsytrf_sp->arenas[DAGUE_zsytrf_sp1dplus_DEFAULT_ARENA], */
    /*                              dague_datatype_double_complex_t, */
    /*                              sopalin_data->solvmtx->gemmmax, 1, 1 ); */

    return (dague_handle_t*)dague_zsytrf_sp;
}

void
dsparse_zsytrf_sp_Destruct( dague_handle_t *o )
{
    (void)o;
    dague_zsytrf_sp1dplus_handle_t *dague_zsytrf_sp = NULL;
    dague_zsytrf_sp = (dague_zsytrf_sp1dplus_handle_t *)o;

    /*dague_matrix_del2arena( dague_zsytrf_sp->arenas[DAGUE_zsytrf_sp1dplus_DEFAULT_ARENA] );*/

    dague_private_memory_fini( dague_zsytrf_sp->p_work1 );
    free( dague_zsytrf_sp->p_work1 );

    dague_private_memory_fini( dague_zsytrf_sp->p_work2 );
    free( dague_zsytrf_sp->p_work2 );

    /* o->destructor(o); */
    /* o = NULL; */
}

int dsparse_zsytrf_sp( dague_context_t *dague,
                       sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    dague_handle_t *dague_zsytrf_sp = NULL;
    int info = 0;

    dague_zsytrf_sp = dsparse_zsytrf_sp_New( A, sopalin_data );

    if ( dague_zsytrf_sp != NULL )
    {
        dague_enqueue( dague, (dague_handle_t*)dague_zsytrf_sp);
        dague_context_wait( dague );
        dsparse_zsytrf_sp_Destruct( dague_zsytrf_sp );
    }
    return info;
}
