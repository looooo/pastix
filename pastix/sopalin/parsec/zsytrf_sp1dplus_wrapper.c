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
#include "common.h"
#include "sopalin_data.h"
#include "parsec/sparse-matrix.h"

typedef struct dague_zsytrf_sp1dplus_handle dague_zsytrf_sp1dplus_handle_t;

extern dague_zsytrf_sp1dplus_handle_t *dague_zsytrf_sp1dplus_new(dague_ddesc_t * dataA /* data dataA */, sopalin_data_t * sopalin_data);

dague_handle_t*
dsparse_zsytrf_sp_New( sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    dague_zsytrf_sp1dplus_handle_t *dague_zsytrf_sp = NULL;

    dague_zsytrf_sp = dague_zsytrf_sp1dplus_new( (dague_ddesc_t*)A, sopalin_data );

    /* dague_zsytrf_sp->p_work = (dague_memory_pool_t*)malloc(sizeof(dague_memory_pool_t)); */
    /* dague_private_memory_init( dague_zsytrf_sp->p_work, (A->pastix_data->solvmatr).coefmax * sizeof(dague_complex64_t) ); */

    /* dsparse_add2arena_tile(((dague_zsytrf_Url_handle_t*)dague_zsytrf)->arenas[DAGUE_zsytrf_Url_DEFAULT_ARENA],  */
    /*                        A->mb*A->nb*sizeof(dague_complex64_t), */
    /*                        DAGUE_ARENA_ALIGNMENT_SSE, */
    /*                        MPI_DOUBLE_COMPLEX, A->mb); */

    return (dague_handle_t*)dague_zsytrf_sp;
}

void
dsparse_zsytrf_sp_Destruct( dague_handle_t *o )
{
    (void)o;
    /* dague_zsytrf_sp1dplus_handle_t *dague_zsytrf_sp = NULL; */
    /* dague_zsytrf_sp = (dague_zsytrf_sp1dplus_handle_t *)o; */

    /*dsparse_datatype_undefine_type( &(osytrf->arenas[DAGUE_zsytrf_Url_DEFAULT_ARENA]->opaque_dtt) );*/

    /* dague_private_memory_fini( dague_zsytrf_sp->p_work ); */
    /* free( dague_zsytrf_sp->p_work ); */

    //DAGUE_INTERNAL_HANDLE_DESTRUCT(o);
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
