/*
 * Copyright (c) 2010      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#include "parsec_internal.h"
#include "dsparse.h"
#include "data_dist/sparse-matrix/pastix_internal/pastix_internal.h"

#include "memory_pool.h"
#include "zpotrs_sp1dplus.h"

parsec_object_t*
dsparse_zpotrs_sp_New(sparse_matrix_desc_t *A, sparse_vector_desc_t *B)
{
    parsec_zpotrs_sp1dplus_object_t *parsec_zpotrs_sp = NULL;

    parsec_zpotrs_sp = parsec_zpotrs_sp1dplus_new(A, B, NULL );

    parsec_zpotrs_sp->p_work = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( parsec_zpotrs_sp->p_work, (A->pastix_data->solvmatr).coefmax * sizeof(parsec_complex64_t) );

    /* dsparse_add2arena_tile(((parsec_zpotrs_Url_object_t*)parsec_zpotrs)->arenas[PARSEC_zpotrs_Url_DEFAULT_ARENA],  */
    /*                        A->mb*A->nb*sizeof(parsec_complex64_t), */
    /*                        PARSEC_ARENA_ALIGNMENT_SSE, */
    /*                        MPI_DOUBLE_COMPLEX, A->mb); */

    return (parsec_object_t*)parsec_zpotrs_sp;
}

void
dsparse_zpotrs_sp_Destruct( parsec_object_t *o )
{
    parsec_zpotrs_sp1dplus_object_t *parsec_zpotrs_sp = NULL;
    parsec_zpotrs_sp = (parsec_zpotrs_sp1dplus_object_t *)o;

    /*dsparse_datatype_undefine_type( &(opotrs->arenas[PARSEC_zpotrs_Url_DEFAULT_ARENA]->opaque_dtt) );*/

    parsec_private_memory_fini( parsec_zpotrs_sp->p_work );
    free( parsec_zpotrs_sp->p_work );

    PARSEC_INTERNAL_OBJECT_DESTRUCT(o);
}

int dsparse_zpotrs_sp( parsec_context_t *parsec, sparse_matrix_desc_t *A, sparse_vector_desc_t *B)
{
    parsec_object_t *parsec_zpotrs_sp = NULL;
    int info = 0;

    parsec_zpotrs_sp = dsparse_zpotrs_sp_New( A, B );

    if ( parsec_zpotrs_sp != NULL )
    {
        parsec_enqueue( parsec, (parsec_object_t*)parsec_zpotrs_sp);
        parsec_progress( parsec );
        dsparse_zpotrs_sp_Destruct( parsec_zpotrs_sp );
    }
    return info;
}
