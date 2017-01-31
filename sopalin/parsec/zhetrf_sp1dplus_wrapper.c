/*
 * Copyright (c) 2010      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> c
 *
 */
#include <parsec.h>
#include <parsec/data_distribution.h>
#include <parsec/private_mempool.h>
#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "parsec/zhetrf_sp1dplus.h"

parsec_handle_t*
dsparse_zhetrf_sp_New( sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    parsec_zhetrf_sp1dplus_handle_t *parsec_zhetrf_sp = NULL;

    parsec_zhetrf_sp = parsec_zhetrf_sp1dplus_new( A, sopalin_data, NULL, NULL );

    parsec_zhetrf_sp->_g_p_work1 = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( parsec_zhetrf_sp->_g_p_work1,
                               pastix_imax(sopalin_data->solvmtx->gemmmax,
                                           sopalin_data->solvmtx->diagmax) * sizeof(pastix_complex64_t) );

    parsec_zhetrf_sp->_g_p_work2 = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( parsec_zhetrf_sp->_g_p_work2, sopalin_data->solvmtx->gemmmax * sizeof(pastix_complex64_t) );

    /* parsec_matrix_add2arena_rect( parsec_zhetrf_sp->arenas[PARSEC_zhetrf_sp1dplus_DEFAULT_ARENA], */
    /*                              parsec_datatype_double_complex_t, */
    /*                              sopalin_data->solvmtx->gemmmax, 1, 1 ); */

    return (parsec_handle_t*)parsec_zhetrf_sp;
}

void
dsparse_zhetrf_sp_Destruct( parsec_handle_t *o )
{
    (void)o;
    parsec_zhetrf_sp1dplus_handle_t *parsec_zhetrf_sp = NULL;
    parsec_zhetrf_sp = (parsec_zhetrf_sp1dplus_handle_t *)o;

    /*parsec_matrix_del2arena( parsec_zhetrf_sp->arenas[PARSEC_zhetrf_sp1dplus_DEFAULT_ARENA] );*/

    parsec_private_memory_fini( parsec_zhetrf_sp->_g_p_work1 );
    free( parsec_zhetrf_sp->_g_p_work1 );

    parsec_private_memory_fini( parsec_zhetrf_sp->_g_p_work2 );
    free( parsec_zhetrf_sp->_g_p_work2 );

    /* o->destructor(o); */
    /* o = NULL; */
}

int dsparse_zhetrf_sp( parsec_context_t *parsec,
                       sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    parsec_handle_t *parsec_zhetrf_sp = NULL;
    int info = 0;

    parsec_zhetrf_sp = dsparse_zhetrf_sp_New( A, sopalin_data );

    if ( parsec_zhetrf_sp != NULL )
    {
        parsec_enqueue( parsec, (parsec_handle_t*)parsec_zhetrf_sp);
        parsec_context_wait( parsec );
        dsparse_zhetrf_sp_Destruct( parsec_zhetrf_sp );
    }
    return info;
}
