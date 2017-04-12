/**
 *
 * @file parsec_zpotrf.c
 *
 *  PaStiX factorization routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @date 2016-01-01
 *
 * @precisions normal z -> s d c
 *
 **/
#include <parsec.h>
#include <parsec/data_distribution.h>
#include <parsec/private_mempool.h>
#include <parsec/arena.h>
#include <data_dist/matrix/matrix.h>
#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "parsec/zpotrf_sp1dplus.h"
#include "parsec/zpotrf_sp2d.h"

parsec_handle_t*
parsec_zpotrf_sp1dplus_New( sparse_matrix_desc_t *A,
                             sopalin_data_t *sopalin_data )
{
    parsec_zpotrf_sp1dplus_handle_t *parsec_zpotrf_sp1dplus = NULL;

    parsec_zpotrf_sp1dplus = parsec_zpotrf_sp1dplus_new( A, sopalin_data, NULL );

    parsec_zpotrf_sp1dplus->_g_p_work = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( parsec_zpotrf_sp1dplus->_g_p_work, sopalin_data->solvmtx->gemmmax * sizeof(pastix_complex64_t) );

    parsec_matrix_add2arena_rect( parsec_zpotrf_sp1dplus->arenas[PARSEC_zpotrf_sp1dplus_DEFAULT_ARENA],
                                  parsec_datatype_double_complex_t,
                                  /*sopalin_data->solvmtx->gemmmax*/ 1, 1, 1 );

    return (parsec_handle_t*)parsec_zpotrf_sp1dplus;
}

void
parsec_zpotrf_sp1dplus_Destruct( parsec_handle_t *handle )
{
    parsec_zpotrf_sp1dplus_handle_t *parsec_zpotrf_sp1dplus = NULL;
    parsec_zpotrf_sp1dplus = (parsec_zpotrf_sp1dplus_handle_t *)handle;

    parsec_matrix_del2arena( parsec_zpotrf_sp1dplus->arenas[PARSEC_zpotrf_sp1dplus_DEFAULT_ARENA] );

    parsec_private_memory_fini( parsec_zpotrf_sp1dplus->_g_p_work );
    free( parsec_zpotrf_sp1dplus->_g_p_work );

    parsec_handle_free( handle );
}

int
parsec_zpotrf_sp1dplus( parsec_context_t *parsec,
                         sparse_matrix_desc_t *A,
                         sopalin_data_t *sopalin_data )
{
    parsec_handle_t *parsec_zpotrf_sp1dplus = NULL;
    int info = 0;

    parsec_zpotrf_sp1dplus = parsec_zpotrf_sp1dplus_New( A, sopalin_data );

    if ( parsec_zpotrf_sp1dplus != NULL )
    {
        parsec_enqueue( parsec, (parsec_handle_t*)parsec_zpotrf_sp1dplus);
        parsec_context_start( parsec );
        parsec_context_wait( parsec );
        parsec_zpotrf_sp1dplus_Destruct( parsec_zpotrf_sp1dplus );
    }
    return info;
}

parsec_handle_t*
parsec_zpotrf_sp2d_New( sparse_matrix_desc_t *A,
                         sopalin_data_t *sopalin_data )
{
    parsec_zpotrf_sp2d_handle_t *parsec_zpotrf_sp2d = NULL;

    parsec_zpotrf_sp2d = parsec_zpotrf_sp2d_new( A, sopalin_data, NULL );

    parsec_zpotrf_sp2d->_g_p_work = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( parsec_zpotrf_sp2d->_g_p_work, sopalin_data->solvmtx->gemmmax * sizeof(pastix_complex64_t) );

    parsec_matrix_add2arena_rect( parsec_zpotrf_sp2d->arenas[PARSEC_zpotrf_sp2d_DEFAULT_ARENA],
                                  parsec_datatype_double_complex_t,
                                  /*sopalin_data->solvmtx->gemmmax*/ 1, 1, 1 );

    return (parsec_handle_t*)parsec_zpotrf_sp2d;
}

void
parsec_zpotrf_sp2d_Destruct( parsec_handle_t *handle )
{
    parsec_zpotrf_sp2d_handle_t *parsec_zpotrf_sp2d = NULL;
    parsec_zpotrf_sp2d = (parsec_zpotrf_sp2d_handle_t *)handle;

    parsec_matrix_del2arena( parsec_zpotrf_sp2d->arenas[PARSEC_zpotrf_sp2d_DEFAULT_ARENA] );

    parsec_private_memory_fini( parsec_zpotrf_sp2d->_g_p_work );
    free( parsec_zpotrf_sp2d->_g_p_work );

    parsec_handle_free( handle );
}

int
parsec_zpotrf_sp2d( parsec_context_t *parsec,
                     sparse_matrix_desc_t *A,
                     sopalin_data_t *sopalin_data )
{
    parsec_handle_t *parsec_zpotrf_sp2d = NULL;
    int info = 0;

    parsec_zpotrf_sp2d = parsec_zpotrf_sp2d_New( A, sopalin_data );

    if ( parsec_zpotrf_sp2d != NULL )
    {
        parsec_enqueue( parsec, (parsec_handle_t*)parsec_zpotrf_sp2d );
        parsec_context_start( parsec );
        parsec_context_wait( parsec );
        parsec_zpotrf_sp2d_Destruct( parsec_zpotrf_sp2d );
    }
    return info;
}

void
parsec_zpotrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    sparse_matrix_desc_t *sdesc = sopalin_data->solvmtx->parsec_desc;
    parsec_context_t *ctx;

    /*
     * Start PaRSEC if not already started
     */
    if (pastix_data->parsec == NULL) {
        int argc = 0;
        pastix_parsec_init( pastix_data, &argc, NULL );
    }
    ctx = pastix_data->parsec;

    if ( sdesc == NULL ) {
        sdesc = (sparse_matrix_desc_t*)malloc(sizeof(sparse_matrix_desc_t));

        /* Create the matrix descriptor */
        sparse_matrix_init( sdesc, sopalin_data->solvmtx,
                            sizeof( pastix_complex64_t ), PastixGeneral,
                            1, 0 );
        sopalin_data->solvmtx->parsec_desc = sdesc;
    }

    /*
     * Select 1D or 2D jdf based on distribution_level
     */
    if ( pastix_data->iparm[IPARM_DISTRIBUTION_LEVEL] >= 0 )
    {
        parsec_zpotrf_sp2d( ctx, sdesc,
                            sopalin_data );
    }
    else {
        parsec_zpotrf_sp1dplus( ctx, sdesc,
                                sopalin_data );
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( sopalin_data->solvmtx, "potrf.txt" );
#endif

    return;
}

