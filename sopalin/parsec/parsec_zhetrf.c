/**
 *
 * @file parsec_zhetrf.c
 *
 * PaStiX zhetrf PaRSEC wrapper.
 *
 * @copyright 2016-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Mathieu Faverge
 * @author Esragul Korkmaz
 * @date 2024-07-05
 * @precisions normal z -> c
 *
 * @addtogroup parsec_hetrf
 * @{
 *
 **/
#include <parsec.h>
#include <parsec/data_distribution.h>
#include <parsec/private_mempool.h>
#include <parsec/arena.h>
#include <parsec/data_dist/matrix/matrix.h>
#include "common.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"
#include "zhetrf_sp1dplus.h"
#include "zhetrf_sp2d.h"
#include "pastix_parsec.h"

/**
 *******************************************************************************
 *
 * @brief Generate the PaRSEC taskpool object for the LDL^h factorization with
 * 1D kernels.
 *
 * The function only return the object that describes the LDL^h factorization
 * of a sparse hermitian matrix A.
 * The factorization has the form
 *
 *    \f[ A = L \times D \times L^h \f]
 *
 * where L is a sparse lower triangular matrix unitary on the diagonal, and D a
 * diagonal matrix.
 *
 * In this object, all operations are panel based operation which might result in
 * lower parallelism for large problems.
 *
 * @warning The computations are not done by this call.
 *
 *******************************************************************************
 *
 * @param[inout] A
 *          Descriptor of the sparse matrix A.
 *          On exit, A is overwritten with the factorized matrix.
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 *******************************************************************************
 *
 * @retval NULL if incorrect parameters are given.
 * @retval The parsec taskpool describing the operation that can be
 *         enqueued in the runtime with parsec_enqueue(). It, then, needs to be
 *         destroy with parsec_zhetrf_sp1dplus_Destruct().
 *
 *******************************************************************************
 *
 * @sa parsec_zhetrf_sp1dplus
 * @sa parsec_zhetrf_sp1dplus_Destruct
 * @sa parsec_zhetrf_sp1dplus
 * @sa parsec_chetrf_sp1dplus_New
 * @sa parsec_dhetrf_sp1dplus_New
 * @sa parsec_shetrf_sp1dplus_New
 *
 ******************************************************************************/
parsec_taskpool_t*
parsec_zhetrf_sp1dplus_New( parsec_sparse_matrix_desc_t *A,
                            sopalin_data_t *sopalin_data )
{
    parsec_zhetrf_sp1dplus_taskpool_t *parsec_zhetrf_sp1dplus = NULL;
    unsigned int gemmmax = sopalin_data->solvmtx->gemmmax;
    unsigned int offdmax = sopalin_data->solvmtx->offdmax;
    unsigned int blokmax = sopalin_data->solvmtx->blokmax;
    pastix_int_t lwork;

    lwork = pastix_imax( gemmmax, blokmax );
    if ( (sopalin_data->solvmtx->lowrank.compress_when != PastixCompressNever) &&
         (sopalin_data->solvmtx->lowrank.ilu_lvl < INT_MAX) )
    {
        lwork = pastix_imax( lwork, 2 * blokmax );
    }

    parsec_zhetrf_sp1dplus = parsec_zhetrf_sp1dplus_new( A, sopalin_data, NULL, lwork );

    parsec_zhetrf_sp1dplus->_g_p_work = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( parsec_zhetrf_sp1dplus->_g_p_work,
                                lwork * sizeof(pastix_complex64_t) );

    parsec_arena_construct( parsec_zhetrf_sp1dplus->arenas[PARSEC_zhetrf_sp1dplus_DEFAULT_ARENA],
                            sizeof(pastix_complex64_t), PARSEC_ARENA_ALIGNMENT_SSE,
                            parsec_datatype_double_complex_t );

    parsec_matrix_add2arena_rect( parsec_zhetrf_sp1dplus->arenas[PARSEC_zhetrf_sp1dplus_CBLK_WS_ARENA],
                                  parsec_datatype_double_complex_t, offdmax, 1, offdmax );

#if defined(PASTIX_GENERATE_MODEL)
    parsec_zhetrf_sp1dplus->_g_forced_pushout = 1;
#endif
    return (parsec_taskpool_t*)parsec_zhetrf_sp1dplus;
}

/**
 *******************************************************************************
 *
 * @brief Free the data structure associated to a taskpool created with
 * parsec_zhetrf_sp1dplus_New().
 *
 *******************************************************************************
 *
 * @param[inout] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 ******************************************************************************/
void
parsec_zhetrf_sp1dplus_Destruct( parsec_taskpool_t *taskpool )
{
    parsec_zhetrf_sp1dplus_taskpool_t *parsec_zhetrf_sp1dplus = NULL;
    parsec_zhetrf_sp1dplus = (parsec_zhetrf_sp1dplus_taskpool_t *)taskpool;

    parsec_matrix_del2arena( parsec_zhetrf_sp1dplus->arenas[PARSEC_zhetrf_sp1dplus_CBLK_WS_ARENA] );

    parsec_private_memory_fini( parsec_zhetrf_sp1dplus->_g_p_work );
    free( parsec_zhetrf_sp1dplus->_g_p_work );

    parsec_taskpool_free( taskpool );
}

/**
 *******************************************************************************
 *
 * @brief Perform a sparse LDL^h factorization with 1D kernels.
 *
 * The function performs the LDL^h factorization of a sparse hermitian matrix A.
 * The factorization has the form
 *
 *    \f[ A = L \times D \times L^h \f]
 *
 * where L is a sparse lower triangular matrix unitary on the diagonal, and D a
 * diagonal matrix.
 *
 *******************************************************************************
 *
 * @param[inout] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[inout] A
 *          Descriptor of the sparse matrix A.
 *          On exit, A is overwritten with the factorized matrix.
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 *******************************************************************************
 *
 * @retval NULL if incorrect parameters are given.
 * @retval The parsec taskpool describing the operation that can be
 *         enqueued in the runtime with parsec_enqueue(). It, then, needs to be
 *         destroy with parsec_zhetrf_sp1dplus_Destruct().
 *
 *******************************************************************************
 *
 * @sa parsec_zhetrf_sp1dplus_New
 * @sa parsec_zhetrf_sp1dplus_Destruct
 *
 ******************************************************************************/
int
parsec_zhetrf_sp1dplus( parsec_context_t *parsec,
                        parsec_sparse_matrix_desc_t *A,
                        sopalin_data_t *sopalin_data )
{
    parsec_taskpool_t *parsec_zhetrf_sp1dplus = NULL;
    int info = 0;

    parsec_zhetrf_sp1dplus = parsec_zhetrf_sp1dplus_New( A, sopalin_data );

    if ( parsec_zhetrf_sp1dplus == NULL ) {
        return -1;
    }

    parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_zhetrf_sp1dplus );
    parsec_context_start( parsec );
    parsec_context_wait( parsec );
    parsec_zhetrf_sp1dplus_Destruct( parsec_zhetrf_sp1dplus );

    return info;
}

/**
 *******************************************************************************
 *
 * @brief Generate the PaRSEC taskpool object for the LDL^h factorization with
 * 1D and 2D kernels.
 *
 * The function only returns the object that describes the LDL^h factorization
 * of a sparse hermitian matrix A.
 * The factorization has the form
 *
 *    \f[ A = L \times D \times L^h \f]
 *
 * where L is a sparse lower triangular matrix unitary on the diagonal, and D a
 * diagonal matrix.
 *
 * In this object, operations are panel based operations for the lower levels of
 * the elimination tree, and the higher levels are taskpool by 2D tasks scheme to
 * create more parallelism and adapt to large architectures.
 *
 * @warning The computations are not done by this call.
 *
 *******************************************************************************
 *
 * @param[inout] A
 *          Descriptor of the sparse matrix A.
 *          On exit, A is overwritten with the factorized matrix.
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 *******************************************************************************
 *
 * @retval NULL if incorrect parameters are given.
 * @retval The parsec taskpool describing the operation that can be
 *         enqueued in the runtime with parsec_enqueue(). It, then, needs to be
 *         destroy with parsec_zhetrf_sp2d_Destruct().
 *
 *******************************************************************************
 *
 * @sa parsec_zhetrf_sp2d
 * @sa parsec_zhetrf_sp2d_Destruct
 * @sa parsec_zhetrf_sp2d
 * @sa parsec_chetrf_sp2d_New
 * @sa parsec_dhetrf_sp2d_New
 * @sa parsec_shetrf_sp2d_New
 *
 ******************************************************************************/
parsec_taskpool_t*
parsec_zhetrf_sp2d_New( parsec_sparse_matrix_desc_t *A,
                        sopalin_data_t *sopalin_data )
{
    parsec_zhetrf_sp2d_taskpool_t *parsec_zhetrf_sp2d = NULL;
    unsigned int gemmmax = sopalin_data->solvmtx->gemmmax;
    unsigned int offdmax = sopalin_data->solvmtx->offdmax;
    unsigned int blokmax = sopalin_data->solvmtx->blokmax;
    pastix_int_t lwork;

    lwork = pastix_imax( gemmmax, blokmax );
    if ( (sopalin_data->solvmtx->lowrank.compress_when != PastixCompressNever) &&
         (sopalin_data->solvmtx->lowrank.ilu_lvl < INT_MAX) )
    {
        lwork = pastix_imax( lwork, 2 * blokmax );
    }

    parsec_zhetrf_sp2d = parsec_zhetrf_sp2d_new( A, sopalin_data, NULL, lwork );

    parsec_zhetrf_sp2d->_g_p_work = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( parsec_zhetrf_sp2d->_g_p_work,
                                lwork * sizeof(pastix_complex64_t) );

    parsec_arena_construct( parsec_zhetrf_sp2d->arenas[PARSEC_zhetrf_sp2d_DEFAULT_ARENA],
                            sizeof(pastix_complex64_t), PARSEC_ARENA_ALIGNMENT_SSE,
                            parsec_datatype_double_complex_t );

    parsec_matrix_add2arena_rect( parsec_zhetrf_sp2d->arenas[PARSEC_zhetrf_sp2d_CBLK_WS_ARENA],
                                  parsec_datatype_double_complex_t,
                                  offdmax, 1, offdmax );

    parsec_matrix_add2arena_rect( parsec_zhetrf_sp2d->arenas[PARSEC_zhetrf_sp2d_BLOK_WS_ARENA],
                                  parsec_datatype_double_complex_t,
                                  blokmax, 1, blokmax );

#if defined(PASTIX_GENERATE_MODEL)
    parsec_zhetrf_sp2d->_g_forced_pushout = 1;
#endif
    return (parsec_taskpool_t*)parsec_zhetrf_sp2d;
}

/**
 *******************************************************************************
 *
 * @brief Free the data structure associated to a taskpool created with
 * parsec_zhetrf_sp2d_New().
 *
 *******************************************************************************
 *
 * @param[inout] taskpool
 *          On entry, the taskpool to destroy.
 *          On exit, the taskpool cannot be used anymore.
 *
 ******************************************************************************/
void
parsec_zhetrf_sp2d_Destruct( parsec_taskpool_t *taskpool )
{
    parsec_zhetrf_sp2d_taskpool_t *parsec_zhetrf_sp2d = NULL;
    parsec_zhetrf_sp2d = (parsec_zhetrf_sp2d_taskpool_t *)taskpool;

    parsec_matrix_del2arena( parsec_zhetrf_sp2d->arenas[PARSEC_zhetrf_sp2d_CBLK_WS_ARENA] );
    if ( parsec_zhetrf_sp2d->arenas[PARSEC_zhetrf_sp2d_BLOK_WS_ARENA]->elem_size > 0 ) {
        parsec_matrix_del2arena( parsec_zhetrf_sp2d->arenas[PARSEC_zhetrf_sp2d_BLOK_WS_ARENA] );
    }

    parsec_private_memory_fini( parsec_zhetrf_sp2d->_g_p_work );
    free( parsec_zhetrf_sp2d->_g_p_work );

    parsec_taskpool_free( taskpool );
}

/**
 *******************************************************************************
 *
 * @brief Perform a sparse LDL^h factorization with 1D and 2D kernels.
 *
 * The function performs the LDL^h factorization of a sparse hermitian matrix A.
 * The factorization has the form
 *
 *    \f[ A = L \times D \times L^h \f]
 *
 * where L is a sparse lower triangular matrix unitary on the diagonal, and D a
 * diagonal matrix.
 *
 *******************************************************************************
 *
 * @param[inout] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[inout] A
 *          Descriptor of the sparse matrix A.
 *          On exit, A is overwritten with the factorized matrix.
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 *******************************************************************************
 *
 * @retval NULL if incorrect parameters are given.
 * @retval The parsec taskpool describing the operation that can be
 *         enqueued in the runtime with parsec_enqueue(). It, then, needs to be
 *         destroy with parsec_zhetrf_sp2d_Destruct().
 *
 *******************************************************************************
 *
 * @sa parsec_zhetrf_sp2d_New
 * @sa parsec_zhetrf_sp2d_Destruct
 *
 ******************************************************************************/
int
parsec_zhetrf_sp2d( parsec_context_t *parsec,
                    parsec_sparse_matrix_desc_t *A,
                    sopalin_data_t *sopalin_data )
{
    parsec_taskpool_t *parsec_zhetrf_sp2d = NULL;
    int info = 0;

    parsec_zhetrf_sp2d = parsec_zhetrf_sp2d_New( A, sopalin_data );

    if ( parsec_zhetrf_sp2d == NULL ) {
        return -1;
    }

    parsec_context_add_taskpool( parsec, (parsec_taskpool_t*)parsec_zhetrf_sp2d );
    parsec_context_start( parsec );
    parsec_context_wait( parsec );
    parsec_zhetrf_sp2d_Destruct( parsec_zhetrf_sp2d );

    return info;
}

/**
 *******************************************************************************
 *
 * @brief Perform a sparse LDL^h factorization using PaRSEC runtime.
 *
 * The function performs the LDL^h factorization of a sparse hermitian matrix A.
 * The factorization has the form
 *
 *    \f[ A = L \times D \times L^h \f]
 *
 * where L is a sparse lower triangular matrix unitary on the diagonal, and D a
 * diagonal matrix.
 *
 * The algorithm is automatically chosen between the 1D and 2D version based on
 * the API parameter IPARM_TASKS2D_LEVEL. If IPARM_TASKS2D_LEVEL != 0
 * the 2D scheme is applied, the 1D otherwise.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 *******************************************************************************
 *
 * @sa parsec_zhetrf_sp1dplus
 * @sa parsec_zhetrf_sp2d
 *
 ******************************************************************************/
void
parsec_zhetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    parsec_sparse_matrix_desc_t *sdesc = sopalin_data->solvmtx->parsec_desc;
    parsec_context_t *ctx;

    /*
     * Start PaRSEC if not already started
     */
    if (pastix_data->parsec == NULL) {
        int argc = 0;
        pastix_parsec_init( pastix_data, &argc, NULL, NULL );
    }
    ctx = pastix_data->parsec;

    if ( sdesc == NULL ) {
        /* Create the matrix descriptor */
        parsec_sparse_matrix_init( sopalin_data->solvmtx,
                                   sizeof( pastix_complex64_t ), PastixHermitian,
                                   pastix_data->inter_node_procnbr,
                                   pastix_data->inter_node_procnum );
        sdesc = sopalin_data->solvmtx->parsec_desc;
    }

    /*
     * Select 1D or 2D jdf based on 2d tasks level
     */
    if ( pastix_data->iparm[IPARM_TASKS2D_LEVEL] != 0 )
    {
        parsec_zhetrf_sp2d( ctx, sdesc,
                            sopalin_data );
    }
    else {
        parsec_zhetrf_sp1dplus( ctx, sdesc,
                                sopalin_data );
    }

    return;
}

/**
 *@}
 */
