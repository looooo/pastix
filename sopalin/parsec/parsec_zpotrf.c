/**
 *
 * @file parsec_zpotrf.c
 *
 * PaStiX zpotrf PaRSEC wrapper.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @date 2013-06-24
 * @precisions normal z -> s d c
 *
 * @addtogroup parsec_potrf
 * @{
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

/**
 *******************************************************************************
 *
 * @brief Generate the PaRSEC handler object for the Cholesky factorization with
 * 1D kernels.
 *
 * The function only return the object that describes the Cholesky factorization
 * of a sparse symmetric positive definite (or Hermitian positive definite in
 * the complex case) matrix A.
 * The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H} \f]
 *
 * where L is a sparse lower triangular matrix.
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
 * @retval The parsec handle describing the operation that can be
 *         enqueued in the runtime with parsec_enqueue(). It, then, needs to be
 *         destroy with parsec_zpotrf_sp1dplus_Destruct().
 *
 *******************************************************************************
 *
 * @sa parsec_zpotrf_sp1dplus
 * @sa parsec_zpotrf_sp1dplus_Destruct
 * @sa parsec_zgetrf_sp1dplus
 * @sa parsec_cpotrf_sp1dplus_New
 * @sa parsec_dpotrf_sp1dplus_New
 * @sa parsec_spotrf_sp1dplus_New
 *
 ******************************************************************************/
parsec_handle_t*
parsec_zpotrf_sp1dplus_New( sparse_matrix_desc_t *A,
                            sopalin_data_t *sopalin_data )
{
    parsec_zpotrf_sp1dplus_handle_t *parsec_zpotrf_sp1dplus = NULL;

    parsec_zpotrf_sp1dplus = parsec_zpotrf_sp1dplus_new( A, sopalin_data, NULL );

    parsec_zpotrf_sp1dplus->_g_p_work = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( parsec_zpotrf_sp1dplus->_g_p_work,
                                sopalin_data->solvmtx->gemmmax * sizeof(pastix_complex64_t) );

    /* This is a default initializer for now, as it is not used in distributed */
    parsec_matrix_add2arena_rect( parsec_zpotrf_sp1dplus->arenas[PARSEC_zpotrf_sp1dplus_DEFAULT_ARENA],
                                  parsec_datatype_double_complex_t,
                                  /*sopalin_data->solvmtx->gemmmax*/ 1, 1, 1 );

    return (parsec_handle_t*)parsec_zpotrf_sp1dplus;
}

/**
 *******************************************************************************
 *
 * @brief Free the data structure associated to an handle created with
 * parsec_zpotrf_sp1dplus_New().
 *
 *******************************************************************************
 *
 * @param[inout] handle
 *          On entry, the handle to destroy.
 *          On exit, the handle cannot be used anymore.
 *
 ******************************************************************************/
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

/**
 *******************************************************************************
 *
 * @brief Perform a sparse Cholesky factorization with 1D kernels.
 *
 * The function performs the Cholesky factorization of a sparse symmetric
 * positive definite (or Hermitian positive definite in the complex case) matrix
 * A.
 * The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H} \f]
 *
 * where L is a sparse lower triangular matrix.
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
 * @retval The parsec handle describing the operation that can be
 *         enqueued in the runtime with parsec_enqueue(). It, then, needs to be
 *         destroy with parsec_zpotrf_sp1dplus_Destruct().
 *
 *******************************************************************************
 *
 * @sa parsec_zpotrf_sp1dplus_New
 * @sa parsec_zpotrf_sp1dplus_Destruct
 *
 ******************************************************************************/
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

/**
 *******************************************************************************
 *
 * @brief Generate the PaRSEC handler object for the Cholesky factorization with
 * 1D and 2D kernels.
 *
 * The function only return the object that describes the Cholesky factorization
 * of a sparse symmetric positive definite (or Hermitian positive definite in
 * the complex case) matrix A.
 * The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H} \f]
 *
 * where L is a sparse lower triangular matrix.
 *
 * In this object, operations are panel based operations for the lower levels of
 * the elimination tree, and the higher levels are handled by 2D tasks scheme to
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
 * @retval The parsec handle describing the operation that can be
 *         enqueued in the runtime with parsec_enqueue(). It, then, needs to be
 *         destroy with parsec_zpotrf_sp2d_Destruct().
 *
 *******************************************************************************
 *
 * @sa parsec_zpotrf_sp2d
 * @sa parsec_zpotrf_sp2d_Destruct
 * @sa parsec_zgetrf_sp2d
 * @sa parsec_cpotrf_sp2d_New
 * @sa parsec_dpotrf_sp2d_New
 * @sa parsec_spotrf_sp2d_New
 *
 ******************************************************************************/
parsec_handle_t*
parsec_zpotrf_sp2d_New( sparse_matrix_desc_t *A,
                        sopalin_data_t *sopalin_data )
{
    parsec_zpotrf_sp2d_handle_t *parsec_zpotrf_sp2d = NULL;

    parsec_zpotrf_sp2d = parsec_zpotrf_sp2d_new( A, sopalin_data, NULL );

    parsec_zpotrf_sp2d->_g_p_work = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( parsec_zpotrf_sp2d->_g_p_work,
                                sopalin_data->solvmtx->gemmmax * sizeof(pastix_complex64_t) );

    parsec_matrix_add2arena_rect( parsec_zpotrf_sp2d->arenas[PARSEC_zpotrf_sp2d_DEFAULT_ARENA],
                                  parsec_datatype_double_complex_t,
                                  /*sopalin_data->solvmtx->gemmmax*/ 1, 1, 1 );

    return (parsec_handle_t*)parsec_zpotrf_sp2d;
}

/**
 *******************************************************************************
 *
 * @brief Free the data structure associated to an handle created with
 * parsec_zpotrf_sp2d_New().
 *
 *******************************************************************************
 *
 * @param[inout] handle
 *          On entry, the handle to destroy.
 *          On exit, the handle cannot be used anymore.
 *
 ******************************************************************************/
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

/**
 *******************************************************************************
 *
 * @brief Perform a sparse Cholesky factorization with 1D and 2D kernels.
 *
 * The function performs the Cholesky factorization of a sparse symmetric
 * positive definite (or Hermitian positive definite in the complex case) matrix
 * A.
 * The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H} \f]
 *
 * where L is a sparse lower triangular matrix.
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
 * @retval The parsec handle describing the operation that can be
 *         enqueued in the runtime with parsec_enqueue(). It, then, needs to be
 *         destroy with parsec_zpotrf_sp2d_Destruct().
 *
 *******************************************************************************
 *
 * @sa parsec_zpotrf_sp2d_New
 * @sa parsec_zpotrf_sp2d_Destruct
 *
 ******************************************************************************/
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

/**
 *******************************************************************************
 *
 * @brief Perform a sparse Cholesky factorization.
 *
 * The function performs the Cholesky factorization of a sparse symmetric
 * positive definite (or Hermitian positive definite in the complex case) matrix
 * A.
 * The factorization has the form
 *
 *    \f[ A = \{_{L\times L^H} \f]
 *
 * where L is a sparse lower triangular matrix.
 *
 * The algorithm is automatically chosen between the 1D and 2D version based on
 * the API parameter IPARM_DISTRIBUTION_LEVEL. If IPARM_DISTRIBUTION_LEVEL >= 0
 * the 2D scheme is applied, the 1D otherwise.
 *
 *******************************************************************************
 *
 * @param[inout] parsec
 *          The parsec context of the application that will run the operation.
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 *******************************************************************************
 *
 * @retval NULL if incorrect parameters are given.
 * @retval The parsec handle describing the operation that can be
 *         enqueued in the runtime with parsec_enqueue(). It, then, needs to be
 *         destroy with parsec_zpotrf_sp2d_Destruct().
 *
 *******************************************************************************
 *
 * @sa parsec_zpotrf_sp2d_New
 * @sa parsec_zpotrf_sp2d_Destruct
 *
 ******************************************************************************/
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

/**
 *@}
 */
