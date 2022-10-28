/**
 *
 * @file bvec_zlapmr.c
 *
 * Functions to compute the permutation of the right hand side.
 *
 * @copyright 2004-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @author Theophile Terraz
 * @author Tony Delarue
 * @author Vincent Bridonneau
 * @date 2022-07-07
 * @precisions normal z -> c d s
 *
 * This file implements the function bvec_zlapmr with the following hiearchy:
 *
 * bvec_zlapmr():
 *    - bvec_zlapmr_shm() for shared memory case
 *    - bvec_zlapmr_rep() for replicated rhs case
 *         - bvec_zlapmr_rep_{vec2bvec,bvec2vec}()
 *    - bvec_zlapmr_dst() for distributed rhs case
 *         - bvec_zlapmr_dst_{vec2bvec,bvec2vec}()
 *
 **/
#include "common.h"
#include <math.h>
#include "lapacke.h"
#include "bcsc/bcsc.h"
#include "bcsc_z.h"
#include "order/order_internal.h"
#include "cblas.h"
#include "blend/solver.h"

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Apply a row permutation to a matrix A (LAPACK xlatmr)
 *
 *******************************************************************************
 *
 * @param[in] thread_safe
 *          Boolean to switch between the thread-safe implementation that
 *          exploits an additional workspace, or the non thread-safe version
 *          that has no memory overhead.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *
 * @param[in] m
 *          The number of rows in the matrix A, and the number of elements in
 *          perm.
 *
 * @param[in] n
 *          The number of columns in the matrix A.
 *
 * @param[inout] A
 *          A matrix of size lda-by-n.
 *          On exit, rowas are permuted and A contains P A.
 *
 * @param[in] lda
 *          The leading dimension of A.
 *
 * @param[inout] perm
 *          The permutation array. Must be 0 based.
 *          If thread_safe is true, then perm array is used only as input, and a
 *          temporary array is allocated to follow the cycles. If thread_safe is
 *          false, then perm array is modified during the swap and restored at
 *          the end of the call.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_zlapmr_dst( __attribute__((unused)) pastix_data_t      *pastix_data,
                 __attribute__((unused)) pastix_dir_t        dir,
                 __attribute__((unused)) pastix_int_t        m,
                 __attribute__((unused)) pastix_int_t        n,
                 __attribute__((unused)) pastix_complex64_t *A,
                 __attribute__((unused)) pastix_int_t        lda,
                 __attribute__((unused)) pastix_rhs_t        PA )
{
    assert( 0 );
    if ( dir == PastixDirForward ) {
        return 0; //bvec_zlapmr_dst_vec2bvec( pastix_data, m, n, A, lda, PA );
    }
    else {
        return 0; //bvec_zlapmr_dst_bvec2vec( pastix_data, m, n, A, lda, PA );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Apply a row permutation to a matrix A (LAPACK xlatmr)
 *
 *******************************************************************************
 *
 * @param[in] thread_safe
 *          Boolean to switch between the thread-safe implementation that
 *          exploits an additional workspace, or the non thread-safe version
 *          that has no memory overhead.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *
 * @param[in] m
 *          The number of rows in the matrix A, and the number of elements in
 *          perm.
 *
 * @param[in] n
 *          The number of columns in the matrix A.
 *
 * @param[inout] A
 *          A matrix of size lda-by-n.
 *          On exit, rowas are permuted and A contains P A.
 *
 * @param[in] lda
 *          The leading dimension of A.
 *
 * @param[inout] perm
 *          The permutation array. Must be 0 based.
 *          If thread_safe is true, then perm array is used only as input, and a
 *          temporary array is allocated to follow the cycles. If thread_safe is
 *          false, then perm array is modified during the swap and restored at
 *          the end of the call.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_zlapmr_rep( __attribute__((unused)) pastix_data_t      *pastix_data,
                 __attribute__((unused)) pastix_dir_t        dir,
                 __attribute__((unused)) pastix_int_t        m,
                 __attribute__((unused)) pastix_int_t        n,
                 __attribute__((unused)) pastix_complex64_t *A,
                 __attribute__((unused)) pastix_int_t        lda,
                 __attribute__((unused)) pastix_rhs_t        PA )
{
    assert( 0 );
    if ( dir == PastixDirForward ) {
        return 0; //bvec_zlapmr_rep_vec2bvec( pastix_data, m, n, A, lda, PA );
    }
    else {
        return 0; //bvec_zlapmr_rep_bvec2vec( pastix_data, m, n, A, lda, PA );
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Apply a row permutation to a matrix A (LAPACK xlatmr)
 *
 *******************************************************************************
 *
 * @param[in] thread_safe
 *          Boolean to switch between the thread-safe implementation that
 *          exploits an additional workspace, or the non thread-safe version
 *          that has no memory overhead.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *
 * @param[in] m
 *          The number of rows in the matrix A, and the number of elements in
 *          perm.
 *
 * @param[in] n
 *          The number of columns in the matrix A.
 *
 * @param[inout] A
 *          A matrix of size lda-by-n.
 *          On exit, rowas are permuted and A contains P A.
 *
 * @param[in] lda
 *          The leading dimension of A.
 *
 * @param[inout] perm
 *          The permutation array. Must be 0 based.
 *          If thread_safe is true, then perm array is used only as input, and a
 *          temporary array is allocated to follow the cycles. If thread_safe is
 *          false, then perm array is modified during the swap and restored at
 *          the end of the call.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_zlapmr_shm( pastix_data_t      *pastix_data,
                 pastix_dir_t        dir,
                 pastix_int_t        m,
                 pastix_int_t        n,
                 pastix_complex64_t *A,
                 pastix_int_t        lda,
                 pastix_rhs_t        PA )
{
    pastix_complex64_t  tmp;
    pastix_int_t        i, j, k, jj;
    pastix_int_t       *perm, *perm_cpy;
    int                 thread_safe = pastix_data->iparm[IPARM_APPLYPERM_WS];

    if ( PA->b != A ) {
        pastix_print_error( "Incorrect definition of the right hand side for in place permutation\n" );
        return PASTIX_ERR_BADPARAMETER;
    }
    assert( PA->allocated == 0 );
    assert( PA->m  == m   );
    assert( PA->n  == n   );
    assert( PA->ld == lda );

    perm = orderGetExpandedPeritab( pastix_data->ordemesh, pastix_data->csc );
    assert( perm != NULL );

    if ( thread_safe ) {
        perm_cpy = malloc( m * sizeof(pastix_int_t) );
        memcpy( perm_cpy, perm, m * sizeof(pastix_int_t) );
    }
    else {
        perm_cpy = perm;
    }

    if ( dir == PastixDirBackward ) {
        for( k = 0; k < m; k++ ) {
            i = k;
            j = perm_cpy[i];

            /* Cycle already seen */
            if ( j < 0 ) {
                continue;
            }

            /* Mark the i^th element as being seen */
            perm_cpy[i] = -j-1;

            while( j != k ) {

                for( jj = 0; jj < n; jj++ ) {
                    tmp             = A[j + jj * lda];
                    A[j + jj * lda] = A[k + jj * lda];
                    A[k + jj * lda] = tmp;
                }

                i = j;
                j = perm_cpy[i];
                perm_cpy[i] = -j-1;

                assert( (j != i) && (j >= 0) );
            }
        }
    }
    else {
        for( k = 0; k < m; k++ ) {
            i = k;
            j = perm_cpy[i];
            perm_cpy[i] = -j-1;

            /* Cycle already seen */
            if ( j < 0 ) {
                continue;
            }

            i = perm_cpy[j];

            /* Mark the i^th element as being seen */
            while( i >= 0 ) {

                for( jj = 0; jj < n; jj++ ) {
                    tmp             = A[j + jj * lda];
                    A[j + jj * lda] = A[i + jj * lda];
                    A[i + jj * lda] = tmp;
                }

                perm_cpy[j] = -i-1;
                j = i;
                i = perm_cpy[j];

                assert( j != i );
            }
        }
    }

    if ( thread_safe ) {
        free( perm_cpy );
    }
    else {
        /* Restore perm array */
        for( k = 0; k < m; k++ ) {
            assert(perm[k] < 0);
            perm[k] = - perm[k] - 1;
        }
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Apply a row permutation to a matrix A (LAPACK xlatmr)
 *
 *******************************************************************************
 *
 * @param[in] thread_safe
 *          Boolean to switch between the thread-safe implementation that
 *          exploits an additional workspace, or the non thread-safe version
 *          that has no memory overhead.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *
 * @param[in] m
 *          The number of rows in the matrix A, and the number of elements in
 *          perm.
 *
 * @param[in] n
 *          The number of columns in the matrix A.
 *
 * @param[inout] A
 *          A matrix of size lda-by-n.
 *          On exit, rowas are permuted and A contains P A.
 *
 * @param[in] lda
 *          The leading dimension of A.
 *
 * @param[inout] perm
 *          The permutation array. Must be 0 based.
 *          If thread_safe is true, then perm array is used only as input, and a
 *          temporary array is allocated to follow the cycles. If thread_safe is
 *          false, then perm array is modified during the swap and restored at
 *          the end of the call.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
bvec_zlapmr( pastix_data_t      *pastix_data,
             pastix_dir_t        dir,
             pastix_int_t        m,
             pastix_int_t        n,
             pastix_complex64_t *A,
             pastix_int_t        lda,
             pastix_rhs_t        PA )
{
    const spmatrix_t *spm = pastix_data->csc;
    int rc;

// #if defined(PASTIX_WITH_MPI)
//     if ( spm->clustnbr > 1 ) {
//         if ( spm->loc2glob != NULL ) {
//             /*
//              * The input vector is distributed, we redispatch it following the
//              * ordering.
//              */
//             rc = bvec_zlapmr_dst( pastix_data, dir, m, n, A, lda, PA );
//         }
//         else {
//             /*
//              * The input vector is replicated, we extract or collect the data
//              * following the ordering.
//              */
//             rc = bvec_zlapmr_rep( pastix_data, dir, m, n, A, lda, PA );
//         }
//     }
//     else
// #endif
    {
        /*
         * We are in the shared memory case, the permutation is applied in place.
         */
        rc = bvec_zlapmr_shm( pastix_data, dir, m, n, A, lda, PA );
        (void)spm;
    }

    return rc;
}
