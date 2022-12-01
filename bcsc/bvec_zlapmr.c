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
 * @date 2022-12-05
 * @precisions normal z -> c d s
 *
 * This file implements the function bvec_zlapmr with the following hierarchy:
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
#include "bcsc/bvec.h"
#include "bcsc_z.h"
#include "order/order_internal.h"
#include "cblas.h"
#include "blend/solver.h"

#if defined( PASTIX_WITH_MPI )
/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Apply a row permutation to a matrix A (LAPACK xlatmr) in the
 * distributed case.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *          If PastixDirForward, A is permuted into PA.
 *          If PastixDirBackward, PA is permuted into A.
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
 *          Referenced as input if dir is PastixDirForward, as output otherwise.
 *
 * @param[in] lda
 *          The leading dimension of A.
 *
 * @param[inout] PA
 *          The rhs structure of the permuted matrix A.
 *          Referenced as inout if dir is PastixDirForward, as input otherwise.
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
 * @ingroup bcsc_internal
 *
 * @brief Applies a row permutation (permtab) to the matrix b. and stores it
 * in Pb.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] m
 *          The number of rows in the matrix b.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[in] b
 *          A matrix of size ldb-by-n.
 *
 * @param[in] ldb
 *          The leading dimension of b >= m.
 *
 * @param[inout] Pb
 *          The rhs structure of the permuted matrix b.
 *          On entry, the structure is initialized. On exit, contains the
 *          permuted matrix b.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_zlapmr_rep_vec2bvec( const pastix_data_t      *pastix_data,
                          pastix_int_t              m,
                          pastix_int_t              nrhs,
                          const pastix_complex64_t *b,
                          pastix_int_t              ldb,
                          pastix_rhs_t              Pb  )
{
    pastix_complex64_t *pb;
    const spmatrix_t   *spm     = pastix_data->csc;
    pastix_int_t        dof     = spm->dof;
    const pastix_int_t *dofs    = spm->dofs;
    pastix_int_t        ldpb    = Pb->ld;
    pastix_int_t        idx_cnt = 0;
    pastix_int_t        val_cnt = 0;
    pastix_int_t        j, ig, ige, ilpe, dofi;
    bvec_data_amount_t *data;
    bvec_handle_comm_t *comm_rhs;

    /* Check on b */
    assert( m == spm->gNexp );
    assert( m <= ldb        );
    //assert( m == spm->nexp  ); // Uncomment with dist2dist
    assert( pastix_data->bcsc->n == Pb->m  );
    assert( pastix_data->bcsc->n == Pb->ld );

    /*
     * Goes through b to fill the data_comm with the data to send and
     * fills pb with the local data.
     */
    pb  = Pb->b;
    ige = 0;
    for ( ig = 0; ig < spm->gN; ig++, ige += dofi ) {
        dofi = ( dof > 0 ) ? dof : dofs[ ig+1 ] - dofs[ ig ];
        ilpe  = bvec_glob2Ploc( pastix_data, ig );

        if ( ilpe < 0 ) {
            continue;
        }

        idx_cnt += 1;
        val_cnt += dofi * nrhs;
        for ( j = 0; j < nrhs; j++ ) {
            memcpy( pb + ilpe + j * ldpb,
                    b  + ige  + j * ldb,
                    dofi * sizeof(pastix_complex64_t) );
        }
    }

    assert( idx_cnt <= val_cnt );

    bvec_handle_comm_init( pastix_data, Pb );

    comm_rhs  = Pb->rhs_comm;
    data = &(comm_rhs->data_comm[comm_rhs->clustnum].send);
    data->idxcnt = idx_cnt;
    data->valcnt = val_cnt;

    bvec_exchange_amount_rep( Pb->rhs_comm );

    (void)m;
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Applies a row permutation (permtab) to the matrix b. and stores it
 * in Pb.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] m
 *          The number of rows in the matrix b.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[inout] b
 *          A matrix of size ldb-by-n.
 *          On entry, the allocated matrix.
 *          On exit, contains the revers permutation of Pb.
 *
 * @param[in] ldb
 *          The leading dimension of b >= m.
 *
 * @param[input] Pb
 *          The rhs structure of the permuted matrix b.
 *
 *******************************************************************************/
static inline int
bvec_zlapmr_rep_bvec2vec( pastix_data_t      *pastix_data,
                          pastix_int_t        m,
                          pastix_int_t        nrhs,
                          pastix_complex64_t *b,
                          pastix_int_t        ldb,
                          pastix_rhs_t        Pb )
{
    pastix_complex64_t   *bp       = Pb->b;
    bvec_handle_comm_t   *comm_rhs = Pb->rhs_comm;
    const pastix_order_t *ord      = pastix_data->ordemesh;
    const spmatrix_t     *spm      = pastix_data->csc;
    pastix_int_t          dof      = spm->dof;
    const pastix_int_t   *dofs     = spm->dofs;
    pastix_int_t         *Ploc2Pglob = Pb->Ploc2Pglob;
    pastix_int_t          clustnum = pastix_data->solvmatr->clustnum;
    pastix_int_t          bcsc_n   = Pb->m;
    pastix_int_t          ilpe, ilp, igp, ig, ige;
    pastix_int_t          j, dofi;
    bvec_data_amount_t   *data_send;
    pastix_int_t         *idxptr;

    assert( Pb->m == pastix_data->bcsc->n );
    // assert( m     == spm->nexp ); Uncomment with dist2dist
    assert( b     != NULL );

    bvec_Ploc2Pglob( pastix_data, Pb );

    Ploc2Pglob = Pb->Ploc2Pglob;
    data_send  = &(comm_rhs->data_comm[clustnum].send);

    /*
     * Allocates the sending indexes buffer.
     */
    if ( ( data_send->idxcnt != 0 ) && ( data_send->idxbuf == NULL ) ) {
        MALLOC_INTERN( data_send->idxbuf, data_send->idxcnt, pastix_int_t );
    }

    /*
     * Sets the counters used to fill indexes and values buffers.
     */
    idxptr = data_send->idxbuf;

    /*
     * Goes through b to fill the data_comm with the data to send and
     * fills pb with the local data.
     */
    ilp  = 0;
    dofi = ( dof > 0 ) ? dof : dofs[ 1 ] - dofs[ 0 ];
    for ( ilpe = 0; ilpe < bcsc_n; ilpe += dofi, ilp ++ ) {
        igp  = Ploc2Pglob[ ilp ];
        ig   = ord->peritab[ igp ];
        ige  = ( dof > 0 ) ? ig * dof : dofs[ ig ];
        dofi = ( dof > 0 ) ? dof : dofs[ ig+1 ] - dofs[ ig ];

        /* Stores the indexes to send to c: (ig, j). */
        *idxptr = ig;
        idxptr++;

        /* Stores the values to send to c. */
        for ( j = 0; j < nrhs; j++ ) {
            memcpy( b + ige + j * ldb, bp + ilpe + j * Pb->ld, dofi * sizeof(pastix_complex64_t) );
        }
    }
    assert( (idxptr - data_send->idxbuf) == data_send->idxcnt );

    bvec_zexchange_data_rep( pastix_data, nrhs, b, ldb, Pb );

    bvec_handle_comm_exit( Pb->rhs_comm );
    memFree_null( Pb->rhs_comm );

    (void)m;
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Apply a row permutation to a matrix A (LAPACK xlatmr) in the
 * replicated case.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *          If PastixDirForward, A is permuted into PA.
 *          If PastixDirBackward, PA is permuted into A.
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
 *          Referenced as input if dir is PastixDirForward, as output otherwise.
 *
 * @param[in] lda
 *          The leading dimension of A.
 *
 * @param[inout] PA
 *          The rhs structure of the permuted matrix A.
 *          Referenced as inout if dir is PastixDirForward, as input otherwise.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_zlapmr_rep( pastix_data_t      *pastix_data,
                 pastix_dir_t        dir,
                 pastix_int_t        m,
                 pastix_int_t        n,
                 pastix_complex64_t *A,
                 pastix_int_t        lda,
                 pastix_rhs_t        PA )
{
    if ( dir == PastixDirForward ) {
        return bvec_zlapmr_rep_vec2bvec( pastix_data, m, n, A, lda, PA );
    }
    else {
        return bvec_zlapmr_rep_bvec2vec( pastix_data, m, n, A, lda, PA );
    }
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Apply a row permutation to a matrix A (LAPACK xlatmr) in the shared
 * memory case.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *          If PastixDirForward, A is permuted into PA.
 *          If PastixDirBackward, PA is permuted into A.
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
 *          Referenced as input if dir is PastixDirForward, as output otherwise.
 *
 * @param[in] lda
 *          The leading dimension of A.
 *
 * @param[inout] PA
 *          The rhs structure of the permuted matrix A.
 *          Referenced as inout if dir is PastixDirForward, as input otherwise.
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
 * @ingroup bcsc_internal
 *
 * @brief Apply a row permutation to a matrix A (LAPACK xlatmr)
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *          If PastixDirForward, A is permuted into PA.
 *          If PastixDirBackward, PA is permuted into A.
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
 *          Referenced as input if dir is PastixDirForward, as output otherwise.
 *
 * @param[in] lda
 *          The leading dimension of A.
 *
 * @param[inout] PA
 *          The rhs structure of the permuted matrix A.
 *          Referenced as inout if dir is PastixDirForward, as input otherwise.
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
    const spmatrix_t    *spm      = pastix_data->csc;
    const pastix_bcsc_t *bcsc     = pastix_data->bcsc;
    const SolverMatrix  *solvmatr = pastix_data->solvmatr;
    int rc;

    assert( lda >= m );
    if ( dir == PastixDirForward ) {
        PA->flttype = PastixComplex64;
        PA->m       = bcsc->n;
        PA->n       = n;

        if ( solvmatr->clustnbr > 1 ) {
            PA->allocated = 1;
            PA->ld        = PA->m;
            PA->b         = malloc( PA->ld * PA->n * pastix_size_of( PA->flttype ) );
        }
        else {
            assert( m == PA->m );
            PA->allocated = 0;
            PA->ld        = lda;
            PA->b         = A;
        }
    }
#if !defined(NDEBUG)
    else {
        assert( PA->allocated >= 0            );
        assert( PA->flttype   == PastixComplex64 );
        assert( PA->m         == bcsc->n      );
        assert( PA->n         == n            );

        if ( PA->allocated == 0 )
        {
            assert( PA->b  == A   );
            assert( PA->ld == lda );
        }
        else {
            assert( PA->b  != A     );
            assert( PA->ld == PA->m );
        }
    }
#endif

#if defined(PASTIX_WITH_MPI)
    if ( spm->clustnbr > 1 ) {
        #if 0
        if ( spm->loc2glob != NULL ) {
            /*
             * The input vector is distributed, we redispatch it following the
             * ordering.
             */
            rc = bvec_zlapmr_dst( pastix_data, dir, m, n, A, lda, PA );
        }
        else
        #endif
        {
            /*
             * The input vector is replicated, we extract or collect the data
             * following the ordering.
             */
            rc = bvec_zlapmr_rep( pastix_data, dir, m, n, A, lda, PA );
        }
    }
    else
#endif
    {
        /*
         * We are in the shared memory case, the permutation is applied in place.
         */
        rc = bvec_zlapmr_shm( pastix_data, dir, m, n, A, lda, PA );
        (void)spm;
    }

    if ( dir == PastixDirBackward ) {
        if ( PA->allocated > 0 ) {
            memFree_null( PA->b );
        }

        PA->allocated = -1;
        PA->flttype   = PastixPattern;
        PA->m         = -1;
        PA->n         = -1;
        PA->ld        = -1;
        PA->b         = NULL;
    }

    return rc;
}
