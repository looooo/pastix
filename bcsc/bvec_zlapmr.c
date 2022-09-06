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
#include "bcsc/bvec.h"
#include "bcsc_z.h"
#include "order/order_internal.h"
#include "cblas.h"
#include "blend/solver.h"

#if defined( PASTIX_WITH_MPI )
/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Applies a row permutation (permtab) to the matrix b. It also
 * sends and receives the part of b according to the block repartition.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *
 * @param[in] replicated
 *          True if the vector b is replicated on all the process, false if the
 *          vector is distributed.
 *
 * @param[in] m
 *          The number of rows in the matrix b.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[inout] b
 *          The matrix b ldb-by-nrhs.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 *******************************************************************************
 *
 * @retval pb which correspond to the vector b permuted and with
 * the correct local data.
 *
 *******************************************************************************/
static inline pastix_int_t
bvec_zglob2Ploc( const pastix_data_t *pastix_data,
                 pastix_int_t         ig )
{
    const spmatrix_t     *spm      = pastix_data->csc;
    const pastix_bcsc_t  *bcsc     = pastix_data->bcsc;
    const pastix_order_t *ord      = pastix_data->ordemesh;
    const SolverMatrix   *solvmatr = pastix_data->solvmatr;
    const pastix_int_t   *col2cblk = bcsc->col2cblk;
    pastix_int_t          basespm  = spm->baseval;
    pastix_int_t          dof      = spm->dof;
    const pastix_int_t   *dofs     = spm->dofs;
    const SolverCblk     *cblk;
    pastix_int_t          igp, igpe, ipe, cblknum;

    assert( ord->baseval == 0 );

    igp     = ord->permtab[ ig ];
    igpe    = ( dof > 0 ) ? igp * dof : dofs[ ig ] - basespm/* vdof incorect */;
    cblknum = col2cblk[ igpe ];
    if ( cblknum >= 0 ) {
        cblk = solvmatr->cblktab + cblknum;
        ipe  = cblk->lcolidx + igpe - cblk->fcolnum;
    }
    else {
        ipe = cblknum;
    }

    return ipe;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Applies a row permutation (permtab) to the matrix b. It also
 * sends and receives the part of b according to the block repartition.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *
 * @param[in] replicated
 *          True if the vector b is replicated on all the process, false if the
 *          vector is distributed.
 *
 * @param[in] m
 *          The number of rows in the matrix b.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[inout] b
 *          The matrix b ldb-by-nrhs.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 *******************************************************************************
 *
 * @retval pb which correspond to the vector b permuted and with
 * the correct local data.
 *
 *******************************************************************************/
static inline pastix_int_t
bvec_zPglob2loc( pastix_data_t *pastix_data,
                 int            replicated,
                 pastix_int_t   idx_Pglob )
{
    const spmatrix_t *spm          = pastix_data->csc;
    pastix_order_t   *ord          = pastix_data->ordemesh;
    pastix_int_t     *perm         = NULL;
    pastix_int_t      baseval_ord  = ord->baseval;
    pastix_int_t      idx_loc;

    perm = orderGetExpandedPeritab( ord, spm );

    idx_loc = perm[ idx_Pglob ] - baseval_ord;

    (void)replicated;
    return idx_loc;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Applies a row permutation (permtab) to the matrix b. It also
 * sends and receives the part of b according to the block repartition.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *
 * @param[in] replicated
 *          True if the vector b is replicated on all the process, false if the
 *          vector is distributed.
 *
 * @param[inout] Pb
 *          The number of columns in the matrix b.
 *
 *******************************************************************************
 *
 * @retval pb which correspond to the vector b permuted and with
 * the correct local data.
 *
 *******************************************************************************/
static inline int
bvec_zcompute_amount_of_data_init( pastix_data_t *pastix_data,
                                   pastix_dir_t   dir,
                                   int            replicated,
                                   pastix_rhs_t   Pb )
{
    const SolverMatrix *solvmtx     = pastix_data->solvmatr;
    bvec_handle_comm_t *rhs_comm    = NULL;
    bvec_proc_comm_t   *data_comm   = NULL;
    pastix_int_t        bcsc_n      = Pb->m;
    pastix_int_t        nrhs        = Pb->n;
    pastix_int_t        size;

    /* Initializes bvec_comm. */
    size = sizeof(bvec_handle_comm_t) + (solvmtx->clustnbr-1)*sizeof(bvec_proc_comm_t);

    Pb->rhs_comm = (bvec_handle_comm_t *)malloc( size );
    rhs_comm = Pb->rhs_comm;

    rhs_comm->flttype  = Pb->flttype;
    rhs_comm->clustnbr = solvmtx->clustnbr;
    rhs_comm->clustnum = solvmtx->clustnum;
    rhs_comm->comm     = solvmtx->solv_comm;

    memset( rhs_comm->data_comm, 0, rhs_comm->clustnbr * sizeof(bvec_proc_comm_t) );
    data_comm = rhs_comm->data_comm;

    /* Case backward + replicated. */
    assert( replicated == 1 );

    /* Sends the same amout of data to all process. */
    data_comm[ rhs_comm->clustnum ].send.idxcnt = bcsc_n;
    data_comm[ rhs_comm->clustnum ].send.valcnt = nrhs * bcsc_n;
    (void)dir;
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Applies a row permutation (permtab) to the matrix b. It also
 * sends and receives the part of b according to the block repartition.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[inout] b
 *          The matrix b lda-by-nrhs.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 *******************************************************************************
 *
 * @retval pb which correspond to the vector b permuted and with
 * the correct local data.
 *
 *******************************************************************************/
static inline int
bvec_zexchange_amount_of_data( bvec_handle_comm_t *rhs_comm,
                               int                 replicated )
{
    bvec_proc_comm_t *data_comm   = rhs_comm->data_comm;
    pastix_int_t      clustnbr    = rhs_comm->clustnbr;
    pastix_int_t      clustnum    = rhs_comm->clustnum;
    pastix_int_t      counter_req = 0;
    pastix_int_t      c;
    MPI_Status        statuses[(clustnbr-1)*2];
    MPI_Request       requests[(clustnbr-1)*2];

    /* Receives the amount of indexes and values. */
    for ( c = 0; c < clustnbr; c++ ) {
        data_comm = rhs_comm->data_comm + c;
        if ( c == clustnum ) {
            continue;
        }

        MPI_Irecv( &(data_comm->recv), 2, PASTIX_MPI_INT,
                   c, PastixTagAmount, rhs_comm->comm, &requests[counter_req++] );

        if ( replicated == 1 ) {
            data_comm = rhs_comm->data_comm + clustnum;
        }
        MPI_Isend( &(data_comm->send), 2, PASTIX_MPI_INT,
                   c, PastixTagAmount, rhs_comm->comm, &requests[counter_req++] );
    }

    MPI_Waitall( counter_req, requests, statuses );

    /* Allocates the indexes and values buffers. */
    for ( c = 0; c < clustnbr; c++ ) {
        bvec_proc_comm_t   *data = rhs_comm->data_comm + c;
        bvec_data_amount_t *send = &( data->send );
        bvec_data_amount_t *recv = &( data->recv );

        if ( ( c == clustnum ) ) {
            if ( ( send->idxcnt != 0 ) && ( send->idxbuf == NULL ) ) {
                MALLOC_INTERN( send->idxbuf, send->idxcnt, pastix_int_t );
            }
            if ( ( send->valcnt != 0 ) && ( send->valbuf == NULL ) ) {
                MALLOC_INTERN( send->valbuf, send->valcnt, pastix_complex64_t );
            }
            continue;
        }

        if ( ( recv->idxcnt != 0 ) && ( recv->idxbuf == NULL ) ) {
            MALLOC_INTERN( recv->idxbuf, recv->idxcnt, pastix_int_t );
        }
        if ( ( recv->valcnt != 0 ) && ( recv->valbuf == NULL ) ) {
            MALLOC_INTERN( recv->valbuf, recv->valcnt, pastix_complex64_t );
        }
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Applies a row permutation (permtab) to the matrix b. It also
 * sends and receives the part of b according to the block repartition.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *
 * @param[in] replicated
 *          True if the vector b is replicated on all the process, false if the
 *          vector is distributed.
 *
 * @param[in] m
 *          The number of rows in the matrix b.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[inout] b
 *          The matrix b ldb-by-nrhs.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 *******************************************************************************
 *
 * @retval pb which correspond to the vector b permuted and with
 * the correct local data.
 *
 *******************************************************************************/
static inline int
bvec_zcompute_and_exchange_amount_of_data( pastix_data_t      *pastix_data,
                                           pastix_dir_t        dir,
                                           int                 replicated,
                                           pastix_rhs_t        Pb )
{
    const spmatrix_t *spm      = pastix_data->csc;
    pastix_int_t     *col2cblk = pastix_data->bcsc->col2cblk;
    pastix_int_t     *cblk2col = NULL;
    pastix_int_t      bcsc_n   = Pb->m;
    pastix_int_t      k, kl;
    if ( ( !( pastix_data->steps & STEP_SOLVE ) || !( pastix_data->steps & STEP_REFINE ) )
       && ( ( ( replicated == 0 ) && ( dir == PastixDirForward ) ) ||
            ( ( replicated == 1 ) && ( dir == PastixDirBackward ) ) ) )
    {
        bvec_zcompute_amount_of_data_init( pastix_data, dir, replicated, Pb );
        bvec_zexchange_amount_of_data( Pb->rhs_comm, replicated );
        /*
         * Creates cblk2col with col2cblk: gives the local index corresponding to the
         * global one (cblk2col[ig] = il).
         */
        if ( bcsc_n > 0 ) {
            MALLOC_INTERN( cblk2col, bcsc_n, pastix_int_t );
            kl = 0;
            for ( k = 0; k < spm->gNexp; k++ ) {
                if ( col2cblk[ k ] >= 0 ) {
                    cblk2col[ kl ] = k;
                    kl ++;
                }
            }
            assert( kl == bcsc_n );
        }
        Pb->cblk2col = cblk2col;
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Applies a row permutation (permtab) to the matrix b. It also
 * sends and receives the part of b according to the block repartition.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[inout] b
 *          The matrix b lda-by-nrhs.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 *******************************************************************************
 *
 * @retval pb which correspond to the vector b permuted and with
 * the correct local data.
 *
 *******************************************************************************/
static inline int
bvec_zexchange_data( bvec_handle_comm_t *rhs_comm,
                     int                 replicated )
{
    bvec_proc_comm_t *data_comm   = rhs_comm->data_comm;
    pastix_int_t      clustnbr    = rhs_comm->clustnbr;
    pastix_int_t      clustnum    = rhs_comm->clustnum;
    pastix_int_t      counter_req = 0;
    pastix_int_t      c;
    MPI_Status        statuses[(clustnbr-1)*4];
    MPI_Request       requests[(clustnbr-1)*4];

    for ( c = 0; c < clustnbr; c++ ) {
        data_comm = rhs_comm->data_comm + c;
        if ( c == clustnum ) {
            continue;
        }

        /* Posts the receptions of the indexes. */
        if ( data_comm->recv.idxcnt != 0 ) {
            MPI_Irecv( data_comm->recv.idxbuf, data_comm->recv.idxcnt,
                       PASTIX_MPI_INT, c, PastixTagIndexes, rhs_comm->comm, &requests[counter_req++] );
        }
        if ( data_comm->recv.valcnt != 0 ) {
            MPI_Irecv( data_comm->recv.valbuf, data_comm->recv.valcnt,
                       PASTIX_MPI_COMPLEX64, c, PastixTagValues, rhs_comm->comm, &requests[counter_req++] );
        }

        data_comm = rhs_comm->data_comm + clustnum;
        /* Posts the emissions of the indexes. */
        if ( data_comm->send.idxcnt != 0 ) {
            MPI_Isend( data_comm->send.idxbuf,  data_comm->send.idxcnt,
                       PASTIX_MPI_INT, c, PastixTagIndexes, rhs_comm->comm, &requests[counter_req++] );
        }
        if ( data_comm->send.valcnt != 0 ) {
            MPI_Isend( data_comm->send.valbuf, data_comm->send.valcnt,
                       PASTIX_MPI_COMPLEX64, c, PastixTagValues, rhs_comm->comm, &requests[counter_req++] );
        }
    }

    MPI_Waitall( counter_req, requests, statuses );

    (void)replicated;
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Applies a row permutation (permtab) to the matrix b. It also
 * sends and receives the part of b according to the block repartition.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *
 * @param[in] replicated
 *          True if the vector b is replicated on all the process, false if the
 *          vector is distributed.
 *
 * @param[in] m
 *          The number of rows in the matrix b.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[inout] b
 *          The matrix b ldb-by-nrhs.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 *******************************************************************************
 *
 * @retval pb which correspond to the vector b permuted and with
 * the correct local data.
 *
 *******************************************************************************/
static inline int
bvec_zhandle_received_data( pastix_data_t      *pastix_data,
                            pastix_dir_t        dir,
                            int                 replicated,
                            pastix_int_t        nrhs,
                            pastix_complex64_t *b,
                            pastix_int_t        ldb,
                            pastix_int_t       *indexes,
                            pastix_complex64_t *values,
                            pastix_int_t        size )
{
    const spmatrix_t *spm        = pastix_data->csc;
    pastix_int_t      dof        = spm->dof;
    pastix_int_t     *dofs       = spm->dofs;
    size_t            size_alloc = sizeof(pastix_complex64_t);
    pastix_int_t      dofi       = 1;
    pastix_int_t      ig, il, idx, j;

    for ( idx = 0; idx < size; idx++, indexes++ ) {
        ig   = indexes[ 0 ];
        if ( dir == PastixDirForward ) {
            il   = bvec_zglob2Ploc( pastix_data, ig );
            dofi = ( dof > 0 ) ? dof : dofs[ ig + 1 ] - dofs[ ig ];
            size_alloc = dofi * sizeof(pastix_complex64_t);
        }
        else {
            il = bvec_zPglob2loc( pastix_data, replicated, ig );
        }
        for ( j = 0; j < nrhs; j++, values += dofi ) {
            memcpy( b + il + j * ldb, values, size_alloc );
        }
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Applies a row permutation (permtab) to the matrix b. It also
 * sends and receives the part of b according to the block repartition.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *
 * @param[in] replicated
 *          True if the vector b is replicated on all the process, false if the
 *          vector is distributed.
 *
 * @param[in] m
 *          The number of rows in the matrix b.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[inout] b
 *          The matrix b ldb-by-nrhs.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 *******************************************************************************
 *
 * @retval pb which correspond to the vector b permuted and with
 * the correct local data.
 *
 *******************************************************************************/
static inline int
bvec_zexchange_and_handle_data( pastix_data_t      *pastix_data,
                                pastix_dir_t        dir,
                                int                 replicated,
                                pastix_complex64_t *b,
                                pastix_int_t        ldb,
                                pastix_rhs_t        Pb )
{
    bvec_handle_comm_t *rhs_comm  = Pb->rhs_comm;
    bvec_proc_comm_t   *data_comm = rhs_comm->data_comm;
    pastix_int_t        clustnbr  = rhs_comm->clustnbr;
    pastix_int_t        clustnum  = rhs_comm->clustnum;
    pastix_complex64_t *vect      = ( dir == PastixDirForward ) ? Pb->b  : b;
    pastix_int_t        ld_vect   = ( dir == PastixDirForward ) ? Pb->ld : ldb;
    pastix_int_t        c;

    bvec_zexchange_data( rhs_comm, replicated );

    for ( c = 0; c < clustnbr; c++ ) {
        data_comm = rhs_comm->data_comm + c;
        if ( c == clustnum ) {
            continue;
        }
        bvec_zhandle_received_data( pastix_data, dir, replicated, Pb->n, vect, ld_vect,
                                    data_comm->recv.idxbuf, data_comm->recv.valbuf,
                                    data_comm->recv.idxcnt );
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
 * @brief Applies a row permutation (permtab) to the matrix b. It also
 * sends and receives the part of b according to the block repartition.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[inout] b
 *          The matrix b lda-by-nrhs.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 *******************************************************************************
 *
 * @retval pb which correspond to the vector b permuted and with
 * the correct local data.
 *
 *******************************************************************************/
static inline int
bvec_zlapmr_rep_vec2bvec( pastix_data_t      *pastix_data,
                          pastix_int_t        m,
                          pastix_int_t        nrhs,
                          pastix_complex64_t *b,
                          pastix_int_t        ldb,
                          pastix_rhs_t        Pb  )
{
    pastix_complex64_t  *pb;
    const spmatrix_t    *spm    = pastix_data->csc;
    const pastix_bcsc_t *bcsc   = pastix_data->bcsc;
    pastix_int_t         dof    = spm->dof;
    pastix_int_t        *dofs   = spm->dofs;
    pastix_int_t         bcsc_n = bcsc->n;
    pastix_int_t         ig, ige, dofi;
    pastix_int_t         j, k, ldpb = Pb->ld;

    /* Check on b */
    assert( m      == spm->gNexp );
    // assert( m      == spm->nexp  ); Uncomment with dist2dist
    assert( ldb    >= m          );
    assert( bcsc_n == Pb->m      );
    assert( bcsc_n == Pb->ld     );

    /*
     * Goes through b to fill the data_comm with the data to send and
     * fills pb with the local data.
     */
    pb = Pb->b;
    for ( ig = 0, ige = 0; ige < m; ige += dofi, ig++ ) {
        k    = bvec_zglob2Ploc( pastix_data, ig );
        dofi = ( dof > 0 ) ? dof : dofs[ ig+1 ] - dofs[ ig ];

        if ( k < 0 ) {
            continue;
        }

        for ( j = 0; j < nrhs; j++ ) {
            memcpy( pb + k   + j * ldpb,
                    b  + ige + j * ldb,
                    dofi * sizeof(pastix_complex64_t) );
        }
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Applies a row permutation (peritab) to the matrix b. It also sends
 * and receives the part of b according to the spm repartition.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[inout] b
 *          The matrix b lda-by-nrhs.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 *******************************************************************************
 *
 * @retval pb which correspond to the vector b expanded, permuted and with
 * the correct local data.
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
    pastix_complex64_t  *bp         = Pb->b;
    bvec_handle_comm_t  *comm_rhs   = NULL;
    pastix_int_t        *cblk2col   = Pb->cblk2col;
    pastix_int_t         clustnum   = pastix_data->solvmatr->clustnum;
    pastix_int_t         bcsc_n     = Pb->m;
    pastix_int_t         replicated = 1;
    pastix_int_t         i, ige, ipe;
    pastix_int_t         j;
    bvec_data_amount_t  *data_send;
    pastix_int_t        *idxptr;
    pastix_complex64_t  *valptr;

    assert( Pb->m == pastix_data->bcsc->n );
    // assert( m     == spm->nexp ); Uncomment with dist2dist
    assert( b     != NULL );

    bvec_zcompute_and_exchange_amount_of_data( pastix_data, PastixDirBackward, replicated, Pb );
    comm_rhs  = Pb->rhs_comm;
    data_send = &(comm_rhs->data_comm[clustnum].send);
    cblk2col  = Pb->cblk2col;

    /*
     * Reset the counters used to fill indexes and values buffers.
     */
    idxptr = data_send->idxbuf;
    valptr = data_send->valbuf;

    /*
     * Goes through b to fill the data_comm with the data to send and
     * fills pb with the local data.
     */
    for ( i = 0; i < bcsc_n; i ++ ) {
        ipe = cblk2col[ i ];
        ige = bvec_zPglob2loc( pastix_data, replicated, ipe );

        /* Stores the indexes to send to c: (ipe, j). */
        *idxptr = ipe;
        idxptr++;

        /* Stores the value to send to c. */
        for ( j = 0; j < nrhs; j++ ) {
            memcpy( valptr, bp + i + j * Pb->ld, sizeof(pastix_complex64_t) );
            valptr++;
            memcpy( b + ige + j * ldb, bp + i + j * Pb->ld, sizeof(pastix_complex64_t) );
        }
    }

    assert( (idxptr - data_send->idxbuf) == data_send->idxcnt );
    assert( (valptr - ((pastix_complex64_t*)(data_send->valbuf))) == data_send->valcnt );

    bvec_zexchange_and_handle_data( pastix_data, PastixDirBackward, replicated, b, ldb, Pb );

    (void)m;
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

    return rc;
}
