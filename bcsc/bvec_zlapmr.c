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
 * @brief Gives the local permuted index corresponding to the global not
 * permuted given in argument or the process to which it belongs to.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] ig
 *          The global index not permuted.
 *
 *******************************************************************************
 *
 * @retval Returns the local permuted index or c = -(p+1) with p the process to
 * which the local permuted column/row belongs to.
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
 * @ingroup bcsc_internal
 *
 * @brief Gives the local index not permuted corresponding to the global
 * index permuted given in argument or the process to which it belongs to.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] replicated
 *          If replicated case then equals to 1 and equals to 0 otherwise.
 *
 * @param[in] igp
 *          The global index permuted.
 *
 *******************************************************************************
 *
 * @retval Returns the local not permuted index or c = -(p+1) with p the process to
 * which the local not permuted column/row belongs to.
 *
 *******************************************************************************/
static inline pastix_int_t
bvec_zPglob2loc( pastix_data_t *pastix_data,
                 int            replicated,
                 pastix_int_t   igp )
{
    const spmatrix_t *spm          = pastix_data->csc;
    pastix_order_t   *ord          = pastix_data->ordemesh;
    pastix_int_t     *perm         = NULL;
    pastix_int_t      baseval_ord  = ord->baseval;
    pastix_int_t      il;

    perm = orderGetExpandedPeritab( ord, spm );

    il = perm[ igp ] - baseval_ord;

    (void)replicated;
    return il;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Computes the amount of data the current process will send to the
 * other process.
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
 *          If replicated case then equals to 1 and equals to 0 otherwise.
 *
 * @param[inout] Pb
 *          The initialized rhs structure that holds the permuted right hand
 *          side matrix.
 *          On exit, the rhs_comm field is initialized.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
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
 * @ingroup bcsc_internal
 *
 * @brief Exchanges the amount of data.
 *
 *******************************************************************************
 *
 * @param[inout] rhs_comm
 *          The rhs_comm of the permuted vector initialized on entry, rhs_comm
 *          with the amount of data exchanged and array allocated at exit.
 *
 * @param[in] replicated
 *          If replicated case then equals to 1 and equals to 0 otherwise.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
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
 * @ingroup bcsc_internal
 *
 * @brief Computes and sends the amount of data the current process will
 * send to the other process and receives the amount of data each processor
 * will send to the current processor. Computes the cblk2col array.
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
 *          If replicated case then equals to 1 and equals to 0 otherwise.
 *
 * @param[inout] Pb
 *          The initialized rhs structure that holds the permuted right hand
 *          side matrix.
 *          On exit, the rhs_comm field is initialized.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_zcompute_and_exchange_amount_of_data( pastix_data_t *pastix_data,
                                           pastix_dir_t   dir,
                                           int            replicated,
                                           pastix_rhs_t   Pb )
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
 * @ingroup bcsc_internal
 *
 * @brief Exchanges the data.
 *
 *******************************************************************************
 *
 * @param[inout] rhs_comm
 *          The rhs_comm of the permuted vector initialized on entry, rhs_comm
 *          with the data exchanged at exit.
 *
 * @param[in] replicated
 *          If replicated case then equals to 1 and equals to 0 otherwise.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
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
 * @ingroup bcsc_internal
 *
 * @brief Copies the received data in the right hand side matrix b.
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
 *          If replicated case then equals to 1 and equals to 0 otherwise.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[inout] b
 *          The right hand side matrix b of size ldb-by-nrhs on entry.
 *          On exit, the matrix integrates the remote values at indexes positions.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 * @param[in] indexes
 *          The indexes array received.
 *
 * @param[in] values
 *          The values array received.
 *
 * @param[in] size
 *          The size of indexes.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
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
 * @ingroup bcsc_internal
 *
 * @brief Exchanges the data and copies the received data in the vector b.
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
 *          If replicated case then equals to 1 and equals to 0 otherwise.
 *
 * @param[in] nrhs
 *          The number of columns in the matrix b.
 *
 * @param[inout] b
 *          The input right hand side matrix b of size ldb-by-nrhs.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 * @param[inout] Pb
 *          The rhs structure that stores the permuted right hand side matrix b.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
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
    const spmatrix_t   *spm  = pastix_data->csc;
    pastix_int_t        dof  = spm->dof;
    const pastix_int_t *dofs = spm->dofs;
    pastix_int_t        ldpb = Pb->ld;
    pastix_int_t        j, ig, ige, ipe, dofi;

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
        ipe  = bvec_zglob2Ploc( pastix_data, ig );

        if ( ipe < 0 ) {
            continue;
        }

        for ( j = 0; j < nrhs; j++ ) {
            memcpy( pb + ipe + j * ldpb,
                    b  + ige + j * ldb,
                    dofi * sizeof(pastix_complex64_t) );
        }
    }

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
            free( PA->b );
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
