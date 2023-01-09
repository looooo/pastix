/**
 *
 * @file bvec_zmpi_comm.c
 *
 * Functions to communicate data between the processes when the right
 * hand side permutation is done on distributed architectures.
 *
 * @copyright 2004-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Alycia Lisito
 * @date 2022-12-01
 * @precisions normal z -> c d s
 *
 * This file implements the communication functions used in bvec_zlapmr
 * when MPI is used in replicated or distributed mode.
 *
 **/
#include "common.h"
#include <math.h>
#include "lapacke.h"
#include "bcsc/bcsc.h"
#include "bcsc/bvec.h"
#include "bcsc/bcsc_z.h"
#include "order/order_internal.h"
#include "cblas.h"
#include "blend/solver.h"

#if defined( PASTIX_WITH_MPI )
/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Copies the received data in the right hand side b in the replicated
 * case.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] nrhs
 *          The number of columns in the right hand side b.
 *
 * @param[inout] b
 *          On entry the right hand side b not fully filled.
 *          At exit b is updated with the remote values.
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
 * @param[in] size_idx
 *          The size of indexes.
 *
 * @param[in] size_val
 *          The size of values.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_zhandle_recv_backward_rep( pastix_data_t      *pastix_data,
                                pastix_int_t        nrhs,
                                pastix_complex64_t *b,
                                pastix_int_t        ldb,
                                pastix_int_t       *indexes,
                                pastix_complex64_t *values,
                                pastix_int_t        size_idx,
                                pastix_int_t        size_val )
{
    const spmatrix_t *spm   = pastix_data->csc;
    pastix_int_t      dof   = spm->dof;
    pastix_int_t     *dofs  = spm->dofs;
    pastix_int_t      ldval = size_val / nrhs ;
    pastix_int_t      ig, ige, idx, j, dofi, cnt;

    /* Checks if ldval is not rounded (there is no error on size_val). */
    assert( nrhs * ldval == size_val );

    cnt = 0;
    for ( idx = 0; idx < size_idx; idx++, indexes++, cnt += dofi ) {
        ig   = indexes[ 0 ];
        ige  = ( dof > 0 ) ? ig * dof : dofs[ig];
        dofi = ( dof > 0 ) ? dof : dofs[ig+1] - dofs[ig];

        for ( j = 0; j < nrhs; j++ ) {
            memcpy( b + ige + j * ldb, values + cnt + j * ldval, dofi * sizeof(pastix_complex64_t) );
        }
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Exchanges the data and copies the received data in the right hand
 * side b in the replicated case.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] nrhs
 *          The number of columns in the right hand side b.
 *
 * @param[inout] b
 *          The right hand side b ldb-by-nrhs.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 * @param[inout] rhs_comm
 *          The rhs_comm of the permuted right hand side initialised on entry,
 *          rhs_comm with the data exchanged at exit.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
bvec_zexchange_data_rep( pastix_data_t      *pastix_data,
                         pastix_int_t        nrhs,
                         pastix_complex64_t *b,
                         pastix_int_t        ldb,
                         pastix_rhs_t        Pb )
{
    bvec_handle_comm_t *rhs_comm  = Pb->rhs_comm;
    bvec_proc_comm_t   *data_comm = rhs_comm->data_comm;
    pastix_int_t        clustnbr  = rhs_comm->clustnbr;
    pastix_int_t        clustnum  = rhs_comm->clustnum;
    pastix_int_t       *idx_buf   = NULL;
    pastix_complex64_t *val_buf   = NULL;
    bvec_data_amount_t *sends;
    bvec_data_amount_t *recvs;
    pastix_int_t        c;

    /* Allocates the receiving indexes and values buffers. */
    if ( rhs_comm->max_idx > 0 ) {
        MALLOC_INTERN( idx_buf, rhs_comm->max_idx, pastix_int_t );
        MALLOC_INTERN( val_buf, rhs_comm->max_val, pastix_complex64_t );
    }

    for ( c = 0; c < clustnbr; c++ ) {
        data_comm = rhs_comm->data_comm + c;
        sends     = &( data_comm->send );
        recvs     = &( data_comm->recv );

        if ( c == clustnum ) {
            /* Posts the emissions of the indexes and values. */
            if ( sends->idxcnt > 0 ) {
                MPI_Bcast( sends->idxbuf, sends->idxcnt, PASTIX_MPI_INT,       c, rhs_comm->comm );
                MPI_Bcast( Pb->b,         sends->valcnt, PASTIX_MPI_COMPLEX64, c, rhs_comm->comm );
            }
            continue;
        }

        /* Posts the receptions of the indexes and values. */
        if ( ( rhs_comm->max_idx > 0 ) && ( recvs->idxcnt > 0 ) ) {
            MPI_Bcast( idx_buf, recvs->idxcnt, PASTIX_MPI_INT,       c, rhs_comm->comm );
            MPI_Bcast( val_buf, recvs->valcnt, PASTIX_MPI_COMPLEX64, c, rhs_comm->comm );

            assert( recvs->idxcnt <= recvs->valcnt );
            bvec_zhandle_recv_backward_rep( pastix_data, nrhs, b, ldb, idx_buf, val_buf,
                                            recvs->idxcnt, recvs->valcnt );
        }
    }

    /* Frees the receiving indexes and values buffers. */
    if ( rhs_comm->max_idx > 0 ) {
        memFree_null( idx_buf );
        memFree_null( val_buf );
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Allocates the sending buffers in rhs_comm->data_comm. These buffer
 * are filled with the sending values.
 *
 *******************************************************************************
 *
 * @param[inout] rhs_comm
 *         On entry the rhs_comm of the permuted right hand side initialized.
 *         At exit the arrays of rhs_comm->data_comm->send are allocated.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
bvec_zallocate_buf_dst( bvec_handle_comm_t *rhs_comm )
{
    bvec_proc_comm_t   *data      = NULL;
    bvec_data_amount_t *data_send = NULL;
    pastix_int_t        clustnbr  = rhs_comm->clustnbr;
    pastix_int_t        clustnum  = rhs_comm->clustnum;
    pastix_int_t        c;

    /* Sends the same amout of data to all process. */
    for ( c = 0; c < clustnbr; c ++ ) {

        data      = rhs_comm->data_comm + c;
        data_send = &(data->send);

        if ( c == clustnum ) {
            continue;
        }

        MALLOC_INTERN( data_send->idxbuf, data_send->idxcnt, pastix_int_t );
        MALLOC_INTERN( data_send->valbuf, data_send->valcnt, pastix_complex64_t );

    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Copies the received data in the right hand side pb in the distributed
 * case.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] nrhs
 *          The number of columns in the right hand side pb.
 *
 * @param[inout] pb
 *          On entry the right hand side pb not fully filled.
 *          At exit pb is updated with the remote values.
 *
 * @param[in] ldpb
 *          The leading dimension of pb.
 *
 * @param[in] indexes
 *          The indexes array received.
 *
 * @param[in] values
 *          The values array received.
 *
 * @param[in] size_idx
 *          The size of indexes.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_zhandle_recv_forward_dst( pastix_data_t      *pastix_data,
                               pastix_int_t        nrhs,
                               pastix_complex64_t *pb,
                               pastix_int_t        ldpb,
                               pastix_int_t       *indexes,
                               pastix_complex64_t *values,
                               pastix_int_t        size_idx )
{
    const spmatrix_t *spm  = pastix_data->csc;
    pastix_int_t      dof  = spm->dof;
    pastix_int_t     *dofs = spm->dofs;
    pastix_int_t      ig, ilpe, idx, j, dofi;

    for ( idx = 0; idx < size_idx; idx++, indexes++ ) {
        ig   = indexes[ 0 ];
        ilpe = bvec_glob2Ploc( pastix_data, ig );
        assert( ilpe >= 0 );
        dofi = ( dof > 0 ) ? dof : dofs[ig+1] - dofs[ig];

        for ( j = 0; j < nrhs; j++, values += dofi ) {
            memcpy( pb + ilpe + j * ldpb, values, dofi * sizeof(pastix_complex64_t) );
        }
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Copies the received data in the right hand side b in the distributed
 * case.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] nrhs
 *          The number of columns in the right hand side b.
 *
 * @param[inout] b
 *          On entry the right hand side b not fully filled.
 *          At exit b is updated with the remote values.
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
 * @param[in] size_idx
 *          The size of indexes.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_zhandle_recv_backward_dst( pastix_data_t      *pastix_data,
                                pastix_int_t        nrhs,
                                pastix_complex64_t *b,
                                pastix_int_t        ldb,
                                pastix_int_t       *indexes,
                                pastix_complex64_t *values,
                                pastix_int_t        size_idx )
{
    const spmatrix_t *spm = pastix_data->csc;
    pastix_int_t      dof = spm->dof;
    pastix_int_t      igp, ile, idx, j, dofi;

    dofi = dof; /* vdof incorrect */
    for ( idx = 0; idx < size_idx; idx++, indexes++ ) {
        igp = indexes[ 0 ];
        ile = bvec_Pglob2loc( pastix_data, igp);
        assert( ile >= 0 );

        for ( j = 0; j < nrhs; j++, values += dofi ) {
            memcpy( b + ile + j * ldb, values, dofi * sizeof(pastix_complex64_t) );
        }
    }
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Copies the received data in the right hand side b or pb in the
 * distributed case.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *          If PastixDirForward, b is permuted into Pb.
 *          If PastixDirBackward, Pb is permuted into b.
 *
 * @param[in] nrhs
 *          The number of columns in the right hand sides b and pb.
 *
 * @param[inout] b
 *          If dir == PastixDirForward:
 *              On entry the right hand side b not fully filled.
 *              At exit b is updated with the remote values.
 *          If dir == PastixDirBackward:
 *              b is not modified.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 * @param[inout] Pb
 *          If dir == PastixDirForward:
 *              On entry the structure of the permuted right hand side pb
 *              not fully filled.
 *              At exit pb is updated with the remote values.
 *          If dir == PastixDirBackward:
 *              Pb is not modified.
 *
 * @param[in] indexes
 *          The indexes array received.
 *
 * @param[in] values
 *          The values array received.
 *
 * @param[in] size_idx
 *          The size of indexes.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_zhandle_recv_dst( pastix_data_t      *pastix_data,
                       pastix_dir_t        dir,
                       pastix_int_t        nrhs,
                       pastix_complex64_t *b,
                       pastix_int_t        ldb,
                       pastix_rhs_t        Pb,
                       pastix_int_t       *indexes,
                       pastix_complex64_t *values,
                       pastix_int_t        size_idx )
{
    if ( dir == PastixDirForward ) {
        bvec_zhandle_recv_forward_dst( pastix_data, nrhs, Pb->b, Pb->ld, indexes, values, size_idx );
    }
    else {
        bvec_zhandle_recv_backward_dst( pastix_data, nrhs, b, ldb, indexes, values, size_idx );
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Exchanges the data and copies the received data in the right hand
 * side b or pb in the distributed case.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] dir
 *          The direction of the permutation.
 *          If PastixDirForward, b is permuted into Pb.
 *          If PastixDirBackward, Pb is permuted into b.
 *
 * @param[in] nrhs
 *          The number of columns in the right hand side b.
 *
 * @param[inout] b
 *          If dir == PastixDirForward:
 *              On entry the right hand side b not fully filled.
 *              At exit b is updated with the remote values.
 *          If dir == PastixDirBackward:
 *              b is not modified.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 * @param[inout] Pb
 *          If dir == PastixDirForward:
 *              On entry the structure of the permuted right hand side pb
 *              not fully filled.
 *              At exit pb is updated with the remote values.
 *          If dir == PastixDirBackward:
 *              Pb is not modified.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
bvec_zexchange_data_dst( pastix_data_t      *pastix_data,
                         pastix_dir_t        dir,
                         pastix_int_t        nrhs,
                         pastix_complex64_t *b,
                         pastix_int_t        ldb,
                         pastix_rhs_t        Pb )
{
    bvec_handle_comm_t *rhs_comm    = Pb->rhs_comm;
    bvec_proc_comm_t   *data_comm   = rhs_comm->data_comm;
    pastix_int_t        clustnbr    = rhs_comm->clustnbr;
    pastix_int_t        clustnum    = rhs_comm->clustnum;
    pastix_int_t       *idx_buf     = NULL;
    pastix_complex64_t *val_buf     = NULL;
    bvec_proc_comm_t   *data_send   = NULL;
    bvec_proc_comm_t   *data_recv   = NULL;
    pastix_int_t        counter_req = 0;
    MPI_Status          statuses[(clustnbr-1)*2];
    MPI_Request         requests[(clustnbr-1)*2];
    bvec_data_amount_t *sends, *recvs;
    pastix_int_t        c_send, c_recv, k;

    /* Allocates the receiving indexes and values buffers. */
    if ( rhs_comm->max_idx > 0 ) {
        MALLOC_INTERN( idx_buf, rhs_comm->max_idx, pastix_int_t );
        MALLOC_INTERN( val_buf, rhs_comm->max_val, pastix_complex64_t );
    }

    c_send = (clustnum+1) % clustnbr;
    for ( k = 0; k < clustnbr-1; k++ ) {
        data_send = data_comm + c_send;
        sends     = &( data_send->send );

        if ( c_send == clustnum ) {
            continue;
        }

        /* Posts the emissions of the indexes. */
        if ( sends->idxcnt > 0 ) {
            MPI_Isend( sends->idxbuf, sends->idxcnt, PASTIX_MPI_INT,       c_send,
                       PastixTagIndexes, rhs_comm->comm, &requests[counter_req++] );
            MPI_Isend( sends->valbuf, sends->valcnt, PASTIX_MPI_COMPLEX64, c_send,
                       PastixTagValues,  rhs_comm->comm, &requests[counter_req++] );
        }
        c_send = (c_send+1) % clustnbr;
    }

    c_recv = (clustnum-1+clustnbr) % clustnbr;
    for ( k = 0; k < clustnbr-1; k++ ) {
        data_recv = data_comm + c_recv;
        recvs     = &( data_recv->recv );
        /* Posts the receptions of the indexes and values. */
        if ( ( rhs_comm->max_idx > 0 ) && ( recvs->idxcnt > 0 ) ) {
            MPI_Recv( idx_buf, recvs->idxcnt, PASTIX_MPI_INT,       c_recv, PastixTagIndexes,
                      rhs_comm->comm, MPI_STATUS_IGNORE );
            MPI_Recv( val_buf, recvs->valcnt, PASTIX_MPI_COMPLEX64, c_recv, PastixTagValues,
                      rhs_comm->comm, MPI_STATUS_IGNORE );

            assert( recvs->idxcnt <= recvs->valcnt );
            bvec_zhandle_recv_dst( pastix_data, dir, nrhs, b, ldb, Pb, idx_buf, val_buf,
                                   recvs->idxcnt );
        }
        c_recv = (c_recv-1+clustnbr) % clustnbr;
    }

    MPI_Waitall( counter_req, requests, statuses );

    /* Frees the receiving indexes and values buffers. */
    if ( rhs_comm->max_idx > 0 ) {
        memFree_null( idx_buf );
        memFree_null( val_buf );
    }

    return PASTIX_SUCCESS;
}
#endif
