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
 *          The number of columns in the matrix b.
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
 * @brief Exchanges the data and copies the received data in the vector b.
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
 *          The matrix b ldb-by-nrhs.
 *
 * @param[in] ldb
 *          The leading dimension of b.
 *
 * @param[inout] rhs_comm
 *          The rhs_comm of the permuted vector initialised on entry, rhs_comm
 *          with the data exchanged at exit.
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
    bvec_data_amount_t *sends;
    bvec_data_amount_t *recvs;
    pastix_int_t       *idx_buf;
    pastix_complex64_t *val_buf;
    pastix_int_t        c;

    /* Allocates the receiving indexes and values buffers. */
    if ( rhs_comm->max_idx != 0 ) {
        MALLOC_INTERN( idx_buf, rhs_comm->max_idx, pastix_int_t );
        MALLOC_INTERN( val_buf, rhs_comm->max_val, pastix_complex64_t );
    }

    for ( c = 0; c < clustnbr; c++ ) {
        data_comm = rhs_comm->data_comm + c;
        sends      = &( data_comm->send );
        recvs      = &( data_comm->recv );

        if ( c == clustnum ) {
            /* Posts the emissions of the indexes and values. */
            if ( sends->idxcnt != 0 ) {
                MPI_Bcast( sends->idxbuf, sends->idxcnt, PASTIX_MPI_INT,       c, rhs_comm->comm );
                MPI_Bcast( Pb->b,         sends->valcnt, PASTIX_MPI_COMPLEX64, c, rhs_comm->comm );
            }
            continue;
        }

        /* Posts the receptions of the indexes and values. */
        if ( recvs->idxcnt != 0 ) {
            MPI_Bcast( idx_buf, recvs->idxcnt, PASTIX_MPI_INT,       c, rhs_comm->comm );
            MPI_Bcast( val_buf, recvs->valcnt, PASTIX_MPI_COMPLEX64, c, rhs_comm->comm );

            assert( recvs->idxcnt <= recvs->valcnt );
            bvec_zhandle_recv_backward_rep( pastix_data, nrhs, b, ldb, idx_buf, val_buf,
                                            recvs->idxcnt, recvs->valcnt );
        }
    }

    /* Frees the receiving indexes and values buffers. */
    if ( idx_buf != NULL ) {
        memFree_null( idx_buf );
        memFree_null( val_buf );
    }

    return PASTIX_SUCCESS;
}
#endif
