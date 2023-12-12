/**
 *
 * @file bvec.c
 *
 * @copyright 2004-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Vincent Bridonneau
 * @author Alycia Lisito
 * @date 2023-07-21
 *
 **/
#include "common.h"
#include <math.h>
#include "lapacke.h"
#include "bcsc/bcsc.h"
#include "bcsc/bvec.h"
#include "order/order_internal.h"
#include "cblas.h"
#include "blend/solver.h"

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Allocate a vector
 *
 *******************************************************************************
 *
 * @param[in] size
 *          The size of the vector
 *
 *******************************************************************************
 *
 * @return The allocated vector
 *
 *******************************************************************************/
void *bvec_malloc( size_t size )
{
    void *x = NULL;
    MALLOC_INTERN(x, size, char);
    return x;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Free a vector
 *
 *******************************************************************************
 *
 * @param[inout] x
 *          The vector to be free
 *
 *******************************************************************************/
void bvec_free( void *x )
{
    memFree_null(x);
}

#if defined( PASTIX_WITH_MPI )
/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initializes the rhs_comm field of an RHS data structure.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[inout] Pb
 *          On entry, the initialized pastix_rhs_t data structure but the rhs_comm
 *          field is NULL.
 *          On exit, rhs_comm is initialized.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit.
 *
 *******************************************************************************/
int
bvec_handle_comm_init( const pastix_data_t *pastix_data,
                       pastix_rhs_t         Pb )
{
    const SolverMatrix *solvmtx  = pastix_data->solvmatr;
    bvec_handle_comm_t *rhs_comm = NULL;

    /* Initializes bvec_comm. */
    if ( Pb->rhs_comm == NULL ) {
        size_t size = sizeof(bvec_handle_comm_t) + (solvmtx->clustnbr-1) * sizeof(bvec_proc_comm_t);

        Pb->rhs_comm = (bvec_handle_comm_t *)malloc( size );
    }

    rhs_comm = Pb->rhs_comm;
    rhs_comm->clustnbr = solvmtx->clustnbr;
    rhs_comm->clustnum = solvmtx->clustnum;
    rhs_comm->comm     = solvmtx->solv_comm;
    rhs_comm->flttype  = Pb->flttype;
    rhs_comm->max_idx  = 0;
    rhs_comm->max_val  = 0;

    memset( rhs_comm->data_comm, 0, rhs_comm->clustnbr * sizeof(bvec_proc_comm_t) );

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Cleans-up a bvec_handle_comm_t structure.
 *
 *******************************************************************************

 * @param[inout] rhs_comm
 *          On entry, the initialized bvec_handle_comm_t data structure.
 *          On exit, the structure is destroyed and should no longer be used.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit.
 *
 *******************************************************************************/
int
bvec_handle_comm_exit( bvec_handle_comm_t *rhs_comm )
{
    int c;
    int clustnbr = rhs_comm->clustnbr;
    bvec_proc_comm_t *data;

    for ( c = 0; c < clustnbr; c++ ) {
        data = rhs_comm->data_comm + c;

        if ( data->send_idxbuf != NULL ) {
            memFree_null( data->send_idxbuf );
        }
        if ( data->send_valbuf != NULL ) {
            memFree_null( data->send_valbuf );
        }
    }
    rhs_comm->max_idx = 0;
    rhs_comm->max_val = 0;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Converts the non permuted global index into the permuted local one.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] ig
 *          Index global not permuted.
 *
 *******************************************************************************
 *
 * @retval Returns the local permuted index or c = -(p+1) with p the process to
 * which the local permuted column/row belongs to.
 *
 *******************************************************************************/
pastix_int_t
bvec_glob2Ploc( const pastix_data_t *pastix_data,
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
    pastix_int_t          igp, igpe, ilpe, cblknum;

    assert( ord->baseval == 0 );

    igp     = ord->permtab[ ig ];
    igpe    = ( dof > 0 ) ? igp * dof : dofs[ ig ] - basespm;/* vdof incorrect */
    cblknum = col2cblk[ igpe ];
    if ( cblknum >= 0 ) {
        cblk = solvmatr->cblktab + cblknum;
        ilpe = cblk->lcolidx + igpe - cblk->fcolnum;
    }
    else {
        ilpe = cblknum;
    }

    return ilpe;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Converts the permuted global index into the non permuted local one.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] glob2loc
 *          The glob2loc array associated to the user spm stored in pastix_data.
 *
 * @param[in] igp
 *          Index global permuted.
 *
 *******************************************************************************
 *
 * @retval Returns the local non permuted index or c = -(p+1) with p the process
 * to which the local permuted column/row belongs to.
 *
 *******************************************************************************/
pastix_int_t
bvec_Pglob2loc( const pastix_data_t *pastix_data,
                const pastix_int_t  *glob2loc,
                pastix_int_t         igp )
{
    const spmatrix_t     *spm      = pastix_data->csc;
    const pastix_order_t *ord      = pastix_data->ordemesh;
    pastix_int_t          basespm  = spm->baseval;
    pastix_int_t          dof      = spm->dof;
    const pastix_int_t   *dofs     = spm->dofs;
    pastix_int_t          ig, il, ile;

    assert( ord->baseval == 0 );

    ig = ord->peritab[ igp ];
    il = glob2loc[ig];
    if ( il >= 0 ) {
        ile = ( dof > 0 ) ? il * dof : dofs[ ig ] - basespm;/* vdof incorrect */
    }
    else {
        ile = il;
    }

    return ile;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Computes the Ploc2Pglob array.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[inout] Pb
 *          On entry, the initialized pastix_rhs_t data structure but the
 *          Ploc2Pglob field is NULL.
 *          On exit, Ploc2Pglob is computed.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
bvec_compute_Ploc2Pglob( pastix_data_t *pastix_data,
                         pastix_rhs_t   Pb )
{
    const spmatrix_t *spm        = pastix_data->csc;
    pastix_int_t     *col2cblk   = pastix_data->bcsc->col2cblk;
    pastix_int_t     *Ploc2Pglob = NULL;
    pastix_int_t      dof        = spm->dof;
    pastix_int_t      bcsc_n     = ( dof > 0 ) ? Pb->m / dof : Pb->m ; // ToDo : mvdof
    pastix_int_t      kgpe, kgp, klp;

    if ( Pb->Ploc2Pglob != NULL ) {
        memFree_null( Pb->Ploc2Pglob );
    }

    /*
     * Creates Ploc2Pglob with col2cblk: gives the local index corresponding to the
     * global one (Ploc2Pglob[ig] = il).
     */
    if ( bcsc_n > 0 ) {
        MALLOC_INTERN( Ploc2Pglob, bcsc_n, pastix_int_t );
        klp = 0;
        kgp = 0;
        for ( kgpe = 0; kgpe < spm->gNexp; kgpe += dof, kgp ++ ) {
            if ( col2cblk[ kgpe ] >= 0 ) {
                Ploc2Pglob[ klp ] = kgp;
                klp ++;
            }
        }
        assert( klp == bcsc_n );
    }
    Pb->Ploc2Pglob = Ploc2Pglob;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Exchanges the amount of data in the replicated case.
 *
 *******************************************************************************
 *
 * @param[inout] rhs_comm
 *          On entry the rhs_comm of the permuted right hand side initialised.
 *          At exit rhs_comm is filled with the amount of data exchanged.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
bvec_exchange_amount_rep( bvec_handle_comm_t *rhs_comm )
{
    bvec_proc_comm_t *data_comm = rhs_comm->data_comm;
    pastix_int_t      clustnbr  = rhs_comm->clustnbr;
    pastix_int_t      clustnum  = rhs_comm->clustnum;
    pastix_int_t      max_idx   = 0;
    pastix_int_t      max_val   = 0;
    pastix_int_t      idx_cnt, val_cnt;
    pastix_int_t      c;

    /* Receives the amount of indexes and values. */
    for ( c = 0; c < clustnbr; c++ ) {
        data_comm = rhs_comm->data_comm + c;
        if ( c == clustnum ) {
            MPI_Bcast( &(data_comm->nsends), 2, PASTIX_MPI_INT, c, rhs_comm->comm );
            continue;
        }

        MPI_Bcast( &(data_comm->nrecvs), 2, PASTIX_MPI_INT, c, rhs_comm->comm );

        idx_cnt = data_comm->nrecvs.idxcnt;
        val_cnt = data_comm->nrecvs.valcnt;

        assert( idx_cnt <= val_cnt );
        if ( idx_cnt == 0 ) {
            assert( val_cnt == 0 );
        }

        if ( max_idx < idx_cnt ) {
            max_idx = idx_cnt;
        }
        if ( max_val < val_cnt ) {
            max_val = val_cnt;
        }
    }

    assert( max_idx <= max_val );

    rhs_comm->max_idx = max_idx;
    rhs_comm->max_val = max_val;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Computes the amount of sending data in the distributed case.
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure.
 *
 * @param[in] m
 *          The number of rows of the right hand side b.
 *
 * @param[in] nrhs
 *          The number of columns in the right hand side b.
 *
 * @param[inout] Pb
 *          On entry the rhs structure initialised.
 *          At exit the field rhs_comm is filled with the amount of sending
 *          data.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_compute_amount_dst( const pastix_data_t *pastix_data,
                         pastix_int_t         m,
                         pastix_int_t         nrhs,
                         pastix_rhs_t         Pb )
{
    const spmatrix_t   *spm         = pastix_data->csc;
    pastix_int_t        baseval_spm = spm->baseval;
    pastix_int_t        dof         = spm->dof;
    pastix_int_t       *dofs        = spm->dofs;
    pastix_int_t       *loc2glob    = spm->loc2glob;
    bvec_handle_comm_t *comm_rhs    = NULL;
    bvec_proc_comm_t   *data        = NULL;
    pastix_int_t        il, ile, ilpe, ig, dofi, c;

    bvec_handle_comm_init( pastix_data, Pb );
    comm_rhs = Pb->rhs_comm;

    /*
     * Goes through b to fill the data_comm with the data to send and
     * fills pb with the local data.
     */
    for ( il = 0, ile = 0; ile < m; ile += dofi, il++ ) {
        ig   = loc2glob[ il ] - baseval_spm;
        ilpe = bvec_glob2Ploc( pastix_data, ig );
        dofi = ( dof > 0 ) ? dof : dofs[ ig+1 ] - dofs[ ig ];

        if ( ilpe < 0 ) {
            c         = - ( ilpe + 1 );
            data = comm_rhs->data_comm + c;

            data->nsends.idxcnt += 1;
            data->nsends.valcnt += dofi * nrhs;
        }
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Computes the maximum size of the sending indexes and values buffers.
 *
 *******************************************************************************
 *
 * @param[inout] rhs_comm
 *         On entry the rhs_comm of the permuted right hand side initialized.
 *         At exit the fields max_idx and max_val of rhs_comm are updated.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_compute_max( bvec_handle_comm_t *rhs_comm )
{
    bvec_proc_comm_t *data     = NULL;
    pastix_int_t      clustnbr = rhs_comm->clustnbr;
    pastix_int_t      clustnum = rhs_comm->clustnum;
    pastix_int_t      max_idx  = 0;
    pastix_int_t      max_val  = 0;
    pastix_int_t      idx_cnt, val_cnt, c;

    /* Receives the amount of indexes and values. */
    for ( c = 0; c < clustnbr; c++ ) {
        data = rhs_comm->data_comm + c;
        if ( c == clustnum ) {
            continue;
        }

        idx_cnt = data->nrecvs.idxcnt;
        val_cnt = data->nrecvs.valcnt;

        if ( max_idx < idx_cnt ) {
            max_idx = idx_cnt;
        }
        if ( max_val < val_cnt ) {
            max_val = val_cnt;
        }
    }

    assert( max_idx <= max_val );

    rhs_comm->max_idx = max_idx;
    rhs_comm->max_val = max_val;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc
 *
 * @brief Switches the amount of the sending and receiving data in the
 * distributed case.
 *
 *******************************************************************************
 *
 * @param[inout] rhs_comm
 *         On entry the rhs_comm of the permuted right hand side initialized.
 *         At exit rhs_comm is updated with the right amount of sending and
 *         receiving  data.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
static inline int
bvec_switch_amount_dst( bvec_handle_comm_t *rhs_comm )
{
    bvec_proc_comm_t *data     = NULL;
    pastix_int_t      clustnbr = rhs_comm->clustnbr;
    pastix_int_t      clustnum = rhs_comm->clustnum;
    pastix_int_t      c, idx_tmp, val_tmp;

    /* Sends the same amout of data to all process. */
    for ( c = 0; c < clustnbr; c ++ ) {

        data = rhs_comm->data_comm + c;

        if ( c == clustnum ) {
            continue;
        }

        idx_tmp             = data->nsends.idxcnt;
        data->nsends.idxcnt = data->nrecvs.idxcnt;
        data->nrecvs.idxcnt = idx_tmp;

        val_tmp             = data->nsends.valcnt;
        data->nsends.valcnt = data->nrecvs.valcnt;
        data->nrecvs.valcnt = val_tmp;

    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Exchanges the amount of data in the distributed case.
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
 * @param[in] m
 *          If dir == PastixDirForward:
 *              m is the number of rows of the right hand side b.
 *          If dir == PastixDirBackward:
 *              m is not used.
 *
 * @param[in] nrhs
 *          The number of columns in the right hand side b.
 *
 * @param[inout] Pb
 *          On entry the rhs structure initialised.
 *          At exit the field rhs_comm is filled with the amount of data
 *          exchanged.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
bvec_exchange_amount_dst( pastix_data_t *pastix_data,
                          pastix_dir_t   dir,
                          pastix_int_t   m,
                          pastix_int_t   nrhs,
                          pastix_rhs_t   Pb )
{
    bvec_handle_comm_t *rhs_comm    = Pb->rhs_comm;
    pastix_int_t        clustnbr    = rhs_comm->clustnbr;
    pastix_int_t        clustnum    = rhs_comm->clustnum;
    bvec_proc_comm_t   *data_comm   = NULL;
    bvec_proc_comm_t   *data_send   = NULL;
    bvec_proc_comm_t   *data_recv   = NULL;
    pastix_int_t        counter_req = 0;
    MPI_Status          statuses[(clustnbr-1)*2];
    MPI_Request         requests[(clustnbr-1)*2];
    bvec_data_amount_t *sends, *recvs;
    pastix_int_t        c_send, c_recv, k;

    if ( dir == PastixDirBackward ) {
        bvec_switch_amount_dst( Pb->rhs_comm );
        bvec_compute_max( Pb->rhs_comm );
        return PASTIX_SUCCESS;
    }

    bvec_compute_amount_dst( pastix_data, m, nrhs, Pb );
    rhs_comm  = Pb->rhs_comm;
    data_comm = rhs_comm->data_comm;

    /* Receives the amount of indexes and values. */
    c_send = (clustnum+1) % clustnbr;
    c_recv = (clustnum-1+clustnbr) % clustnbr;
    for ( k = 0; k < clustnbr-1; k++ ) {
        data_send = data_comm + c_send;
        data_recv = data_comm + c_recv;
        sends     = &( data_send->nsends );
        recvs     = &( data_recv->nrecvs );
        if ( c_send == clustnum ) {
            continue;
        }

        MPI_Irecv( recvs, 2, PASTIX_MPI_INT, c_recv, PastixTagAmount,
                   rhs_comm->comm, &requests[counter_req++] );

        MPI_Isend( sends, 2, PASTIX_MPI_INT, c_send, PastixTagAmount,
                   rhs_comm->comm, &requests[counter_req++] );

        c_send = (c_send+1) % clustnbr;
        c_recv = (c_recv-1+clustnbr) % clustnbr;
    }

    MPI_Waitall( counter_req, requests, statuses );

    bvec_compute_max( rhs_comm );

    return PASTIX_SUCCESS;
}
#endif
