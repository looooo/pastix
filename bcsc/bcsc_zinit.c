/**
 *
 * @file bcsc_zinit.c
 *
 * @copyright 2004-2022 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Theophile Terraz
 * @author Tony Delarue
 * @author Vincent Bridonneau
 * @author Alycia Lisito
 * @date 2022-10-11
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "pastix/order.h"
#include <spm.h>
#include "blend/solver.h"
#include "bcsc/bcsc.h"
#include "bcsc_z.h"

static inline pastix_complex64_t
__fct_id( pastix_complex64_t val ) {
    return val;
}

static inline pastix_complex64_t
__fct_conj( pastix_complex64_t val ) {
#if defined(PRECISION_c) || defined(PRECISION_z)
    return conj( val );
#else
    /* This function should not be called in this case. */
    assert(0);
    return val;
#endif
}

#if defined(PASTIX_WITH_MPI)
/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Stores the data the current processor has to send to the other processors
 *        in the bcsc_comm structure.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] col2cblk
 *          The array of matching columns with cblk indexes.
 *
 * @param[inout] bcsc_comm
 *          The structure in which the data is stored.
 *
 *******************************************************************************/
void
bcsc_zstore_data( const spmatrix_t     *spm,
                  const pastix_order_t *ord,
                  const pastix_int_t   *col2cblk,
                  bcsc_handle_comm_t   *bcsc_comm )
{
    pastix_complex64_t *values_c;
    pastix_complex64_t *values    = (pastix_complex64_t*)(spm->values);
    pastix_int_t       *colptr    = spm->colptr;
    pastix_int_t       *rowptr    = spm->rowptr;
    pastix_int_t       *loc2glob  = spm->loc2glob;
    pastix_int_t       *dofs      = spm->dofs;
    pastix_int_t        dof       = spm->dof;
    bcsc_proc_comm_t   *data_comm = bcsc_comm->data_comm;
    pastix_int_t        clustnbr  = bcsc_comm->clustnbr;
    pastix_int_t        k, j, ig, jg, ip, jp, ipe, jpe;
    pastix_int_t        itercblk_A, itercblk_At, baseval;
    pastix_int_t        dofi, dofj, c;
    bcsc_send_proc_t   *data_send;
    bcsc_data_amount_t *data_cntA, *data_cntAt;

    /* Allocates the indexes and values buffers. */
    bcsc_allocate_buf( bcsc_comm, PastixTagMemSend );

    /*
     * Allocates and initialises the counters used to fill bcsc_comm->values
     * and bcsc_comm->indexes.
     */
    MALLOC_INTERN( data_cntA, 2 * clustnbr, bcsc_data_amount_t );
    memset( data_cntA, 0, 2 * clustnbr * sizeof(bcsc_data_amount_t) );
    data_cntAt = data_cntA + clustnbr;

    baseval = spm->baseval;
    /*
     * Initializes the values of the matrix A in the blocked csc format. This
     * applies the permutation to the values array.
     *
     * Goes through the columns of the spm.
     * j is the initial local column index.
     * jg is the initial global column index.
     * jp is the "new" jg index -> the one resulting from the permutation.
     * jpe is the final index in the expanded matrix after permutation.
     */
    for ( j = 0; j < spm->n; j++, colptr++, loc2glob++ ) {
        jg   = *loc2glob - baseval;
        jp   = ord->permtab[jg];
        jpe  = ( dof > 0 ) ? jp * dof : dofs[ jg ] - baseval;
        dofj = ( dof > 0 ) ? dof : dofs[ jg+1 ] - dofs[ jg ];

        itercblk_A = col2cblk[ jpe ];

        /*
         * Goes through the row of the column j.
         * ig is the initial global row index (for the row index local = global).
         * ip is the "new" ig index -> the one resulting from the permutation.
         * ipe is the final index in the expanded matrix after permutation.
         */
        for ( k = colptr[0]; k < colptr[1]; k++, rowptr++ ) {
            ig   = *rowptr - baseval;
            ip   = ord->permtab[ig];
            ipe  = ( dof > 0 ) ? ip * dof : dofs[ ig ] - baseval;
            dofi = ( dof > 0 ) ? dof : dofs[ ig+1 ] - dofs[ ig ];

            /* The block of the column jp does not belong to the current processor. */
            if ( itercblk_A < 0 ) {
                /* Gets the owner of the block. */
                c         = - ( itercblk_A + 1 );
                data_comm = bcsc_comm->data_comm + c;
                data_send = &( data_comm->sendA );

                /* Stores the indexes to send to c: (ip, jp). */
                data_send->idxbuf[data_cntA[c].idxcnt  ] = ip;
                data_send->idxbuf[data_cntA[c].idxcnt+1] = jp;
                data_cntA[c].idxcnt += 2;

                /*
                 * Stores the values to send to c:
                 * dofi values in the row ip.
                 * dofj values in the column jp.
                 */
                values_c = ((pastix_complex64_t*)(data_send->valbuf)) + data_cntA[c].valcnt;

                memcpy( values_c, values, dofi * dofj * sizeof(pastix_complex64_t) );
                data_cntA[c].valcnt += dofi * dofj;
            }

            itercblk_At = col2cblk[ ipe ];
            /* The block of the row ip does not belong to the current processor. */
            if ( (itercblk_At < 0) && (ip != jp) ) {
                /* Gets the owner of the block. */
                c         = - ( itercblk_At + 1 );
                data_comm = bcsc_comm->data_comm + c;
                data_send = &( data_comm->sendAt );

                /* Stores the indexes to send to c: (jp, ip). */
                data_send->idxbuf[data_cntAt[c].idxcnt  ] = ip;
                data_send->idxbuf[data_cntAt[c].idxcnt+1] = jp;
                data_cntAt[c].idxcnt += 2;

                /*
                 * Stores the values to send to c:
                 * dofi values in the row ip.
                 * dofj values in the column jp.
                 */
                values_c = ((pastix_complex64_t*)(data_send->valbuf)) + data_cntAt[c].valcnt;

                memcpy( values_c, values, dofi * dofj * sizeof(pastix_complex64_t) );
                data_cntAt[c].valcnt += dofi * dofj;
            }
            values += dofi * dofj;
        }
    }
    memFree_null( data_cntA );
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Adds the remote data of A to the bcsc structure.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure which describes the data distribution.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc fields are updated.
 *
 * @param[in] values_buf
 *          The array containing the remote values.
 *
 * @param[in] idx_cnt
 *          The index for the begining of the indexes array buffer
 *          corresponding to the values_buf.
 *
 * @param[in] idx_size
 *          The number of indexes corresponding to the values_buf.
 *
 *******************************************************************************/
static inline pastix_int_t
bcsc_zhandle_recv_A( const spmatrix_t     *spm,
                     const pastix_order_t *ord,
                     const SolverMatrix   *solvmtx,
                     pastix_bcsc_t        *bcsc,
                     pastix_complex64_t   *values_buf,
                     pastix_int_t          idx_cnt,
                     pastix_int_t          idx_size )
{
    bcsc_handle_comm_t *bcsc_comm = bcsc->bcsc_comm;
    pastix_int_t       *dofs = spm->dofs;
    pastix_int_t        dof  = spm->dof;
    SolverCblk         *cblk;
    pastix_int_t       *coltab;
    pastix_int_t        k, ip, jp, ig, jg, ipe, jpe;
    pastix_int_t        itercblk;
    pastix_int_t        idofi, idofj, dofi, dofj;
    pastix_int_t        colidx, rowidx, pos;
    pastix_int_t        clustnum = bcsc_comm->clustnum;
    pastix_complex64_t *Values   = bcsc->Lvalues;
    pastix_int_t       *indexes  = bcsc_comm->data_comm[clustnum].sendA.idxbuf;
    pastix_int_t        nbelt    = 0;

    /*
     * Goes through the received indexes in order to add the data
     * received.
     */
    indexes += idx_cnt;
    for ( k = 0; k < idx_size; k+=2, indexes+=2 ) {
        /* Gets received indexes jp and ip. */
        ip = indexes[0];
        jp = indexes[1];

        /* Gets dofi and dofj with ig and jg. */
        ig   = ord->peritab[ip];
        jg   = ord->peritab[jp];
        ipe  = (dof > 0) ? ip * dof : dofs[ ig ] - spm->baseval;
        jpe  = (dof > 0) ? jp * dof : dofs[ jg ] - spm->baseval;
        dofi = (dof > 0) ? dof : dofs[ig+1] - dofs[ig];
        dofj = (dof > 0) ? dof : dofs[jg+1] - dofs[jg];

        itercblk = bcsc->col2cblk[ jpe ];
        assert( itercblk >= 0 );

        /* Gets the block on which the data received will be added. */
        cblk   = solvmtx->cblktab + itercblk;
        coltab = bcsc->cscftab[cblk->bcscnum].coltab;

        /* Goes through the values at (ip, jp). */
        colidx = jpe - cblk->fcolnum;
        for ( idofj = 0; idofj < dofj; idofj++, colidx++ ) {
            rowidx = ipe;
            pos    = coltab[ colidx ];
            for ( idofi = 0; idofi < dofi; idofi++, rowidx++, pos++, values_buf++ ) {
                /* Adds the row ip. */
                assert( rowidx >= 0 );
                assert( rowidx < spm->gNexp );
                bcsc->rowtab[ pos ] = rowidx;
                /* Adds the data at row ip and column jp. */
                Values[ pos ]       = *values_buf;
                nbelt++;
            }
            coltab[ colidx ] += dofi;
            assert( coltab[colidx] <= coltab[colidx+1] );
        }
    }

    return nbelt;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Adds the remote data of At to the bcsc structure.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure which describes the data distribution.
 *
 * @param[inout] rowtab
 *          The row tab of the bcsc or the row tab associated to At.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc fields are updated.
 *
 * @param[in] values_buf
 *          The array containing the remote values.
 *
 * @param[in] idx_cnt
 *          The index for the begining of the indexes array buffer
 *          corresponding to the values_buf.
 *
 * @param[in] idx_size
 *          The number of indexes corresponding to the values_buf.
 *
 *******************************************************************************/
static inline pastix_int_t
bcsc_zhandle_recv_At( const spmatrix_t     *spm,
                      const pastix_order_t *ord,
                      const SolverMatrix   *solvmtx,
                      pastix_int_t         *rowtab,
                      pastix_bcsc_t        *bcsc,
                      pastix_complex64_t   *values_buf,
                      pastix_int_t          idx_cnt,
                      pastix_int_t          idx_size )
{
    bcsc_handle_comm_t *bcsc_comm = bcsc->bcsc_comm;
    pastix_int_t       *dofs = spm->dofs;
    pastix_int_t        dof  = spm->dof;
    SolverCblk         *cblk;
    pastix_int_t       *coltab;
    pastix_int_t        k, ip, jp, ig, jg, ipe, jpe;
    pastix_int_t        itercblk;
    pastix_int_t        idofi, idofj, dofi, dofj;
    pastix_int_t        colidx, rowidx, pos;
    pastix_int_t        clustnum = bcsc_comm->clustnum;
    pastix_int_t       *indexes  = bcsc_comm->data_comm[clustnum].sendAt.idxbuf;
    pastix_complex64_t *Values;
    pastix_complex64_t (*_bcsc_conj)(pastix_complex64_t) = __fct_id;
    pastix_int_t        nbelt     = 0;

    /* Gets the right values in which the data will be added. */
    if ( spm->mtxtype == SpmGeneral ) {
        _bcsc_conj = __fct_id;
        Values = (pastix_complex64_t*)(bcsc->Uvalues);
    }
    else {
        if( spm->mtxtype == SpmHermitian ) {
            _bcsc_conj = __fct_conj;
        }
        if( spm->mtxtype == SpmSymmetric ) {
            _bcsc_conj = __fct_id;
        }
        Values = (pastix_complex64_t*)(bcsc->Lvalues);
    }

    /*
     * Goes through the received indexes in order to add the data
     * received.
     */
    indexes += idx_cnt;
    for ( k = 0; k < idx_size; k+=2, indexes+=2 ) {
        /* Gets received indexes ip and jp. */
        ip = indexes[0];
        jp = indexes[1];

        /* Gets dofi and dofj with ig and jg. */
        ig = ord->peritab[ip];
        jg = ord->peritab[jp];
        ipe  = (dof > 0) ? ip * dof : dofs[ ig ] - spm->baseval;
        jpe  = (dof > 0) ? jp * dof : dofs[ jg ] - spm->baseval;
        dofi = (dof > 0) ? dof : dofs[ig+1] - dofs[ig];
        dofj = (dof > 0) ? dof : dofs[jg+1] - dofs[jg];

        itercblk = bcsc->col2cblk[ ipe ];
        assert( itercblk >= 0 );

        /* Gets the block on which the data received will be added. */
        cblk   = solvmtx->cblktab + itercblk;
        coltab = bcsc->cscftab[cblk->bcscnum].coltab;

        /* Goes through the values at (ip, jp). */
        for ( idofj = 0; idofj < dofj; idofj++ ) {
            /* Gets the column ip. */
            rowidx = jpe + idofj;
            colidx = ipe - cblk->fcolnum;
            for ( idofi = 0; idofi < dofi; idofi++, colidx++, values_buf++ ) {
                pos = coltab[ colidx ];

                /* Adds the row jp. */
                assert( rowidx >= 0 );
                assert( rowidx < spm->gNexp );
                rowtab[ pos ] = rowidx;

                /* Adds the data at row jp and column ip. */
                Values[ pos ] = _bcsc_conj( *values_buf );
                coltab[ colidx ] ++;
                assert( coltab[colidx] <= coltab[colidx+1] );
                nbelt++;
            }
        }
    }
    return nbelt;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Exchanges and stores the remote values of A.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering which needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure which describes the data distribution.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc fields are updated with the remote values of A.
 *
 *******************************************************************************/
static inline pastix_int_t
bcsc_zexchange_values_A( const spmatrix_t     *spm,
                         const pastix_order_t *ord,
                         const SolverMatrix   *solvmtx,
                         pastix_bcsc_t        *bcsc )
{
    bcsc_handle_comm_t *bcsc_comm   = bcsc->bcsc_comm;
    pastix_int_t        clustnbr    = bcsc_comm->clustnbr;
    pastix_int_t        clustnum    = bcsc_comm->clustnum;
    bcsc_proc_comm_t   *data_comm   = bcsc_comm->data_comm;
    bcsc_proc_comm_t   *data_send   = NULL;
    bcsc_proc_comm_t   *data_recv   = NULL;
    pastix_complex64_t *val_buf     = NULL;
    pastix_int_t        idx_size    = 0;
    pastix_int_t        counter_req = 0;
    pastix_int_t        nbelt_recv  = 0;
    pastix_int_t        cnt         = 0;
    pastix_int_t        idx_cnt[clustnbr];
    MPI_Status          statuses[clustnbr-1];
    MPI_Request         requests[clustnbr-1];
    bcsc_data_amount_t *sends, *recvs;
    pastix_int_t        c_send, c_recv, k;

    /* Allocates the receiving indexes and values buffers. */
    if ( bcsc_comm->max_idx != 0 ) {
        MALLOC_INTERN( val_buf, bcsc_comm->max_val, pastix_complex64_t );
    }

    for ( k = 0; k < clustnbr; k++ ) {
        if ( k == clustnum ) {
            idx_cnt[k] = 0;
            continue;
        }
        idx_cnt[ k ] = cnt;
        cnt += data_comm[k].recvA.idxcnt;
    }

    /* Sends the values. */
    c_send = (clustnum+1) % clustnbr;
    for ( k = 0; k < clustnbr-1; k++ ) {
        data_send = data_comm + c_send;
        sends     = &( data_send->sendA.size );
        if ( c_send == clustnum ) {
            continue;
        }

        /* Posts the emissions of the values. */
        if ( sends->valcnt != 0 ) {
            MPI_Isend( data_send->sendA.valbuf, sends->valcnt, PASTIX_MPI_COMPLEX64,
                       c_send, PastixTagValuesA, bcsc_comm->comm, &requests[counter_req++] );
        }
	c_send = (c_send+1) % clustnbr;
    }

    /* Receives the values. */
    c_recv = (clustnum-1+clustnbr) % clustnbr;
    for ( k = 0; k < clustnbr-1; k++ ) {
        data_recv = data_comm + c_recv;
        recvs     = &( data_recv->recvA );
        if ( c_recv == clustnum ) {
            continue;
        }

        /* Posts the receptions of the values. */
        if ( recvs->valcnt != 0 ) {
            MPI_Recv( val_buf, recvs->valcnt, PASTIX_MPI_COMPLEX64,
                      c_recv, PastixTagValuesA, bcsc_comm->comm, MPI_STATUS_IGNORE );
            idx_size = recvs->idxcnt;
            nbelt_recv += bcsc_zhandle_recv_A( spm, ord, solvmtx, bcsc,
                                               val_buf, idx_cnt[c_recv], idx_size );
        }
        c_recv = (c_recv-1+clustnbr) % clustnbr;
    }

    MPI_Waitall( counter_req, requests, statuses );
    free( val_buf );
    assert( nbelt_recv == bcsc_comm->data_comm[clustnum].recvA.valcnt );
    bcsc_free_buf( bcsc_comm, PastixTagMemSendValA );
    bcsc_free_buf( bcsc_comm, PastixTagMemRecvIdxA );
    return nbelt_recv;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Exchanges and stores the remote values of At.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering which needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure which describes the data distribution.
 *
 * @param[inout] rowtab
 *          The row tab of the bcsc or the row tab associated to At.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc fields are updated with the remote values of At.
 *
 *******************************************************************************/
static inline pastix_int_t
bcsc_zexchange_values_At( const spmatrix_t     *spm,
                          const pastix_order_t *ord,
                          const SolverMatrix   *solvmtx,
                          pastix_int_t         *rowtab,
                          pastix_bcsc_t        *bcsc )
{
    bcsc_handle_comm_t *bcsc_comm   = bcsc->bcsc_comm;
    pastix_int_t        clustnbr    = bcsc_comm->clustnbr;
    pastix_int_t        clustnum    = bcsc_comm->clustnum;
    bcsc_proc_comm_t   *data_comm   = bcsc_comm->data_comm;
    bcsc_proc_comm_t   *data_send   = NULL;
    bcsc_proc_comm_t   *data_recv   = NULL;
    pastix_complex64_t *val_buf     = NULL;
    pastix_int_t        idx_size    = 0;
    pastix_int_t        counter_req = 0;
    pastix_int_t        nbelt_recv  = 0;
    pastix_int_t        cnt         = 0;
    pastix_int_t        idx_cnt[clustnbr];
    MPI_Status          statuses[clustnbr-1];
    MPI_Request         requests[clustnbr-1];
    bcsc_data_amount_t *sends, *recvs;
    pastix_int_t        c_send, c_recv, k;

    /* Allocates the receiving indexes and values buffers. */
    if ( bcsc_comm->max_idx != 0 ) {
        MALLOC_INTERN( val_buf, bcsc_comm->max_val, pastix_complex64_t );
    }

    for ( k = 0; k < clustnbr; k++ ) {
        if ( k == clustnum ) {
            idx_cnt[k] = 0;
            continue;
        }
        idx_cnt[ k ] = cnt;
        cnt += data_comm[k].recvAt.idxcnt;
    }

    /* Sends the values. */
    c_send = (clustnum+1) % clustnbr;
    for ( k = 0; k < clustnbr-1; k++ ) {
        data_send = data_comm + c_send;
        sends     = &( data_send->sendAt.size );
        if ( c_send == clustnum ) {
            continue;
        }

        /* Posts the emissions of the values. */
        if (  sends->valcnt != 0 ) {
            MPI_Isend( data_send->sendAt.valbuf, sends->valcnt, PASTIX_MPI_COMPLEX64,
                       c_send, PastixTagValuesAt, bcsc_comm->comm, &requests[counter_req++] );
        }
        c_send = (c_send+1) % clustnbr;
    }

    /* Receives the values. */
    c_recv = (clustnum-1+clustnbr) % clustnbr;
    for ( k = 0; k < clustnbr-1; k++ ) {
        data_recv = data_comm + c_recv;
        recvs     = &( data_recv->recvAt );
        if ( c_recv == clustnum ) {
            continue;
        }

        /* Posts the receptions of the values. */
        if ( recvs->valcnt != 0 ) {
            MPI_Recv( val_buf, recvs->valcnt, PASTIX_MPI_COMPLEX64,
                      c_recv, PastixTagValuesAt, bcsc_comm->comm, MPI_STATUS_IGNORE );
            idx_size = recvs->idxcnt;
            nbelt_recv += bcsc_zhandle_recv_At( spm, ord, solvmtx, rowtab, bcsc,
                                                val_buf, idx_cnt[c_recv], idx_size );
        }
        c_recv = (c_recv-1+clustnbr) % clustnbr;
    }

    MPI_Waitall( counter_req, requests, statuses );
    free( val_buf );
    assert( nbelt_recv == bcsc_comm->data_comm[clustnum].recvAt.valcnt );
    bcsc_free_buf( bcsc_comm, PastixTagMemSendValAt );
    bcsc_free_buf( bcsc_comm, PastixTagMemRecvIdxAt );
    return nbelt_recv;
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initializes the A values of the block csc stored in the given spm.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering which needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure which describes the data distribution.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc fields are updated.
 *
 *******************************************************************************/
static inline pastix_int_t
bcsc_zinit_A( const spmatrix_t     *spm,
              const pastix_order_t *ord,
              const SolverMatrix   *solvmtx,
              pastix_bcsc_t        *bcsc )
{
    pastix_complex64_t *values   = (pastix_complex64_t*)(spm->values);
    pastix_complex64_t *Lvalues  = (pastix_complex64_t*)(bcsc->Lvalues);
    pastix_int_t       *loc2glob = spm->loc2glob;
    pastix_int_t       *colptr   = spm->colptr;
    pastix_int_t       *rowptr   = spm->rowptr;
    pastix_int_t       *dofs     = spm->dofs;
    pastix_int_t        dof      = spm->dof;
    pastix_int_t        nbelt    = 0;
    SolverCblk         *cblk;
    pastix_int_t       *coltab;
    pastix_int_t        k, j, ig, jg, ip, jp, ipe, jpe;
    pastix_int_t        itercblk, baseval;
    pastix_int_t        idofi, idofj, dofi, dofj;
    pastix_int_t        colidx, rowidx, pos;

    baseval = spm->baseval;
    /*
     * Initializes the values of the matrix A in the blocked csc format. This
     * applies the permutation to the values array.
     *
     * Goes through the columns of the spm.
     * j is the initial local column index.
     * jg is the initial global column index.
     * jp is the "new" jg index -> the one resulting from the permutation.
     */
    for ( j = 0; j < spm->n; j++, colptr++, loc2glob++ ) {
        jg   = ( spm->loc2glob == NULL ) ? j : *loc2glob - baseval;
        jp   = ord->permtab[jg];
        jpe  = ( dof > 0 ) ? jp * dof : dofs[ jg ] - baseval;
        dofj = ( dof > 0 ) ? dof : dofs[ jg+1 ] - dofs[ jg ];

        itercblk = bcsc->col2cblk[ jpe ];

        /*
         * If MPI is used in shared memory, the block can belong to another
         * processor, in this case the column is skiped.
         */
        /* The block of the column jp belongs to the current processor. */
        if ( itercblk >= 0 ) {
            /* Gets the block in which the data will be added. */
            cblk   = solvmtx->cblktab + itercblk;
            coltab = bcsc->cscftab[cblk->bcscnum].coltab;
            /*
             * Goes through the row of the column j.
             * ig is the initial global row index (for the row index local = global).
             * ip is the "new" ig index -> the one resulting from the permutation.
             */
            for ( k = colptr[0]; k < colptr[1]; k++, rowptr++ ) {
                ig   = *rowptr - baseval;
                ip   = ord->permtab[ig];
                ipe  = ( dof > 0 ) ? ip * dof : dofs[ ig ] - baseval;
                dofi = ( dof > 0 ) ? dof : dofs[ ig+1 ] - dofs[ ig ];

                /* Copies the values from element (ig, jg) to position (ipe, jpe) into expanded bcsc. */
                colidx = jpe - cblk->fcolnum;
                for ( idofj = 0; idofj < dofj; idofj++, colidx++ ) {
                    rowidx = ipe;
                    pos    = coltab[ colidx ];
                    for ( idofi = 0; idofi < dofi; idofi++, rowidx++, pos++, values++ ) {
                        assert( rowidx >= 0 );
                        assert( rowidx < spm->gNexp );
                        bcsc->rowtab[ pos ] = rowidx;
                        Lvalues[ pos ]      = *values;
                        nbelt++;
                    }
                    coltab[ colidx ] += dofi;
                    assert( coltab[ colidx ] <= coltab[ colidx+1 ] );
                }
            }
        }
        /* The block of the column jp belongs to another processor. */
        else {
            for ( k = colptr[0]; k < colptr[1]; k++, rowptr++ ) {
                ig   = *rowptr - baseval;
                dofi = ( dof > 0 ) ? dof : dofs[ ig+1 ] - dofs[ ig ];

                values += dofi * dofj;
            }
        }
    }

    return nbelt;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initializes the At values in the block cscstored in the given spm.
 *
 * This routine initializes either :
 *      The symmetric upper part (L^t)
 *      The hermitian upper part (L^h)
 *      The transpose part of A (A^t -> U)
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering which needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure which describes the data distribution.
 *
 * @param[inout] rowtab
 *          The row tab of the bcsc or the row tab associated to At.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc fields are updated.
 *
 *******************************************************************************/
static inline pastix_int_t
bcsc_zinit_At( const spmatrix_t     *spm,
               const pastix_order_t *ord,
               const SolverMatrix   *solvmtx,
               pastix_int_t         *rowtab,
               pastix_bcsc_t        *bcsc )
{
    pastix_complex64_t *values  = (pastix_complex64_t*)(spm->values);
    pastix_complex64_t *Uvalues;
    pastix_int_t       *loc2glob = spm->loc2glob;
    pastix_int_t       *colptr   = spm->colptr;
    pastix_int_t       *rowptr   = spm->rowptr;
    pastix_int_t       *dofs     = spm->dofs;
    pastix_int_t        dof      = spm->dof;
    pastix_int_t        nbelt    = 0;
    SolverCblk         *cblk;
    pastix_int_t       *coltab;
    pastix_int_t        k, j, ig, jg, ip, jp, ipe, jpe;
    pastix_int_t        itercblk, baseval;
    pastix_int_t        idofi, idofj, dofi, dofj;
    pastix_int_t        colidx, rowidx, pos;
    pastix_complex64_t (*_bcsc_conj)(pastix_complex64_t) = NULL;

    /* We're working on U. */
    if ( spm->mtxtype == SpmGeneral ) {
        _bcsc_conj = __fct_id;
        Uvalues = (pastix_complex64_t*)(bcsc->Uvalues);
    }
    /* L^[t|h] part of the matrix. */
    else {
        /*
         * precision_generator/sub.py will change SpmHermitian to SpmSymmetric
         * Don't use else or else if.
         */
        if( spm->mtxtype == SpmHermitian ) {
            _bcsc_conj = __fct_conj;
        }
        if( spm->mtxtype == SpmSymmetric ) {
            _bcsc_conj = __fct_id;
        }
        Uvalues = (pastix_complex64_t*)(bcsc->Lvalues);
    }
    assert( _bcsc_conj != NULL );

    baseval = spm->baseval;
    /*
     * Initializes the values of the matrix A^t in the blocked csc format. This
     * applies the permutation to the values array.
     *
     * Goes through the columns of the spm.
     * j is the initial local column index.
     * jg is the initial global column index.
     * jp is the "new" jg index -> the one resulting from the permutation.
     */
    for ( j = 0; j < spm->n; j++, colptr++, loc2glob++ ) {
        jg   = ( spm->loc2glob == NULL ) ? j : *loc2glob - baseval;
        jp   = ord->permtab[jg];
        jpe  = ( dof > 0 ) ? jp * dof : dofs[ jg ] - baseval;
        dofj = ( dof > 0 ) ? dof : dofs[ jg+1 ] - dofs[ jg ];

        /*
         * Goes through the row of the column j.
         * ig is the initial global row index (for the row index local = global).
         * ip is the "new" ig index -> the one resulting from the permutation.
         */
        for ( k = colptr[0]; k < colptr[1]; k++, rowptr++ ) {
            ig   = *rowptr - baseval;
            ip   = ord->permtab[ig];
            ipe  = ( dof > 0 ) ? ip * dof : dofs[ ig ] - baseval;
            dofi = ( dof > 0 ) ? dof : dofs[ ig+1 ] - dofs[ ig ];

            /* Diagonal block of a symmetric matrix. */
            if ( ( ig == jg ) && ( spm->mtxtype != SpmGeneral ) ) {
                values += dofi * dofj;
                continue;
            }
            itercblk = bcsc->col2cblk[ ipe ];

            /*
             * If MPI is used in shared memory, the block can belong to another
             * processor, in this case the row is skiped.
             */
            /* The block of the row ip belongs to the current processor. */
            if ( itercblk >= 0 ) {
                /* Gets the block in which the data will be added. */
                cblk   = solvmtx->cblktab + itercblk;
                coltab = bcsc->cscftab[cblk->bcscnum].coltab;

                /* Copies the values from element (ig, jg) to position (ipe, jpe) into expanded bcsc. */
                for ( idofj = 0; idofj < dofj; idofj++ ) {
                    rowidx = jpe + idofj;
                    colidx = ipe - cblk->fcolnum;
                    for ( idofi = 0; idofi < dofi; idofi++, colidx++, values++ ) {
                        pos = coltab[ colidx ];

                        /* Copies the index/value */
                        assert( rowidx >= 0 );
                        assert( rowidx < spm->gNexp );
                        rowtab[ pos ]  = rowidx;
                        Uvalues[ pos ] = _bcsc_conj( *values );

                        coltab[ colidx ]++;
                        nbelt++;
                    }
                }
            }
            /* The block of the row ip belongs to another processor. */
            else {
                values += dofi * dofj;
                continue;
            }
        }
    }
    return nbelt;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Sorts the block csc subarray associated to each column block.
 *
 *******************************************************************************
 *
 * @param[in] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *
 * @param[in] rowtab
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] valtab
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 *******************************************************************************/
static inline void
bcsc_zsort( const pastix_bcsc_t *bcsc,
            pastix_int_t        *rowtab,
            pastix_complex64_t  *valtab )
{
    bcsc_cblk_t *block;
    pastix_int_t block_num, col, len_col_block, block_nbr, col_nbr;
    void        *sortptr[2];

    block_nbr = bcsc->cscfnbr;
    block     = bcsc->cscftab;
    for ( block_num = 0; block_num < block_nbr; block_num++, block++ ) {

        col_nbr = block->colnbr;
        for ( col = 0; col < col_nbr; col++ ) {

            sortptr[0] = (void*)(rowtab + block->coltab[col]);
            sortptr[1] = (void*)(valtab + block->coltab[col]);

            len_col_block = block->coltab[col+1] - block->coltab[col];
#if !defined(NDEBUG)
            {
                pastix_int_t gN = bcsc->gN;
                pastix_int_t i;

                for ( i = 0; i < len_col_block; i++ ) {
                    assert( rowtab[ block->coltab[col] + i ] >= 0 );
                    assert( rowtab[ block->coltab[col] + i ] < gN );
                }
            }
#endif

            z_qsortIntFloatAsc( sortptr, len_col_block );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Initializes a centralize pastix_complex64_t block csc.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] initAt
 *          A flag to enable/disable the initialization of A'.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped according to the distribution described in solvmtx.
 *
 * @param[in] valuesize
 *         The number of non zero unknowns in the matrix.
 *
 *******************************************************************************/
void
bcsc_zinit( const spmatrix_t     *spm,
            const pastix_order_t *ord,
            const SolverMatrix   *solvmtx,
            int                   initAt,
            pastix_bcsc_t        *bcsc,
            pastix_int_t          valuesize )
{
    pastix_int_t nbelt = 0;
#if defined(PASTIX_WITH_MPI)
    pastix_int_t nbelt_recv = 0;
    bcsc_handle_comm_t *bcsc_comm = bcsc->bcsc_comm;
#endif

    /* Exchanges the data if the spm is distributed and adds the received data of A. */
#if defined(PASTIX_WITH_MPI)
    if ( bcsc_comm != NULL ) {
        nbelt_recv = bcsc_zexchange_values_A( spm, ord, solvmtx, bcsc );
        nbelt = nbelt_recv;
    }
#endif

    /* Initializes the blocked structure of the matrix A. */
    nbelt += bcsc_zinit_A( spm, ord, solvmtx, bcsc );

    /* Initializes the blocked structure of the matrix At if the matrix is not general. */
    if ( spm->mtxtype != SpmGeneral ) {
        nbelt += bcsc_zinit_At( spm, ord, solvmtx, bcsc->rowtab, bcsc );

        /* Exchanges the data if the spm is distributed and adds the received data of A. */
#if defined(PASTIX_WITH_MPI)
        if ( bcsc_comm != NULL ) {
            nbelt_recv = bcsc_zexchange_values_At( spm, ord, solvmtx, bcsc->rowtab, bcsc );
            nbelt += nbelt_recv;
        }
#endif
    }

    /* Restores the correct coltab arrays. */
    bcsc_restore_coltab( bcsc );

    /* Sorts the csc. */
    bcsc_zsort( bcsc, bcsc->rowtab, bcsc->Lvalues );

    /* Initializes the blocked structure of the matrix At if the matrix is general. */
    if ( spm->mtxtype == SpmGeneral ) {
        /* If only the refinement is performed A^t is not required. */
        if ( initAt ) {
            pastix_int_t *trowtab, i;
            MALLOC_INTERN( bcsc->Uvalues, valuesize,  pastix_complex64_t );
            MALLOC_INTERN( trowtab, valuesize, pastix_int_t );

            for (i=0; i<valuesize; i++) {
                trowtab[i] = -1;
            }
            nbelt = bcsc_zinit_At( spm, ord, solvmtx, trowtab, bcsc );
#if defined(PASTIX_WITH_MPI)
            if ( bcsc_comm != NULL ) {
                nbelt_recv = bcsc_zexchange_values_At( spm, ord, solvmtx, trowtab, bcsc );
                nbelt += nbelt_recv;
            }
#endif
            /* Restores the correct coltab arrays */
            bcsc_restore_coltab( bcsc );

            /* Sorts the transposed csc */
            bcsc_zsort( bcsc, trowtab, bcsc->Uvalues );
            memFree( trowtab );
        }
    }
    /* In the case of SpmHermitian, conj is applied when used to save memory space */
    else {
        bcsc->Uvalues = bcsc->Lvalues;
    }

    (void) nbelt;
#if defined(PASTIX_WITH_MPI)
    (void) nbelt_recv;
#endif
}
