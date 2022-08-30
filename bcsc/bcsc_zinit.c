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
 * @date 2022-07-07
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
    pastix_complex64_t *values   = (pastix_complex64_t*)(spm->values);
    pastix_complex64_t *values_c;
    pastix_int_t       *colptr   = spm->colptr;
    pastix_int_t       *rowptr   = spm->rowptr;
    pastix_int_t       *loc2glob = spm->loc2glob;
    pastix_int_t       *dofs     = spm->dofs;
    pastix_int_t        dof      = spm->dof;
    pastix_int_t        k, j, ig, jg, ip, jp, ipe, jpe;
    pastix_int_t        itercblk_A, itercblk_At, baseval;
    pastix_int_t        dofi, dofj;
    pastix_int_t        c;
    bcsc_data_amount_t *data_cnt;
    bcsc_proc_comm_t   *data_comm   = bcsc_comm->data_comm;
    pastix_int_t        clustnbr    = bcsc_comm->clustnbr;

    /*
     * Allocates and initialises the counters used to fill bcsc_comm->values
     * and bcsc_comm->indexes.
     */
    MALLOC_INTERN( data_cnt, clustnbr, bcsc_data_amount_t );
    memset( data_cnt, 0, clustnbr * sizeof(bcsc_data_amount_t) );

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
                c = - ( itercblk_A + 1 );
                data_comm = bcsc_comm->data_comm + c;

                /* Stores the indexes to send to c: (ip, jp). */
                data_comm->indexes_A[data_cnt[c].idx_A  ] = ip;
                data_comm->indexes_A[data_cnt[c].idx_A+1] = jp;
                data_cnt[c].idx_A += 2;

                /*
                 * Stores the values to send to c:
                 * dofi values in the row ip.
                 * dofj values in the column jp.
                 */
                values_c = ((pastix_complex64_t*)(data_comm->values_A)) + data_cnt[c].val_A;

                memcpy( values_c, values, dofi * dofj * sizeof(pastix_complex64_t) );
                data_cnt[c].val_A += dofi * dofj;
            }

            itercblk_At = col2cblk[ ipe ];
            /* The block of the row ip does not belong to the current processor. */
            if ( (itercblk_At < 0) && (ip != jp) ) {
                /* Gets the owner of the block. */
                c = - ( itercblk_At + 1 );
                data_comm = bcsc_comm->data_comm + c;

                /* Stores the indexes to send to c: (jp, ip). */
                data_comm->indexes_At[data_cnt[c].idx_At  ] = ip;
                data_comm->indexes_At[data_cnt[c].idx_At+1] = jp;
                data_cnt[c].idx_At += 2;

                /*
                 * Stores the values to send to c:
                 * dofi values in the row ip.
                 * dofj values in the column jp.
                 */
                values_c = ((pastix_complex64_t*)(data_comm->values_At)) + data_cnt[c].val_At;

                memcpy( values_c, values, dofi * dofj * sizeof(pastix_complex64_t) );
                data_cnt[c].val_At += dofi * dofj;
            }
            values += dofi * dofj;
        }
    }
    memFree_null( data_cnt );
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Exchanges the values stored in the bcsc_comm structure.
 *
 *******************************************************************************
 *
 * @param[inout] bcsc_comm
 *          The structure in which the sending and receiving data are stored.
 *
 *******************************************************************************/
static inline void
bcsc_zexchange_values( bcsc_handle_comm_t *bcsc_comm )
{
    pastix_int_t      c;
    pastix_int_t      clustnbr    = bcsc_comm->clustnbr;
    pastix_int_t      clustnum     = bcsc_comm->clustnum;
    bcsc_proc_comm_t *data_comm   = bcsc_comm->data_comm;
    bcsc_proc_comm_t *data_local  = bcsc_comm->data_comm + clustnum;
    pastix_int_t      val_A_cnt   = 0;
    pastix_int_t      val_At_cnt  = 0;
    pastix_int_t      counter_req = 0;
    MPI_Status        statuses[(clustnbr-1)*4];
    MPI_Request       requests[(clustnbr-1)*4];

    /* Receives the values. */
    for ( c = 0; c < clustnbr; c++ ) {
        data_comm = bcsc_comm->data_comm + c;
        if ( c == clustnum ) {
            continue;
        }

        /* Posts the receptions of the values. */
        if ( data_comm->nrecvs.val_A != 0 ) {
            MPI_Irecv( ((pastix_complex64_t*)(data_local->values_A)) + val_A_cnt, data_comm->nrecvs.val_A,
                       PASTIX_MPI_COMPLEX64, c, PastixTagValuesA, bcsc_comm->comm, &requests[counter_req++] );
            val_A_cnt += data_comm->nrecvs.val_A;
        }
        if ( data_comm->nrecvs.val_At != 0 ) {
            MPI_Irecv( ((pastix_complex64_t*)(data_local->values_At)) + val_At_cnt, data_comm->nrecvs.val_At,
                       PASTIX_MPI_COMPLEX64, c, PastixTagValuesAt, bcsc_comm->comm, &requests[counter_req++] );
            val_At_cnt += data_comm->nrecvs.val_At;
        }

        /* Posts the emissions of the values. */
        if ( data_comm->nsends.val_A != 0 ) {
            MPI_Isend( data_comm->values_A, data_comm->nsends.val_A,
                       PASTIX_MPI_COMPLEX64, c, PastixTagValuesA, bcsc_comm->comm, &requests[counter_req++] );
        }
        if ( data_comm->nsends.val_At != 0 ) {
            MPI_Isend( data_comm->values_At, data_comm->nsends.val_At,
                       PASTIX_MPI_COMPLEX64, c, PastixTagValuesAt, bcsc_comm->comm, &requests[counter_req++] );
        }
    }

    MPI_Waitall( counter_req, requests, statuses );

    /* Checks the total amount of values and values received. */
    assert( data_local->nrecvs.val_A  == val_A_cnt );
    assert( data_local->nrecvs.val_At == val_At_cnt );
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
 * @param[in] col2cblk
 *          The array of matching columns with cblk indexes.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc fields are updated.
 *
 * @param[in] bcsc_comm
 *          The structure which stores the exchanged data.
 *
 *******************************************************************************/
static inline pastix_int_t
bcsc_zadd_remote_A( const spmatrix_t         *spm,
                    const pastix_order_t     *ord,
                    const SolverMatrix       *solvmtx,
                    const pastix_int_t       *col2cblk,
                    pastix_bcsc_t            *bcsc,
                    const bcsc_handle_comm_t *bcsc_comm)
{
    pastix_int_t       *dofs = spm->dofs;
    pastix_int_t        dof  = spm->dof;
    SolverCblk         *cblk;
    pastix_int_t       *coltab;
    pastix_int_t        k, ip, jp, ig, jg, ipe, jpe;
    pastix_int_t        itercblk;
    pastix_int_t        idofi, idofj, dofi, dofj;
    pastix_int_t        colidx, rowidx, pos;
    pastix_int_t        clustnum  = bcsc_comm->clustnum;
    pastix_complex64_t *Values    = bcsc->Lvalues;
    pastix_int_t        size      = bcsc_comm->data_comm[clustnum].nrecvs.idx_A;
    pastix_int_t       *indexes   = bcsc_comm->data_comm[clustnum].indexes_A;
    pastix_complex64_t *rvalues   = bcsc_comm->data_comm[clustnum].values_A;
    pastix_int_t        nbelt     = 0;

    /*
     * Goes through the received indexes in order to add the data
     * received.
     */
    for ( k = 0; k < size; k+=2, indexes+=2 ) {
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

        itercblk = col2cblk[ jpe ];
        assert( itercblk >= 0 );

        /* Gets the block on which the data received will be added. */
        cblk   = solvmtx->cblktab + itercblk;
        coltab = bcsc->cscftab[cblk->bcscnum].coltab;

        /* Goes through the values at (ip, jp). */
        colidx = jpe - cblk->fcolnum;
        for ( idofj = 0; idofj < dofj; idofj++, colidx++ ) {
            rowidx = ipe;
            pos    = coltab[ colidx ];
            for ( idofi = 0; idofi < dofi; idofi++, rowidx++, pos++, rvalues++ ) {
                /* Adds the row ip. */
                assert( rowidx >= 0 );
                assert( rowidx < spm->gNexp );
                bcsc->rowtab[ pos ] = rowidx;
                /* Adds the data at row ip and column jp. */
                Values[ pos ]       = *rvalues;
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
 * @param[in] col2cblk
 *          The array of matching columns with cblk indexes.
 *
 * @param[inout] rowtab
 *          The row tab of the bcsc or the row tab associated to At.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc fields are updated.
 *
 * @param[in] bcsc_comm
 *          The structure which stores the exchanged data.
 *
 *******************************************************************************/
static inline pastix_int_t
bcsc_zadd_remote_At( const spmatrix_t         *spm,
                     const pastix_order_t     *ord,
                     const SolverMatrix       *solvmtx,
                     const pastix_int_t       *col2cblk,
                     pastix_int_t             *rowtab,
                     pastix_bcsc_t            *bcsc,
                     const bcsc_handle_comm_t *bcsc_comm )
{
    pastix_int_t       *dofs = spm->dofs;
    pastix_int_t        dof  = spm->dof;
    SolverCblk         *cblk;
    pastix_int_t       *coltab;
    pastix_int_t        k, ip, jp, ig, jg, ipe, jpe;
    pastix_int_t        itercblk;
    pastix_int_t        idofi, idofj, dofi, dofj;
    pastix_int_t        colidx, rowidx, pos;
    pastix_int_t        clustnum  = bcsc_comm->clustnum;
    pastix_complex64_t *Values;
    pastix_int_t        size      = bcsc_comm->data_comm[clustnum].nrecvs.idx_At;
    pastix_int_t       *indexes   = bcsc_comm->data_comm[clustnum].indexes_At;
    pastix_complex64_t *rvalues   = bcsc_comm->data_comm[clustnum].values_At;
    pastix_complex64_t (*_bcsc_conj)(pastix_complex64_t) = NULL;
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
    for ( k = 0; k < size; k+=2, indexes+=2 ) {
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

        itercblk = col2cblk[ ipe ];
        assert( itercblk >= 0 );

        /* Gets the block on which the data received will be added. */
        cblk   = solvmtx->cblktab + itercblk;
        coltab = bcsc->cscftab[cblk->bcscnum].coltab;

        /* Goes through the values at (ip, jp). */
        for ( idofj = 0; idofj < dofj; idofj++ ) {
            /* Gets the column ip. */
            rowidx = jpe + idofj;
            colidx = ipe - cblk->fcolnum;
            for ( idofi = 0; idofi < dofi; idofi++, colidx++, rvalues++ ) {
                pos = coltab[ colidx ];

                /* Adds the row jp. */
                assert( rowidx >= 0 );
                assert( rowidx < spm->gNexp );
                rowtab[ pos ] = rowidx;

                /* Adds the data at row jp and column ip. */
                Values[ pos ] = _bcsc_conj( *rvalues );
                coltab[ colidx ] ++;
                assert( coltab[colidx] <= coltab[colidx+1] );
                nbelt++;
            }
        }
    }
    return nbelt;
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
 * @param[in] col2cblk
 *          The array of matching columns with local cblk indexes. If the column
 *          k belongs to a remote block then col2cblk[k] = -(owner_proc + 1).
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc fields are updated.
 *
 * @param[inout] bcsc_comm
 *          On entry, an empty structure.
 *          On exit, stores the amount of data to send and allocates the
 *          reception buffer.
 *
 *******************************************************************************/
static inline pastix_int_t
bcsc_zinit_A( const spmatrix_t     *spm,
              const pastix_order_t *ord,
              const SolverMatrix   *solvmtx,
              const pastix_int_t   *col2cblk,
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

        itercblk = col2cblk[ jpe ];

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
 * @param[in] col2cblk
 *          The array of matching columns with local cblk indexes. If the column
 *          k belongs to a remote block then col2cblk[k] = -(owner_proc + 1).
 *
 * @param[inout] rowtab
 *          The row tab of the bcsc or the row tab associated to At.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc fields are updated.
 *
 * @param[inout] bcsc_comm
 *          On entry, contains the amount of data to send.
 *          On exit, the data have been exchanged and stored in the bcsc.
 *
 *******************************************************************************/
static inline pastix_int_t
bcsc_zinit_At( const spmatrix_t     *spm,
               const pastix_order_t *ord,
               const SolverMatrix   *solvmtx,
               const pastix_int_t   *col2cblk,
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
            itercblk = col2cblk[ ipe ];

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
 * @param[in] col2cblk
 *          Array of matching column with cblk indexes.
 *
 * @param[in] initAt
 *          A flag to enable/disable the initialization of A'.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped according to the distribution described in solvmtx.
 *
 *******************************************************************************/
void
bcsc_zinit( const spmatrix_t     *spm,
            const pastix_order_t *ord,
            const SolverMatrix   *solvmtx,
            const pastix_int_t   *col2cblk,
            int                   initAt,
            pastix_bcsc_t        *bcsc,
            bcsc_handle_comm_t   *bcsc_comm,
            pastix_int_t          valuesize )
{
    pastix_int_t nbelt;
#if defined(PASTIX_WITH_MPI)
    pastix_int_t nbelt_recv;
#endif

    /* Exchanges the data if the spm is distributed */
#if defined(PASTIX_WITH_MPI)
    if ( bcsc_comm != NULL ) {
        bcsc_zexchange_values( bcsc_comm );
    }
#endif

    /* Initializes the blocked structure of the matrix A. */
    nbelt = bcsc_zinit_A( spm, ord, solvmtx, col2cblk, bcsc );

    /* Adds the received data of A if the spm is distributed */
#if defined(PASTIX_WITH_MPI)
    if ( bcsc_comm != NULL ) {
        nbelt_recv = bcsc_zadd_remote_A( spm, ord, solvmtx, col2cblk, bcsc, bcsc_comm );
        assert( nbelt_recv == bcsc_comm->data_comm[bcsc_comm->clustnum].nrecvs.val_A );
        nbelt += nbelt_recv;
    }
#endif

    /* Initializes the blocked structure of the matrix At if the matrix is not general. */
    if ( spm->mtxtype != SpmGeneral ) {
        nbelt += bcsc_zinit_At( spm, ord, solvmtx, col2cblk, bcsc->rowtab, bcsc );

        /* Adds the received data of At if the spm is distributed */
#if defined(PASTIX_WITH_MPI)
        if ( bcsc_comm != NULL ) {
            nbelt_recv = bcsc_zadd_remote_At( spm, ord, solvmtx, col2cblk, bcsc->rowtab, bcsc, bcsc_comm );
            assert( nbelt_recv == bcsc_comm->data_comm[bcsc_comm->clustnum].nrecvs.val_At );
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
            nbelt = bcsc_zinit_At( spm, ord, solvmtx, col2cblk, trowtab, bcsc );
#if defined(PASTIX_WITH_MPI)
            if ( bcsc_comm != NULL ) {
                nbelt += bcsc_zadd_remote_At( spm, ord, solvmtx, col2cblk, trowtab, bcsc, bcsc_comm );
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

    (void) bcsc_comm;
    (void) nbelt;
#if defined(PASTIX_WITH_MPI)
    (void) nbelt_recv;
#endif
}
