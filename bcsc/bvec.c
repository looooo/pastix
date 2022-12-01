/**
 *
 * @file bvec.c
 *
 * @copyright 2004-2021 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.2.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Vincent Bridonneau
 * @author Alycia Lisito
 * @date 2022-12-05
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
    pastix_int_t        size;

    /* Initializes bvec_comm. */
    size = sizeof(bvec_handle_comm_t) + (solvmtx->clustnbr-1)*sizeof(bvec_proc_comm_t);

    Pb->rhs_comm = (bvec_handle_comm_t *)malloc( size );
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

        if ( data->send.idxbuf != NULL ) {
            memFree_null( data->send.idxbuf );
        }
        if ( data->send.valbuf != NULL ) {
            memFree_null( data->send.valbuf );
        }
        if ( data->recv.idxbuf != NULL ) {
            memFree_null( data->recv.idxbuf );
        }
        if ( data->recv.valbuf != NULL ) {
            memFree_null( data->recv.valbuf );
        }
    }
    rhs_comm->max_idx  = 0;
    rhs_comm->max_val  = 0;

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
    igpe    = ( dof > 0 ) ? igp * dof : dofs[ ig ] - basespm/* vdof incorect */;
    cblknum = col2cblk[ igpe ];
    if ( cblknum >= 0 ) {
        cblk = solvmatr->cblktab + cblknum;
        ilpe  = cblk->lcolidx + igpe - cblk->fcolnum;
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
bvec_Ploc2Pglob( pastix_data_t *pastix_data,
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
 *          On entry the rhs_comm of the permuted vector initialised.
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
            MPI_Bcast( &(data_comm->send), 2, PASTIX_MPI_INT, c, rhs_comm->comm );
            continue;
        }

        MPI_Bcast( &(data_comm->recv), 2, PASTIX_MPI_INT, c, rhs_comm->comm );

        idx_cnt = data_comm->recv.idxcnt;
        val_cnt = data_comm->recv.valcnt;

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
