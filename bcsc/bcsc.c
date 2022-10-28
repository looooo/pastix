/**
 *
 * @file bcsc.c
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
 * @author Alycia Lisito
 * @date 2022-10-11
 *
 **/
#include "common.h"
#include "pastix/order.h"
#include <spm.h>
#include "blend/solver.h"
#include "bcsc/bcsc.h"

#include "bcsc/bcsc_z.h"
#include "bcsc/bcsc_c.h"
#include "bcsc/bcsc_d.h"
#include "bcsc/bcsc_s.h"

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initializes the bcsc_handle_comm_t structure.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix structure which describes the data distribution.
 *
 * @param[out] bcsc_comm
 *          The bcsc_handle_comm initialised.
 *
 *******************************************************************************/
void
bcsc_init_handle_comm( const SolverMatrix *solvmtx,
                       pastix_bcsc_t      *bcsc )
{
    pastix_int_t        size = sizeof(bcsc_handle_comm_t) + (solvmtx->clustnbr-1)*sizeof(bcsc_proc_comm_t);
    bcsc_handle_comm_t *bcsc_comm;

    bcsc->bcsc_comm = (bcsc_handle_comm_t *)malloc( size );
    bcsc_comm = bcsc->bcsc_comm;

    bcsc_comm->flttype  = bcsc->flttype;
    bcsc_comm->clustnbr = solvmtx->clustnbr;
    bcsc_comm->clustnum = solvmtx->clustnum;
    bcsc_comm->comm     = solvmtx->solv_comm;

    memset( bcsc_comm->data_comm, 0, bcsc_comm->clustnbr * sizeof(bcsc_proc_comm_t) );
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Frees the bcsc_handle_comm pointers.
 *
 *******************************************************************************
 *
 * @param[inout] bcsc_comm
 *          The bcsc_handle_comm_t structure.
 *
 *******************************************************************************/
void
bcsc_exit_handle_comm( bcsc_handle_comm_t *bcsc_comm )
{
    int c;
    int clustnbr = bcsc_comm->clustnbr;

    for ( c = 0; c < clustnbr; c++ ) {
        if( bcsc_comm->data_comm[c].indexes_A != NULL ) {
            memFree_null(bcsc_comm->data_comm[c].indexes_A);
        }
        if( bcsc_comm->data_comm[c].values_A != NULL ) {
            memFree_null(bcsc_comm->data_comm[c].values_A);
        }
        if( bcsc_comm->data_comm[c].indexes_At != NULL ) {
            memFree_null(bcsc_comm->data_comm[c].indexes_At);
        }
        if( bcsc_comm->data_comm[c].values_At != NULL ) {
            memFree_null(bcsc_comm->data_comm[c].values_At);
        }
    }
}

#if defined(PASTIX_WITH_MPI)
/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Exchanges the amount of data the current processor will send to and
 *        receive from each processor.
 *
 *******************************************************************************
 *
 * @param[in] bcsc_comm
 *          The bcsc_handle_comm_t structure.
 *
 *******************************************************************************/
void
bcsc_exchange_amount_of_data( bcsc_handle_comm_t *bcsc_comm )
{
    int               c;
    int               clustnbr    = bcsc_comm->clustnbr;
    int               clustnum    = bcsc_comm->clustnum;
    pastix_int_t      idx_A_cnt   = 0;
    pastix_int_t      val_A_cnt   = 0;
    pastix_int_t      idx_At_cnt  = 0;
    pastix_int_t      val_At_cnt  = 0;
    pastix_int_t      counter_req = 0;
    bcsc_proc_comm_t *data_comm   = bcsc_comm->data_comm;
    MPI_Status        statuses[(clustnbr-1)*2];
    MPI_Request       requests[(clustnbr-1)*2];
    size_t            size;

    /* Receives the amount of indexes and values. */
    for ( c = 0; c < clustnbr; c++ ) {
        if ( c == clustnum ) {
            continue;
        }

        MPI_Irecv( &(data_comm[c].nrecvs), 4, PASTIX_MPI_INT,
                   c, PastixTagCount, bcsc_comm->comm, &requests[counter_req++] );

        MPI_Isend( &(data_comm[c].nsends), 4, PASTIX_MPI_INT,
                   c, PastixTagCount, bcsc_comm->comm, &requests[counter_req++] );
    }

    MPI_Waitall( counter_req, requests, statuses );

    /* Saves the total amount of indexes and values received. */
    for ( c = 0; c < clustnbr; c++ ) {
        if ( c == clustnum ) {
            continue;
        }

        idx_A_cnt  += data_comm[c].nrecvs.idx_A;
        val_A_cnt  += data_comm[c].nrecvs.val_A;
        idx_At_cnt += data_comm[c].nrecvs.idx_At;
        val_At_cnt += data_comm[c].nrecvs.val_At;
    }
    data_comm[clustnum].nrecvs.idx_A  = idx_A_cnt;
    data_comm[clustnum].nrecvs.val_A  = val_A_cnt;
    data_comm[clustnum].nrecvs.idx_At = idx_At_cnt;
    data_comm[clustnum].nrecvs.val_At = val_At_cnt;

    /* Allocates the indexes and values buffers. */
    for ( c = 0; c < clustnbr; c++ ) {
        bcsc_proc_comm_t   *data   = bcsc_comm->data_comm + c;
        bcsc_data_amount_t *amount = ( c == clustnum ) ? &(data->nrecvs) : &(data->nsends);

        if ( ( amount->idx_A != 0 ) && ( data->indexes_A == NULL ) ) {
            MALLOC_INTERN( data->indexes_A,  amount->idx_A , pastix_int_t );
        }
        if ( ( amount->idx_At != 0 ) && ( data->indexes_At == NULL ) ) {
            MALLOC_INTERN( data->indexes_At, amount->idx_At, pastix_int_t );
        }
        if ( ( amount->val_A != 0 ) && ( data->values_A == NULL ) ) {
            size = amount->val_A * pastix_size_of( bcsc_comm->flttype );
            MALLOC_INTERN( data->values_A,  size, char );
        }
        if ( ( amount->val_At != 0 ) && ( data->values_At == NULL ) ) {
            size = amount->val_At * pastix_size_of( bcsc_comm->flttype );
            MALLOC_INTERN( data->values_At, size, char );
        }
    }

    return;
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Creates the array which represents the repartition of each column
 *        in the block structure. The array size is spm->gNexp where:
 *            - col2cblk[k] = cblknum, with cblknum the index of the block column
 *              where the column k is stored.
 *        This routine is called when the matrix is in shared memory.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solvmtx structure associated to the problem.
 *
 * @param[in,out] bcsc
 *           The internal block CSC structure.
 *           The number of local columns is updated.
 *
 *******************************************************************************
 *
 * @return The col2cblk array which gives the repartition of the solvmtx columns
 *         into the block structure.
 *
 *******************************************************************************/
pastix_int_t *
bcsc_init_col2cblk_shm( const SolverMatrix  *solvmtx,
                        const pastix_bcsc_t *bcsc )
{
    pastix_int_t      j;
    pastix_int_t      cblknum;
    pastix_int_t     *col2cblk;

    /* Allocates the col2cblk. */
    MALLOC_INTERN( col2cblk, bcsc->gN, pastix_int_t );
    memset( col2cblk, 0xff, bcsc->gN * sizeof(pastix_int_t) );

    const SolverCblk *cblk    = solvmtx->cblktab;
    pastix_int_t      cblknbr = solvmtx->cblknbr;
    /* Goes through the blocks. */
    for ( cblknum = 0; cblknum < cblknbr; cblknum++, cblk++ ) {
        if ( cblk->cblktype & (CBLK_FANIN|CBLK_RECV) ) {
            continue;
        }
        /*
         * Goes through the columns of the block and adds the number of
         * the block in col2cblk at the corresponding index.
         */
        for ( j = cblk->fcolnum; j <= cblk->lcolnum; j++ ) {
            col2cblk[j] = cblknum;
        }
    }

    return col2cblk;
}

#if defined(PASTIX_WITH_MPI)
/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Creates the array which represents the repartition of each column
 *        in the block structure. The array size is spm->gNexp where:
 *            - col2cblk[k] = - (owner + 1) if the column is not stored in a local block
 *            - col2cblk[k] = cblknum, if the column k is stored in a local block, with
 *              cblknum the index of this block column.
 *       This routine is called when the matrix is in distributed memory.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solvmtx structure associated to the problem.
 *
 * @param[in,out] bcsc
 *           The internal block CSC structure.
 *           The number of local columns is updated.
 *
 *******************************************************************************
 *
 * @return The col2cblk array which gives the repartition of the solvmtx columns
 *         into the block structure.
 *
 *******************************************************************************/
pastix_int_t *
bcsc_init_col2cblk_dst( const SolverMatrix  *solvmtx,
                        const pastix_bcsc_t *bcsc )
{
    pastix_int_t  n, nr = 0;
    pastix_int_t  k, j, c;
    pastix_int_t  clustnum = solvmtx->clustnum;
    pastix_int_t  clustnbr = solvmtx->clustnbr;
    pastix_int_t  fcolnum, lcolnum, cblknum;
    pastix_int_t *col2cblk;
    pastix_int_t *col2cblk_bcast = NULL;

    /* Allocates the col2cblk. */
    MALLOC_INTERN( col2cblk, bcsc->gN, pastix_int_t );
    memset( col2cblk, 0xff, bcsc->gN * sizeof(pastix_int_t) );

    for( c = 0; c < clustnbr; c++ ) {
        if ( c == clustnum ) {
            const SolverCblk *cblk    = solvmtx->cblktab;
            pastix_int_t      cblknbr = solvmtx->cblknbr;
            pastix_int_t      colcount;

            n = (solvmtx->cblknbr - solvmtx->faninnbr - solvmtx->recvnbr) * 2;

            /* Sends the size of data. */
            MPI_Bcast( &n, 1, PASTIX_MPI_INT, c, solvmtx->solv_comm );

            if ( n > nr ) {
                nr = n;
                col2cblk_bcast = (pastix_int_t *)realloc( col2cblk_bcast, nr * sizeof(pastix_int_t) );
            }

            colcount = 0;
            k = 0;
            /* Goes through the blocks. */
            for ( cblknum = 0; cblknum < cblknbr; cblknum++, cblk++ ) {
                if ( cblk->cblktype & (CBLK_FANIN|CBLK_RECV) ) {
                    continue;
                }
                /* Adds the first and last columns of the block in col2cblk_bcast. */
                col2cblk_bcast[k]   = cblk->fcolnum;
                col2cblk_bcast[k+1] = cblk->lcolnum;
                k += 2;
                /*
                 * Goes through the columns of the block and adds the
                 * block number in col2cblk.
                 */
                for ( j = cblk->fcolnum; j <= cblk->lcolnum; j++ ) {
                    colcount++;
                    col2cblk[j] = cblknum;
                }
            }
            assert( colcount == bcsc->n );

            /* Sends the col2cblk_bcast. */
            MPI_Bcast( col2cblk_bcast, n, PASTIX_MPI_INT, c, solvmtx->solv_comm );
        }
        else {
            /* Receives the size of data from c. */
            MPI_Bcast( &n, 1, PASTIX_MPI_INT, c, solvmtx->solv_comm );

            if ( n > nr ) {
                nr = n;
                col2cblk_bcast = (pastix_int_t *)realloc( col2cblk_bcast, nr * sizeof(pastix_int_t) );
            }

            if ( n == 0 ) {
                continue;
            }

            /* Receives the col2cblk_bcast from c. */
            MPI_Bcast( col2cblk_bcast, n, PASTIX_MPI_INT, c, solvmtx->solv_comm );
            /*
             * Goes through the columns in col2cblk_bcast and adds the processor
             * number in col2cblk.
             */
            for ( k = 0; k < n; k += 2 ) {
                fcolnum = col2cblk_bcast[k];
                lcolnum = col2cblk_bcast[k+1];
                for ( j = fcolnum; j <= lcolnum; j++ ) {
                    col2cblk[j] = - c - 1;
                }
            }
        }
    }

    free( col2cblk_bcast );

    return col2cblk;
}
#endif

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Creates the array which represents the repartition of each column
 *        in the block structure. This routine calls bcsc_init_col2cblk_shm or
 *        bcsc_init_col2cblk_dst according to the way the matrix is stored in the
 *        memory.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solvmtx structure associated to the problem.
 *
 * @param[in] bcsc
 *           The internal block CSC structure.
 *           The number of local columns is updated.
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 *******************************************************************************
 *
 * @return The col2cblk array which gives the repartition of the solvmtx columns
 *         into the block structure.
 *
 *******************************************************************************/
pastix_int_t *
bcsc_init_col2cblk( const SolverMatrix  *solvmtx,
                    const pastix_bcsc_t *bcsc,
                    const spmatrix_t    *spm )
{
    pastix_int_t *col2cblk;
    /* Tests if the spm is in shared or distributed memory. */
    if ( spm->loc2glob == NULL ) {
        col2cblk = bcsc_init_col2cblk_shm( solvmtx, bcsc );
    }
#if defined(PASTIX_WITH_MPI)
    else {
        col2cblk = bcsc_init_col2cblk_dst( solvmtx, bcsc );
    }
#endif
    return col2cblk;
}

/**
 *******************************************************************************
 *
 * @brief Initializes the dofshit array of size gNexp which gives
 *        dofshift[index_permuted] = index. This corresponds to the inverse of
 *        the permutation given in ord->permtab.
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
 *******************************************************************************
 *
 * @return The dofshift array.
 *
 *******************************************************************************/
static inline pastix_int_t*
bcsc_init_dofshift( const spmatrix_t     *spm,
                    const pastix_order_t *ord )
{
    pastix_int_t *dofshift, *ptr;
    pastix_int_t *dofs;
    pastix_int_t  dof;
    pastix_int_t  idof, dofj, dofidx;
    pastix_int_t  jg, jp;

    /* Allocates the dofshift array. */
    MALLOC_INTERN( dofshift, spm->gNexp, pastix_int_t );

    dofs = spm->dofs;
    dof  = spm->dof;
    ptr  = dofshift;
    for ( jg = 0; jg < spm->gN; jg++ ) {
        jp     = ord->permtab[jg];
        dofidx = (dof > 0) ? jp * dof : dofs[jg];
        ptr    = dofshift + dofidx;
        dofj   = (dof > 0) ? dof : dofs[jg+1] - dofs[jg];
        for ( idof = 0; idof < dofj; idof++, ptr++ ) {
            *ptr = jp;
        }
    }
    return dofshift;
}

/**
 *******************************************************************************
 *
 * @brief Initializes the coltab of a block csc matrix. The coltab corresponds to
 *        the number of rows (expended) per column (non expended). This rountine
 *        is called when the matrix is stored in shared memory or the matrix is
 *        replicated on the processors and the matrix's degree of liberty is
 *        constant.
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
 * @param[out] globcol
 *          The array which contains, for each column, its beginning in the
 *          smp->colptr.
 *
 *******************************************************************************/
static inline void
bcsc_init_global_coltab_shm_cdof( const spmatrix_t     *spm,
                                  const pastix_order_t *ord,
                                  pastix_int_t         *globcol )
{
    pastix_int_t *colptr   = spm->colptr;
    pastix_int_t *rowptr   = spm->rowptr;
    pastix_int_t  dof      = spm->dof;
    pastix_int_t  baseval  = spm->baseval;
    pastix_int_t  frow, lrow;
    pastix_int_t  k, j, ig, jg, ip, jp;
    int           sym = (spm->mtxtype == SpmSymmetric) || (spm->mtxtype == SpmHermitian);

    assert( dof > 0 );
    assert( spm->loc2glob == NULL );

    /* Goes through the column of the spm. */
    for ( j = 0; j < spm->n; j++, colptr++ ) {
        jg   = j;
        jp   = ord->permtab[jg];
        frow = colptr[0] - baseval;
        lrow = colptr[1] - baseval;
        assert( (lrow - frow) >= 0 );
        /* Adds the number of values in the column jg. */
        globcol[jp] += (lrow - frow) * dof;

        /*
         * Adds for At the number of values in the row ig and column jg. This
         * is not required for the general case as the spm has a symmetric
         * pattern.
         */
        if ( !sym ) {
            continue;
        }

        for ( k = frow; k < lrow; k++ ) {
            ig = rowptr[k] - baseval;
            if ( ig != jg ) {
                ip = ord->permtab[ig];
                globcol[ip] += dof;
            }
        }
    }

    return;
}

/**
 *******************************************************************************
 *
 * @brief Initializes the coltab of a block csc matrix. The coltab corresponds to
 *        the number of rows (expended) per column (non expended). This rountine
 *        is called when the matrix is stored in shared memory or the matrix is
 *        replicated on the processors and the matrix's degree of liberty is
 *        variadic.
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
 * @param[out] globcol
 *          The array which contains, for each column, its begining in the
 *          smp->colptr.
 *
 *******************************************************************************/
static inline void
bcsc_init_global_coltab_shm_vdof( const spmatrix_t     *spm,
                                  const pastix_order_t *ord,
                                  pastix_int_t         *globcol )
{
    pastix_int_t *colptr   = spm->colptr;
    pastix_int_t *rowptr   = spm->rowptr;
    pastix_int_t *dofs     = spm->dofs;
    pastix_int_t  baseval  = spm->baseval;
    pastix_int_t  frow, lrow;
    pastix_int_t  k, j, ig, jg, ip, jp;
    pastix_int_t  dofj, dofi;
    int           sym = (spm->mtxtype == SpmSymmetric) || (spm->mtxtype == SpmHermitian);

    assert( spm->dof <= 0 );
    assert( spm->loc2glob == NULL );

    /* Goes through the column of the spm. */
    for ( j=0; j<spm->n; j++, colptr++ ) {
        jg   = j;
        dofj = dofs[jg+1] - dofs[jg];
        jp   = ord->permtab[jg];
        frow = colptr[0] - baseval;
        lrow = colptr[1] - baseval;
        assert( (lrow - frow) >= 0 );

        for ( k=frow; k<lrow; k++ ) {
            ig   = rowptr[k] - baseval;
            dofi = dofs[ig+1] - dofs[ig];
            /* Adds the number of values in the row ig and column jg. */
            globcol[jp] += dofi;

            /* Adds for At the number of values in the row ig and column jg. */
            if ( sym && (ig != jg) ) {
                ip = ord->permtab[ig];
                globcol[ip] += dofj;
            }
        }
    }

    return;
}

#if defined(PASTIX_WITH_MPI)
/**
 *******************************************************************************
 *
 * @brief Initializes the coltab of a block csc matrix. The coltab corresponds to
 *        the number of rows (expended) per column (non expended). This rountine
 *        is called when the matrix is distributed in the memory and the matrix's
 *        degree of liberty is constant.
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
 * @param[in] col2cblk
 *          The array which contains the repartition of the matrix columns
 *          into the block structure.
 *
 * @param[out] globcol
 *          The array which contains, for each column, its begining in the
 *          smp->colptr.
 *
 * @param[in,out] bcsc_comm
 *          On entry, the initialised bcsc_comm structure. On exit, the
 *          bcsc_handle_comm structure which contains the amount of data to
 *          send to the other processors.
 *
 *******************************************************************************/
static inline void
bcsc_init_global_coltab_dst_cdof( const spmatrix_t     *spm,
                                  const pastix_order_t *ord,
                                  const pastix_int_t   *col2cblk,
                                  pastix_int_t         *globcol,
                                  bcsc_handle_comm_t   *bcsc_comm )
{
    pastix_int_t     *colptr    = spm->colptr;
    pastix_int_t     *rowptr    = spm->rowptr;
    pastix_int_t     *loc2glob  = spm->loc2glob;
    pastix_int_t      dof       = spm->dof;
    pastix_int_t      baseval   = spm->baseval;
    bcsc_proc_comm_t *data_comm = bcsc_comm->data_comm;
    pastix_int_t      frow, lrow;
    pastix_int_t      k, j, ig, jg, ip, jp;
    int               sym = (spm->mtxtype == SpmSymmetric) || (spm->mtxtype == SpmHermitian);
    pastix_int_t      owner;

    assert( dof > 0 );

    /* Goes through the columns of spm. */
    for ( j = 0; j < spm->n; j++, colptr++, loc2glob++ ) {
        jg = *loc2glob - baseval;
        jp = ord->permtab[jg];

        frow = colptr[0] - baseval;
        lrow = colptr[1] - baseval;
        assert( (lrow - frow) >= 0 );

        owner = col2cblk[jp * dof];

        /* The column is in a block which does not belong to the current processor. */
        if ( owner < 0 ) {
            owner = - owner - 1;
            /* Adds the number of indexes to send to the owner. */
            data_comm[owner].nsends.idx_A += ( lrow - frow ) * 2;
            /* Adds the number of values to send to the owner. */
            data_comm[owner].nsends.val_A += ( lrow - frow ) * dof * dof;
        }
        else {
            /* Adds number of values in column jp. */
            globcol[jp] += (lrow - frow) * dof;
        }

        /*
         * Adds for At the number of values in the row ig and column jg. This
         * is not required for the general case as the spm has a symmetric
         * pattern.
         */
        if ( !sym ) {
            continue;
        }

        for ( k = frow; k < lrow; k++ ) {
            ig = rowptr[k] - baseval;
            if ( ig == jg ) {
                continue;
            }
            ip    = ord->permtab[ig];
            owner = col2cblk[ip * dof];

            if ( owner < 0 ) {
                owner = - owner - 1;
                /* Adds the number of indexes to send to owner. */
                data_comm[owner].nsends.idx_At += 2;
                /* Adds the number of values to send to owner. */
                data_comm[owner].nsends.val_At += dof * dof;
            }
            else {
                /* Adds for At the number of values in column jg and row ip. */
                globcol[ip] += dof;
            }
        }
    }

    return;
}

/**
 *******************************************************************************
 *
 * @brief Initializes the coltab of a block csc matrix. The coltab corresponds to
 *        the number of rows (expended) per column (non expended). This rountine
 *        is called when the matrix is distributed in the memory and the matrix's
 *        degree of liberty is variadic.
 *
 * DO NOT CURRENTLY WORKS
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
 * @param[in] col2cblk
 *          The array which contains the repartition of the matrix columns
 *          into the block structure.
 *
 * @param[out] globcol
 *          The array which contains, for each column, its begining in the
 *          smp->colptr.
 *
 * @param[out] bcsc_comm
 *          The bcsc_handle_comm structure which contains the amount of
 *          data to send to the other processors.
 *
 *******************************************************************************/
static inline void
bcsc_init_global_coltab_dst_vdof( __attribute__((unused)) const spmatrix_t     *spm,
                                  __attribute__((unused)) const pastix_order_t *ord,
                                  __attribute__((unused)) const pastix_int_t   *col2cblk,
                                  __attribute__((unused)) pastix_int_t         *globcol,
                                  __attribute__((unused)) bcsc_handle_comm_t   *bcsc_comm )
{
    // pastix_int_t *colptr   = spm->colptr;
    // pastix_int_t *rowptr   = spm->rowptr;
    // pastix_int_t *loc2glob = spm->loc2glob;
    // pastix_int_t *dofs     = spm->dofs;
    // pastix_int_t  dof      = spm->dof;
    // pastix_int_t  baseval  = spm->baseval;
    // pastix_int_t  frow, lrow;
    // pastix_int_t  k, j, ig, jg, ip, jp;
    // pastix_int_t  dofj, dofi;
    // int           sym = (spm->mtxtype == SpmSymmetric) || (spm->mtxtype == SpmHermitian);

    // assert( dof <= 0 );

    // for ( j=0; j<spm->n; j++, colptr++, loc2glob++ )
    // {
    //     jg   = *loc2glob - baseval;
    //     jp   = ord->permtab[jg];
    //     dofj = dofs[jg+1] - dofs[jg];

    //     frow = colptr[0] - baseval;
    //     lrow = colptr[1] - baseval;
    //     assert( (lrow - frow) >= 0 );

    //     jpe    = ...;
    //     ownerj = col2cblk[jpe]; // FAUX
    //     localj = ( ownerj >= 0 );
    //     ownerj = - ownerj - 1;

    //     for ( k=frow; k<lrow; k++ )
    //     {
    //         ig   = rowptr[k] - baseval;
    //         dofi = dofs[ig+1] - dofs[ig];

    //         if ( localj ) {
    //             /* The column is local */
    //             globcol[jp] += dofi;
    //         }
    //         else {
    //             /* The column is remote */
    //             //update_counter_tosend( ownerj, 1 /* Nbr Elt */, dofi /* Nbr values */ );
    //         }

    //         if ( sym && (ig != jg) ) {
    //             ip     = ord->permtab[ig];
    //             ipe    = ...;
    //             owneri = col2cblk[ipe]; // FAUX

    //             if ( owneri >= 0 ) {
    //                 globcol[ip] += dofj;
    //             }
    //             else {
    //                 owneri = - owneri - 1;
    //                 //update_counter_tosend( owneri, 1 /* Nbr Elt */, dofj /* Nbr values */ );
    //             }
    //         }
    //     }
    // }

    return;
}

/**
 *******************************************************************************
 *
 * @brief Exchanges the indexes with the other processors.
 *
 *******************************************************************************
 *
 * @param[in,out] bcsc_comm
 *          The bcsc_handle_comm structure which contains the data the current
 *          processor has to send to the other processors on entry. On exit,
 *          the structure is updated with the received data from the other
 *          processors.
 *
 *******************************************************************************/
void
bcsc_exchange_indexes( bcsc_handle_comm_t *bcsc_comm )
{
    pastix_int_t      c;
    pastix_int_t      clustnbr    = bcsc_comm->clustnbr;
    pastix_int_t      clustnum    = bcsc_comm->clustnum;
    bcsc_proc_comm_t *data_comm   = bcsc_comm->data_comm;
    bcsc_proc_comm_t *data_local  = bcsc_comm->data_comm + clustnum;
    pastix_int_t      idx_A_cnt   = 0;
    pastix_int_t      idx_At_cnt  = 0;
    pastix_int_t      counter_req = 0;
    MPI_Status        statuses[(clustnbr-1)*4];
    MPI_Request       requests[(clustnbr-1)*4];

    for ( c = 0; c < clustnbr; c++ ) {
        data_comm = bcsc_comm->data_comm + c;
        if ( c == clustnum ) {
            continue;
        }

        /* Posts the receptions of the indexes. */
        if ( data_comm->nrecvs.idx_A != 0 ) {
            MPI_Irecv( data_local->indexes_A + idx_A_cnt, data_comm->nrecvs.idx_A,
                       PASTIX_MPI_INT, c, PastixTagIndexesA, bcsc_comm->comm, &requests[counter_req++] );
            idx_A_cnt += data_comm->nrecvs.idx_A;
        }
        if ( data_comm->nrecvs.idx_At != 0 ) {
            MPI_Irecv( data_local->indexes_At + idx_At_cnt, data_comm->nrecvs.idx_At,
                       PASTIX_MPI_INT, c, PastixTagIndexesAt, bcsc_comm->comm, &requests[counter_req++] );
            idx_At_cnt += data_comm->nrecvs.idx_At;
        }

        /* Posts the emissions of the indexes. */
        if ( data_comm->nsends.idx_A != 0 ) {
            MPI_Isend( data_comm->indexes_A,  data_comm->nsends.idx_A,
                       PASTIX_MPI_INT, c, PastixTagIndexesA, bcsc_comm->comm, &requests[counter_req++] );
        }
        if ( data_comm->nsends.idx_At != 0 ) {
            MPI_Isend( data_comm->indexes_At, data_comm->nsends.idx_At,
                       PASTIX_MPI_INT, c, PastixTagIndexesAt, bcsc_comm->comm, &requests[counter_req++] );
        }
    }

    MPI_Waitall( counter_req, requests, statuses );

    /* Checks the total amount of indexes and values received. */
    assert( data_local->nrecvs.idx_A  == idx_A_cnt  );
    assert( data_local->nrecvs.idx_At == idx_At_cnt );
}

/**
 *******************************************************************************
 *
 * @brief Updates globcol with the received indexes.
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
 * @param[out] globcol
 *          The array which contains, for each column, its begining in the
 *          smp->colptr. This array is updated with the data received from the
 *          other processors.
 *
 * @param[in] bcsc_comm
 *          The bcsc_handle_comm structure which contains the received data
 *          from the other processors.
 *
 *******************************************************************************/
static inline void
bcsc_update_globcol( const spmatrix_t     *spm,
                     const pastix_order_t *ord,
                     pastix_int_t         *globcol,
                     bcsc_handle_comm_t   *bcsc_comm )
{
    pastix_int_t     *dofs     = spm->dofs;
    pastix_int_t      dof      = spm->dof;
    pastix_int_t      k, ip, jp, jg, ig, baseval;
    pastix_int_t      clustnum    = bcsc_comm->clustnum;
    bcsc_proc_comm_t *data_local  = bcsc_comm->data_comm + clustnum;
    pastix_int_t     *indexes_A;
    pastix_int_t     *indexes_At;

    assert( ord->baseval == 0 );
    baseval = ord->baseval;

    /* Updates globcol. */
    indexes_A  = data_local->indexes_A;
    indexes_At = data_local->indexes_At;

    /* Goes through data_local->indexes_A. */
    for ( k = 0; k < data_local->nrecvs.idx_A; k += 2, indexes_A += 2 ) {
        /* Adds the element (ip, jp) received to column jp. */
        ip = indexes_A[0];
        jp = indexes_A[1];
        ig = ord->peritab[ip] - baseval;
        globcol[jp] += ( dof < 0 ) ? dofs[ ig+1 ] - dofs[ig] : dof;
    }

    /* Goes through data_local->indexes_At. */
    if ( spm->mtxtype != SpmGeneral ) {
        for ( k = 0; k < data_local->nrecvs.idx_At; k += 2, indexes_At += 2 ) {
            /* Adds the element (ip, jp) received to column ip (transpose). */
            ip = indexes_At[0];
            jp = indexes_At[1];
            jg = ord->peritab[jp] - baseval;
            globcol[ip] += ( dof < 0 ) ? dofs[ jg+1 ] - dofs[jg] : dof;
        }
    }
}
#endif

/**
 *******************************************************************************
 *
 * @brief Initializes the coltab of a block csc matrix. The coltab corresponds to
 *        the number of rows (expended) per column (non expended). This routine
 *        is calls bcsc_init_global_coltab_[shm,dst]_[c,v]dof according to the way
 *        the matrix is stored and if the degree of liberty of the matrix is
 *        constant or variadic. If the matrix is distributed in the memory, this
 *        function also calls the routines which exchange the amount of data for
 *        the communication, store the indexes and values to send and exchange
 *        the indexes.
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
 *          The array which contains the repartition of the matrix columns
 *          into the block structure.
 *
 * @param[in,out] bcsc_comm
 *          The handle_comm_structure updated with the amount of data the current
 *          processor has to send to the other processor if PASTIX_WITH_MPI = ON
 *          and the matrix is distributed in memory. If it is not the case,
 *          bcsc_comm = NULL.
 *
 *******************************************************************************
 *
 * @returns The array which contains, for each column, its begining in the
 *          smp->colptr.
 *
 *******************************************************************************/
static inline pastix_int_t*
bcsc_init_global_coltab( const spmatrix_t     *spm,
                         const pastix_order_t *ord,
                         const SolverMatrix   *solvmtx,
                         const pastix_int_t   *col2cblk,
                         bcsc_handle_comm_t   *bcsc_comm )
{
    spm_int_t *globcol;

    /*
     * Allocates and initializes globcol which contains the number of elements in
     * each column of the input matrix.
     * Globcol is equivalent to the classic colptr for the internal blocked
     * csc. The blocked csc integrates the permutation computed within order
     * structure.
     */
    MALLOC_INTERN( globcol, spm->gN+1, pastix_int_t );
    memset( globcol, 0, (spm->gN+1) * sizeof(pastix_int_t) );

    if ( bcsc_comm == NULL ) {
        if ( spm->dof > 0  ) {
            bcsc_init_global_coltab_shm_cdof( spm, ord, globcol );
        }
        else {
            bcsc_init_global_coltab_shm_vdof( spm, ord, globcol );
        }
    }
#if defined(PASTIX_WITH_MPI)
    else {
        if ( spm->dof > 0 ) {
            bcsc_init_global_coltab_dst_cdof( spm, ord, col2cblk, globcol, bcsc_comm );
        }
        else {
            bcsc_init_global_coltab_dst_vdof( spm, ord, col2cblk, globcol, bcsc_comm );
        }

        /* Exchanges the amount of data which will be sent and received. */
        bcsc_exchange_amount_of_data( bcsc_comm );

        /* Stores the indexes and values the current processor has to send to the others. */
        switch( spm->flttype ) {
            case SpmFloat:
                bcsc_sstore_data( spm, ord, col2cblk, bcsc_comm );
                break;
            case SpmDouble:
                bcsc_dstore_data( spm, ord, col2cblk, bcsc_comm );
                break;
            case SpmComplex32:
                bcsc_cstore_data( spm, ord, col2cblk, bcsc_comm );
                break;
            case SpmComplex64:
                bcsc_zstore_data( spm, ord, col2cblk, bcsc_comm );
                break;
            case SpmPattern:
            default:
                fprintf(stderr, "bcsc_init: Error unknown floating type for input spm\n");
        }

        /* Exchanges the indexes and updates globcol with the received indexes. */
        bcsc_exchange_indexes( bcsc_comm );
        bcsc_update_globcol( spm, ord, globcol, bcsc_comm );

#if !defined(NDEBUG)
        if ( spm->dof > 0 ) {
            pastix_int_t ig, ip, ipe, dofi;
            pastix_int_t nnzl = 0;
            pastix_int_t nnzg = 0;
            pastix_int_t nnz;
            for( ig=0; ig<spm->gN; ig++ ) {
                ip  = ord->permtab[ig];
                ipe = ( spm->dof > 0 ) ? ip * spm->dof : spm->dofs[ ig ] - spm->baseval;
                if ( col2cblk[ipe] < 0 ) {
                    continue;
                }

                dofi = ( spm->dof > 0 ) ? spm->dof: spm->dofs[ig+1] - spm->dofs[ig];
                nnzl += globcol[ip] * dofi;
            }
            MPI_Allreduce( &nnzl, &nnzg, 1, PASTIX_MPI_INT, MPI_SUM, bcsc_comm->comm );

            if ( spm->mtxtype != SpmGeneral ) {
                nnz = spm->gnnzexp * 2 - (spm->gN * spm->dof * spm->dof);
            }
            else {
                nnz = spm->gnnzexp;
            }
            assert( nnzg == nnz );
        }
#endif
    }

#endif

    (void)solvmtx;
    (void)col2cblk;
    return globcol;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initializes the coltab of a block csc matrix. The coltab corresponds
 *        to the number of rows (expended) per column (non expended). If the
 *        matrix is distributed in the memory, this function also calls the
 *        routines which exchange the amount of data for the communication,
 *        store the indexes and values to send and exchange the indexes.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The spm structure that stores the dofs.
 *
 * @param[in] ord
 *          The ordering which needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure which describes the data distribution.

 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the initialized coltab split per block
 *          corresponding to the input spm with the permutation applied
 *          and grouped accordingly to the distribution described in solvmtx.
 *
 *******************************************************************************
 *
 * @return The number of non zero unknowns in the matrix.
 *
 *******************************************************************************/
pastix_int_t
bcsc_init_coltab( const spmatrix_t     *spm,
                  const pastix_order_t *ord,
                  const SolverMatrix   *solvmtx,
                  pastix_bcsc_t        *bcsc )
{
    SolverCblk   *cblk;
    bcsc_cblk_t  *blockcol;
    pastix_int_t *dofshift = NULL;
    pastix_int_t *globcol  = NULL;
    pastix_int_t  cblknum, bcscnum, iter, idxcol, nodeidx, colsize;

    bcsc->cscfnbr = solvmtx->cblknbr - solvmtx->faninnbr - solvmtx->recvnbr;
    MALLOC_INTERN( bcsc->cscftab, bcsc->cscfnbr, bcsc_cblk_t );

    /* Creates an array to convert expanded indexes to not expanded indexes. */
    dofshift = bcsc_init_dofshift( spm, ord );

    /* Computes the number of rows (expanded) per column (not expanded). */
    globcol  = bcsc_init_global_coltab( spm, ord, solvmtx, bcsc->col2cblk, bcsc->bcsc_comm );

    idxcol   = 0;
    bcscnum  = 0;
    cblk     = solvmtx->cblktab;
    blockcol = bcsc->cscftab;
    for ( cblknum = 0; cblknum < solvmtx->cblknbr; cblknum++, cblk++ ) {
        if ( cblk->cblktype & (CBLK_FANIN|CBLK_RECV) ) {
            continue;
        }

        blockcol->cblknum = cblknum;
        blockcol->colnbr  = cblk_colnbr( cblk );
        assert( cblk->bcscnum == bcscnum );
        MALLOC_INTERN( blockcol->coltab, blockcol->colnbr + 1, pastix_int_t );

        blockcol->coltab[0] = idxcol;
        for ( iter = 0; iter < blockcol->colnbr; iter++ ) {
            nodeidx = dofshift[ cblk->fcolnum + iter ];
            colsize = globcol[nodeidx];
            //jpe = cblk->fcolnum + iter;
            //jp  = dofshift[ jpe ];
            //colsize = globcol[jp];
            blockcol->coltab[iter+1] = blockcol->coltab[iter] + colsize;
        }
        idxcol = blockcol->coltab[blockcol->colnbr];

        blockcol++;
        bcscnum++;
    }
    assert( (blockcol - bcsc->cscftab) == bcsc->cscfnbr );
    assert( bcscnum == bcsc->cscfnbr );

    memFree_null( globcol );
    memFree_null( dofshift );

    if ( idxcol > 0 ) {
        MALLOC_INTERN( bcsc->rowtab,  idxcol, pastix_int_t);
        MALLOC_INTERN( bcsc->Lvalues, idxcol * pastix_size_of( bcsc->flttype ), char );
    }
    else {
        bcsc->rowtab  = NULL;
        bcsc->Lvalues = NULL;
    }
    bcsc->Uvalues = NULL;

    return idxcol;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Restores the coltab array when it has been modified to initialize
 *        the row and values arrays.
 *
 *******************************************************************************
 *
 * @param[inout] bcsc
 *          On entry, the bcsc to restore.
 *          On exit, the coltab array of the bcsc is restored to the correct
 *          indexes.
 *
 *******************************************************************************/
void
bcsc_restore_coltab( pastix_bcsc_t *bcsc )
{
    bcsc_cblk_t *blockcol;
    pastix_int_t index, iter, idxcol, idxcoltmp;

    idxcol   = 0;
    blockcol = bcsc->cscftab;
    for ( index=0; index<bcsc->cscfnbr; index++, blockcol++ )
    {
        for ( iter=0; iter <= blockcol->colnbr; iter++ )
        {
            idxcoltmp = blockcol->coltab[iter];
            blockcol->coltab[iter] = idxcol;
            idxcol = idxcoltmp;
        }
    }
    return;
}

/**
 *******************************************************************************
 *
 * @brief Initializes a block csc.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] solvmtx
 *          The solver matrix structure which describes the data distribution.
 *
 * @param[in] initAt
 *          The test to know if At has to be initialized:
 *          - if initAt = 0 then the matrix is symmetric of hermitian which
 *            means that Lvalues = Uvalues so At does not need to be
 *            initialised.
 *          - if initAt = 1 then the matrix is general and which means that
 *            At needs to be initialised and computed.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped accordingly to the distribution described in solvmtx.
 *
 *******************************************************************************/
void
bcsc_init_struct( const spmatrix_t   *spm,
                  const SolverMatrix *solvmtx,
                  pastix_bcsc_t      *bcsc )
{
    pastix_int_t       *col2cblk  = NULL;

    bcsc->mtxtype = spm->mtxtype;
    bcsc->flttype = spm->flttype;
    bcsc->gN      = spm->gNexp;
    bcsc->n       = solvmtx->nodenbr;

    /*
     * Creates the col2cblk array which associates each column to a cblk
     * (expanded).
     */
    col2cblk = bcsc_init_col2cblk( solvmtx, bcsc, spm );
    bcsc->col2cblk = col2cblk;

    /*
     * Initializes the coltab array of the bcsc and allocates the rowtab and
     * Lvalues arrays.
     */
    bcsc->bcsc_comm = NULL;
    if ( spm->loc2glob != NULL ) {
        bcsc_init_handle_comm( solvmtx, bcsc );
    }
}

/**
 *******************************************************************************
 *
 * @brief Cleanup the bcsc struct. (symmetric of bcsc_init_struct)
 *
 *******************************************************************************
 *
 * @param[inout] bcsc
 *          On entry, the pointer to the initialized bcsc.
 *          On exit, the bcsc freed from the informations initialized by
 *          bcsc_init_struct().
 *
 *******************************************************************************/
void
bcsc_exit_struct( pastix_bcsc_t *bcsc )
{
    if ( bcsc->col2cblk != NULL ) {
        memFree_null( bcsc->col2cblk );
    }

    if ( bcsc->bcsc_comm != NULL ) {
        bcsc_exit_handle_comm( bcsc->bcsc_comm );
        memFree_null( bcsc->bcsc_comm );
    }
}

/**
 *******************************************************************************
 *
 * @brief Initializes a block csc.
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
 * @param[in] initAt
 *          The test to know if At has to be initialized:
 *          - if initAt = 0 then the matrix is symmetric of hermitian which
 *            means that Lvalues = Uvalues so At does not need to be
 *            initialised.
 *          - if initAt = 1 then the matrix is general and which means that
 *            At needs to be initialised and computed.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped accordingly to the distribution described in solvmtx.
 *
 *******************************************************************************/
static inline void
bcsc_init( const spmatrix_t     *spm,
           const pastix_order_t *ord,
           const SolverMatrix   *solvmtx,
           pastix_int_t          initAt,
           pastix_bcsc_t        *bcsc )
{
    pastix_int_t valuesize;

    bcsc_init_struct( spm, solvmtx, bcsc );
    valuesize = bcsc_init_coltab( spm, ord, solvmtx, bcsc );

    /*
     * Fills in the lower triangular part of the blocked csc with values and
     * rows. The upper triangular part is done later if required through LU
     * factorization.
     */
    switch( spm->flttype ) {
    case SpmFloat:
        bcsc_sinit( spm, ord, solvmtx, initAt, bcsc, valuesize );
        break;
    case SpmDouble:
        bcsc_dinit( spm, ord, solvmtx, initAt, bcsc, valuesize );
        break;
    case SpmComplex32:
        bcsc_cinit( spm, ord, solvmtx, initAt, bcsc, valuesize );
        break;
    case SpmComplex64:
        bcsc_zinit( spm, ord, solvmtx, initAt, bcsc, valuesize );
        break;
    case SpmPattern:
    default:
        fprintf(stderr, "bcsc_init: Error unknown floating type for input spm\n");
    }
}

/**
 *******************************************************************************
 *
 * @brief Initializes the block csc matrix.
 *
 * The block csc matrix is used to initialize the factorized matrix, and to
 * perform the matvec operations in refinement.
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
 * @param[in] initAt
 *          The test to know if At has to be initialized:
 *          - if initAt = 0 then the matrix is symmetric of hermitian which
 *            means that Lvalues = Uvalues so At does not need to be
 *            initialised.
 *          - if initAt = 1 then the matrix is general which means that
 *            At needs to be initialised and computed.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped accordingly to the distribution described in solvmtx.
 *
 *******************************************************************************
 *
 * @return The time spent to initialize the bcsc structure.
 *
 *******************************************************************************/
double
bcscInit( const spmatrix_t     *spm,
          const pastix_order_t *ord,
          const SolverMatrix   *solvmtx,
          pastix_int_t          initAt,
          pastix_bcsc_t        *bcsc )
{
    double time = 0.;

    assert( ord->baseval == 0 );
    assert( ord->vertnbr == spm->gN );

    clockStart(time);
    bcsc_init( spm, ord, solvmtx, initAt, bcsc );
    clockStop(time);

    return time;
}

/**
 *******************************************************************************
 *
 * @brief Frees the block csc structure but do not free the bcsc pointer.
 *
 *******************************************************************************
 *
 * @param[inout] bcsc
 *          The block csc matrix to free.
 *
 *******************************************************************************/
void
bcscExit( pastix_bcsc_t *bcsc )
{
    bcsc_cblk_t *cblk;
    pastix_int_t i;

    if ( bcsc->cscftab == NULL ) {
        return;
    }

    for ( i=0, cblk=bcsc->cscftab; i < bcsc->cscfnbr; i++, cblk++ ) {
        memFree_null( cblk->coltab );
    }

    memFree_null( bcsc->cscftab );
    memFree_null( bcsc->rowtab );

    if ( (bcsc->Uvalues != NULL) &&
         (bcsc->Uvalues != bcsc->Lvalues) ) {
        memFree_null( bcsc->Uvalues );
    }

    memFree_null( bcsc->Lvalues );

    bcsc_exit_struct( bcsc );
}
