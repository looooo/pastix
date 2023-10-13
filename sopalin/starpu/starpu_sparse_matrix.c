/**
 *
 * @file starpu_sparse_matrix.c
 *
 * PaStiX sparse matrix descriptor for StarPU.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Tony Delarue
 * @author Alycia Lisito
 * @author Nolan Bredel
 * @date 2023-12-01
 *
 * @ingroup pastix_starpu
 * @{
 *
 **/
#include "common.h"
#include "blend/solver.h"
#include "pastix_starpu.h"
#include "pastix_zcores.h"
#include <starpu_data.h>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
static inline void
pastix_starpu_filter_interface( void                      *father_interface,
                                void                      *child_interface,
                                struct starpu_data_filter *f,
                                unsigned                   id,
                                unsigned                   nchunks )
{
    pastix_starpu_interface_t *father   = (pastix_starpu_interface_t *)father_interface;
    pastix_starpu_interface_t *child    = (pastix_starpu_interface_t *)child_interface;
    size_t                    *sizetab  = (size_t *)f->filter_arg_ptr;
    size_t                     childoff = 0;

    assert( father->id == PASTIX_STARPU_INTERFACE_ID );
    assert( father->offset == -1 );
    assert( father->cblk->cblktype & CBLK_LAYOUT_2D );

    child->id        = father->id;
    child->flttype   = father->flttype;
    child->offset    = sizetab[id];
    child->nbblok    = sizetab[id+1] - sizetab[id];
    child->allocsize = 0;
    child->cblk      = father->cblk;
    child->dataptr   = NULL;

    assert( child->offset >= 0 );

    if ( father->dataptr == NULL ) {
        return;
    }

    if ( father->cblk->cblktype & CBLK_COMPRESSED ) {
        childoff         = sizetab[id]   * sizeof( pastix_lrblock_t );
        child->allocsize = child->nbblok * sizeof( pastix_lrblock_t );
    }
    else {
        SolverBlok *blok = father->cblk->fblokptr + sizetab[id];
        SolverBlok *lblk = father->cblk->fblokptr + sizetab[id+1];
        childoff         = pastix_size_of( father->flttype ) * blok->coefind;

        if ( lblk->coefind == 0 ) {
            int i;
            int nbrow = 0;
            for ( i=0; i<child->nbblok; i++, blok++) {
                nbrow += blok_rownbr( blok );
            }
            child->allocsize = pastix_size_of( father->flttype ) * nbrow * cblk_colnbr( father->cblk );
        }
        else {
            child->allocsize = pastix_size_of( father->flttype ) * (lblk->coefind - blok->coefind);
        }
    }

#if defined(PASTIX_STARPU_INTERFACE_DEBUG)
    fprintf( stderr,
             "blok (%9s, size=%8zu, nbblok=%2ld )\n",
             child->cblk->cblktype & CBLK_COMPRESSED ? "Low-rank" : "Full-rank",
             child->allocsize, (long)(child->nbblok) );
#endif

    assert( child->allocsize > 0 );

    child->dataptr = ((char*)father->dataptr) + childoff;

    (void)nchunks;
}

static inline void
pastix_starpu_register_interface( const starpu_sparse_matrix_desc_t *spmtx,
                                  SolverCblk                        *cblk,
                                  int                                side,
                                  pastix_coeftype_t                  flttype )
{
    starpu_data_handle_t *handler = ( (starpu_data_handle_t *)( cblk->handler ) ) + side;

    pastix_starpu_register( handler, cblk, side, flttype );
#if defined( PASTIX_WITH_MPI )
    int64_t tag_cblk;
    if ( cblk->cblktype & CBLK_FANIN ) {
        tag_cblk = 2 * cblk->gfaninnum + side;
        starpu_mpi_data_register( *handler, spmtx->mpitag + tag_cblk, spmtx->solvmtx->clustnum );
        starpu_data_set_reduction_methods( *handler, NULL, &cl_fanin_init_cpu );
    }
    else {
        if ( cblk->cblktype & CBLK_RECV ) {
            tag_cblk = 2 * cblk->gfaninnum + side;
        }
        else {
            tag_cblk = 2 * cblk->gcblknum + side;
        }
        starpu_mpi_data_register( *handler, spmtx->mpitag + tag_cblk, cblk->ownerid );
    }
#if defined(PASTIX_DEBUG_STARPU)
    fprintf( stderr, "[%2d][pastix][%s] Matrix cblk=%d, owner=%d, tag=%ld, size=%ld\n",
             spmtx->solvmtx->clustnum, __func__, cblk->gcblknum, cblk->ownerid,
             spmtx->mpitag + tag_cblk,
             cblk->stride * cblk_colnbr( cblk ) * pastix_size_of( spmtx->solvmtx->flttype ) );
#endif
#endif
    (void)spmtx;
}

static inline void
pastix_starpu_register_cblk( const starpu_sparse_matrix_desc_t *spmtx,
                             SolverCblk                        *cblk,
                             pastix_coeftype_t                  flttype )
{
    pastix_starpu_register_interface( spmtx, cblk, PastixLCoef, flttype );
    if ( spmtx->mtxtype == PastixGeneral ) {
        pastix_starpu_register_interface( spmtx, cblk, PastixUCoef, flttype );
    }
}

#if defined( PASTIX_WITH_MPI )
static inline void
pastix_starpu_mpi_register_blok( const starpu_sparse_matrix_desc_t *spmtx,
                                 const SolverCblk                  *cblk,
                                 SolverBlok                        *blok,
                                 int64_t                            tag_desc )
{
    int64_t tag;
    int     ownerid = cblk->ownerid;

    if ( cblk->cblktype & CBLK_FANIN ) {
        tag = tag_desc + 2 * blok->gfaninnm;
        ownerid = spmtx->solvmtx->clustnum;
    }
    else if ( cblk->cblktype & CBLK_RECV ) {
        tag = tag_desc + 2 * blok->gfaninnm;
    }
    else {
        tag = tag_desc + 2 * blok->gbloknm;
    }

    starpu_mpi_data_register( blok->handler[0], tag, ownerid );
    if ( spmtx->mtxtype == PastixGeneral ) {
        tag = tag + 1;
        starpu_mpi_data_register( blok->handler[1], tag, ownerid );
    }
}
#else
static inline void
pastix_starpu_mpi_register_blok( const starpu_sparse_matrix_desc_t *spmtx,
                                 const SolverCblk                  *cblk,
                                 SolverBlok                        *blok,
                                 int64_t                            tag_desc )
{
    (void)spmtx;
    (void)cblk;
    (void)blok;
    (void)tag_desc;
}
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 *******************************************************************************
 *
 * @brief Generate the StarPU descriptor of the sparse matrix.
 *
 * This function creates the StarPU descriptor that will provide tha data
 * mapping and memory location to StarPU for the computation.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix structure that describes the sparse matrix for
 *          PaStiX.
 *
 * @param[in] mtxtype
 *          The type of sparse matrix to describe.
 *          @arg PastixGeneral:   The sparse matrix is general.
 *          @arg PastixSymmetric: The sparse matrix is lower triangular symmetric.
 *          @arg PastixHermitian: The sparse matrix is lower triangular hermitian.
 *
 * @param[in] nodes
 *          The number of processes used to solve the problem.
 *
 * @param[in] myrank
 *          The rank of the calling process.
 *
 * @param[in] flttype
 *          The memory size of the arithmetic used to store the matrix
 *          coefficients.
 *
 ******************************************************************************/
void
starpu_sparse_matrix_init( SolverMatrix     *solvmtx,
                           pastix_mtxtype_t  mtxtype,
                           int               nodes,
                           int               myrank,
                           pastix_coeftype_t flttype )
{
    pastix_int_t cblknbr, cblkmin2d;
    size_t       key1, key2;
    SolverCblk  *cblk;
    SolverBlok  *blok, *lblok;
    pastix_int_t n = 0, cblknum;
    pastix_int_t nbrow;
    size_t       size;
    int64_t      tag_desc;

    starpu_sparse_matrix_desc_t *spmtx = solvmtx->starpu_desc;
    if ( spmtx != NULL ) {
        starpu_sparse_matrix_destroy( spmtx );
    }
    else {
        spmtx = (starpu_sparse_matrix_desc_t *)malloc( sizeof( starpu_sparse_matrix_desc_t ) );
    }

    tag_desc              = ( (int64_t) (solvmtx->gcblknbr + solvmtx->gbloknbr + solvmtx->gfanincblknbr + solvmtx->gfaninbloknbr) ) * 2;
    spmtx->mpitag         = pastix_starpu_tag_book( tag_desc );
    tag_desc              = spmtx->mpitag + 2 * ( solvmtx->gcblknbr + solvmtx->gfanincblknbr );
    spmtx->typesze        = pastix_size_of( flttype );
    spmtx->mtxtype        = mtxtype;
    spmtx->solvmtx        = solvmtx;
    spmtx->cblktab_handle = NULL;
    spmtx->gpu_blocktab   = NULL;

    cblknbr   = solvmtx->cblknbr;
    cblkmin2d = solvmtx->cblkmin2d;
    key1      = 2 * cblknbr;

    /* Initialize 1D cblk handlers */
    cblk = spmtx->solvmtx->cblktab;
    for ( cblknum = 0; cblknum < cblkmin2d; cblknum++, n++, cblk++ ) {
        pastix_starpu_register_cblk( spmtx, cblk, flttype );
    }

    /* Initialize 2D cblk handlers */
    if ( cblkmin2d < cblknbr ) {
        struct starpu_data_filter filter = { .filter_func = pastix_starpu_filter_interface };
        starpu_cblk_t            *cblkhandle;
        size_t                   *sizetab = NULL;
        pastix_int_t              nchildren, sizenbr = 0;

        spmtx->cblktab_handle =
            (starpu_cblk_t *)malloc( ( cblknbr - cblkmin2d ) * sizeof( starpu_cblk_t ) );

        cblk       = spmtx->solvmtx->cblktab + cblkmin2d;
        cblkhandle = spmtx->cblktab_handle;

        sizenbr = ( cblk[1].fblokptr - cblk[0].fblokptr ) + 1;
        sizetab = malloc( sizenbr * sizeof( size_t ) );
        assert( sizenbr >= 1 );

        for ( cblknum = cblkmin2d, n = 0; cblknum < cblknbr;
              cblknum++, n++, cblk++, cblkhandle++ ) {
            pastix_starpu_register_cblk( spmtx, cblk, flttype );

            if ( !(cblk->cblktype & CBLK_TASKS_2D) ) {
                continue;
            }

            /* Let's build the sizetab array */
            blok  = cblk[0].fblokptr;
            lblok = cblk[1].fblokptr;

            if ( ( lblok - blok ) >= sizenbr ) {
                sizenbr = ( lblok - blok ) + 1;
                free( sizetab );
                sizetab = malloc( sizenbr * sizeof( size_t ) );
            }
            nchildren  = 0;
            sizetab[0] = 0;

            /*
             * Diagonal block
             */
            sizetab[nchildren + 1] = 1;
            nchildren++;

            /*
             * Off-diagonal blocks
             */
            blok++;
            for ( ; blok < lblok; blok++ ) {
                nbrow = 1;

                while ( ( blok + 1 < lblok ) &&
                        ( blok[0].fcblknm == blok[1].fcblknm ) &&
                        ( blok[0].lcblknm == blok[1].lcblknm ) )
                {
                    blok++;
                    nbrow++;
                }
                size = nbrow;

                sizetab[nchildren + 1] = sizetab[nchildren] + size;
                nchildren++;
            }
            filter.nchildren      = nchildren;
            filter.filter_arg_ptr = sizetab;

            cblkhandle->handlenbr = nchildren;
            if ( mtxtype == PastixGeneral ) {
                cblkhandle->handletab = (starpu_data_handle_t *)malloc(
                    2 * nchildren * sizeof( starpu_data_handle_t ) );

                starpu_data_partition_plan( cblk->handler[0], &filter, cblkhandle->handletab );

                starpu_data_partition_plan(
                    cblk->handler[1], &filter, cblkhandle->handletab + nchildren );
            }
            else {
                cblkhandle->handletab =
                    (starpu_data_handle_t *)malloc( nchildren * sizeof( starpu_data_handle_t ) );

                starpu_data_partition_plan( cblk->handler[0], &filter, cblkhandle->handletab );
            }

            nchildren = 0;
            blok      = cblk[0].fblokptr;
            lblok     = cblk[1].fblokptr;

            /*
             * Diagonal block
             */
            blok->handler[0] = cblkhandle->handletab[nchildren];
            if ( mtxtype == PastixGeneral ) {
                blok->handler[1] = cblkhandle->handletab[cblkhandle->handlenbr + nchildren];
            }
            else {
                blok->handler[1] = NULL;
            }
            pastix_starpu_mpi_register_blok( spmtx, cblk, blok, tag_desc );
            nchildren++;

            /*
             * Off-diagonal blocks
             */
            blok++;
            for ( ; blok < lblok; blok++ ) {
                blok->handler[0] = cblkhandle->handletab[nchildren];
                if ( mtxtype == PastixGeneral ) {
                    blok->handler[1] = cblkhandle->handletab[cblkhandle->handlenbr + nchildren];
                }
                else {
                    blok->handler[1] = NULL;
                }
                pastix_starpu_mpi_register_blok( spmtx, cblk, blok, tag_desc );
                nchildren++;

                while ( ( blok < lblok ) && ( blok[0].fcblknm == blok[1].fcblknm ) &&
                        ( blok[0].lcblknm == blok[1].lcblknm ) ) {
                    blok++;
                    blok->handler[0] = NULL;
                    blok->handler[1] = NULL;
                }
            }
        }

        if ( sizetab != NULL ) {
            free( sizetab );
        }
    }
    solvmtx->starpu_desc = spmtx;

    (void)key1;
    (void)key2;
    (void)nodes;
    (void)myrank;
    (void)tag_desc;
}

/**
 *******************************************************************************
 *
 * @brief Submit asynchronous calls to retrieve the data on main memory.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          The cblk for which the data needs to be flushed.
 *
 * @param[inout] cblk
 *          The cblk for which the data needs to be flushed.
 *
 ******************************************************************************/
void
starpu_sparse_cblk_wont_use( pastix_coefside_t side, SolverCblk *cblk )
{
    if ( side != PastixUCoef ) {
        assert( cblk->handler[0] );
        starpu_data_wont_use( cblk->handler[0] );
    }
    if ( side != PastixLCoef ) {
        assert( cblk->handler[1] );
        starpu_data_wont_use( cblk->handler[1] );
    }
}

/**
 *******************************************************************************
 *
 * @brief Submit asynchronous calls to retrieve the data on main memory.
 *
 *******************************************************************************
 *
 * @param[inout] spmtx
 *          The sparse matrix descriptor to retrieve on main memory.
 *
 ******************************************************************************/
void
starpu_sparse_matrix_getoncpu( starpu_sparse_matrix_desc_t *spmtx )
{
    SolverCblk  *cblk;
    pastix_int_t i;

    cblk = spmtx->solvmtx->cblktab;
    for ( i = 0; i < spmtx->solvmtx->cblknbr; i++, cblk++ ) {
        assert( cblk->handler[0] );

        if ( cblk->ownerid != spmtx->solvmtx->clustnum ) {
            continue;
        }

        starpu_data_acquire_cb( cblk->handler[0],
                                STARPU_R,
                                ( void( * )( void * ) ) & starpu_data_release,
                                cblk->handler[0] );

        if ( cblk->ucoeftab ) {
            starpu_data_acquire_cb( cblk->handler[1],
                                    STARPU_R,
                                    ( void( * )( void * ) ) & starpu_data_release,
                                    cblk->handler[1] );
        }
    }
}

#if !defined(DOXYGEN_SHOULD_SKIP_THIS)
#if !defined(HAVE_STARPU_DATA_PARTITION_CLEAN_NODE)
/**
 *******************************************************************************
 *
 * @brief StarPU function only avalaible with StarPU >= 1.4.0
 *
 ******************************************************************************/
static inline void
pastix_starpu_data_partition_clean_node( int                   unpartition,
                                         starpu_data_handle_t  root_handle,
                                         unsigned              nparts,
                                         starpu_data_handle_t *children,
                                         int                   gather_node )
{
    if ( unpartition ) {
        starpu_data_unpartition_submit( root_handle, nparts, children, gather_node );
    }
    starpu_data_partition_clean( root_handle, nparts, children );
}
#else
static inline void
pastix_starpu_data_partition_clean_node( int                   unpartition __attribute__((unused)),
                                         starpu_data_handle_t  root_handle,
                                         unsigned              nparts,
                                         starpu_data_handle_t *children,
                                         int                   gather_node )
{
    starpu_data_partition_clean_node( root_handle, nparts, children, gather_node );
}
#endif
#endif

/**
 *******************************************************************************
 *
 * @brief Destroy a single cblk StarPU data structure of the sparse matrix
 *
 *******************************************************************************
 *
 * @param[in] is_owner
 *          Boolean to specify if the calling process owns the cblk or not
 *
 * @param[in,out] cblk
 *          Pointer to the cblk data structure
 *
 * @param[in,out] starpu_cblk
 *          Pointer to the starpu data structure associated to the cblk
 *
 ******************************************************************************/
void
pastix_starpu_cblk_destroy( int            is_owner,
                            SolverCblk    *cblk,
                            starpu_cblk_t *cblkhandle )
{
    if ( cblk->cblktype & CBLK_TASKS_2D ) {
        int gather_node = -1;

        if ( is_owner && !(cblk->cblktype & CBLK_FANIN) ) {
            gather_node = STARPU_MAIN_RAM;
        }

        /*
         * First, let's unpartition ourself as long as StarPU does not use the
         * correct gather_node in starpu_data_partition_clean()
         */
        if ( cblk->handler[0] ) {
            pastix_starpu_data_partition_clean_node( !is_owner && (cblk->partitioned & (PastixLCoef + 1)),
                                                     cblk->handler[0],
                                                     cblkhandle->handlenbr,
                                                     cblkhandle->handletab,
                                                     gather_node );
        }

        if ( cblk->handler[1] ) {
            pastix_starpu_data_partition_clean_node( !is_owner && (cblk->partitioned & (PastixUCoef + 1)),
                                                     cblk->handler[1],
                                                     cblkhandle->handlenbr,
                                                     cblkhandle->handletab + cblkhandle->handlenbr,
                                                     gather_node );
        }

        free( cblkhandle->handletab );
        cblkhandle->handletab = NULL;
    }

    if ( cblk->handler[0] ) {
        starpu_data_unregister( cblk->handler[0] );
    }
    if ( cblk->handler[1] ) {
        starpu_data_unregister( cblk->handler[1] );
    }
    cblk->handler[0] = NULL;
    cblk->handler[1] = NULL;
}

/**
 *******************************************************************************
 *
 * @brief Free the StarPU descriptor of the sparse matrix.
 *
 * This function destroys the StarPU descriptor, but do not free the matrix data
 * that are managed by PaStiX.
 *
 *******************************************************************************
 *
 * @param[inout] spmtx
 *          The descriptor to free.
 *
 ******************************************************************************/
void
starpu_sparse_matrix_destroy( starpu_sparse_matrix_desc_t *spmtx )
{
    starpu_cblk_t *cblkhandle;
    SolverCblk    *cblk;
    pastix_int_t   i, cblkmin2d;
    pastix_int_t   rank = spmtx->solvmtx->clustnum;

    cblkmin2d = spmtx->solvmtx->cblkmin2d;
    cblk      = spmtx->solvmtx->cblktab;
    for ( i = 0; i < cblkmin2d; i++, cblk++ ) {
        pastix_starpu_cblk_destroy( cblk->ownerid == rank, cblk, NULL );
    }

    cblkhandle = spmtx->cblktab_handle;
    for ( i = cblkmin2d; i < spmtx->solvmtx->cblknbr; i++, cblk++, cblkhandle++ ) {
        pastix_starpu_cblk_destroy( cblk->ownerid == rank, cblk, cblkhandle );
    }

    if ( spmtx->cblktab_handle != NULL ) {
        free( spmtx->cblktab_handle );
    }

    pastix_starpu_tag_release( spmtx->mpitag );
}

/**
 * @}
 */
