/**
 *
 * @file starpu_rhs.c
 *
 * PaStiX dense matrix descriptor for StarPU.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2023-07-21
 *
 * @ingroup pastix_starpu
 * @{
 *
 **/
#include "common.h"
#include "blend/solver.h"
#include "pastix_starpu.h"
#include <starpu_data.h>

/**
 *******************************************************************************
 *
 * @brief Generate the StarPU descriptor of the dense matrix.
 *
 * This function creates the StarPU descriptor that will provide tha data
 * mapping and memory location to StarPU for the computation.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver matrix structure that describes the dense matrix for
 *          PaStiX.
 *
 * @param[inout] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 * @param[in] typesize
 *          The memory size of the arithmetic used to store the matrix
 *          coefficients.
 *
 * @param[in] nodes
 *          The number of processes used to solve the problem.
 *
 * @param[in] myrank
 *          The rank of the calling process.
 *
 ******************************************************************************/
void
starpu_rhs_init( SolverMatrix *solvmtx,
                 pastix_rhs_t  rhsb,
                 int           typesize,
                 int           nodes,
                 int           myrank )
{
    starpu_data_handle_t *handler;
    SolverCblk           *cblk;
    pastix_int_t          cblknbr, cblknum, nrow;
    pastix_int_t          ncol = rhsb->n;

    starpu_rhs_desc_t *rhsdesc = rhsb->starpu_desc;
    if ( rhsdesc != NULL ) {
        if ( ( ncol    == rhsdesc->ncol ) &&
             ( rhsb->b == rhsdesc->dataptr ) ) {
            return;
        }
        starpu_rhs_destroy( rhsdesc );
    }
    else {
        rhsdesc = (starpu_rhs_desc_t*)malloc(sizeof(starpu_rhs_desc_t));
    }

    cblknbr = solvmtx->cblknbr;

    rhsdesc->ncol      = ncol;
    rhsdesc->typesze   = pastix_size_of( typesize );
    rhsdesc->solvmtx   = solvmtx;
    rhsdesc->handletab = malloc( cblknbr * sizeof(starpu_data_handle_t) );
    rhsdesc->dataptr   = rhsb->b;

    /* Initialize 1D cblk handlers */
    cblk    = rhsdesc->solvmtx->cblktab;
    handler = rhsdesc->handletab;
    for( cblknum = 0;
         cblknum < cblknbr;
         cblknum++, cblk++, handler++ )
    {
        nrow = cblk_colnbr( cblk );

        if( cblk->ownerid == myrank ) {
            starpu_matrix_data_register( handler, STARPU_MAIN_RAM,
                                         (uintptr_t)(rhsb->b + (cblk->lcolidx * rhsdesc->typesze)),
                                         rhsb->ld, nrow, ncol, rhsdesc->typesze );
        }
        else {
            starpu_matrix_data_register( handler, -1, 0,
                                         rhsb->ld, nrow, ncol, rhsdesc->typesze );
        }
#if defined(PASTIX_WITH_MPI)
        rhsdesc->mpitag = pastix_starpu_tag_book( solvmtx->gcblknbr );
        starpu_mpi_data_register( *handler, rhsdesc->mpitag + ((int64_t)cblknum), cblk->ownerid );
#endif
    }

    rhsb->starpu_desc = rhsdesc;

    (void)nodes;
    (void)myrank;
}

/**
 *******************************************************************************
 *
 * @brief Submit asynchronous calls to retrieve the data on main memory.
 *
 *******************************************************************************
 *
 * @param[inout] rhsdesc
 *          The dense matrix descriptor to retrieve on main memory.
 *
 ******************************************************************************/
void
starpu_rhs_getoncpu( starpu_rhs_desc_t *rhsdesc )
{
    starpu_data_handle_t *handler = rhsdesc->handletab;
    SolverCblk *cblk;
    pastix_int_t cblknbr, cblknum;

    cblk    = rhsdesc->solvmtx->cblktab;
    cblknbr = rhsdesc->solvmtx->cblknbr;
    for(cblknum=0; cblknum<cblknbr; cblknum++, cblk++, handler++)
    {
        assert( handler );

#if defined(PASTIX_WITH_MPI)
        starpu_mpi_cache_flush( rhsdesc->solvmtx->solv_comm, *handler );
#endif
        if ( cblk->ownerid == rhsdesc->solvmtx->clustnum ) {
            starpu_data_acquire_cb( *handler, STARPU_R,
                                    (void (*)(void*))&starpu_data_release,
                                    *handler );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Free the StarPU descriptor of the dense matrix.
 *
 * This function destroys the StarPU descriptor, but do not free the matrix data
 * that are managed by PaStiX.
 *
 *******************************************************************************
 *
 * @param[inout] rhsdesc
 *          The descriptor to free.
 *
 ******************************************************************************/
void
starpu_rhs_destroy( starpu_rhs_desc_t *rhsdesc )
{
    starpu_data_handle_t *handler = rhsdesc->handletab;
    SolverCblk *cblk;
    pastix_int_t cblknbr, cblknum;

    cblk    = rhsdesc->solvmtx->cblktab;
    cblknbr = rhsdesc->solvmtx->cblknbr;
    for( cblknum = 0; cblknum < cblknbr; cblknum++, cblk++, handler++ ) {
        assert( handler );
        starpu_data_unregister( *handler );
    }

    free( rhsdesc->handletab );
    rhsdesc->handletab = NULL;

    pastix_starpu_tag_release( rhsdesc->mpitag );
}

/**
 * @}
 */
