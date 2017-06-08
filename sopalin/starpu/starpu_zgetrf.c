/**
 *
 * @file starpu_zgetrf.c
 *
 * PaStiX zgetrf StarPU wrapper.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2017-06-24
 * @precisions normal z -> s d c
 *
 * @addtogroup starpu_getrf
 * @{
 *
 **/
#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_zstarpu.h"

void
starpu_zgetrf_sp1dplus( sopalin_data_t              *sopalin_data,
                        starpu_sparse_matrix_desc_t *desc )
{
    const SolverMatrix *solvmtx = sopalin_data->solvmtx;
    SolverCblk         *cblk, *fcblk;
    SolverBlok         *blok, *lblk;
    pastix_int_t  i;

    cblk = solvmtx->cblktab;
    for (i=0; i<solvmtx->cblknbr; i++, cblk++){

        if ( cblk->cblktype & CBLK_IN_SCHUR )
            break;

        starpu_task_cblk_zgetrfsp1d_panel( sopalin_data, cblk );

        blok = cblk->fblokptr + 1; /* this diagonal block */
        lblk = cblk[1].fblokptr;   /* the next diagonal block */

        /* if there are off-diagonal supernodes in the column */
        for( ; blok < lblk; blok++ )
        {
            fcblk = (solvmtx->cblktab + blok->fcblknm);

            /* Update on L */
            starpu_task_cblk_zgemmsp( PastixLCoef, PastixUCoef, PastixTrans,
                                      cblk, blok, fcblk, sopalin_data );

            /* Update on U */
            if ( blok+1 < lblk ) {
                starpu_task_cblk_zgemmsp( PastixUCoef, PastixLCoef, PastixTrans,
                                          cblk, blok, fcblk, sopalin_data );
            }
        }
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( datacode, "getrf_L.txt" );
#endif
    (void)desc;
}

void
starpu_zgetrf_sp2d( sopalin_data_t              *sopalin_data,
                    starpu_sparse_matrix_desc_t *desc )
{
    const SolverMatrix *solvmtx = sopalin_data->solvmtx;
    SolverCblk         *cblk, *fcblk;
    SolverBlok         *blok, *lblk;
    pastix_int_t  i;

    cblk = solvmtx->cblktab;
    for (i=0; i<solvmtx->cblknbr; i++, cblk++){

        if ( cblk->cblktype & CBLK_IN_SCHUR )
            break;

        starpu_task_cblk_zgetrfsp1d_panel( sopalin_data, cblk );

        blok = cblk->fblokptr + 1; /* this diagonal block */
        lblk = cblk[1].fblokptr;   /* the next diagonal block */

        /* if there are off-diagonal supernodes in the column */
        for( ; blok < lblk; blok++ )
        {
            fcblk = (solvmtx->cblktab + blok->fcblknm);

            /* Update on L */
            starpu_task_cblk_zgemmsp( PastixLCoef, PastixUCoef, PastixTrans,
                                      cblk, blok, fcblk, sopalin_data );

            /* Update on U */
            if ( blok+1 < lblk ) {
                starpu_task_cblk_zgemmsp( PastixUCoef, PastixLCoef, PastixTrans,
                                          cblk, blok, fcblk, sopalin_data );
            }
        }
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( datacode, "getrf_L.txt" );
#endif
    (void)desc;
}

void
starpu_zgetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    starpu_sparse_matrix_desc_t *sdesc = sopalin_data->solvmtx->starpu_desc;

    /*
     * Start StarPU if not already started
     */
    if (pastix_data->starpu == NULL) {
        int argc = 0;
        pastix_starpu_init( pastix_data, &argc, NULL );
    }

    if ( sdesc == NULL ) {
        /* Create the matrix descriptor */
        starpu_sparse_matrix_init( sopalin_data->solvmtx,
                                   sizeof( pastix_complex64_t ), PastixGeneral,
                                   1, 0 );
        sdesc = sopalin_data->solvmtx->starpu_desc;
    }

    /*
     * Select 1D or 2D jdf based on distribution_level
     */
    if ( pastix_data->iparm[IPARM_DISTRIBUTION_LEVEL] >= 0 )
    {
        starpu_zgetrf_sp2d( sopalin_data, sdesc );
    }
    else
    {
        starpu_zgetrf_sp1dplus( sopalin_data, sdesc );
    }

    starpu_sparse_matrix_getoncpu( sdesc );
    starpu_task_wait_for_all();
#if defined(PASTIX_WITH_MPI)
    starpu_mpi_barrier(MPI_COMM_WORLD);
#endif

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( sopalin_data->solvmtx, "getrf.txt" );
#endif

    return;
}

/**
 *@}
 */
