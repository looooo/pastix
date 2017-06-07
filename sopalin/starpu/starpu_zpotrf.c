/**
 *
 * @file starpu_zpotrf.c
 *
 *  PaStiX factorization routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"

void
starpu_task_cblk_zpotrfsp1d_panel( sopalin_data_t *sopalin_data,
                                   SolverCblk     *cblk );

void
starpu_task_cblk_zgemmsp( pastix_coefside_t sideA,
                          pastix_coefside_t sideB,
                          pastix_trans_t    trans,
                          const SolverCblk *cblk,
                          const SolverBlok *blok,
                          SolverCblk       *fcblk,
                          sopalin_data_t   *sopalin_data );


void
starpu_zpotrf_sp1dplus( sopalin_data_t              *sopalin_data,
                        starpu_sparse_matrix_desc_t *desc )
{
    const SolverMatrix *solvmtx = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    SolverBlok         *blok;
    pastix_int_t  i;

    cblk = solvmtx->cblktab;
    for (i=0; i<solvmtx->cblknbr; i++, cblk++){

        if ( cblk->cblktype & CBLK_IN_SCHUR )
            break;

        starpu_task_cblk_zpotrfsp1d_panel( sopalin_data, cblk );

        for(blok=cblk->fblokptr + 1; blok<cblk[1].fblokptr; blok++) {

            starpu_task_cblk_zgemmsp( PastixLCoef, PastixLCoef, PastixConjTrans,
                                      cblk, blok,
                                      solvmtx->cblktab + blok->fcblknm, sopalin_data );

        }
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( datacode, "potrf_L.txt" );
#endif
    (void)desc;
}

void
starpu_zpotrf_sp2d( sopalin_data_t              *sopalin_data,
                    starpu_sparse_matrix_desc_t *desc )
{
    const SolverMatrix *solvmtx = sopalin_data->solvmtx;
    SolverCblk         *cblk;
    SolverBlok         *blok;
    pastix_int_t  i;

    cblk = solvmtx->cblktab;
    for (i=0; i<solvmtx->cblknbr; i++, cblk++){

        if ( cblk->cblktype & CBLK_IN_SCHUR )
            break;

        starpu_task_cblk_zpotrfsp1d_panel( sopalin_data, cblk );

        for(blok=cblk->fblokptr + 1; blok<cblk[1].fblokptr; blok++) {

            starpu_task_cblk_zgemmsp( PastixLCoef, PastixLCoef, PastixConjTrans,
                                      cblk, blok,
                                      solvmtx->cblktab + blok->fcblknm, sopalin_data );

        }
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( datacode, "potrf_L.txt" );
#endif

    (void)desc;
}

void
starpu_zpotrf( pastix_data_t  *pastix_data,
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
        sdesc = (starpu_sparse_matrix_desc_t*)malloc(sizeof(starpu_sparse_matrix_desc_t));

        /* Create the matrix descriptor */
        starpu_sparse_matrix_init( sdesc, sopalin_data->solvmtx,
                                   sizeof( pastix_complex64_t ), PastixGeneral,
                                   1, 0 );
        sopalin_data->solvmtx->starpu_desc = sdesc;
    }

    /*
     * Select 1D or 2D jdf based on distribution_level
     */
    if ( pastix_data->iparm[IPARM_DISTRIBUTION_LEVEL] >= 0 )
    {
        starpu_zpotrf_sp2d( sopalin_data, sdesc );
    }
    else
    {
        starpu_zpotrf_sp1dplus( sopalin_data, sdesc );
    }

#if defined(PASTIX_DEBUG_FACTO)
    coeftab_zdump( sopalin_data->solvmtx, "potrf.txt" );
#endif

    return;
}
