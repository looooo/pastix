/**
 *
 * @file starpu_zpxtrf.c
 *
 * PaStiX zpxtrf StarPU wrapper.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.1
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2023-07-21
 * @precisions normal z -> z c
 *
 * @addtogroup starpu_pxtrf
 * @{
 *
 **/
#include "common.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_starpu.h"
#include "pastix_zstarpu.h"

/**
 *******************************************************************************
 *
 * @brief Perform a sparse LL^t factorization with 1D kernels.
 *
 * The function performs the LL^t factorization of a sparse symmetric complex
 * matrix A.
 * The factorization has the form
 *
 *    \f[ A = L\times L^t \f]
 *
 * where L is a sparse lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 * @param[inout] desc
 *          StarPU descriptor of the sparse matrix.
 *
 ******************************************************************************/
void
starpu_zpxtrf_sp1dplus_rl( sopalin_data_t              *sopalin_data,
                           starpu_sparse_matrix_desc_t *desc )
{
    const SolverMatrix *solvmtx = sopalin_data->solvmtx;
    SolverCblk         *cblk, *fcblk;
    SolverBlok         *blok, *lblk;
    pastix_int_t        k, m, cblknbr, cblk_n;

    cblknbr = solvmtx->cblknbr;
    cblk = solvmtx->cblktab;
    for (k=0; k<solvmtx->cblknbr; k++, cblk++){

        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            break;
        }

        starpu_task_cblk_zpxtrfsp( sopalin_data, cblk,
                                   cblknbr - k );

        blok = cblk->fblokptr + 1; /* this diagonal block */
        lblk = cblk[1].fblokptr;   /* the next diagonal block */

        /* if there are off-diagonal supernodes in the column */
        for(m=0; blok < lblk; blok++, m++ )
        {
            fcblk = (solvmtx->cblktab + blok->fcblknm);
            cblk_n = fcblk - solvmtx->cblktab;

            starpu_task_cblk_zgemmsp( sopalin_data, PastixLCoef, PastixLCoef, PastixConjTrans,
                                      cblk, blok, fcblk,
                                      cblknbr - pastix_imin( k + m, cblk_n ) );
        }
        starpu_sparse_cblk_wont_use( PastixLCoef, cblk );
    }
    (void)desc;
}

/**
 *******************************************************************************
 *
 * @brief Perform a sparse LL^t factorization with 1D kernels.
 *
 * The function performs the LL^t factorization of a sparse symmetric complex
 * matrix A.
 * The factorization has the form
 *
 *    \f[ A = L\times L^t \f]
 *
 * where L is a sparse lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 * @param[inout] desc
 *          StarPU descriptor of the sparse matrix.
 *
 ******************************************************************************/
void
starpu_zpxtrf_sp1dplus_ll( sopalin_data_t              *sopalin_data,
                           starpu_sparse_matrix_desc_t *desc )
{
    const SolverMatrix *solvmtx = sopalin_data->solvmtx;
    SolverCblk         *cblk, *fcblk, *lcblk;
    SolverBlok         *blok;
    pastix_int_t        k, m, cblknbr, cblk_n;

    cblknbr = solvmtx->cblknbr;
    cblk = solvmtx->cblktab;
    for ( k = 0; k < solvmtx->cblknbr; k++, cblk++ ) {

        for ( m = cblk[0].brownum; m < cblk[1].brownum; m++ ) {
            blok  = solvmtx->bloktab + solvmtx->browtab[m];
            lcblk = solvmtx->cblktab + blok->lcblknm;

            if ( lcblk->cblktype & CBLK_IN_SCHUR ) {
                break;
            }

            fcblk  = solvmtx->cblktab + blok->fcblknm;
            cblk_n = fcblk - solvmtx->cblktab;

            assert( fcblk == cblk );

            starpu_task_cblk_zgemmsp( sopalin_data, PastixLCoef, PastixLCoef, PastixConjTrans,
                                      lcblk, blok, cblk,
                                      cblknbr - pastix_imin( k + m, cblk_n ) );
        }

        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            continue;
        }

        starpu_task_cblk_zpxtrfsp( sopalin_data, cblk,
                                   cblknbr - k );
    }

    cblk = solvmtx->cblktab;
    for ( k = 0; k < solvmtx->cblknbr; k++, cblk++ ) {
        starpu_sparse_cblk_wont_use( PastixLCoef, cblk );
    }
    (void)desc;
}

/**
 *******************************************************************************
 *
 * @brief Perform a sparse LL^t factorization with 1D and 2D kernels.
 *
 * The function performs the LL^t factorization of a sparse symmetric complex
 * matrix A.
 * The factorization has the form
 *
 *    \f[ A = L\times L^t \f]
 *
 * where L is a sparse lower triangular matrix.
 *
 *******************************************************************************
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 * @param[inout] desc
 *          StarPU descriptor of the sparse matrix.
 *
 ******************************************************************************/
void
starpu_zpxtrf_sp2d( sopalin_data_t              *sopalin_data,
                    starpu_sparse_matrix_desc_t *desc )
{
    const SolverMatrix *solvmtx = sopalin_data->solvmtx;
    SolverCblk         *cblk, *fcblk;
    SolverBlok         *blok, *lblk, *blokA, *blokB;
    starpu_cblk_t      *cblkhandle;
    pastix_int_t        k, m, cblknbr, cblk_n;

    cblknbr = solvmtx->cblknbr;

    /* Let's submit all 1D tasks first */
    cblk = solvmtx->cblktab;
    for (k=0; k<=solvmtx->cblkmax1d; k++, cblk++){

        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            break;
        }

        if ( cblk->cblktype & CBLK_TASKS_2D ) {
            continue;
        }

        starpu_task_cblk_zpxtrfsp( sopalin_data, cblk,
                                   cblknbr - k );

        blok = cblk->fblokptr + 1; /* this diagonal block     */
        lblk = cblk[1].fblokptr;   /* the next diagonal block */

        /* if there are off-diagonal supernodes in the column */
        for(m=0; blok < lblk; blok++, m++ )
        {
            fcblk = (solvmtx->cblktab + blok->fcblknm);
            cblk_n = fcblk - solvmtx->cblktab;

            starpu_task_cblk_zgemmsp( sopalin_data, PastixLCoef, PastixLCoef, PastixConjTrans,
                                      cblk, blok, fcblk,
                                      cblknbr - pastix_imin( k + m, cblk_n ) );
        }
        starpu_sparse_cblk_wont_use( PastixLCoef, cblk );
    }

    /* Now we submit all 2D tasks */
    cblk       = solvmtx->cblktab + solvmtx->cblkmin2d;
    cblkhandle = desc->cblktab_handle;
    for (k=solvmtx->cblkmin2d; k<solvmtx->cblknbr; k++, cblk++, cblkhandle++){

        if ( !(cblk->cblktype & CBLK_TASKS_2D) ) {
            continue; /* skip 1D cblk */
        }

        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            continue;
        }

        starpu_task_blok_zpxtrf( sopalin_data, cblk,
                                 cblknbr - k );

        lblk = cblk[1].fblokptr;
        for(blokA=cblk->fblokptr + 1, m=0; blokA<lblk; blokA++, m++) {

            cblk_n = blokA->fcblknm;

            starpu_task_blok_ztrsmsp( sopalin_data, PastixLCoef, PastixRight, PastixLower,
                                      PastixConjTrans, PastixNonUnit,
                                      cblk, blokA,
                                      cblknbr - k );

            for(blokB=cblk->fblokptr + 1; blokB<=blokA; blokB++) {

                starpu_task_blok_zgemmsp( sopalin_data, PastixLCoef, PastixLCoef, PastixConjTrans,
                                          cblk, solvmtx->cblktab + blokB->fcblknm,
                                          blokA, blokB,
                                          cblknbr - pastix_imin( k + m, cblk_n ) );

                /* Skip B blocks facing the same cblk */
                while( (blokB < blokA) &&
                       (blokB[0].fcblknm == blokB[1].fcblknm) &&
                       (blokB[0].lcblknm == blokB[1].lcblknm) )
                {
                    blokB++;
                }
            }

            /* Skip A blocks facing the same cblk */
            while( (blokA < lblk) &&
                   (blokA[0].fcblknm == blokA[1].fcblknm) &&
                   (blokA[0].lcblknm == blokA[1].lcblknm) )
            {
                blokA++;
            }
        }
        starpu_sparse_cblk_wont_use( PastixLCoef, cblk );
    }
    (void)desc;
}

/**
 *******************************************************************************
 *
 * @brief Perform a sparse LL^t factorization using StarPU runtime.
 *
 * The function performs the LL^t factorization of a sparse symmetric complex
 * matrix A.
 * The factorization has the form
 *
 *    \f[ A = L\times L^t \f]
 *
 * where L is a sparse lower triangular matrix.
 *
 * The algorithm is automatically chosen between the 1D and 2D version based on
 * the API parameter IPARM_TASKS2D_LEVEL. If IPARM_TASKS2D_LEVEL != 0
 * the 2D scheme is applied, the 1D otherwise.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *
 * @param[inout] sopalin_data
 *          Solver matrix information structure that will guide the algorithm.
 *
 ******************************************************************************/
void
starpu_zpxtrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data )
{
    starpu_sparse_matrix_desc_t *sdesc = sopalin_data->solvmtx->starpu_desc;
    double sub = 0.;
    double com = 0.;

    /*
     * Start StarPU if not already started
     */
    if (pastix_data->starpu == NULL) {
        int argc = 0;
        pastix_starpu_init( pastix_data, &argc, NULL, NULL );
    }

    if ( sdesc == NULL ) {
        /* Create the matrix descriptor */
        starpu_sparse_matrix_init( sopalin_data->solvmtx,
                                   PastixSymmetric,
                                   pastix_data->inter_node_procnbr,
                                   pastix_data->inter_node_procnum,
                                   PastixComplex64 );
        sdesc = sopalin_data->solvmtx->starpu_desc;
    }

    starpu_profiling_status_set(STARPU_PROFILING_ENABLE);
#if defined(STARPU_USE_FXT)
    if (pastix_data->iparm[IPARM_TRACE] & PastixTraceNumfact) {
        starpu_fxt_start_profiling();
    }
#endif
#if defined(PASTIX_STARPU_STATS)
    clockStart( sub );
#else
    starpu_resume();
#endif
    /*
     * Select 1D or 2D algorithm based on 2d tasks level
     */
    if ( pastix_data->iparm[IPARM_TASKS2D_LEVEL] != 0 )
    {
        starpu_zpxtrf_sp2d( sopalin_data, sdesc );
    }
    else
    {
        if ( pastix_data->iparm[IPARM_FACTO_LOOK_SIDE] == PastixFactLeftLooking ) {
            starpu_zpxtrf_sp1dplus_ll( sopalin_data, sdesc );
        }
        else {
            starpu_zpxtrf_sp1dplus_rl( sopalin_data, sdesc );
        }
    }

    starpu_sparse_matrix_getoncpu( sdesc );
#if defined(PASTIX_STARPU_STATS)
    clockStop( sub );
    clockStart( com );
    starpu_resume();
#endif
    starpu_task_wait_for_all();
#if defined(PASTIX_WITH_MPI)
    starpu_mpi_wait_for_all( pastix_data->pastix_comm );
    starpu_mpi_barrier( pastix_data->inter_node_comm );
#endif
    starpu_pause();
#if defined(STARPU_USE_FXT)
    if (pastix_data->iparm[IPARM_TRACE] & PastixTraceNumfact) {
        starpu_fxt_stop_profiling();
    }
#endif
    starpu_profiling_status_set(STARPU_PROFILING_DISABLE);
#if defined(PASTIX_STARPU_STATS)
    clockStop( com );
    print_stats( sub, com, pastix_data->solvmatr );
#endif

    (void)com;
    (void)sub;
    return;
}

/**
 *@}
 */
