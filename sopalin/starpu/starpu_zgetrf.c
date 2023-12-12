/**
 *
 * @file starpu_zgetrf.c
 *
 * PaStiX zgetrf StarPU wrapper.
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Alycia Lisito
 * @author Nolan Bredel
 * @author Tom Moenne-Loccoz
 * @date 2023-11-07
 * @precisions normal z -> s d c
 *
 * @addtogroup starpu_getrf
 * @{
 *
 **/
#include "common.h"
#include "blend/solver.h"
#include "sopalin/sopalin_data.h"
#include "pastix_zcores.h"
#include "pastix_starpu.h"
#include "pastix_zstarpu.h"
#include "timing.h"

/**
 *******************************************************************************
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] sopalin_data
 *          TODO
 *
 * @param[in] cblk
 *          TODO
 *
 * @param[in] fcblk
 *          TODO
 *
 * @param[in] blokA
 *          TODO
 *
 * @param[in] prio
 *          TODO
 *
 *******************************************************************************/
void
starpu_task_getrf_zgemmsp( sopalin_data_t   *sopalin_data,
                           SolverCblk       *cblk,
                           const SolverBlok *blokB,
                           SolverCblk       *fcblk,
                           int               prio )
{
    const SolverBlok *blokA, *lblk;
    lblk = cblk[1].fblokptr;

    if ( cblk->cblktype & CBLK_TASKS_2D ) {
        for ( blokA = blokB; blokA < lblk; blokA++ ) {
            if ( blokA == blokB ) {
                starpu_task_blok_zgemmsp( sopalin_data, PastixLCoef, PastixUCoef, PastixTrans,
                            cblk, fcblk, blokB, blokB, prio );

                /* Skip A blocks facing the same cblk */
                while ( ( blokA < lblk ) &&
                        ( blokA[0].fcblknm == blokA[1].fcblknm ) &&
                        ( blokA[0].lcblknm == blokA[1].lcblknm ) )
                {
                    blokA++;
                }
                continue;
            }
            starpu_task_blok_zgemmsp( sopalin_data, PastixLCoef, PastixUCoef, PastixTrans,
                                      cblk, fcblk, blokA, blokB, prio );

            starpu_task_blok_zgemmsp( sopalin_data, PastixUCoef, PastixLCoef, PastixTrans,
                                      cblk, fcblk, blokA, blokB, prio );

            /* Skip A blocks facing the same cblk */
            while ( ( blokA < lblk ) &&
                    ( blokA[0].fcblknm == blokA[1].fcblknm ) &&
                    ( blokA[0].lcblknm == blokA[1].lcblknm ) )
            {
                blokA++;
            }
        }
    }
    else {
        /* Update on L */
        starpu_task_cblk_zgemmsp( sopalin_data, PastixLCoef, PastixUCoef, PastixTrans,
                                    cblk, blokB, fcblk, prio );

        /* Update on U */
        if ( blokB+1 < lblk ) {
            starpu_task_cblk_zgemmsp( sopalin_data, PastixUCoef, PastixLCoef, PastixTrans,
                                        cblk, blokB, fcblk, prio );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief TODO
 *
 *******************************************************************************
 *
 * @param[in] sopalin_data
 *          TODO
 *
 * @param[in] cblk
 *          TODO
 *
 * @param[in] prio
 *          TODO
 *
 *******************************************************************************/
void
starpu_task_zgetrfsp( sopalin_data_t *sopalin_data,
                      SolverCblk     *cblk,
                      int             prio )
{
    SolverBlok   *lblk, *blok;
    pastix_int_t  m;

    if ( cblk->cblktype & CBLK_TASKS_2D ) {
        starpu_task_blok_zgetrf( sopalin_data, cblk, prio );

        lblk = cblk[1].fblokptr;
        for ( blok = cblk->fblokptr + 1, m = 0; blok < lblk; blok++, m++ ) {

            starpu_task_blok_ztrsmsp( sopalin_data, PastixLCoef, PastixRight, PastixUpper,
                                      PastixNoTrans, PastixNonUnit,
                                      cblk, blok, prio );

            starpu_task_blok_ztrsmsp( sopalin_data, PastixUCoef, PastixRight, PastixUpper,
                                      PastixNoTrans, PastixUnit,
                                      cblk, blok, prio );

            /* Skip blocks facing the same cblk */
            while ( ( blok < lblk ) &&
                    ( blok[0].fcblknm == blok[1].fcblknm ) &&
                    ( blok[0].lcblknm == blok[1].lcblknm ) )
            {
                blok++;
            }
        }
    }
    else {
        starpu_task_cblk_zgetrfsp( sopalin_data, cblk, prio );
    }
}

/**
 *******************************************************************************
 *
 * @brief Perform a sparse LU factorization with 1D kernels.
 *
 * The function performs the LU factorization of a sparse general matrix A.
 * The factorization has the form
 *
 *    \f[ A = L\times U \f]
 *
 * where L is a sparse lower triangular matrix, and U a sparse upper triangular
 * with the same pattern as L^t.
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
starpu_zgetrf_sp1dplus_rl( sopalin_data_t              *sopalin_data,
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

        starpu_task_cblk_zgetrfsp( sopalin_data, cblk,
                                   cblknbr - k );

        blok = cblk->fblokptr + 1; /* this diagonal block */
        lblk = cblk[1].fblokptr;   /* the next diagonal block */

        /* if there are off-diagonal supernodes in the column */
        for(m=0; blok < lblk; blok++, m++ )
        {
            fcblk = (solvmtx->cblktab + blok->fcblknm);
            cblk_n = fcblk - solvmtx->cblktab;

            /* Update on L */
            starpu_task_cblk_zgemmsp( sopalin_data, PastixLCoef, PastixUCoef, PastixTrans,
                                      cblk, blok, fcblk,
                                      cblknbr - pastix_imin( k + m, cblk_n ) );

            /* Update on U */
            if ( blok+1 < lblk ) {
                starpu_task_cblk_zgemmsp( sopalin_data, PastixUCoef, PastixLCoef, PastixTrans,
                                          cblk, blok, fcblk,
                                          cblknbr - pastix_imin( k + m, cblk_n ) );
            }
        }
        starpu_sparse_cblk_wont_use( PastixLUCoef, cblk );
    }
    (void)desc;
}

/**
 *******************************************************************************
 *
 * @brief Perform a sparse LU factorization with 1D kernels.
 *
 * The function performs the LU factorization of a sparse general matrix A.
 * The factorization has the form
 *
 *    \f[ A = L\times U \f]
 *
 * where L is a sparse lower triangular matrix, and U a sparse upper triangular
 * with the same pattern as L^t.
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
starpu_zgetrf_sp1dplus_ll( sopalin_data_t              *sopalin_data,
                           starpu_sparse_matrix_desc_t *desc )
{
    const SolverMatrix *solvmtx = sopalin_data->solvmtx;
    SolverCblk         *cblk, *fcblk, *lcblk;
    SolverBlok         *blok, *lblk;
    pastix_int_t        k, m, cblknbr, cblk_n;

    cblknbr = solvmtx->cblknbr;
    cblk = solvmtx->cblktab;
    for (k=0; k<solvmtx->cblknbr; k++, cblk++){

        for ( m = cblk[0].brownum; m < cblk[1].brownum; m++ ) {
            blok  = solvmtx->bloktab + solvmtx->browtab[m];
            lcblk = solvmtx->cblktab + blok->lcblknm;

            if ( lcblk->cblktype & CBLK_IN_SCHUR ) {
                break;
            }

            fcblk  = solvmtx->cblktab + blok->fcblknm;
            cblk_n = fcblk - solvmtx->cblktab;

            assert( fcblk == cblk );

            /* Update on L */
            starpu_task_cblk_zgemmsp( sopalin_data, PastixLCoef, PastixUCoef, PastixTrans,
                                      lcblk, blok, cblk,
                                      cblknbr - pastix_imin( k + m, cblk_n ) );

            lblk = fcblk[1].fblokptr;

            /* Update on U */
            if ( blok+1 < lblk ) {
                starpu_task_cblk_zgemmsp( sopalin_data, PastixUCoef, PastixLCoef, PastixTrans,
                                          lcblk, blok, cblk,
                                          cblknbr - pastix_imin( k + m, cblk_n ) );
            }
        }

        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            continue;
        }

        starpu_task_cblk_zgetrfsp( sopalin_data, cblk,
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
 * @brief Perform a sparse LU factorization with 1D and 2D kernels.
 *
 * The function performs the LU factorization of a sparse general matrix A.
 * The factorization has the form
 *
 *    \f[ A = L\times U \f]
 *
 * where L is a sparse lower triangular matrix, and U a sparse upper triangular
 * with the same pattern as L^t.
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
starpu_zgetrf_sp2d_rl( sopalin_data_t              *sopalin_data,
                       starpu_sparse_matrix_desc_t *desc )
{
    const SolverMatrix *solvmtx = sopalin_data->solvmtx;
    SolverCblk         *cblk, *fcblk;
    SolverBlok         *blok, *lblk;
    pastix_int_t        k, m, cblknbr, cblk_n;

    cblknbr = solvmtx->cblknbr;

    /* Let's submit all 1D tasks first */
    cblk = solvmtx->cblktab;
    for ( k = 0; k <= solvmtx->cblkmax1d; k++, cblk++ ) {

        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            break;
        }

        if ( cblk->cblktype & CBLK_TASKS_2D ) {
            continue;
        }

        starpu_task_zgetrfsp( sopalin_data, cblk, cblknbr - k );

        blok = cblk->fblokptr + 1; /* this diagonal block     */
        lblk = cblk[1].fblokptr;   /* the next diagonal block */

        /* if there are off-diagonal supernodes in the column */
        for(m=0; blok < lblk; blok++, m++ )
        {
            fcblk = (solvmtx->cblktab + blok->fcblknm);
            cblk_n = fcblk - solvmtx->cblktab;

            starpu_task_getrf_zgemmsp( sopalin_data, cblk, blok, fcblk,
                                       cblknbr - pastix_imin( k + m, cblk_n ) );
        }
        starpu_sparse_cblk_wont_use( PastixLUCoef, cblk );
    }

    /* Now we submit all 2D tasks */
    cblk = solvmtx->cblktab + solvmtx->cblkmin2d;
    for ( k = solvmtx->cblkmin2d; k < solvmtx->cblknbr; k++, cblk++ ) {

        if ( cblk->cblktype & CBLK_IN_SCHUR ) {
            continue;
        }

        if ( ! ( cblk->cblktype & CBLK_TASKS_2D ) ) {
            continue; /* skip 1D cblk */
        }

        starpu_task_zgetrfsp( sopalin_data, cblk, cblknbr - k );

        lblk = cblk[1].fblokptr;

        for ( blok = cblk->fblokptr + 1, m = 0; blok < lblk; blok++, m++ ) {
            fcblk  = ( solvmtx->cblktab + blok->fcblknm );
            cblk_n = ( cblk->cblktype & CBLK_TASKS_2D ) ? blok->fcblknm : fcblk - solvmtx->cblktab;

            starpu_task_getrf_zgemmsp( sopalin_data, cblk, blok, fcblk,
                                       cblknbr - pastix_imin( k + m, cblk_n ) );

            /* Skip blocks facing the same cblk */
            while( ( blok < lblk ) &&
                    ( blok[0].fcblknm == blok[1].fcblknm ) &&
                    ( blok[0].lcblknm == blok[1].lcblknm ) )
            {
                blok++;
            }
        }
        starpu_sparse_cblk_wont_use( PastixLUCoef, cblk );
    }

    (void)desc;
}

/**
 *******************************************************************************
 *
 * @brief Perform a sparse LU factorization with 1D and 2D kernels.
 *
 * The function performs the LU factorization of a sparse general matrix A.
 * The factorization has the form
 *
 *    \f[ A = L\times U \f]
 *
 * where L is a sparse lower triangular matrix, and U a sparse upper triangular
 * with the same pattern as L^t.
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
starpu_zgetrf_sp2d_ll( sopalin_data_t              *sopalin_data,
                       starpu_sparse_matrix_desc_t *desc )
{
    const SolverMatrix *solvmtx = sopalin_data->solvmtx;
    SolverCblk         *cblk, *fcblk;
    SolverBlok         *blok = NULL;
    SolverBlok         *blok_prev;
    pastix_int_t        k, m, cblknbr;

    cblknbr = solvmtx->cblknbr;
    fcblk   = solvmtx->cblktab;

    for ( k = 0; k < cblknbr; k++, fcblk++ ) {

        for ( m = fcblk[0].brownum; m < fcblk[1].brownum; m++ ) {
            blok_prev = blok;
            blok = solvmtx->bloktab + solvmtx->browtab[m];
            cblk = solvmtx->cblktab + blok->lcblknm;

            if ( cblk->cblktype & CBLK_IN_SCHUR ) {
                continue;
            }

            if( ( cblk->cblktype & CBLK_TASKS_2D )      &&
                ( blok_prev          != NULL          ) &&
                ( blok_prev->fcblknm == blok->fcblknm ) &&
                ( blok_prev->lcblknm == blok->lcblknm ) )
            {
                continue;
            }

            starpu_task_getrf_zgemmsp( sopalin_data, cblk, blok, fcblk,
                                       cblknbr - k );
        }

        if ( fcblk->cblktype & CBLK_IN_SCHUR ) {
            continue;
        }

        starpu_task_zgetrfsp( sopalin_data, fcblk, cblknbr - k );
    }

    cblk = solvmtx->cblktab;
    for ( k = 0; k < solvmtx->cblknbr; k++, cblk++ ) {
        starpu_sparse_cblk_wont_use( PastixLUCoef, cblk );
    }

    (void)desc;
}

/**
 *******************************************************************************
 *
 * @brief Perform a sparse LU factorization using StarPU runtime.
 *
 * The function performs the LU factorization of a sparse general matrix A.
 * The factorization has the form
 *
 *    \f[ A = L\times U \f]
 *
 * where L is a sparse lower triangular matrix, and U a sparse upper triangular
 * with the same pattern as L^t.
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
starpu_zgetrf( pastix_data_t  *pastix_data,
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
                                   PastixGeneral,
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
        if ( pastix_data->iparm[IPARM_FACTO_LOOK_SIDE] == PastixFactLeftLooking ) {
            starpu_zgetrf_sp2d_ll( sopalin_data, sdesc );
        }
        else {
            starpu_zgetrf_sp2d_rl( sopalin_data, sdesc );
        }
    }
    else
    {
        if ( pastix_data->iparm[IPARM_FACTO_LOOK_SIDE] == PastixFactLeftLooking ) {
            starpu_zgetrf_sp1dplus_ll( sopalin_data, sdesc );
        }
        else {
            starpu_zgetrf_sp1dplus_rl( sopalin_data, sdesc );
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
