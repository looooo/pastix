/**
 * @file solver.c
 *
 * PaStiX solver structure basic functions.
 *
 * @copyright 2004-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.4.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Tony Delarue
 * @date 2024-07-05
 *
 **/
#include "common.h"
#include "blend/solver.h"
#include "blend/solver_comm_matrix.h"
#include "sopalin/coeftab.h"

#if defined(PASTIX_WITH_PARSEC)
#include "sopalin/parsec/pastix_parsec.h"
#endif

#if defined(PASTIX_WITH_STARPU)
#include "sopalin/starpu/pastix_starpu.h"
#endif

/**
 *******************************************************************************
 *
 * @ingroup blend_dev_solver_null
 *
 * @brief Compute the memory size used by the solver sturcture itself.
 *
 * This function doesn't count the memory space of the numerical information,
 * but only the sapce of the data structure that describes the matrix.
 *
 *******************************************************************************
 *
 * @param[in] solvptr
 *          The pointer to the solver matrix structure.
 *
 *******************************************************************************
 *
 * @return the memory size in octet of the solver structure.
 *
 *******************************************************************************/
static inline size_t
solver_size( const SolverMatrix *solvptr )
{
    size_t mem = sizeof(SolverMatrix);

    /* cblk and blocks arrays */
    if ( solvptr->cblktab ) {
        mem += solvptr->cblknbr * sizeof( SolverCblk );
    }
    if ( solvptr->bloktab ) {
        mem += solvptr->bloknbr * sizeof( SolverBlok );
    }
    if ( solvptr->browtab ) {
        mem += solvptr->brownbr * sizeof( pastix_int_t );
    }
#if defined(PASTIX_WITH_PARSEC)
    if ( solvptr->parsec_desc ) {
        mem += sizeof( parsec_sparse_matrix_desc_t );
    }
#endif
#if defined(PASTIX_WITH_STARPU)
    if ( solvptr->starpu_desc ) {
        mem += sizeof( starpu_sparse_matrix_desc_t );
    }
#endif

    /* /\* BubbleTree *\/ */
    /* if ( solvptr->btree ) { */
    /*     mem += solvptr->bublnbr * sizeof( BubbleTree ); */
    /*     mem += solvptr->btree->nodemax * sizeof( BubbleTreeNode ); */
    /*     mem += solvptr->btree->nodemax * sizeof( int ); */
    /* } */

    /* Tasks */
    if ( solvptr->tasktab ) {
        mem += solvptr->tasknbr * sizeof(Task);
    }
    if ( solvptr->ttsknbr ) {
        int i;
        mem += solvptr->thrdnbr * sizeof(pastix_int_t);
        mem += solvptr->thrdnbr * sizeof(pastix_int_t*);

        for( i=0; i<solvptr->thrdnbr; i++ ) {
            mem += solvptr->ttsknbr[i] * sizeof(pastix_int_t);
        }
    }

    return mem;
}

/**
 * @addtogroup blend_dev_solver
 * @{
 *
 */

/**
 *******************************************************************************
 *
 * @brief Initialize the solver structure.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The solver structure to initialize.
 *
 *******************************************************************************/
void
solverInit( SolverMatrix *solvmtx )
{
    memset(solvmtx, 0, sizeof (SolverMatrix));
    solvmtx->cblkmax1d  = -1;
    solvmtx->cblkmaxblk = 1;
    return;
}

/**
 *******************************************************************************
 *
 * @brief Free the content of the solver matrix structure.
 *
 * All the arrays from the structure are freed and the structure is memset to 0
 * at exit, but the solver itself is not freed. It will require a new call to
 * solverInit if the memory space area needs to be reused for a new solver
 * matrix.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The pointer to the structure to free.
 *
 *******************************************************************************/
void
solverExit(SolverMatrix *solvmtx)
{
    pastix_int_t i;

    coeftabExit( solvmtx );

    /* Free arrays of solvmtx */
    if(solvmtx->cblktab) {
        memFree_null(solvmtx->cblktab);
    }
    if(solvmtx->bloktab) {
        memFree_null(solvmtx->bloktab);
    }
    if(solvmtx->browtab) {
        memFree_null(solvmtx->browtab);
    }
    if(solvmtx->gcbl2loc) {
        memFree_null(solvmtx->gcbl2loc);
    }
    if(solvmtx->tasktab) {
        memFree_null(solvmtx->tasktab);
    }
    memFree_null(solvmtx->ttsknbr);
    for (i=0;i<solvmtx->bublnbr;i++)
    {
        if (solvmtx->ttsktab[i] != NULL) {
            memFree_null(solvmtx->ttsktab[i]);
        }
    }
    memFree_null(solvmtx->ttsktab);
}

/**
 *******************************************************************************
 *
 * @brief Print statistical information about the solver matrix structure.
 *
 *******************************************************************************
 *
 * @param[in] solvptr
 *          The pointer to the solver matrix structure.
 *
 *******************************************************************************/
void
solverPrintStats( const SolverMatrix *solvptr )
{
    SolverCblk *cblk;
    SolverBlok *blok;
    size_t memstruct, memcoef;
    pastix_int_t itercblk;
    int64_t      cblknbr;

    /* 0: Total, 1: 1d, 2: 2d */
    int64_t fcol[3], lcol[3];
    /* Average width of cblks, and height of off-diagonal blocks */
    int64_t width[3]  = { 0, 0, 0 };
    int64_t height[3] = { 0, 0, 0 };
    /* Number of cblks, and of blocks without diagonal blocks */
    int64_t nbcblk[3] = { 0, 0, 0 };
    int64_t nbblok[3] = { 0, 0, 0 };
    int64_t fblok[3], lblok[3];
    /* Number of blocks with teh regular diagonal partion of the matrix that includes subblocks */
    int64_t nbpartblok[3] = { 0, 0, 0 };

    /* Compute the number of GEMM tasks in the multiple cases */
    int64_t gemm_dense  = 0;
    int64_t gemm_nopart_full2  = 0;
    int64_t gemm_nopart_hybrid = 0;
    int64_t gemm_parsec_full2  = 0;
    int64_t gemm_parsec_hybrid = 0;
    int64_t gemm_starpu_full2  = 0;
    int64_t gemm_starpu_hybrid = 0;
    int64_t gemm_full1  = 0;

    cblknbr = solvptr->cblknbr;
    cblk    = solvptr->cblktab;
    memcoef = 0;
    for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
    {
        pastix_int_t colnbr = cblk->lcolnum - cblk->fcolnum + 1;
        pastix_int_t rownbr = cblk->stride;
        int64_t      bcol_size  = cblk[1].fblokptr - cblk[0].fblokptr;
        int64_t      brow_size[3];
        int64_t      nbpblok = 0;

        brow_size[0] = cblk[1].brownum - cblk[0].brownum;
        brow_size[1] = cblk[0].brown2d - cblk[0].brownum;
        brow_size[2] = cblk[1].brownum - cblk[0].brown2d;
        assert( brow_size[0] == (brow_size[1] + brow_size[2]) );

        memcoef += colnbr * rownbr;

        /* Update counters when no diagonal partitionning is considered and all blocks are considers */
        gemm_nopart_full2  += brow_size[0] * bcol_size;
        gemm_nopart_hybrid += brow_size[1] + (brow_size[2] * bcol_size);

        /* Compute the compressed version of the brow size */
#if !defined(NDEBUG)
        {
            pastix_int_t  b;
            pastix_int_t  lcblk        = -1;
            pastix_int_t *browptr      = solvptr->browtab + cblk[0].brownum;
            int64_t       brow_csze[3] = { 0, 0, 0 };

            for ( b = cblk[0].brownum; b < cblk[1].brownum; b++, browptr++ ) {
                blok = solvptr->bloktab + (*browptr);
                if ( blok->lcblknm != lcblk ) {
                    lcblk = blok->lcblknm;

                    brow_csze[0]++;
                    if ( (solvptr->cblktab + lcblk)->cblktype & CBLK_TASKS_2D ) {
                        brow_csze[2]++;
                    }
                    else {
                        brow_csze[1]++;
                    }
                }
            }
            assert( brow_csze[0] == (brow_csze[1] + brow_csze[2]) );
            assert( brow_csze[0] <= brow_size[0] );
            assert( brow_csze[1] <= brow_size[1] );
            assert( brow_csze[2] <= brow_size[2] );
        }
#endif

        /*
         * Compute the compressed version of the bcol size
         * The algorithm goes backward to avoid a double pass to compute StarPU
         * number of tasks.
         */
        blok = cblk->fblokptr + 1;
        while( blok < cblk[1].fblokptr ) {
            while( (blok < cblk[1].fblokptr-1) &&
                   (blok[0].fcblknm == blok[1].fcblknm) &&
                   (blok[0].lcblknm == blok[1].lcblknm) )
            {
                blok++;
            }
            nbpblok++;
            blok++;
        }

        /* Compute the number of PaRSEC GEMM tasks in a left looking manner */
        gemm_parsec_full2  +=  (nbpblok+1) * brow_size[0];
        gemm_parsec_hybrid += ((nbpblok+1) * brow_size[2]) + brow_size[1];

        /* Compute the number of StarPU GEMM tasks in a right looking manner */
        gemm_starpu_full2 += (nbpblok * (nbpblok+1)) / 2;

        if (cblk->cblktype & CBLK_TASKS_2D) {
            gemm_starpu_hybrid += (nbpblok * (nbpblok+1)) / 2;

            nbpartblok[2] += nbpblok;
            width[2]      += colnbr;
            height[2]     += rownbr - colnbr;
        }
        else {
            gemm_starpu_hybrid += bcol_size - 1;

            nbpartblok[1] += nbpblok;
            width[1]      += colnbr;
            height[1]     += rownbr - colnbr;
        }

        nbpartblok[0] += nbpblok;
        width[0]      += colnbr;
        height[0]     += rownbr - colnbr;
    }

    assert( (width[1] + width[2]) == solvptr->nodenbr );
    assert( (width[1] + width[2]) == width[0] );
    assert( (height[1] + height[2]) == height[0] );

    memstruct = solver_size( solvptr );

    gemm_dense = (cblknbr * ( cblknbr * cblknbr - 1 )) / 6;
    gemm_full1 = solvptr->bloknbr - solvptr->cblknbr;

    fcol[0] = 0;
    fcol[1] = 0;
    fcol[2] = (solvptr->cblktab + solvptr->cblkmin2d)->fcolnum;

    lcol[0] = (solvptr->cblktab + solvptr->cblknbr  )->fcolnum;
    lcol[1] = (solvptr->cblktab + solvptr->cblkmax1d)->lcolnum + 1;
    lcol[2] = (solvptr->cblktab + solvptr->cblknbr  )->fcolnum;

    nbcblk[0] = cblknbr;
    nbcblk[1] = (cblknbr - solvptr->nb2dcblk);
    nbcblk[2] = solvptr->nb2dcblk;

    nbblok[0] = solvptr->bloknbr - cblknbr;
    nbblok[1] = (solvptr->bloknbr - cblknbr) - (solvptr->nb2dblok - solvptr->nb2dcblk);
    nbblok[2] = solvptr->nb2dblok - solvptr->nb2dcblk;

    fblok[0] = 0;
    fblok[1] = 0;
    fblok[2] = ((solvptr->cblktab + solvptr->cblkmin2d)->fblokptr - solvptr->bloktab);

    lblok[0] = solvptr->bloknbr;
    lblok[1] = (solvptr->cblktab + solvptr->cblkmax1d + 1)->fblokptr - solvptr->bloktab;
    lblok[2] = solvptr->bloknbr;

    fprintf( stdout,
             "    Solver Matrix statistics:         | %-12s | %-12s | %-12s |\n"
             "    --------------------------------------------------------------------------------\n"
             "      Number of cblk                  | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " |\n"
             "      Number of block                 | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " |\n"
             "      Number of block (diag part.)    | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " |\n"
             "      Cblk:   first                   | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " |\n"
             "              last                    | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " |\n"
             "      Block:  first                   | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " |\n"
             "              last                    | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " |\n"
             "      rownum: first                   | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " |\n"
             "              last                    | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " |\n"
             "      Average width                   | %12.2lf | %12.2lf | %12.2lf |\n"
             "      Average height                  | %12.2lf | %12.2lf | %12.2lf |\n"
             "      Structure memory space           %11.2lf %co\n"
             "      Number of coeficients stored      %10ld\n",
             /* Header */
             "All", "1d", "2d",
             /* Number of cblk */
             nbcblk[0], nbcblk[1], nbcblk[2],
             /* Number of blok without diagonal blocks (Add cblknbr to get the total) */
             nbblok[0], nbblok[1], nbblok[2],
             /*
              * Number of block in the compressed version (count 1 per block in
              * the matrix partition following the diagonal blocks
              */
             nbpartblok[0], nbpartblok[1], nbpartblok[2],
             /* Cblk */
             (int64_t)0,         (int64_t)0,                        (int64_t)(solvptr->cblkmin2d),
             (int64_t)(cblknbr), (int64_t)(solvptr->cblkmax1d + 1), (int64_t)(cblknbr),
             /* Blok */
             fblok[0], fblok[1], fblok[2],
             lblok[0], lblok[1], lblok[2],
             /* Rownum/Colnum */
             fcol[0], fcol[1], fcol[2],
             lcol[0], lcol[1], lcol[2],
             /* Average width */
             (double)(width[0]) / (double)(nbcblk[0]),
             (double)(width[1]) / (double)(nbcblk[1]),
             (double)(width[2]) / (double)(nbcblk[2]),
             /* Average height */
             (double)(height[0]) / (double)(nbblok[0]),
             (double)(height[1]) / (double)(nbblok[1]),
             (double)(height[2]) / (double)(nbblok[2]),
             /* Memory space */
             pastix_print_value( memstruct ),
             pastix_print_unit( memstruct ),
             /* Memory coefficient */
             (long)memcoef );

    fprintf( stdout,
             "      Number of GEMM tasks:           | %-12s | %-12s | %-12s | %-12s |\n"
             "        - All blocks                  | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " |\n"
             "        - PaRSEC                      | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " |\n"
             "        - StarPU                      | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " | %12" PRIi64 " |\n",
             "Dense", "Full2d", "Hybrid", "Full1d",
             gemm_dense, gemm_nopart_full2, gemm_nopart_hybrid, gemm_full1,
             gemm_dense, gemm_parsec_full2, gemm_parsec_hybrid, gemm_full1,
             gemm_dense, gemm_starpu_full2, gemm_starpu_hybrid, gemm_full1 );
}

/**
 *******************************************************************************
 *
 * @brief Instanciate the arrays for the requests according to the scheduler.
 *
 *******************************************************************************
 *
 * @param[in] solve_step
 *          Define which step of the solve is concerned.
 *          @arg PastixSolveForward
 *          @arg PastixSolveBackward
 *          @arg PastixFacto
 *
 * @param[inout] solvmtx
 *          The pointer to the solver matrix structure.
 *
 *******************************************************************************/
void
solverRequestInit( solve_step_t  solve_step,
                   SolverMatrix *solvmtx )
{
    MPI_Request  *request;
    pastix_int_t *reqindx;
    pastix_int_t  i, reqnbr;

    if ( solve_step == PastixSolveBackward ) {
        reqnbr = solvmtx->recvnbr + 1;
    }
    else {
        reqnbr = solvmtx->faninnbr + 1;
    }

    /* Restore recv and fanin nbr */
    solvmtx->fanincnt = solvmtx->faninnbr;
    solvmtx->recvcnt  = solvmtx->recvnbr;

    solvmtx->reqnbr  = reqnbr;
    solvmtx->reqlock = PASTIX_ATOMIC_UNLOCKED;

    MALLOC_INTERN( solvmtx->reqtab, reqnbr, MPI_Request  );
    MALLOC_INTERN( solvmtx->reqidx, reqnbr, pastix_int_t );

    request = solvmtx->reqtab;
    reqindx = solvmtx->reqidx;
    for ( i = 0; i < reqnbr; i++, request++, reqindx++ )
    {
        *request = MPI_REQUEST_NULL;
        *reqindx = -1;
    }
    solverComMatrixInit( solvmtx );

    return;
}

/**
 *******************************************************************************
 *
 * @brief Free the arrays related to the requests
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The pointer to the solver matrix structure.
 *
 *******************************************************************************/
void
solverRequestExit( SolverMatrix *solvmtx )
{
    assert( solvmtx->reqnum == 0 );
    assert( solvmtx->reqlock == PASTIX_ATOMIC_UNLOCKED );

    if( solvmtx->reqtab ) {
        memFree_null( solvmtx->reqtab );
    }
    if( solvmtx->reqidx ) {
        memFree_null( solvmtx->reqidx );
    }
    solverComMatrixGather( solvmtx );
    solverComMatrixExit( solvmtx );
}

/**
 *******************************************************************************
 *
 * @brief Allocate the reception buffer, and initiate the first persistant
 * reception
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the cblk must be tested.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] solvmtx
 *          The pointer to the solver matrix structure.
 *
 * @param[in] flttype
 *          Define which type are the coefficients.
 *          @arg PastixFloat
 *          @arg PastixDouble
 *          @arg PastixComplex32
 *          @arg PastixComplex64
 *
 *******************************************************************************/
void
solverRecvInit( pastix_coefside_t side,
                SolverMatrix     *solvmtx,
                pastix_coeftype_t flttype )
{
    /* Compute the max size (in bytes) for the communication buffer */
    pastix_int_t size = pastix_size_of(flttype) * solvmtx->maxrecv;
    size *= (side == PastixLUCoef) ? 2 : 1;

    if( solvmtx->recvnbr == 0 ) {
        return;
    }

    assert( solvmtx->maxrecv > 0 );

    /* Init communication */
    MALLOC_INTERN( solvmtx->rcoeftab, size, char );
    MPI_Recv_init( solvmtx->rcoeftab, size,
                   MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
                   solvmtx->solv_comm, solvmtx->reqtab );
    MPI_Start( solvmtx->reqtab );

    assert( solvmtx->reqnum == 0 );
    solvmtx->reqnum++;
#if defined(PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] Start persistant recv from any source\n",
             solvmtx->clustnum );
#endif
}

/**
 *******************************************************************************
 *
 * @brief Free the array linked to pending reception.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The pointer to the solver matrix structure.
 *
 *******************************************************************************/
void
solverRecvExit( SolverMatrix *solvmtx )
{
    /* In fact, the pointer should never been freed by this call */
    assert( solvmtx->reqtab == NULL );
    if( solvmtx->rcoeftab ) {
        memFree_null( solvmtx->rcoeftab );
    }
}

/**
 *******************************************************************************
 *
 * @brief Computes the max size of recv cblk.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The pointer to the solver matrix structure.
 *
 *******************************************************************************
 *
 * @return maximun recv size
 *
 *******************************************************************************/
static inline pastix_int_t
solverRhsRecvMax( SolverMatrix *solvmtx )
{
    pastix_int_t      cblknbr = solvmtx->cblknbr;
    const SolverCblk *cblk;
    pastix_int_t      k, max = 0;

    cblk = solvmtx->cblktab;
    for ( k = 0; k < cblknbr; k++, cblk++ ) {
        if ( cblk->cblktype & (CBLK_RECV | CBLK_FANIN) ) {
            max = pastix_imax( max, cblk_colnbr( cblk ) );
        }
    }

    return max;
}

/**
 *******************************************************************************
 *
 * @brief Allocates the reception buffer, and initiate the first persistant
 * reception
 *
 *******************************************************************************
 *
 * @param[in] solve_step
 *          Define which step of the solve is concerned.
 *          @arg PastixSolveForward
 *          @arg PastixSolveBackward
 *
 * @param[inout] solvmtx
 *          The pointer to the solver matrix structure.
 *
 * @param[in] flttype
 *          Define which type are the coefficients.
 *          @arg PastixFloat
 *          @arg PastixDouble
 *          @arg PastixComplex32
 *          @arg PastixComplex64
 *
 * @param[inout] rhsb
 *          The pointer to the rhs data structure that holds the vectors of the
 *          right hand side.
 *
 *******************************************************************************/
void
solverRhsRecvInit( solve_step_t      solve_step,
                   SolverMatrix     *solvmtx,
                   pastix_coeftype_t flttype,
                   pastix_rhs_t      rhsb  )
{
    /* Computes the max size (in bytes) for the communication buffer */
    pastix_int_t size;

    if ( ( ( solve_step == PastixSolveForward  ) && ( solvmtx->recvnbr  == 0 ) ) ||
         ( ( solve_step == PastixSolveBackward ) && ( solvmtx->faninnbr == 0 ) ) )
    {
        return;
    }

    size = pastix_size_of(flttype) * solverRhsRecvMax( solvmtx ) * rhsb->n;

    /* Init communication */
    MALLOC_INTERN( solvmtx->rcoeftab, size, char );
    MPI_Recv_init( solvmtx->rcoeftab, size,
                   MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
                   solvmtx->solv_comm, solvmtx->reqtab );
    MPI_Start( solvmtx->reqtab );

    assert( solvmtx->reqnum == 0 );
    solvmtx->reqnum++;
#if defined(PASTIX_DEBUG_MPI)
    fprintf( stderr, "[%2d] Start persistant recv from any source (max = %ld B)\n",
             solvmtx->clustnum, (long)size );
#endif
}

/**
 *******************************************************************************
 *
 * @brief Frees the array linked to pending reception.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          The pointer to the solver matrix structure.
 *
 *******************************************************************************/
void
solverRhsRecvExit( SolverMatrix *solvmtx )
{
    solverRecvExit( solvmtx );
}

/**
 *@}
 */
