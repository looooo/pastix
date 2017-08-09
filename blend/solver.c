/**
 * @file solver.c
 *
 * PaStiX solver structure basic functions.
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 *
 **/
#include "common.h"
#include "solver.h"
#include "coeftab.h"

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

    /* Fanins */
    if ( solvptr->ftgttab ) {
        mem += solvptr->ftgtnbr * sizeof( solver_ftgt_t );
    }
    if ( solvptr->indtab ) {
        mem += solvptr->indnbr  * sizeof( pastix_int_t );
    }

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

    /* proc2clust */
    if ( solvptr->proc2clust ) {
        mem += solvptr->procnbr * sizeof(pastix_int_t);
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
solverInit(SolverMatrix *solvmtx)
{
    memset(solvmtx, 0, sizeof (SolverMatrix));
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
    if(solvmtx->ftgttab) {
        memFree_null(solvmtx->ftgttab);
    }
    if(solvmtx->tasktab) {
        memFree_null(solvmtx->tasktab);
    }
    if(solvmtx->indtab) {
        memFree_null(solvmtx->indtab);
    }
    memFree_null(solvmtx->ttsknbr);
    for (i=0;i<solvmtx->bublnbr;i++)
    {
        if (solvmtx->ttsktab[i] != NULL) {
            memFree_null(solvmtx->ttsktab[i]);
        }
    }
    memFree_null(solvmtx->ttsktab);
    memFree_null(solvmtx->proc2clust);
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
    size_t memstruct, memcoef;
    pastix_int_t fcol2d, lcol2d;
    pastix_int_t itercblk, cblknbr;
    pastix_int_t gemm2d = 0;
    pastix_int_t gemm2dtot = 0;
    double avg2d;

    fcol2d = (solvptr->cblktab + solvptr->cblkmin2d )->fcolnum;
    lcol2d = (solvptr->cblktab + solvptr->cblknbr   )->fcolnum;
    avg2d  = 0.;
    if ( (lcol2d - fcol2d) > 0 ) {
        avg2d  = (double)(lcol2d - fcol2d) / (double)( solvptr->cblknbr - solvptr->cblkmin2d);
    }
    cblknbr = solvptr->cblknbr;
    cblk    = solvptr->cblktab;
    memcoef = 0;
    for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
    {
        pastix_int_t colnbr = cblk->lcolnum - cblk->fcolnum + 1;
        pastix_int_t rownbr = cblk->stride;

        memcoef += colnbr * rownbr;

        if (cblk->cblktype & CBLK_TASKS_2D) {
            pastix_int_t contrib = cblk[1].brownum - cblk[0].brown2d;
            gemm2d += contrib;
            gemm2dtot += contrib * (cblk[1].fblokptr - cblk[0].fblokptr);
        }
    }

    memstruct = solver_size( solvptr );

    fprintf( stdout,
             "    Solver Matrix statistics:\n"
             "      Number of cblk (1D|2D)            %10ld (%10ld | %10ld )\n"
             "      Number of block (1D|2D)           %10ld (%10ld | %10ld )\n"
             "      Cblk: last1D/first2d              %10ld | %10ld\n"
             "      First row handled in 2D           %10ld\n"
             "      Average 2D cblk size             %11.2lf\n"
             "      Maximum block 1D                  %10ld\n"
             "      Structure memory space           %11.2lf %s\n"
             "      Number of coeficients stored      %10ld\n"
             "      Number of 2D brow                 %10ld\n"
             "      Number of 2D Gemm tasks           %10ld\n",
             (long)solvptr->cblknbr, (long)(solvptr->cblknbr - solvptr->nb2dcblk), (long)solvptr->nb2dcblk,
             (long)solvptr->bloknbr, (long)(solvptr->bloknbr - solvptr->nb2dblok), (long)solvptr->nb2dblok,
             (long)solvptr->cblkmax1d, (long)solvptr->cblkmin2d,
             (long)fcol2d, (double)avg2d,
             (long)(((solvptr->cblktab + solvptr->cblkmax1d + 1)->fblokptr - solvptr->bloktab) - 1),
             MEMORY_WRITE( memstruct ), MEMORY_UNIT_WRITE( memstruct ),
             (long)memcoef, (long)gemm2d, (long)gemm2dtot );
}

/**
 *@}
 */
