/**
 *
 * @file coeftab_zinit.c
 *
 * Precision dependent coeficient array initialization routines.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> s d c
 *
 **/
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "common.h"
#include "solver.h"
#include "bcsc.h"
#include "sopalin/coeftab_z.h"
#include "pastix_zcores.h"

/*
 * Function: z_CoefMatrix_Allocate
 *
 * Allocate matrix coefficients in coeftab and ucoeftab.
 *
 * Should be first called with me = -1 to allocated coeftab.
 * Then, should be called with me set to thread ID
 * to allocate column blocks coefficients arrays.
 *
 * Parameters
 *
 *    datacode  - solverMatrix
 *    factotype - factorization type (LU, LLT ou LDLT)
 *    me        - thread number. (-1 for first call,
 *                from main thread. >=0 to allocate column blocks
 *     assigned to each thread.)
 */
void
coeftab_zcblkfill( const SolverMatrix  *solvmtx,
                   const pastix_bcsc_t *bcsc,
                   pastix_int_t itercblk,
                   int factoLU )
{
    pastix_coefside_t side = factoLU ? PastixLUCoef : PastixLCoef;
    SolverCblk       *cblk = solvmtx->cblktab + itercblk;

    cpucblk_zalloc( side, cblk );
    cpucblk_zffbcsc( factoLU ? PastixLUCoef : PastixLCoef,
                     solvmtx, bcsc, itercblk );
}

/*
 * Function: z_CoefMatrix_Allocate
 *
 * Allocate matrix coefficients in coeftab and ucoeftab.
 *
 * Should be first called with me = -1 to allocated coeftab.
 * Then, should be called with me set to thread ID
 * to allocate column blocks coefficients arrays.
 *
 * Parameters
 *
 *    datacode  - solverMatrix
 *    factotype - factorization type (LU, LLT ou LDLT)
 *    me        - thread number. (-1 for first call,
 *                from main thread. >=0 to allocate column blocks
 *     assigned to each thread.)
 */
void
coeftab_zcblkinit( const SolverMatrix  *solvmtx,
                   const pastix_bcsc_t *bcsc,
                   pastix_int_t         itercblk,
                   int                  factoLU,
                   char               **directory )
{
    SolverCblk *cblk           = solvmtx->cblktab + itercblk;
    pastix_int_t compress_when = solvmtx->lowrank.compress_when;

    coeftab_zcblkfill( solvmtx, bcsc, itercblk, factoLU );

#if defined(PASTIX_DEBUG_DUMP_COEFTAB)
    {
        FILE *f = NULL;
        char *filename;
        int rc;

        rc = asprintf( &filename, "Lcblk%05ld_init.txt", itercblk );
        f  = pastix_fopenw( directory, filename, "w" );
        if ( f != NULL ) {
            coeftab_zcblkdump( cblk, PastixLower, f );
            fclose( f );
        }
        free( filename );

        if ( cblk->ucoeftab ) {
            rc = asprintf( &filename, "Ucblk%05ld_init.txt", itercblk );
            f  = pastix_fopenw( directory, filename, "w" );
            if ( f != NULL ) {
                coeftab_zcblkdump( cblk, PastixUpper, f );
                fclose( f );
            }
            free( filename );
        }
        (void)rc;
    }
#endif /* defined(PASTIX_DEBUG_DUMP_COEFTAB) */

    /**
     * Try to compress the cblk if needs to be compressed
     * TODO: change the criteria based on the level in the tree
     */
    if ( cblk->cblktype & CBLK_LAYOUT_2D )
    {
        if ( compress_when == PastixCompressWhenBegin ){
            coeftab_zcompress_one( cblk, solvmtx->lowrank );
        }
        else if ( compress_when == PastixCompressWhenEnd ){
            coeftab_zalloc_one( cblk );
        }
    }

    (void)directory;
}
