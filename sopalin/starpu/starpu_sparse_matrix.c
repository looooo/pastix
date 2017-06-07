/**
 *
 * @file starpu_sparse_matrix.c
 *
 * PaStiX sparse matrix descriptor for StarPU.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2013-06-24
 *
 * @addtogroup pastix_starpu
 * @{
 *
 **/
#include "common.h"
#include "solver.h"
#include <starpu_data.h>

void
starpu_sparse_matrix_init( SolverMatrix *solvmtx,
                           int typesize, int mtxtype,
                           int nodes, int myrank )
{
    pastix_int_t   cblknbr, cblkmin2d;
    parsec_data_key_t key1, key2;
    SolverCblk *cblk;
    SolverBlok *blok, *fblok, *lblok;
    pastix_int_t m=0, n=0, cblknum;
    pastix_int_t nbcol, nbrow, ld;
    size_t offset;
    char *ptrL, *ptrU;

    starpu_sparse_matrix_desc_t *spmtx = solvmtx->starpu_desc;
    if ( spmtx != NULL ) {
        starpu_sparse_matrix_destroy( spmtx );
    }
    else {
        spmtx = (starpu_sparse_matrix_desc_t*)malloc(sizeof(starpu_sparse_matrix_desc_t));
    }

    spmtx->typesze = typesize;
    spmtx->mtxtype = mtxtype;
    spmtx->solvmtx = solvmtx;

    cblknbr   = solvmtx->cblknbr;
    cblkmin2d = solvmtx->cblkmin2d;
    ld        = solvmtx->cblkmaxblk * 2;
    key1      = 2 * cblknbr;

    /* Initialize 1D cblk handlers */
    cblk = spmtx->solvmtx->cblktab;
    for(cblknum = 0;
        cblknum < cblkmin2d;
        cblknum++, n++, cblk++ )
    {
        starpu_data_handle_t *handler = (starpu_data_handle_t*)(cblk->handler);
        nbrow = cblk->stride;
        nbcol = cblk_colnbr( cblk );
        ld = nbrow;

        starpu_matrix_data_register( handler, STARPU_MAIN_RAM,
                                     (uintptr_t)(cblk->lcoeftab), ld, nbrow, nbcol, spmtx->typesze );

        if ( mtxtype == PastixGeneral ) {
            starpu_matrix_data_register( handler + 1, STARPU_MAIN_RAM,
                                         (uintptr_t)(cblk->ucoeftab), ld, nbrow, nbcol, spmtx->typesze );
        }
    }

    /* Initialize 2D cblk handlers */
    cblk = spmtx->solvmtx->cblktab + cblkmin2d;
    for(cblknum = cblkmin2d, n = 0;
        cblknum < cblknbr;
        cblknum++, n++, cblk++ )
    {
        starpu_data_handle_t *handler = (starpu_data_handle_t*)(cblk->handler);
        nbrow = cblk->stride;
        nbcol = cblk_colnbr( cblk );
        ld = nbrow;

        starpu_matrix_data_register( handler, STARPU_MAIN_RAM,
                                     (uintptr_t)(cblk->lcoeftab), ld, nbrow, nbcol, spmtx->typesze );

        if ( mtxtype == PastixGeneral ) {
            starpu_matrix_data_register( handler + 1, STARPU_MAIN_RAM,
                                         (uintptr_t)(cblk->ucoeftab), ld, nbrow, nbcol, spmtx->typesze );
        }

        if ( !(cblk->cblktype & CBLK_TASKS_2D) )
            continue;

        /**
         * Diagonal block
         */
        ptrL   = cblk->lcoeftab;
        ptrU   = cblk->ucoeftab;
        blok   = cblk->fblokptr;
        nbrow  = blok_rownbr( blok );
        ld     = nbrow;
        offset = blok->coefind * (size_t)spmtx->typesze;
        key2   = n * ld;

        assert(offset == 0);
        starpu_matrix_data_register( (starpu_data_handle_t*)(&(blok->handler[0])), STARPU_MAIN_RAM,
                                     (uintptr_t)(ptrL + offset), ld, nbrow, nbcol, spmtx->typesze );

        if ( mtxtype == PastixGeneral ) {
            starpu_matrix_data_register( (starpu_data_handle_t*)(&(blok->handler[1])), STARPU_MAIN_RAM,
                                         (uintptr_t)(ptrU + offset), ld, nbrow, nbcol, spmtx->typesze );
        }
        else {
            blok->handler[1] = NULL;
        }

        /**
         * Lower Part
         */
        blok++; key2 += 2;
        lblok = cblk[1].fblokptr;
        for( ; blok < lblok; blok++, key2+=2 )
        {
            fblok = blok;
            m = 0;
            nbrow  = blok_rownbr( blok );
            offset = blok->coefind * (size_t)spmtx->typesze;

            while( (blok < lblok) &&
                   (blok[0].fcblknm == blok[1].fcblknm) &&
                   (blok[0].lcblknm == blok[1].lcblknm) )
            {
                blok++; m++;
                nbrow += blok_rownbr( blok );
            }

            ld = nbrow;
            starpu_matrix_data_register( (starpu_data_handle_t*)(&(fblok->handler[0])), STARPU_MAIN_RAM,
                                         (uintptr_t)(ptrL + offset), ld, nbrow, nbcol, spmtx->typesze );

            if ( mtxtype == PastixGeneral ) {
                starpu_matrix_data_register( (starpu_data_handle_t*)(&(fblok->handler[1])), STARPU_MAIN_RAM,
                                             (uintptr_t)(ptrU + offset), ld, nbrow, nbcol, spmtx->typesze );
            }
            else {
                fblok->handler[1] = NULL;
            }

            key2 += m * 2;
        }
    }

    solvmtx->starpu_desc = spmtx;

    (void)key1; (void)key2;
    (void)nodes; (void)myrank;
}

void
starpu_sparse_matrix_getoncpu( starpu_sparse_matrix_desc_t *spmtx )
{
    SolverCblk *cblk;
    SolverBlok *blok, *lblok;
    pastix_int_t i, cblkmin2d;

    cblkmin2d = spmtx->solvmtx->cblkmin2d;
    cblk = spmtx->solvmtx->cblktab;
    for(i=0; i<cblkmin2d; i++, cblk++)
    {
        assert( cblk->handler[0] );
        starpu_data_acquire( cblk->handler[0], STARPU_R );
        starpu_data_release( cblk->handler[0] );

        if ( spmtx->mtxtype == PastixGeneral ) {
            assert( cblk->handler[1] );
            starpu_data_acquire( cblk->handler[1], STARPU_R );
            starpu_data_release( cblk->handler[1] );
        }
    }

    for(i=cblkmin2d; i<spmtx->solvmtx->cblknbr; i++, cblk++)
    {
        assert( cblk->handler[0] );
        starpu_data_acquire( cblk->handler[0], STARPU_R );
        starpu_data_release( cblk->handler[0] );

        if ( spmtx->mtxtype == PastixGeneral ) {
            assert( cblk->handler[1] );
            starpu_data_acquire( cblk->handler[1], STARPU_R );
            starpu_data_release( cblk->handler[1] );
        }

        blok  = cblk->fblokptr;
        lblok = cblk[1].fblokptr;
        while( blok < cblk[1].fblokptr )
        {
            assert( blok->handler[0] );
            starpu_data_acquire( blok->handler[0], STARPU_R );
            starpu_data_release( blok->handler[0] );

            if ( spmtx->mtxtype == PastixGeneral ) {
                assert( blok->handler[1] );
                starpu_data_acquire( blok->handler[1], STARPU_R );
                starpu_data_release( blok->handler[1] );
            }
            while( (blok < lblok) &&
                   (blok[0].fcblknm == blok[1].fcblknm) &&
                   (blok[0].lcblknm == blok[1].lcblknm) )
            {
                blok++;
            }
            blok++;
        }
    }
}

void
starpu_sparse_matrix_destroy( starpu_sparse_matrix_desc_t *spmtx )
{
    SolverCblk *cblk;
    SolverBlok *blok;
    pastix_int_t i, cblkmin2d;

    cblkmin2d = spmtx->solvmtx->cblkmin2d;
    cblk = spmtx->solvmtx->cblktab;
    for(i=0; i<cblkmin2d; i++, cblk++)
    {
        if ( cblk->handler[0] ) {
            starpu_data_unregister( cblk->handler[0] );

            if ( spmtx->mtxtype == PastixGeneral ) {
                starpu_data_unregister( cblk->handler[1] );
            }
        }

        cblk->handler[0] = NULL;
        cblk->handler[1] = NULL;
    }

    for(i=cblkmin2d; i<spmtx->solvmtx->cblknbr; i++, cblk++)
    {
        if ( cblk->handler[0] ) {
            starpu_data_unregister( cblk->handler[0] );
            if ( spmtx->mtxtype == PastixGeneral ) {
                starpu_data_unregister( cblk->handler[1] );
            }
        }

        cblk->handler[0] = NULL;
        cblk->handler[1] = NULL;

        blok = cblk->fblokptr;
        while( blok < cblk[1].fblokptr )
        {
            if ( blok->handler[0] ) {
                starpu_data_unregister( blok->handler[0] );
                if ( spmtx->mtxtype == PastixGeneral ) {
                    starpu_data_unregister( blok->handler[1] );
                }
            }

            blok->handler[0] = NULL;
            blok->handler[1] = NULL;

            blok++;
        }
    }
}

/**
 *@}
 */
