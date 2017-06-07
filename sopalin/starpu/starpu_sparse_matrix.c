/**
 *
 * @file starpu_sparse_matrix.c
 *
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2011-10-17
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "solver.h"

void
starpu_sparse_matrix_init( starpu_sparse_matrix_desc_t *spmtx,
                           SolverMatrix *solvmtx,
                           int typesize, int mtxtype,
                           int nodes, int myrank )
{
    pastix_int_t   cblknbr, cblkmin2d, ld;
    parsec_data_key_t key1, key2;
    SolverCblk *cblk;
    SolverBlok *blok, *lblok;
    pastix_int_t m, n, cblknum;
    pastix_int_t nbcol, nbrow, ld;
    size_t size, offset;
    char *ptrL, *ptrU;
    int ratio = ( mtxtype == PastixGeneral ) ? 2 : 1;

    spmtx->typesze = typesize;
    spmtx->mtxtype = mtxtype;
    spmtx->solvmtx = solvmtx;

    cblknbr   = solvmtx->cblknbr;
    cblkmin2d = solvmtx->cblkmin2d;
    ld        = solvmtx->cblkmaxblk * ratio;
    key1      = ratio * cblknbr;

    /* Initialize 1D cblk handlers */
    cblk = spmtx->solvmtx->cblktab;
    for(cblknum = 0;
        cblknum < cblkmin2d;
        cblknum++, n++, cblk++ )
    {
        nbrow = cblk->stride;
        nbcol = cblk_colnbr( cblk );
        ld = nbrow;

        starpu_matrix_data_register( (starpu_data_handle_t*)(cblk->handler[0]), STARPU_MAIN_RAM,
                                     cblk->lcoeftab, ld, nbrow, nbcol, spmtx->typesze );

        if ( ratio == 2 ) {
            starpu_matrix_data_register( (starpu_data_handle_t*)(cblk->handler[1]), STARPU_MAIN_RAM,
                                         cblk->ucoeftab, ld, nbrow, nbcol, spmtx->typesze );
        }
    }

    /* Initialize 2D cblk handlers */
    cblk = spmtx->solvmtx->cblktab + cblkmin2d;
    for(cblknum = cblkmin2d, n = 0;
        cblknum < cblknbr;
        cblknum++, n++, cblk++ )
    {
        nbrow = cblk->stride;
        nbcol = cblk_colnbr( cblk );
        ld = nbrow;

        starpu_matrix_data_register( (starpu_data_handle_t*)(cblk->handler[0]), STARPU_MAIN_RAM,
                                     cblk->lcoeftab, ld, nbrow, nbcol, spmtx->typesze );

        if ( ratio == 2 ) {
            starpu_matrix_data_register( (starpu_data_handle_t*)(cblk->handler[1]), STARPU_MAIN_RAM,
                                         cblk->ucoeftab, ld, nbrow, nbcol, spmtx->typesze );
        }

        if ( !(cblk->cblktype & CBLK_TASK_2D) )
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
        starpu_matrix_data_register( (starpu_data_handle_t*)(blok->handler[0]), STARPU_MAIN_RAM,
                                     ptrL + offset, ld, nbrow, nbcol, spmtx->typesze );

        if ( ratio == 2 ) {
            starpu_matrix_data_register( (starpu_data_handle_t*)(blok->handler[1]), STARPU_MAIN_RAM,
                                         ptrU + offset, ld, nbrow, nbcol, spmtx->typesze );
        }
        else {
            blok->handler[1] = NULL;
        }

        /**
         * Lower Part
         */
        blok++; key2 += ratio;
        lblok = cblk[1].fblokptr;
        for( ; blok < lblok; blok++, key2+=ratio )
        {
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
            starpu_matrix_data_register( (starpu_data_handle_t*)(blok->handler[0]), STARPU_MAIN_RAM,
                                         ptrL + offset, ld, nbrow, nbcol, spmtx->typesze );

            if ( ratio == 2 ) {
                starpu_matrix_data_register( (starpu_data_handle_t*)(blok->handler[1]), STARPU_MAIN_RAM,
                                             ptrU + offset, ld, nbrow, nbcol, spmtx->typesze );
            }
            else {
                blok->handler[1] = NULL;
            }

            key2 += m * ratio;
        }
    }
}

void
starpu_sparse_matrix_destroy( starpu_sparse_matrix_desc_t *spmtx )
{
    SolverCblk *cblk;
    SolverBlok *blok;
    pastix_int_t i, cblkmin2d;
    int ratio = ( spmtx->mtxtype == PastixGeneral ) ? 2 : 1;

    cblkmin2d = spmtx->solvmtx->cblkmin2d;
    cblk = spmtx->solvmtx->cblktab;
    for(i=0; i<cblkmin2d; i++, cblk++)
    {
        if ( cblk->handler[0] ) {
            parsec_data_destroy( cblk->handler[0] );

            if (ratio == 2) {
                parsec_data_destroy( cblk->handler[1] );
            }
        }
    }

    for(i=cblkmin2d; i<spmtx->solvmtx->cblknbr; i++, cblk++)
    {
        if ( cblk->handler[0] ) {
            parsec_data_destroy( cblk->handler[0] );
            if (ratio == 2) {
                parsec_data_destroy( cblk->handler[1] );
            }
        }

        blok = cblk->fblokptr;
        while( blok < cblk[1].fblokptr )
        {
            if ( blok->handler[0] ) {
                parsec_data_destroy( blok->handler[0] );
                if (ratio == 2) {
                    parsec_data_destroy( blok->handler[1] );
                }
            }
            blok++;
        }
    }

    parsec_ddesc_destroy( (parsec_ddesc_t*)spmtx );
}
