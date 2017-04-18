/**
 *
 * @file sparse-matrix.c
 *
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2011-10-17
 * @precisions normal z -> c d s
 *
 **/
#define _GNU_SOURCE
#include "common.h"
#include "solver.h"
#include "sopalin/parsec/pastix_parsec.h"

#include <parsec.h>
#include <parsec/data.h>
#include <parsec/data_distribution.h>

#include "common.h"

static inline void
spm_data_key_to_value( parsec_data_key_t key,
                       int ratio, const SolverMatrix *solvmtx,
                       int *uplo,
                       pastix_int_t *cblknum,
                       pastix_int_t *bloknum)
{
    parsec_data_key_t key2;
    pastix_int_t cblkmin2d, cblknbr;

    cblknbr   = solvmtx->cblknbr;
    cblkmin2d = solvmtx->cblkmin2d;
    key2 = ratio * cblknbr;

    /* This is a block */
    if ( key >= key2 ) {
        pastix_int_t m, n, ld;

        key2 = key - key2;
        ld   = solvmtx->cblkmaxblk * ratio;

        m = key2 % ld;
        n = key2 / ld;

        *uplo    = m % ratio;
        *bloknum = m / ratio;
        *cblknum = cblkmin2d + n;
    }
    /* This is a cblk */
    else {
        *uplo    = key % ratio;
        *cblknum = key / ratio;
        *bloknum = -1;
    }
}

static uint32_t
sparse_matrix_data_key(parsec_ddesc_t *mat, ... )
{
    va_list ap;
    sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat;
    int uplo, ratio;
    pastix_int_t cblknum, bloknum;

    va_start(ap, mat);
    uplo    = va_arg(ap, int);
    cblknum = va_arg(ap, int);
    bloknum = va_arg(ap, int) - 1;
    va_end(ap);

    ratio = spmtx->mtxtype == PastixGeneral ? 2 : 1;
    uplo = uplo ? 1 : 0;
    assert( ratio == 2 || uplo == 0 );

    if ( bloknum == -1 ) {
        return cblknum * ratio + uplo;
    }
    else {
        pastix_int_t offset, ld, cblknbr;
        pastix_int_t cblkmin2d, n;

        cblknbr   = spmtx->solvmtx->cblknbr;
        cblkmin2d = spmtx->solvmtx->cblkmin2d;
        ld        = spmtx->solvmtx->cblkmaxblk * ratio;
        offset    = cblknbr * ratio;
        n         = cblknum - cblkmin2d;

        return offset + n * ld + bloknum * ratio + uplo;
    }
}

static uint32_t
sparse_matrix_rank_of(parsec_ddesc_t *mat, ... )
{
    (void)mat;
    return 0;
}

static uint32_t
sparse_matrix_rank_of_key(parsec_ddesc_t *mat, parsec_data_key_t key )
{
    (void)mat; (void)key;
    return 0;
}

static int32_t
sparse_matrix_vpid_of(parsec_ddesc_t *mat, ... )
{
    (void)mat;
    return 0;
}

static int32_t
sparse_matrix_vpid_of_key(parsec_ddesc_t *mat, parsec_data_key_t key )
{
    (void)mat; (void)key;
    return 0;
}

static parsec_data_t *
sparse_matrix_data_of(parsec_ddesc_t *mat, ... )
{
    sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat;
    SolverCblk *cblk;
    va_list ap;
    int uplo;
    pastix_int_t cblknum, bloknum;

    va_start(ap, mat);
    uplo    = va_arg(ap, int);
    cblknum = va_arg(ap, int);
    bloknum = va_arg(ap, int) - 1;
    va_end(ap);

    uplo = uplo ? 1 : 0;

    cblk = spmtx->solvmtx->cblktab + cblknum;

    /* This is a cblk */
    if ( bloknum == -1 ) {
        assert( cblk->handler[uplo] );
        return (parsec_data_t*)(cblk->handler[uplo]);
    }
    /* This is a blok */
    else {
        SolverBlok *blok = cblk->fblokptr + bloknum;

        assert( blok->handler[uplo] );
        return (parsec_data_t*)(blok->handler[uplo]);
    }
}

static parsec_data_t *
sparse_matrix_data_of_key(parsec_ddesc_t *mat, parsec_data_key_t key )
{
    sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat;
    SolverMatrix *solvmtx = spmtx->solvmtx;
    SolverCblk *cblk;
    int uplo, ratio;
    pastix_int_t cblknum, bloknum;

    ratio = (spmtx->mtxtype == PastixGeneral) ? 2 : 1;
    spm_data_key_to_value( key, ratio, solvmtx,
                           &uplo, &cblknum, &bloknum );

    cblk = solvmtx->cblktab + cblknum;

    /* This is a cblk */
    if ( bloknum == -1 ) {
        assert( cblk->handler[uplo] );
        return (parsec_data_t*)(cblk->handler[uplo]);
    }
    /* This is a blok */
    else {
        SolverBlok *blok = cblk->fblokptr + bloknum;

        assert( blok->handler[uplo] );
        return (parsec_data_t*)(blok->handler[uplo]);
    }
}

#ifdef PARSEC_PROF_TRACE
static int
sparse_matrix_key_to_string( parsec_ddesc_t *mat,
                             uint32_t datakey,
                             char *buffer, uint32_t buffer_size )
{
    sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat;
    int uplo;
    pastix_int_t cblknum, bloknum;
    int res, ratio;

    ratio = (spmtx->mtxtype == PastixGeneral) ? 2 : 1;
    spm_data_key_to_value( datakey, ratio, spmtx->solvmtx,
                           &uplo, &cblknum, &bloknum );

    res = snprintf(buffer, buffer_size, "(%d, %ld, %ld)",
                   uplo, (long int)cblknum, (long int)bloknum);
    if (res < 0)
    {
        printf("error in key_to_string for tile (%d, %ld, %ld) key: %u\n",
               uplo, (long int)cblknum, (long int)bloknum, datakey);
    }
    return res;
}
#endif

void
sparse_matrix_init( sparse_matrix_desc_t *spmtx,
                    SolverMatrix *solvmtx,
                    int typesize, int mtxtype,
                    int nodes, int myrank )
{
    parsec_ddesc_t *o = (parsec_ddesc_t*)spmtx;
    pastix_int_t   cblknbr, cblkmin2d, ld;
    parsec_data_key_t key1, key2;
    SolverCblk *cblk;
    SolverBlok *blok, *lblok;
    pastix_int_t m, n, cblknum;
    size_t size, offset;
    char *ptrL, *ptrU;
    int ratio = ( mtxtype == PastixGeneral ) ? 2 : 1;

    parsec_ddesc_init( o, nodes, myrank );

    o->data_key      = sparse_matrix_data_key;
#if defined(PARSEC_PROF_TRACE)
    o->key_to_string = sparse_matrix_key_to_string;
#endif

    o->rank_of     = sparse_matrix_rank_of;
    o->rank_of_key = sparse_matrix_rank_of_key;
    o->vpid_of     = sparse_matrix_vpid_of;
    o->vpid_of_key = sparse_matrix_vpid_of_key;
    o->data_of     = sparse_matrix_data_of;
    o->data_of_key = sparse_matrix_data_of_key;

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
        parsec_data_t **handler = (parsec_data_t**)(cblk->handler);
        size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;

        parsec_data_create( handler,
                            o, cblknum * ratio, cblk->lcoeftab, size );

        if ( ratio == 2 ) {
            parsec_data_create( handler+1,
                                o, cblknum * ratio + 1, cblk->ucoeftab, size );
        }
    }

    /* Initialize 2D cblk handlers */
    cblk = spmtx->solvmtx->cblktab + cblkmin2d;
    for(cblknum = cblkmin2d, n = 0;
        cblknum < cblknbr;
        cblknum++, n++, cblk++ )
    {
        parsec_data_t **handler = (parsec_data_t**)(cblk->handler);
        size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;

        parsec_data_create( handler,
                            o, cblknum * ratio, cblk->lcoeftab, size );

        if ( ratio == 2 ) {
            parsec_data_create( handler+1,
                                o, cblknum * ratio + 1, cblk->ucoeftab, size );
        }

        if ( !(cblk->cblktype & CBLK_LAYOUT_2D) )
            continue;

        /**
         * Diagonal block
         */
        ptrL   = cblk->lcoeftab;
        ptrU   = cblk->ucoeftab;
        blok   = cblk->fblokptr;
        size   = blok_rownbr( blok ) * cblk_colnbr( cblk ) * (size_t)spmtx->typesze;
        offset = blok->coefind * (size_t)spmtx->typesze;
        key2   = n * ld;

        assert(offset == 0);
        parsec_data_create( (parsec_data_t**)&(blok->handler[0]),
                            o, key1 + key2,
                            ptrL + offset, size );

        if ( ratio == 2 ) {
            parsec_data_create( (parsec_data_t**)&(blok->handler[1]),
                                o, key1 + key2 + 1,
                                ptrU + offset, size );
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
            size   = blok_rownbr( blok );
            offset = blok->coefind * (size_t)spmtx->typesze;

            while( (blok < lblok) &&
                   (blok[0].fcblknm == blok[1].fcblknm) &&
                   (blok[0].lcblknm == blok[1].lcblknm) )
            {
                blok++; m++;
                size += blok_rownbr( blok );
            }
            size *= cblk_colnbr( cblk )
                *  (size_t)spmtx->typesze;

            parsec_data_create( (parsec_data_t**)&(blok->handler[0]),
                                &spmtx->super, key1 + key2,
                                ptrL + offset, size );

            if ( ratio == 2 ) {
                parsec_data_create( (parsec_data_t**)&(blok->handler[1]),
                                    &spmtx->super, key1 + key2 + 1,
                                    ptrU + offset, size );
            }
            else {
                blok->handler[1] = NULL;
            }

            key2 += m * ratio;
        }
    }
}

void
sparse_matrix_destroy( sparse_matrix_desc_t *spmtx )
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

        cblk->handler[0] = NULL;
        cblk->handler[1] = NULL;
    }

    for(i=cblkmin2d; i<spmtx->solvmtx->cblknbr; i++, cblk++)
    {
        if ( cblk->handler[0] ) {
            parsec_data_destroy( cblk->handler[0] );
            if (ratio == 2) {
                parsec_data_destroy( cblk->handler[1] );
            }
        }

        cblk->handler[0] = NULL;
        cblk->handler[1] = NULL;

        blok = cblk->fblokptr;
        while( blok < cblk[1].fblokptr )
        {
            if ( blok->handler[0] ) {
                parsec_data_destroy( blok->handler[0] );
                if (ratio == 2) {
                    parsec_data_destroy( blok->handler[1] );
                }
            }

            blok->handler[0] = NULL;
            blok->handler[1] = NULL;

            blok++;
        }
    }

    parsec_ddesc_destroy( (parsec_ddesc_t*)spmtx );
}
