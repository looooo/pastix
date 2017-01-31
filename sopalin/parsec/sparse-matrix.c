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
    char *ptr;
    va_list ap;
    int uplo, ratio;
    pastix_int_t cblknum, bloknum;
    parsec_data_key_t key1, key2;
    size_t size;

    va_start(ap, mat);
    uplo    = va_arg(ap, int);
    cblknum = va_arg(ap, int);
    bloknum = va_arg(ap, int) - 1;
    va_end(ap);

    ratio = spmtx->mtxtype == PastixGeneral ? 2 : 1;
    uplo = uplo ? 1 : 0;
    assert( ratio == 2 || uplo == 0 );

    cblk = spmtx->solvmtx->cblktab + cblknum;
    ptr  = uplo ? cblk->ucoeftab : cblk->lcoeftab;

    /* This is a cblk */
    if ( bloknum == -1 ) {
        key1 = ratio * cblknum + uplo;
        size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;
        return parsec_data_create( (parsec_data_t**)(&(cblk->handler[uplo])),
                                   mat, key1, ptr, size );
    }
    /* This is a blok */
    else {
        pastix_int_t n, cblkmin2d, cblknbr, ld;
        SolverBlok *blok = cblk->fblokptr + bloknum;

        cblknbr   = spmtx->solvmtx->cblknbr;
        cblkmin2d = spmtx->solvmtx->cblkmin2d;
        ld        = spmtx->solvmtx->cblkmaxblk * ratio;
        n         = cblknum - cblkmin2d;

        key1 = ratio * cblknbr;
        key2 = n * ld + bloknum * ratio + uplo;

        assert( spmtx->datamap_blok[key2] != NULL );
        return parsec_data_create( (parsec_data_t**)(&(blok->handler[uplo])),
                                   mat, key1+key2, NULL, 0 );
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
    size_t size;
    char *ptr;

    ratio = (spmtx->mtxtype == PastixGeneral) ? 2 : 1;
    spm_data_key_to_value( key, ratio, solvmtx,
                           &uplo, &cblknum, &bloknum );

    cblk = solvmtx->cblktab + cblknum;
    ptr  = uplo ? cblk->ucoeftab : cblk->lcoeftab;

    /* This is a cblk */
    if ( bloknum == -1 ) {
        size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;
        return parsec_data_create( (parsec_data_t**)(&(cblk->handler[uplo])),
                                   mat, key, ptr, size );
    }
    /* This is a blok */
    else {
        SolverBlok *blok = cblk->fblokptr + bloknum;

        return parsec_data_create( (parsec_data_t**)(&(blok->handler[uplo])),
                                   mat, key, NULL, 0 );
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

    spmtx->datamap_cblk = NULL;
    spmtx->datamap_blok = NULL;

    if ( cblkmin2d < cblknbr ) {
        SolverCblk *cblkN;
        SolverBlok *blok, *lblok;
        pastix_int_t m, n, cblknumN;
        size_t size, offset;
        char *ptrL, *ptrU;

        cblkN = spmtx->solvmtx->cblktab + cblkmin2d;
        for(cblknumN = cblkmin2d, n = 0;
            cblknumN < cblknbr;
            cblknumN++, n++, cblkN++ )
        {
            if ( !(cblkN->cblktype & CBLK_LAYOUT_2D) )
                continue;

            /**
             * Diagonal block
             */
            ptrL   = cblkN->lcoeftab;
            ptrU   = cblkN->ucoeftab;
            blok   = cblkN->fblokptr;
            size   = blok_rownbr( blok ) * cblk_colnbr( cblkN ) * (size_t)spmtx->typesze;
            offset = blok->coefind * (size_t)spmtx->typesze;
            key2   = n * ld;

            assert(offset == 0);
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

            /**
             * Lower Part
             */
            blok++; key2+=ratio;
            lblok = cblkN[1].fblokptr;
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
                size *= cblk_colnbr( cblkN )
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
