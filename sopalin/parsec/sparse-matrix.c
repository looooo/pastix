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

#include <dague.h>
#include <dague/data.h>
#include <dague/data_distribution.h>

#include "common.h"

static inline pastix_int_t
spm_data_key( int ratio, int uplo,
              pastix_int_t cblknum,
              pastix_int_t bloknum )
{
    /**
     * Key is strictly negative for a cblk, postive or null for a blok, so
     * either blokum == 0, or cblknum == cblknbr toa ply a shift in the postive
     * range
     */
    return  ratio * cblknum + bloknum + uplo;
}

static inline void
spm_data_key_to_value( dague_data_key_t key,
                       int ratio, const SolverMatrix *solvmtx,
                       int *uplo,
                       pastix_int_t *cblknum,
                       pastix_int_t *bloknum)
{
    dague_data_key_t key2;
    pastix_int_t cblkmin2d, cblknbr;

    /* Refers to a block */
    cblknbr   = solvmtx->cblknbr;
    cblkmin2d = solvmtx->cblkmin2d;
    key2 = ratio * cblknbr;
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
    else {
        *uplo    = key % ratio;
        *cblknum = key / ratio;
        *bloknum = -1;
    }
}

static uint32_t
sparse_matrix_data_key(dague_ddesc_t *mat, ... )
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
sparse_matrix_rank_of(dague_ddesc_t *mat, ... )
{
    (void)mat;
    return 0;
}

static uint32_t
sparse_matrix_rank_of_key(dague_ddesc_t *mat, dague_data_key_t key )
{
    (void)mat; (void)key;
    return 0;
}

static int32_t
sparse_matrix_vpid_of(dague_ddesc_t *mat, ... )
{
    (void)mat;
    return 0;
}

static int32_t
sparse_matrix_vpid_of_key(dague_ddesc_t *mat, dague_data_key_t key )
{
    (void)mat; (void)key;
    return 0;
}

static dague_data_t *
sparse_matrix_data_of(dague_ddesc_t *mat, ... )
{
    sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat;
    SolverCblk *cblk;
    char *ptr;
    va_list ap;
    int uplo, ratio;
    pastix_int_t cblknum, bloknum;
    dague_data_key_t key1, key2;
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

    if ( bloknum == -1 ) {
        key1 = ratio * cblknum + uplo;
        size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;
        return dague_data_create( spmtx->datamap_cblk + key1,
                                  mat, key1, ptr, size );
    }
    else {
        pastix_int_t n, cblkmin2d, cblknbr, ld;

        cblknbr   = spmtx->solvmtx->cblknbr;
        cblkmin2d = spmtx->solvmtx->cblkmin2d;
        ld        = spmtx->solvmtx->cblkmaxblk * ratio;
        n         = cblknum - cblkmin2d;

        key1 = ratio * cblknbr;
        key2 = n * ld + bloknum * ratio + uplo;

        assert( spmtx->datamap_blok[key2] != NULL );
        return dague_data_create( spmtx->datamap_blok + key2, mat, key1+key2, NULL, 0 );
    }
}

static dague_data_t *
sparse_matrix_data_of_key(dague_ddesc_t *mat, dague_data_key_t key )
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

    if ( bloknum == -1 ) {
        size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;
        return dague_data_create( spmtx->datamap_cblk + key,
                                  mat, key, ptr, size );
    }
    else {
        dague_data_key_t key2;
        pastix_int_t n, cblkmin2d, ld;

        cblkmin2d = spmtx->solvmtx->cblkmin2d;
        ld        = spmtx->solvmtx->cblkmin2d * ratio;
        n         = cblknum - cblkmin2d;
        key2 = n * ld + bloknum * ratio + uplo;

        assert( spmtx->datamap_blok[key2] != NULL );
        return dague_data_create( spmtx->datamap_blok + key2, mat, key, NULL, 0 );
    }
}

#ifdef DAGUE_PROF_TRACE
static int sparse_matrix_key_to_string(dague_ddesc_t *mat, uint32_t datakey, char *buffer, uint32_t buffer_size)
{
    sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat;
    int uplo;
    pastix_int_t cblknum, bloknum;
    int res;

    spm_data_key_to_value( datakey, ratio, spm->solvmtx,
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

void sparse_matrix_init( sparse_matrix_desc_t *spmtx,
                         SolverMatrix *solvmtx,
                         int typesize, int mtxtype,
                         int nodes, int myrank)
{
    dague_ddesc_t *o    = (dague_ddesc_t*)spmtx;
    dague_data_t **datamap;
    pastix_int_t   cblknbr, cblkmin2d, ld;
    dague_data_key_t offkey, key2;
    int ratio = ( mtxtype == PastixGeneral ) ? 2 : 1;

    dague_ddesc_init( o, nodes, myrank );

    o->data_key      = sparse_matrix_data_key;
#if defined(DAGUE_PROF_TRACE)
    o->key_to_string = sparse_matrix_key_to_string;
#endif

    o->rank_of     = sparse_matrix_rank_of;
    o->rank_of_key = sparse_matrix_rank_of_key;
    o->vpid_of     = sparse_matrix_vpid_of;
    o->vpid_of_key = sparse_matrix_vpid_of_key;
    o->data_of     = sparse_matrix_data_of;
    o->data_of_key = sparse_matrix_data_of_key;

    spmtx->typesze   = typesize;
    spmtx->mtxtype   = mtxtype;
    spmtx->solvmtx   = solvmtx;

    cblknbr   = solvmtx->cblknbr;
    cblkmin2d = solvmtx->cblkmin2d;
    ld        = solvmtx->cblkmaxblk * ratio;
    offkey    = ratio * cblknbr;

    spmtx->datamap_cblk = (dague_data_t**)calloc( ratio * cblknbr,
                                                  sizeof(dague_data_t*) );
    if ( cblkmin2d < cblknbr ) {
        SolverCblk *cblkN;
        SolverBlok *blok, *lblok;
        pastix_int_t m, n, cblknumN;
        size_t size, offset;
        char *ptr;

        spmtx->datamap_blok = (dague_data_t**)calloc( ld * (cblknbr - cblkmin2d),
                                                      sizeof(dague_data_t*) );
        datamap = spmtx->datamap_blok;

        cblkN = spmtx->solvmtx->cblktab + cblkmin2d;
        for(cblknumN = cblkmin2d, n = 0;
            cblknumN < cblknbr;
            cblknumN++, n++, cblkN++ )
        {
            if ( !(cblkN->cblktype & CBLK_SPLIT) )
                continue;

            /**
             * Upper Part
             */
            /* if ( ratio == 2 ) { */
            /*     cblkM = spmtx->solvmtx->cblktab + cblkmin2d; */
            /*     for(cblknumM = cblkmin2d, m = 0; */
            /*         cblknumM < cblknumN; */
            /*         cblknumM++, m++, cblkM++ ) */
            /*     { */
            /*     } */
            /* } */

            /**
             * Diagonal block
             */
            ptr    = cblkN->lcoeftab;
            blok   = cblkN->fblokptr;
            size   = blok_rownbr( blok ) * cblk_colnbr( cblkN ) * (size_t)spmtx->typesze;
            offset = blok->coefind * (size_t)spmtx->typesze;
            key2   = n * ld;

            assert(offset == 0);
            dague_data_create( datamap + key2, &spmtx->super, offkey + key2,
                               ptr + offset, size );


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

                while( (blok < lblok) && (blok[0].fcblknm == blok[1].fcblknm) ) {
                    blok++; m++;
                    size += blok_rownbr( blok );
                }
                size *= cblk_colnbr( cblkN )
                    *  (size_t)spmtx->typesze;

                dague_data_create( datamap + key2, &spmtx->super,
                                   offkey + key2,
                                   ptr + offset, size );
                key2 += m * ratio;
            }
        }
    }
    else {
        spmtx->datamap_blok = NULL;
    }
}

void sparse_matrix_destroy( sparse_matrix_desc_t *spmtx )
{
    pastix_int_t i;
    dague_data_t **data;
    int ratio = ( spmtx->mtxtype == PastixGeneral ) ? 2 : 1;

    if ( spmtx->datamap_cblk != NULL ) {
        data = spmtx->datamap_cblk;

        for(i=0; i < spmtx->solvmtx->cblknbr; i++, data+=ratio)
        {
            dague_data_destroy( data[0] );
            if (ratio == 2)
                dague_data_destroy( data[1] );
        }
        free( spmtx->datamap_cblk );
        spmtx->datamap_cblk = NULL;
    }

    if ( spmtx->datamap_blok != NULL ) {
        /* pastix_int_t nbblock2d = (spmtx->solvmtx->cblknbr - spmtx->solvmtx->cblkmin2d); */

        /* nbblock2d = nbblock2d * nbblock2d; */
        /* data = spmtx->datamap_blok; */

        /* for(i=0; i < nbblock2d; i++, data++) */
        /* { */
        /*     dague_data_destroy( *data ); */
        /* } */
        free( spmtx->datamap_blok );
        spmtx->datamap_blok = NULL;
    }

    dague_ddesc_destroy( (dague_ddesc_t*)spmtx );
}
