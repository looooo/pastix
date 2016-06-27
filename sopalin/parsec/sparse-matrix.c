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
    return  ratio * (cblknum + bloknum) + uplo;
}

static inline void
spm_data_key_to_value( dague_data_key_t key,
                       int ratio, const SolverMatrix *solvmtx,
                       int *uplo,
                       pastix_int_t *cblknum,
                       pastix_int_t *bloknum)
{
    dague_data_key_t key2;
    const SolverBlok *blok;

    /* Refers to a block */
    key2 = ratio * solvmtx->cblknbr;
    if ( key >= key2 ) {
        key2 = key - key2;

        *uplo    = key2 % ratio;
        *bloknum = key2 / ratio;
        blok     = solvmtx->bloktab + (*bloknum);
        *cblknum = blok->lcblknm;
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
    pastix_int_t cblknum, bloknum, fbloknum;
    SolverCblk *cblk;

    va_start(ap, mat);
    uplo    = va_arg(ap, int);
    cblknum = va_arg(ap, int);
    bloknum = va_arg(ap, int);
    va_end(ap);

    ratio = spmtx->mtxtype == PastixGeneral ? 2 : 1;
    uplo = uplo ? 1 : 0;
    assert( ratio == 2 || uplo == 0 );

    cblk = spmtx->solvmtx->cblktab + cblknum;
    fbloknum = cblk->fblokptr - spmtx->solvmtx->bloktab;

    if ( (cblk->cblktype & CBLK_SPLIT) && (bloknum != -1) ) {
        return spm_data_key( ratio, uplo, spmtx->solvmtx->cblknbr, fbloknum + bloknum );
    }
    else {
        return spm_data_key( ratio, uplo, cblknum, 0 );
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
    pastix_int_t cblknum, bloknum, fbloknum;
    dague_data_key_t key;
    size_t size, offset, pos;
    dague_data_t **dataptr;

    va_start(ap, mat);
    uplo    = va_arg(ap, int);
    cblknum = va_arg(ap, int);
    bloknum = va_arg(ap, int) - 1;
    va_end(ap);

    ratio = spmtx->mtxtype == PastixGeneral ? 2 : 1;
    uplo = uplo ? 1 : 0;
    assert( ratio == 2 || uplo == 0 );

    pos  = ratio * cblknum + uplo;
    cblk = spmtx->solvmtx->cblktab + cblknum;
    dataptr = spmtx->data_map + pos;
    ptr  = uplo ? cblk->ucoeftab : cblk->lcoeftab;

    if ( cblk->cblktype & CBLK_SPLIT ) {
        /* Return the data for all the cblk to process as 1d */
        if ( bloknum == -1 ) {
            key  = spm_data_key( ratio, uplo, cblknum, 0 );
            size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;
            offset = 0;
        }
        else {
            fbloknum = cblk->fblokptr - spmtx->solvmtx->bloktab;

            /* Return the data for one of the block in the cblk to process as 2d */
            key  = spm_data_key( ratio, uplo, spmtx->solvmtx->cblknbr, fbloknum + bloknum );
            size = blok_rownbr( cblk->fblokptr + bloknum ) * cblk_colnbr( cblk )  * (size_t)spmtx->typesze;
            offset = (cblk->fblokptr + bloknum)->coefind * (size_t)spmtx->typesze;
        }

        /* Extra level of indirection */
        dataptr = ((dague_data_t**)(*dataptr)) + bloknum + 1;
    }
    else {
        /* Return the data for all the cblk to process as 1d */
        key = spm_data_key( ratio, uplo, cblknum, 0 );
        size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;
        offset = 0;
    }
    return dague_data_create( dataptr, mat, key, ptr + offset, size );
}

static dague_data_t *
sparse_matrix_data_of_key(dague_ddesc_t *mat, dague_data_key_t key )
{
    sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat;
    SolverMatrix *solvmtx = spmtx->solvmtx;
    SolverCblk *cblk;
    int uplo, ratio;
    pastix_int_t cblknum, bloknum;
    size_t size, pos, offset;
    dague_data_t **dataptr;
    char *ptr;

    ratio = (spmtx->mtxtype == PastixGeneral) ? 2 : 1;
    spm_data_key_to_value( key, ratio, solvmtx,
                           &uplo, &cblknum, &bloknum );

    cblk = solvmtx->cblktab + cblknum;
    pos  = ratio * cblknum + uplo;
    ptr  = uplo ? cblk->ucoeftab : cblk->lcoeftab;
    dataptr = spmtx->data_map + pos;

    if ( cblk->cblktype & CBLK_SPLIT ) {
        /* Return the data for all the cblk to process as 1d */
        if ( bloknum == -1 ) {
            size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;
            offset = 0;
        }
        else {
            /* Return the data for one of the block in the cblk to process as 2d */
            size = blok_rownbr( cblk->fblokptr + bloknum ) * cblk_colnbr( cblk )  * (size_t)spmtx->typesze;
            offset = (cblk->fblokptr + bloknum)->coefind * (size_t)spmtx->typesze;
        }
        /* Extra level of indirection */
        dataptr = ((dague_data_t**)(*dataptr)) + bloknum + 1;
    }
    else {
        /* Return the data for all the cblk to process as 1d */
        size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;
        offset = 0;
    }
    return dague_data_create( dataptr, mat, key,
                              ptr + offset, size );
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

void sparse_matrix_init( sparse_matrix_desc_t *desc,
                         SolverMatrix *solvmtx,
                         int typesize, int mtxtype,
                         int nodes, int myrank)
{
    dague_ddesc_t *o    = (dague_ddesc_t*)desc;
    SolverCblk    *cblk = solvmtx->cblktab;
    dague_data_t **data;
    pastix_int_t   cblknum, cblknbr;
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

    desc->typesze   = typesize;
    desc->mtxtype   = mtxtype;
    desc->solvmtx   = solvmtx;
    desc->data_map  = (dague_data_t**)calloc( ratio * solvmtx->cblknbr, sizeof(dague_data_t*) );

    cblknbr = solvmtx->cblknbr;
    data = desc->data_map;
    for(cblknum=0; cblknum<cblknbr; cblknum++, cblk++, data += ratio)
    {
        if ( cblk->cblktype & CBLK_SPLIT ) {
            data[0] = calloc( (cblk[1].fblokptr - cblk[0].fblokptr + 1), sizeof(dague_data_t*) );
            if (ratio == 2) {
                data[1] = calloc( (cblk[1].fblokptr - cblk[0].fblokptr + 1), sizeof(dague_data_t*) );
            }
        }
    }
}

void sparse_matrix_destroy( sparse_matrix_desc_t *desc )
{
    if ( desc->data_map != NULL ) {
        dague_data_t **data = desc->data_map;
        dague_data_t **dataptrL;
        dague_data_t **dataptrU;
        SolverCblk    *cblk = desc->solvmtx->cblktab;
        int ratio = ( desc->mtxtype == PastixGeneral ) ? 2 : 1;
        pastix_int_t i, j;

        for(i=0; i < desc->solvmtx->cblknbr; i++, cblk++, data+=ratio)
        {
            if ( cblk->cblktype & CBLK_SPLIT ) {
                dataptrL = (dague_data_t**)data[0];
                dataptrU = (dague_data_t**)data[1];
                for( j=0; j<(cblk[1].fblokptr - cblk[0].fblokptr + 1); j++, dataptrL++, dataptrU++ )
                {
                    dague_data_destroy( *dataptrL );
                    if (ratio == 2)
                        dague_data_destroy( *dataptrU );
                }
                free( data[0] );
                if (ratio == 2)
                    free( data[1] );
            }
            else {
                dague_data_destroy( *data );
                if (ratio == 2)
                    dague_data_destroy( data[1] );
            }
        }
        free( desc->data_map );
        desc->data_map = NULL;
    }
    dague_ddesc_destroy( (dague_ddesc_t*)desc );
}
