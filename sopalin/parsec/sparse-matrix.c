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

static inline uint32_t
spm_data_key( const SolverMatrix *solvmtx,
              int cblknum, int uplo )
{
    return  uplo * solvmtx->cblknbr + cblknum;
}


static uint32_t
sparse_matrix_data_key(dague_ddesc_t *mat, ... )
{
    va_list ap;
    sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat;
    int uplo, cblknum;

    va_start(ap, mat);
    uplo    = va_arg(ap, int);
    cblknum = va_arg(ap, int);
    va_end(ap);

    return spm_data_key( spmtx->solvmtx, cblknum, uplo );
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
    va_list ap;
    int uplo, cblknum, bloknum, pos;
    dague_data_key_t key;
    size_t size;

    va_start(ap, mat);
    uplo    = va_arg(ap, int);
    cblknum = va_arg(ap, int);
    bloknum = va_arg(ap, int) - 1;
    va_end(ap);

    cblk = spmtx->solvmtx->cblktab + cblknum;
    key  = spm_data_key( spmtx->solvmtx, cblknum, (uplo ? 1 : 0) );
    pos  = key;
    size = (size_t)cblk->stride * (size_t)cblk_colnbr( cblk ) * (size_t)spmtx->typesze;

    assert(bloknum == -1);
    return dague_data_create( spmtx->data_map + pos, mat, key,
                              (uplo == 1) ? cblk->ucoeftab : cblk->lcoeftab,
                              size );
}

static dague_data_t *
sparse_matrix_data_of_key(dague_ddesc_t *mat, dague_data_key_t key )
{
    sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat;
    SolverMatrix *solvmtx = spmtx->solvmtx;
    SolverCblk *cblk;
    int cblknbr = solvmtx->cblknbr;
    int uplo, cblknum, pos;
    size_t size;

    uplo = ( key >= (dague_data_key_t)cblknbr ) ? PastixUpper : PastixLower;
    cblknum = key % cblknbr;

    cblk = solvmtx->cblktab + cblknum;

    pos  = key;
    size = cblk->stride * (cblk->lcolnum - cblk->fcolnum + 1) * spmtx->typesze;

    return dague_data_create( spmtx->data_map + pos, mat, key,
                              (uplo == 1) ? cblk->ucoeftab : cblk->lcoeftab,
                              size );
}

#ifdef DAGUE_PROF_TRACE
static int sparse_matrix_key_to_string(dague_ddesc_t *mat, uint32_t datakey, char *buffer, uint32_t buffer_size)
{
    sparse_matrix_desc_t *spmtx = (sparse_matrix_desc_t*)mat;
    pastix_int_t uplo, cblknum, cblknbr;
    int res;

    cblknbr = spmtx->solvmtx->cblknbr;
    cblknum = (pastix_int_t)datakey % cblknbr;
    uplo    = (pastix_int_t)datakey / cblknbr;

    res = snprintf(buffer, buffer_size, "(%ld, %ld)",
                   (long int)uplo,
                   (long int)cblknum);
    if (res < 0)
    {
        printf("error in key_to_string for tile (%ld, %ld) key: %u\n",
               (long int)uplo, (long int)cblknum, datakey);
    }
    return res;
}
#endif

void sparse_matrix_init( sparse_matrix_desc_t *desc,
                         SolverMatrix *solvmtx,
                         int typesize, int mtxtype,
                         int nodes, int myrank)
{
    dague_ddesc_t *o = (dague_ddesc_t*)desc;
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
}

void sparse_matrix_destroy( sparse_matrix_desc_t *desc )
{
    if ( desc->data_map != NULL ) {
        dague_data_t **data = desc->data_map;
        int ratio = ( desc->mtxtype == PastixGeneral ) ? 2 : 1;
        int i;

        for(i=0; i<ratio*desc->solvmtx->cblknbr; i++, data++)
        {
            dague_data_destroy( *data );
        }

        free( desc->data_map );
        desc->data_map = NULL;
    }
    dague_ddesc_destroy( (dague_ddesc_t*)desc );
}
