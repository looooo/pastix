/**
 *
 * @file coeftab.c
 *
 * PaStiX coefficient array initialization and free routines.
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
 **/
#include "common.h"
#include "bcsc.h"
#include "isched.h"
#include "solver.h"
#include "coeftab.h"
#include "pastix_zcores.h"
#include "pastix_ccores.h"
#include "pastix_dcores.h"
#include "pastix_scores.h"

#if defined(PASTIX_WITH_PARSEC)
#include "sopalin/parsec/pastix_parsec.h"
#endif

#if defined(PASTIX_WITH_STARPU)
#include "sopalin/starpu/pastix_starpu.h"
#endif

coeftab_fct_diff_t coeftabDiff[4] =
{
    coeftab_sdiff, coeftab_ddiff, coeftab_cdiff, coeftab_zdiff
};

coeftab_fct_memory_t coeftabMemory[4] =
{
    coeftab_smemory, coeftab_dmemory, coeftab_cmemory, coeftab_zmemory
};

coeftab_fct_uncompress_t coeftabUncompress[4] =
{
    coeftab_suncompress, coeftab_duncompress, coeftab_cuncompress, coeftab_zuncompress
};

coeftab_fct_compress_t coeftabCompress[4] =
{
    coeftab_scompress, coeftab_dcompress, coeftab_ccompress, coeftab_zcompress
};

struct coeftabinit_s {
    const SolverMatrix  *datacode;
    const pastix_bcsc_t *bcsc;
    char               **dirtemp;
    pastix_coefside_t    side;
};

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
pcoeftabInit( isched_thread_t *ctx, void *args )
{
    struct coeftabinit_s *ciargs   = (struct coeftabinit_s*)args;
    const SolverMatrix   *datacode = ciargs->datacode;
    const pastix_bcsc_t  *bcsc     = ciargs->bcsc;
    char                **dirtemp  = ciargs->dirtemp;
    pastix_coefside_t     side     = ciargs->side;
    pastix_int_t i, itercblk;
    pastix_int_t task;
    int rank = ctx->rank;

    void (*initfunc)( pastix_coefside_t, const SolverMatrix*,
                      const pastix_bcsc_t*, pastix_int_t, char **) = NULL;

    switch( bcsc->flttype ) {
    case PastixComplex32:
        initfunc = cpucblk_cinit;
        break;
    case PastixComplex64:
        initfunc = cpucblk_zinit;
        break;
    case PastixFloat:
        initfunc = cpucblk_sinit;
        break;
    case PastixDouble:
    case PastixPattern:
    default:
        initfunc = cpucblk_dinit;
    }

    for (i=0; i < datacode->ttsknbr[rank]; i++)
    {
        task = datacode->ttsktab[rank][i];
        itercblk = datacode->tasktab[task].cblknum;

        /* Init as full rank */
        initfunc( side, datacode, bcsc, itercblk, dirtemp );
    }
}

void
coeftabInit( pastix_data_t    *pastix_data,
             pastix_coefside_t side )
{
    struct coeftabinit_s args;

    args.datacode   = pastix_data->solvmatr;
    args.bcsc       = pastix_data->bcsc;
    args.side       = side;
    args.dirtemp    = &(pastix_data->dirtemp);

#if defined(PASTIX_DEBUG_DUMP_COEFTAB)
    /* Make sure dirtemp is initialized before calling it with multiple threads */
    pastix_gendirtemp( args.dirtemp );
#endif

    isched_parallel_call( pastix_data->isched, pcoeftabInit, &args );
}

void
coeftabExit( SolverMatrix *solvmtx )
{
    pastix_int_t i;

#if defined(PASTIX_WITH_PARSEC)
    {
        if ( solvmtx->parsec_desc != NULL ) {
            parsec_sparse_matrix_destroy( solvmtx->parsec_desc );
            free( solvmtx->parsec_desc );
        }
        solvmtx->parsec_desc = NULL;
    }
#endif
#if defined(PASTIX_WITH_STARPU)
    {
        if ( solvmtx->starpu_desc != NULL ) {
            starpu_sparse_matrix_destroy( solvmtx->starpu_desc );
            free( solvmtx->starpu_desc );
        }
        solvmtx->starpu_desc = NULL;
    }
#endif

    /* Free arrays of solvmtx */
    if(solvmtx->cblktab)
    {
        for (i = 0; i < solvmtx->cblknbr; i++)
        {
            if (solvmtx->cblktab[i].lcoeftab)
                memFree_null(solvmtx->cblktab[i].lcoeftab);

            if (solvmtx->cblktab[i].ucoeftab)
                memFree_null(solvmtx->cblktab[i].ucoeftab);

            if (solvmtx->cblktab[i].cblktype & CBLK_COMPRESSED) {
                SolverBlok *blok  = solvmtx->cblktab[i].fblokptr;
                SolverBlok *lblok = solvmtx->cblktab[i+1].fblokptr;

                for (; blok<lblok; blok++) {
                    core_zlrfree(blok->LRblock);
                    if (solvmtx->factotype == PastixFactLU)
                        core_zlrfree(blok->LRblock+1);
                }
                free(solvmtx->cblktab[i].fblokptr->LRblock);
            }
        }
    }
}
