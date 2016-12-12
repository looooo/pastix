/**
 *
 * @file coeftab.c
 *
 *  PaStiX factorization routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
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

int (*coeftabDiff[4])(const SolverMatrix*, SolverMatrix*) =
{
    coeftab_sdiff, coeftab_ddiff, coeftab_cdiff, coeftab_zdiff
};

void (*coeftabMemory[4])(SolverMatrix*) =
{
    coeftab_smemory, coeftab_dmemory, coeftab_cmemory, coeftab_zmemory
};

void (*coeftabUncompress[4])(SolverMatrix*) =
{
    coeftab_suncompress, coeftab_duncompress, coeftab_cuncompress, coeftab_zuncompress
};

void (*coeftabCompress[4])(SolverMatrix*) =
{
    coeftab_scompress, coeftab_dcompress, coeftab_ccompress, coeftab_zcompress
};

struct coeftabinit_s {
    const SolverMatrix  *datacode;
    const pastix_bcsc_t *bcsc;
    int fakefillin;
    int factoLU;
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
pcoeftabInit( int rank, void *args )
{
    struct coeftabinit_s *ciargs = (struct coeftabinit_s*)args;
    const SolverMatrix  *datacode = ciargs->datacode;
    const pastix_bcsc_t *bcsc     = ciargs->bcsc;
    int fakefillin = ciargs->fakefillin;
    int factoLU    = ciargs->factoLU;
    pastix_int_t i, itercblk;
    pastix_int_t task;
    void (*initfunc)(const SolverMatrix*,
                     const pastix_bcsc_t*,
                     pastix_int_t,
                     int, int) = NULL;
    void (*dumpfunc)(const SolverMatrix*,
                     const char *) = NULL;

    switch( bcsc->flttype ) {
    case PastixComplex32:
        initfunc = coeftab_cinitcblk;
        dumpfunc = coeftab_cdump;
        break;
    case PastixComplex64:
        initfunc = coeftab_zinitcblk;
        dumpfunc = coeftab_zdump;
        break;
    case PastixFloat:
        initfunc = coeftab_sinitcblk;
        dumpfunc = coeftab_sdump;
        break;
    case PastixDouble:
    case PastixPattern:
    default:
        initfunc = coeftab_dinitcblk;
        dumpfunc = coeftab_ddump;
    }

    for (i=0; i < datacode->ttsknbr[rank]; i++)
    {
        task = datacode->ttsktab[rank][i];
        itercblk = datacode->tasktab[task].cblknum;

        initfunc( datacode, bcsc, itercblk, fakefillin, factoLU );
    }

    (void)dumpfunc;
    //dumpfunc( datacode, "lcoeftab.txt" );
}

void
coeftabInit( const pastix_data_t *pastix_data,
             int fakefillin, int factoLU )
{
    struct coeftabinit_s args;

    args.datacode   = pastix_data->solvmatr;
    args.bcsc       = pastix_data->bcsc;
    args.fakefillin = fakefillin;
    args.factoLU    = factoLU;

    isched_parallel_call( pastix_data->isched, pcoeftabInit, &args );

    //dumpfunc( datacode, "lcoeftab.txt" );
}

void
coeftabExit( SolverMatrix *solvmtx )
{
    pastix_int_t i;

#if defined(PASTIX_WITH_PARSEC)
    {
        if ( solvmtx->parsec_desc != NULL ) {
            sparse_matrix_destroy( solvmtx->parsec_desc );
            free( solvmtx->parsec_desc );
        }
        solvmtx->parsec_desc = NULL;
    }
#endif

    /** Free arrays of solvmtx **/
    if(solvmtx->cblktab)
    {
        for (i = 0; i < solvmtx->cblknbr; i++)
        {
            if (solvmtx->cblktab[i].lcoeftab)
                memFree_null(solvmtx->cblktab[i].lcoeftab);

            if (solvmtx->cblktab[i].ucoeftab)
                memFree_null(solvmtx->cblktab[i].ucoeftab);

            if (! (solvmtx->cblktab[i].cblktype & CBLK_DENSE)) {
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
