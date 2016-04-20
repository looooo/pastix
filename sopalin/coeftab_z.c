/**
 *
 * @file coeftab_z.c
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
 * @precisions normal z -> s d c
 *
 **/
#define _GNU_SOURCE
#include "common.h"
#include "solver.h"
#include "bcsc.h"
#include <lapacke.h>
#include "pastix_zcores.h"

void
coeftab_zdumpcblk( const SolverCblk *cblk,
                   void *array,
                   FILE *stream );

/* Section: Functions */
void
coeftab_zcompress_one( SolverCblk *cblk,
                       double      tol )
{
    pastix_lrblock_t   *LRblocks;
    SolverBlok         *blok     = cblk[0].fblokptr;
    SolverBlok         *lblok    = cblk[1].fblokptr;
    pastix_complex64_t *lcoeftab = cblk->lcoeftab;
    pastix_complex64_t *ucoeftab = cblk->ucoeftab;
    pastix_int_t        ncols    = cblk_colnbr( cblk );
    pastix_int_t        gainL    = ncols * cblk->stride;
    pastix_int_t        gainU    = ncols * cblk->stride;
    int factoLU = (cblk->ucoeftab == NULL) ? 0 : 1;

    /* One allocation per cblk */
    LRblocks = malloc( (factoLU+1) * (lblok - blok) * sizeof(pastix_lrblock_t) );

    /**
     * Diagonal block (Not compressed)
     */
    LRblocks->rk    = -1;
    LRblocks->rkmax = -1;
    LRblocks->u     = malloc( ncols * ncols * sizeof(pastix_complex64_t) );
    LRblocks->v     = NULL;
    LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', ncols, ncols,
                         lcoeftab, cblk->stride, LRblocks->u, ncols );
    blok->coefL_u_LR = LRblocks->u;
    blok->coefL_v_LR = LRblocks->v;
    blok->rankL      = LRblocks->rk;
    blok->LRblock = LRblocks; LRblocks++;

    gainL -= ncols * ncols;

    if (factoLU) {
        LRblocks->rk    = -1;
        LRblocks->rkmax = -1;
        LRblocks->u     = malloc( ncols * ncols * sizeof(pastix_complex64_t) );
        LRblocks->v     = NULL;
        LAPACKE_zlacpy_work( LAPACK_COL_MAJOR, 'A', ncols, ncols,
                             ucoeftab, cblk->stride, LRblocks->u, ncols );
        blok->coefU_u_LR = LRblocks->u;
        blok->coefU_v_LR = LRblocks->v;
        blok->rankU      = LRblocks->rk;
        LRblocks++;

        gainU -= ncols * ncols;
    }

    for (blok++; blok<lblok; blok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );

        blok->LRblock = LRblocks;
        core_zge2lr( tol, nrows, ncols,
                     lcoeftab + blok->coefind, cblk->stride,
                     blok->LRblock );
        blok->coefL_u_LR = LRblocks->u;
        blok->coefL_v_LR = LRblocks->v;
        blok->rankL      = LRblocks->rk;
        LRblocks++;

        gainL -= ((nrows+ncols) * blok->LRblock[0].rk);

        if (factoLU) {
            core_zge2lr( tol, nrows, ncols,
                         ucoeftab + blok->coefind, cblk->stride,
                         blok->LRblock+1 );
            blok->coefU_u_LR = LRblocks->u;
            blok->coefU_v_LR = LRblocks->v;
            blok->rankU      = LRblocks->rk;
            LRblocks++;

            gainU -= ((nrows+ncols) * blok->LRblock[1].rk);
        }
    }

    /**
     * Free the dense version
     */
    free(cblk->lcoeftab); cblk->lcoeftab = NULL;
    //gain_L += gainL * sizeof(pastix_complex64_t) * 1.e-6;
    if (cblk->ucoeftab) {
        free(cblk->ucoeftab); cblk->ucoeftab = NULL;
        //gain_U += gainU * sizeof(pastix_complex64_t) * 1.e-6;
    }
}

void
coeftab_zuncompress_one( SolverCblk *cblk, int factoLU )
{
    SolverBlok *blok  = cblk[0].fblokptr;
    SolverBlok *lblok = cblk[1].fblokptr;

    pastix_int_t ncols = cblk_colnbr( cblk );
    pastix_complex64_t *lcoeftab = NULL;
    pastix_complex64_t *ucoeftab = NULL;

    /* One allocation per cblk */
    assert( cblk->lcoeftab == NULL );
    lcoeftab = malloc( cblk->stride * ncols * sizeof(pastix_complex64_t) );

    if ( factoLU ) {
        assert( cblk->ucoeftab == NULL );
        ucoeftab = malloc( cblk->stride * ncols * sizeof(pastix_complex64_t) );
    }

    for (; blok<lblok; blok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );

        core_zlr2ge( nrows, ncols,
                     blok->LRblock,
                     lcoeftab + blok->coefind, cblk->stride );

        if (factoLU) {
            core_zlr2ge( nrows, ncols,
                         blok->LRblock+1,
                         ucoeftab + blok->coefind, cblk->stride );
        }
    }

    cblk->lcoeftab = lcoeftab;
    cblk->ucoeftab = ucoeftab;

    /**
     * Free all the LRblock structures associated to the cblk
     */
    free(cblk->fblokptr->LRblock);
}


void
coeftab_zuncompress( SolverMatrix *solvmtx )
{
    SolverCblk *cblk  = solvmtx->cblktab;
    pastix_int_t cblknum;

    for(cblknum=0; cblknum<solvmtx->cblknbr; cblknum++, cblk++) {
        if (!(cblk->cblktype & CBLK_DENSE)) {
            coeftab_zuncompress_one( cblk, 1 );
            cblk->cblktype |= CBLK_DENSE;
        }
    }
}

void
coeftab_zffbcsc( const SolverMatrix  *solvmtx,
                 const pastix_bcsc_t *bcsc,
                 pastix_int_t         itercblk )
{
    const bcsc_format_t *csccblk = bcsc->cscftab + itercblk;
    SolverCblk *solvcblk = solvmtx->cblktab + itercblk;
    SolverBlok *solvblok;
    SolverBlok *solvblok2 = (solvcblk+1)->fblokptr;
    pastix_complex64_t *lcoeftab = solvcblk->lcoeftab;
    //pastix_complex64_t *dcoeftab = solvcblk->dcoeftab;
    pastix_complex64_t *ucoeftab = solvcblk->ucoeftab;
    pastix_complex64_t *Lvalues = bcsc->Lvalues;
    pastix_complex64_t *Uvalues = bcsc->Uvalues;
    pastix_int_t itercoltab, iterval, coefindx;
    //    pastix_int_t ncols = solvcblk->lcolnum - solvcblk->fcolnum + 1;

    for (itercoltab=0; itercoltab<csccblk->colnbr; itercoltab++)
    {
        pastix_int_t frow = csccblk->coltab[itercoltab];
        pastix_int_t lrow = csccblk->coltab[itercoltab+1];
        solvblok = solvcblk->fblokptr;

        for (iterval=frow; iterval<lrow; iterval++)
        {
            pastix_int_t rownum = bcsc->rowtab[iterval];

            /* If values in the lower part of the matrix */
            if (rownum >= (solvcblk->fcolnum+itercoltab))
            {
                while ((solvblok < solvblok2) &&
                       ((solvblok->lrownum < rownum) ||
                        (solvblok->frownum > rownum)))
                {
                    solvblok++;
                }

                if ( solvblok < solvblok2 )
                {
                    coefindx  = solvblok->coefind;
                    coefindx += rownum - solvblok->frownum;
                    coefindx += solvcblk->stride * itercoltab;

                    lcoeftab[coefindx] = Lvalues[iterval];

                    if ( (ucoeftab != NULL) &&
                         (rownum > (solvcblk->fcolnum + itercoltab)) )
                    {
#if defined(PRECISION_z) || defined(PRECISION_c)
                        if (bcsc->mtxtype == PastixHermitian)
                            ucoeftab[coefindx] = conj(Uvalues[iterval]);
                        else
#endif
                            ucoeftab[coefindx] = Uvalues[iterval];
                    }
/*                     pastix_int_t i, j; */
/*                     i = solvblok->coefind + rownum - solvblok->frownum; */
/*                     j = itercoltab; */

/*                     if (i < ncols){ */
/*                         dcoeftab[i*ncols + j] = Uvalues[iterval]; */
/*                         dcoeftab[j*ncols + i] = Lvalues[iterval]; */
/*                     } */
/*                     else{ */
/*                         coefindx  = solvblok->coefind; */
/*                         coefindx += rownum - solvblok->frownum; */
/*                         coefindx += solvcblk->stride * itercoltab; */
/*                         lcoeftab[coefindx] = Lvalues[iterval]; */

/*                         if ( (ucoeftab != NULL) && */
/*                              (rownum > (solvcblk->fcolnum + itercoltab)) ) */
/*                         { */
/* #if defined(PRECISION_z) || defined(PRECISION_c) */
/*                             if (bcsc->mtxtype == PastixHermitian) */
/*                                 ucoeftab[coefindx] = conj(Uvalues[iterval]); */
/*                             else */
/* #endif */
/*                                 ucoeftab[coefindx] = Uvalues[iterval]; */
/*                         } */
/*                     } */
                }
                else {
                    /* printf("ILU: csc2solv drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n", */
                    /*        (long)datacode->cblktab[itercblk].fcolnum+ */
                    /*        (long)itercoltab,(long)itercoltab, */
                    /*        (long)CSC_ROW(cscmtx,iterval),(long)iterval, */
                    /*        (long)itercblk, */
                    /*        (long)datacode->cblktab[itercblk].fcolnum, */
                    /*        (long)datacode->cblktab[itercblk].lcolnum); */
                }
            }
        }
    }
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
coeftab_zinitcblk( const SolverMatrix  *solvmtx,
                   const pastix_bcsc_t *bcsc,
                   pastix_int_t itercblk,
                   int fakefillin, int factoLU )
{
    SolverCblk *cblk = solvmtx->cblktab + itercblk;
    pastix_int_t coefnbr = cblk->stride * cblk_colnbr( cblk );
    pastix_int_t j;

    /* If not NULL, allocated to store the shur complement for exemple */
    assert( cblk->lcoeftab == NULL );

    MALLOC_INTERN( cblk->lcoeftab, coefnbr, pastix_complex64_t );
    memset( cblk->lcoeftab, 0, coefnbr * sizeof(pastix_complex64_t) );

    if ( factoLU ) {
        MALLOC_INTERN( cblk->dcoeftab, cblk_colnbr( cblk ) * cblk_colnbr( cblk ), pastix_complex64_t );
        memset( cblk->dcoeftab, 0, cblk_colnbr( cblk ) * cblk_colnbr( cblk ) * sizeof(pastix_complex64_t) );

        /* Extra diagonal block for low-rank updates */
        MALLOC_INTERN( cblk->ucoeftab, coefnbr, pastix_complex64_t );
        memset( cblk->ucoeftab, 0, coefnbr * sizeof(pastix_complex64_t) );

    }
    else {
        cblk->dcoeftab = NULL;
        cblk->ucoeftab = NULL;
    }

    /**
     * Fake initialization of the bloc column such that:
     *       - L is filled with 1.
     *       - U is filled with 2.
     *       - The diagonal is made dominant
     */
    if ( fakefillin ) {
        pastix_complex64_t *L = cblk->lcoeftab;
        pastix_complex64_t *U = cblk->ucoeftab;

        for (j=0; j<coefnbr; j++, L++)
        {
            *L = (pastix_complex64_t)1.;
        }

        if ( factoLU ) {
            for (j=0; j<coefnbr; j++, U++)
            {
                *U = (pastix_complex64_t)2.;
            }
        }

        /* for (j=0; j<size; itercol++) */
        /* { */
        /*     /\* On s'assure que la matrice est diagonale dominante *\/ */
        /*     SOLV_COEFTAB(itercblk)[index+itercol*stride+itercol] = (pastix_complex64_t) (UPDOWN_GNODENBR*UPDOWN_GNODENBR); */
        /* } */
    }
    else
    {
        coeftab_zffbcsc( solvmtx, bcsc, itercblk );
    }

#if defined(PASTIX_DUMP_COEFTAB)
    {
        FILE *f;
        char *filename;

        asprintf( &filename, "Lcblk%05ld.txt", itercblk );
        f = fopen( filename, "w" );
        coeftab_zdumpcblk( cblk, cblk->lcoeftab, f );
        fclose( f );
        free( filename );

        if ( cblk->ucoeftab ) {
            asprintf( &filename, "Ucblk%05ld.txt", itercblk );
            f = fopen( filename, "w" );
            coeftab_zdumpcblk( cblk, cblk->ucoeftab, f );
            fclose( f );
            free( filename );
        }
    }
#endif /* defined(PASTIX_DUMP_COEFTAB) */

    /**
     * Try to compress the cblk if needs to be compressed
     * TODO: change the criteria based on the level in the tree
     */
    {
        /* TODO: cleanup to pass that as arguments */
        int compress_size = 80;
        char  *tolerance = getenv("TOLERANCE");
        double tol = atof(tolerance);

        if (cblk_colnbr( cblk ) >= compress_size)
        {
            fprintf(stderr, "Try to compress a block\n");
            coeftab_zcompress_one( cblk, tol );
            cblk->cblktype &= ~(CBLK_DENSE);
        }
    }
}

/*
 * Function: z_dump3
 *
 * Prints z_solver matrix informations, in (i,j,v) format, in a file,
 * for LLt or LDLt decomposition.
 *
 * Parameters:
 *   datacode - z_SolverMatrix.
 *   stream   - FILE * opened in write mode.
 */
void
coeftab_zdumpcblk( const SolverCblk *cblk,
                   void *array,
                   FILE *stream )
{
    pastix_complex64_t *coeftab = (pastix_complex64_t*)array;
    SolverBlok *blok;
    pastix_int_t itercol;
    pastix_int_t iterrow;
    pastix_int_t coefindx;

    for (itercol  = cblk->fcolnum;
         itercol <= cblk->lcolnum;
         itercol++)
    {
        /* Diagonal Block */
        blok     = cblk->fblokptr;
        coefindx = blok->coefind;
        coefindx += (itercol - cblk->fcolnum) * cblk->stride;

        for (iterrow  = blok->frownum;
             iterrow <= blok->lrownum;
             iterrow++, coefindx++)
        {
            if ((cabs( coeftab[coefindx] ) > 0.) &&
                (itercol <= iterrow))
            {
#if defined(PRECISION_z) || defined(PRECISION_c)
                fprintf(stream, "%ld %ld (%13e,%13e)\n",
                        (long)itercol, (long)iterrow,
                        creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                fprintf(stream, "%ld %ld %13e\n",
                        (long)itercol, (long)iterrow,
                        coeftab[coefindx]);
#endif
            }
        }

        /* Off diagonal blocks */
        blok++;
        while( blok < (cblk+1)->fblokptr )
        {
            coefindx  = blok->coefind;
            coefindx += (itercol - cblk->fcolnum) * cblk->stride;

            for (iterrow  = blok->frownum;
                 iterrow <= blok->lrownum;
                 iterrow++, coefindx++)
            {
                if (cabs( coeftab[coefindx]) > 0.)
                {
#if defined(PRECISION_z) || defined(PRECISION_c)
                    fprintf(stream, "%ld %ld (%13e,%13e)\n",
                            (long)itercol, (long)iterrow,
                            creal(coeftab[coefindx]), cimag(coeftab[coefindx]));
#else
                    fprintf(stream, "%ld %ld %13e\n",
                            (long)itercol, (long)iterrow,
                            coeftab[coefindx]);
#endif
                }
            }
            blok++;
        }
    }
}

void
coeftab_zdump( const SolverMatrix *solvmtx,
               const char   *filename )
{
    SolverCblk *cblk = solvmtx->cblktab;
    pastix_int_t itercblk;
    FILE *stream = fopen( filename, "w" );

    for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++, cblk++)
    {
        coeftab_zdumpcblk( cblk, cblk->lcoeftab, stream );
        if ( NULL != cblk->ucoeftab )
            coeftab_zdumpcblk( cblk, cblk->lcoeftab, stream );
    }

    fclose( stream );
}
