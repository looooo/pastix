/**
 *
 * @file cpucblk_zinit.c
 *
 * Precision dependent coeficient array initialization routines.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Esragul Korkmaz
 * @author Tony Delarue
 * @author Alycia Lisito
 * @author Nolan Bredel
 * @date 2023-10-25
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common/common.h"
#include "blend/solver.h"
#include "bcsc/bcsc.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"

/**
 *******************************************************************************
 *
 * @brief Initialize a lrblock structure from a workspace from a specific block to the end of all
 * blocks.
 *
 * The lrblock structure must be allocated before.
 *
 *******************************************************************************
 *
 * @param[inout] cblk
 *          The column block associated to the initialization.
 *
 * @param[in] blok
 *          The block representation associated to the initialization.
 *
 * @param[in] lrblok
 *          The structure blok to initialize. Must be allocated before.
 *
 * @param[in] ws
 *          The workspace associated with the data that will be used for initialize lrblok.
 *
 *******************************************************************************/
void
cpublok_zalloc_lrws( const SolverCblk   *cblk,
                     const SolverBlok   *blok,
                     pastix_lrblock_t   *lrblok,
                     pastix_complex64_t *ws )
{
    SolverBlok  *lblok   = cblk[1].fblokptr;
    pastix_int_t fcblknm = blok->fcblknm;
    pastix_int_t ncols   = cblk_colnbr( cblk );

    /* H then split */
    assert( cblk->cblktype & CBLK_LAYOUT_2D );
    assert( lrblok != NULL );

    for (; (blok < lblok) && (blok->fcblknm == fcblknm); blok++, lrblok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );
        lrblok->rk    = -1;
        lrblok->rkmax = nrows;
        lrblok->u     = ws;
        lrblok->v     = NULL;

        ws += nrows * ncols;
    }
}

/**
 *******************************************************************************
 *
 * @brief Initialize lrblock structure from a workspace for all blocks of the cblk associated.
 *
 * The lrblock structure must be allocated before.
 *
 *******************************************************************************
 *
 * @param[inout] cblk
 *          The column block associated to the initialization.
 *
 * @param[in] lrblok
 *          The structure blok to initialize. Must be allocated before.
 *
 * @param[in] ws
 *          The workspace associated with the data that will be used for initialize lrblok.
 *
 *******************************************************************************/
void
cpucblk_zalloc_lrws( const SolverCblk   *cblk,
                     pastix_lrblock_t   *lrblok,
                     pastix_complex64_t *ws )
{
    pastix_int_t      ncols = cblk_colnbr( cblk );
    const SolverBlok *blok  = cblk[0].fblokptr;
    const SolverBlok *lblok = cblk[1].fblokptr;

    /* H then split */
    assert( cblk->cblktype & CBLK_LAYOUT_2D );
    assert( lrblok != NULL );

    for (; blok<lblok; blok++, lrblok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );
        lrblok->rk    = -1;
        lrblok->rkmax = nrows;
        lrblok->u     = ws;
        lrblok->v     = NULL;

        ws += nrows * ncols;
    }
}

/**
 *******************************************************************************
 *
 * @brief Allocate the cblk structure to store the coefficient
 *
 * When stored in low-rank format, the data pointer in the low-rank structure of
 * each block must be initialized.
 * This routines performs only the allocation and is thread-safe if called in
 * parallel on the Lower and upper part.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] cblk
 *          The column block to allocate.
 *
 * @param[in] rkmax
 *          TODO
 *
 *******************************************************************************/
void
cpucblk_zalloc_lr( pastix_coefside_t  side,
                   SolverCblk        *cblk,
                   int                rkmax )
{
    pastix_int_t      ncols    = cblk_colnbr( cblk );
    pastix_lrblock_t *LRblocks = NULL;
    SolverBlok       *blok     = cblk[0].fblokptr;
    SolverBlok       *lblok    = cblk[1].fblokptr;
    size_t            size     = lblok - blok;

    /* H then split */
    assert( cblk->cblktype & CBLK_LAYOUT_2D );

    LRblocks = blok->LRblock[0];

    if ( LRblocks == NULL ) {
        /* One allocation per cblk */
        LRblocks = malloc( 2 * size * sizeof(pastix_lrblock_t) );
        memset( LRblocks, 0, 2 * size * sizeof(pastix_lrblock_t) );
        if (!pastix_atomic_cas_xxb( &(blok->LRblock[0]), (uint64_t)NULL, (uint64_t)LRblocks, sizeof(void*) )) {
            free( LRblocks );
            LRblocks = blok->LRblock[0];
        }
    }
    assert( LRblocks != NULL );

    for (; blok<lblok; blok++)
    {
        pastix_int_t nrows = blok_rownbr( blok );
        blok->LRblock[0] = LRblocks;
        blok->LRblock[1] = LRblocks + size;

        if ( side != PastixUCoef ) {
            core_zlralloc( nrows, ncols, rkmax, blok->LRblock[0] );
        }

        if ( side != PastixLCoef ) {
            core_zlralloc( nrows, ncols, rkmax, blok->LRblock[1] );
        }
        LRblocks++;
    }

    /* Backup the fact that the cblk has been initialized */
    if ( side != PastixUCoef ) {
        cblk->lcoeftab = (void*)-1;
    }
    if ( side != PastixLCoef ) {
        cblk->ucoeftab = (void*)-1;
    }
}

/**
 *******************************************************************************
 *
 * @brief Allocate the cblk structure to store the coefficient
 *
 * When stored in low-rank format, the data pointer in the low-rank structure of
 * each block must be initialized.
 * This routines performs only the allocation and is thread-safe if called in
 * parallel on the Lower and upper part.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] cblk
 *          The column block to allocate.
 *
 *******************************************************************************/
void
cpucblk_zalloc_fr( pastix_coefside_t  side,
                   SolverCblk        *cblk )
{
    size_t ncols   = cblk_colnbr( cblk );
    size_t coefnbr = cblk->stride * ncols;

    if ( side == PastixLCoef ) {
        assert( cblk->lcoeftab == NULL );
        MALLOC_INTERN( cblk->lcoeftab, coefnbr, pastix_complex64_t );
        memset( cblk->lcoeftab, 0, coefnbr * sizeof(pastix_complex64_t) );
    }
    else {
        assert( cblk->lcoeftab == NULL );
        assert( cblk->ucoeftab == NULL );

        MALLOC_INTERN( cblk->lcoeftab, 2 * coefnbr, pastix_complex64_t );
        memset( cblk->lcoeftab, 0, 2 * coefnbr * sizeof(pastix_complex64_t) );

        cblk->ucoeftab = (pastix_complex64_t *)cblk->lcoeftab + coefnbr;
        assert( cblk->ucoeftab );
    }
    assert( cblk->lcoeftab );
}

/**
 *******************************************************************************
 *
 * @brief Allocate the cblk structure to store the coefficient
 *
 * When stored in low-rank format, the data pointer in the low-rank structure of
 * each block must be initialized.
 * This routines performs only the allocation and is thread-safe if called in
 * parallel on the Lower and upper part.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] cblk
 *          The column block to allocate.
 *
 *******************************************************************************/
void
cpucblk_zalloc( pastix_coefside_t  side,
                SolverCblk        *cblk )
{
    /* Make sure they have the correct values */
    assert( PastixLCoef == 0 );
    assert( PastixUCoef == 1 );

    pastix_cblk_lock( cblk );

    /* Shift to play with bitmasks */
    side += 1;
    if ( cblk->lcoeftab != NULL ) {
        side &= ( ~(PastixLCoef+1) );
    }
    if ( cblk->ucoeftab != NULL ) {
        side &= ( ~(PastixUCoef+1) );
    }
    if ( !side ) {
        pastix_cblk_unlock( cblk );
        return;
    }
    side -= 1;

    if ( cblk->cblktype & CBLK_COMPRESSED ) {
        cpucblk_zalloc_lr( side, cblk, cblk->cblktype & CBLK_FANIN ? 0 : -1 );
    }
    else {
        cpucblk_zalloc_fr( side, cblk );
    }
    pastix_cblk_unlock( cblk );
}

/**
 *******************************************************************************
 *
 * @brief Free the cblk structure that store the coefficient
 *
 * This routines performs the free and is thread-safe if called in
 * parallel on the Lower and upper part.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[inout] cblk
 *          The column block to free.
 *
 *******************************************************************************/
void
cpucblk_zfree( pastix_coefside_t  side,
               SolverCblk        *cblk )
{
    pastix_cblk_lock( cblk );
    if ( (side != PastixUCoef) && (cblk->lcoeftab != NULL) ) {

        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            SolverBlok *blok  = cblk[0].fblokptr;
            SolverBlok *lblok = cblk[1].fblokptr;

            assert( blok->LRblock[0] != NULL );
            for (; blok<lblok; blok++) {
                core_zlrfree(blok->LRblock[0]);
            }

            if ( cblk->lcoeftab != (void*)-1 ) {
                memFree_null( cblk->lcoeftab );
            }
        }
        else {
            memFree_null( cblk->lcoeftab );
        }
        cblk->lcoeftab = NULL;
    }
    if ( (side != PastixLCoef) && (cblk->ucoeftab != NULL) ) {

        if ( cblk->cblktype & CBLK_COMPRESSED ) {
            SolverBlok *blok  = cblk[0].fblokptr;
            SolverBlok *lblok = cblk[1].fblokptr;

            assert( blok->LRblock[1] != NULL );
            for (; blok<lblok; blok++) {
                core_zlrfree(blok->LRblock[1]);
            }
        }
        cblk->ucoeftab = NULL;
    }
    if ( (cblk->cblktype & CBLK_COMPRESSED) &&
         (cblk->lcoeftab == NULL)           &&
         (cblk->ucoeftab == NULL) )
    {
        free( cblk->fblokptr->LRblock[0] );
        cblk->fblokptr->LRblock[0] = NULL;
        cblk->fblokptr->LRblock[1] = NULL;
    }
    pastix_cblk_unlock( cblk );
}

/**
 *******************************************************************************
 *
 * @brief Initialize the full-rank coeftab structure from the internat bcsc.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] solvmtx
 *          PaStiX structure to store numerical data and flags
 *
 * @param[in] bcsc
 *          The internal bcsc structure that hold the graph with permutation
 *          stored by cblk.
 *
 * @param[in] itercblk
 *          The index of the cblk to fill in both bcsc and solvmtx structures.
 *
 *******************************************************************************/
static inline void
cpucblk_zfillin_fr( pastix_coefside_t    side,
                    const SolverMatrix  *solvmtx,
                    const pastix_bcsc_t *bcsc,
                    pastix_int_t         itercblk )
{
    SolverCblk *solvcblk       = solvmtx->cblktab + itercblk;
    const bcsc_cblk_t *csccblk = bcsc->cscftab + solvcblk->bcscnum;
    SolverBlok *solvblok;
    SolverBlok *lsolvblok = (solvcblk+1)->fblokptr;
    pastix_complex64_t *lcoeftab = solvcblk->lcoeftab;
    pastix_complex64_t *ucoeftab = solvcblk->ucoeftab;
    pastix_complex64_t *Lvalues  = bcsc->Lvalues;
    pastix_complex64_t *Uvalues  = bcsc->Uvalues;
    pastix_int_t ldd = solvcblk->stride;
    pastix_int_t itercoltab, iterval, coefindx;
    int is2d = solvcblk->cblktype & CBLK_LAYOUT_2D;

    assert( (side != PastixUCoef) || (ucoeftab != NULL) );

    for (itercoltab=0; itercoltab<csccblk->colnbr; itercoltab++)
    {
        pastix_int_t frow = csccblk->coltab[itercoltab];
        pastix_int_t lrow = csccblk->coltab[itercoltab+1];
        solvblok = solvcblk->fblokptr;
        if ( is2d ) {
            ldd = blok_rownbr( solvblok );
        }

        for (iterval=frow; iterval<lrow; iterval++)
        {
            pastix_int_t rownum = bcsc->rowtab[iterval];
#if !defined(NDEBUG) && defined(PASTIX_DEBUG_DUMP_COEFTAB)
            if ( isnan( (double)Lvalues[iterval] ) || isinf( (double)Lvalues[iterval] ) ) {
                printf( "cpucblk_zfillin_fr: Lvalues not initialised correctly.\n" );
                assert( 0 );
            }
            if ( isnan( (double)Uvalues[iterval] ) || isinf( (double)Uvalues[iterval] ) ) {
                printf( "cpucblk_zfillin_fr: Uvalues not initialised correctly.\n" );
                assert( 0 );
            }
#endif
            /* If values in the lower part of the matrix */
            if (rownum >= (solvcblk->fcolnum+itercoltab))
            {
                while ((solvblok < lsolvblok) &&
                       ((solvblok->lrownum < rownum) ||
                        (solvblok->frownum > rownum)))
                {
                    solvblok++;
                    if ( is2d ) {
                        ldd = blok_rownbr( solvblok );
                    }
                }

                if ( solvblok < lsolvblok )
                {
                    coefindx  = solvblok->coefind;
                    coefindx += rownum - solvblok->frownum; /* Row shift    */
                    coefindx += itercoltab * ldd;           /* Column shift */
                    pastix_cblk_lock( solvcblk );
                    solvblok->iluklvl = 0;
                    pastix_cblk_unlock( solvcblk );

                    if ( side != PastixUCoef ) {
                        lcoeftab[coefindx] = Lvalues[iterval];
                    }

                    if ( (side != PastixLCoef) &&
                         (rownum > (solvcblk->fcolnum + itercoltab)) )
                    {
#if defined(PRECISION_z) || defined(PRECISION_c)
                        if (bcsc->mtxtype == PastixHermitian) {
                            ucoeftab[coefindx] = conj(Uvalues[iterval]);
                        }
                        else
#endif
                        {
                            ucoeftab[coefindx] = Uvalues[iterval];
                        }
                    }
                }
                else {
#if defined(PASTIX_DEBUG_COEFTAB)
                    fprintf(stderr, "cpucblk_zfillin: drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n",
                            (long)solvcblk->fcolnum + itercoltab, (long)itercoltab,
                            (long)rownum, (long)iterval, (long)itercblk,
                            (long)solvcblk->fcolnum, (long)solvcblk->lcolnum );
#endif
                }
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Initialize the low-rank coeftab structure from the internal bcsc.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] solvmtx
 *          PaStiX structure to store numerical data and flags
 *
 * @param[in] bcsc
 *          The internal bcsc structure that hold the graph with permutation
 *          stored by cblk.
 *
 * @param[in] itercblk
 *          The index of the cblk to fill in both bcsc and solvmtx structures.
 *
 *******************************************************************************/
static inline void
cpucblk_zfillin_lr( pastix_coefside_t    side,
                    const SolverMatrix  *solvmtx,
                    const pastix_bcsc_t *bcsc,
                    pastix_int_t         itercblk )
{
    SolverCblk *solvcblk       = solvmtx->cblktab + itercblk;
    const bcsc_cblk_t *csccblk = bcsc->cscftab + solvcblk->bcscnum;
    SolverBlok *solvblok;
    SolverBlok *lsolvblok = (solvcblk+1)->fblokptr;
    pastix_complex64_t *lcoeftab, *ucoeftab;
    pastix_complex64_t *Lvalues  = bcsc->Lvalues;
    pastix_complex64_t *Uvalues  = bcsc->Uvalues;
    pastix_int_t itercoltab, iterval, coefindx, ldd;

    assert( solvcblk->cblktype & CBLK_LAYOUT_2D );

    for (itercoltab=0; itercoltab<csccblk->colnbr; itercoltab++)
    {
        pastix_int_t frow = csccblk->coltab[itercoltab];
        pastix_int_t lrow = csccblk->coltab[itercoltab+1];

        solvblok = solvcblk->fblokptr;
        ldd = blok_rownbr( solvblok );
        lcoeftab = (pastix_complex64_t*)(solvblok->LRblock[0]->u);
        ucoeftab = (pastix_complex64_t*)(solvblok->LRblock[1]->u);

        for (iterval=frow; iterval<lrow; iterval++)
        {
            pastix_int_t rownum = bcsc->rowtab[iterval];

#if !defined(NDEBUG)
            if ( isnan( (double)Lvalues[iterval] ) || isinf( (double)Lvalues[iterval] ) ) {
                printf( "cpucblk_zfillin_lr: Lvalues not initialised correctly.\n" );
                assert( 0 );
            }
            if ( isnan( (double)Uvalues[iterval] ) || isinf( (double)Uvalues[iterval] ) ) {
                printf( "cpucblk_zfillin_lr: Uvalues not initialised correctly.\n" );
                assert( 0 );
            }
#endif
            /* If values in the lower part of the matrix */
            if (rownum >= (solvcblk->fcolnum+itercoltab))
            {
                while ((solvblok < lsolvblok) &&
                       ((solvblok->lrownum < rownum) ||
                        (solvblok->frownum > rownum)))
                {
                    solvblok++;
                    ldd = blok_rownbr( solvblok );
                    lcoeftab = (pastix_complex64_t*)(solvblok->LRblock[0]->u);
                    ucoeftab = (pastix_complex64_t*)(solvblok->LRblock[1]->u);
                }

                if ( solvblok < lsolvblok )
                {
                    coefindx  = rownum - solvblok->frownum; /* Row shift    */
                    coefindx += itercoltab * ldd;           /* Column shift */
                    pastix_cblk_lock( solvcblk );
                    solvblok->iluklvl = 0;
                    pastix_cblk_unlock( solvcblk );

                    if ( side != PastixUCoef ) {
                        lcoeftab[coefindx] = Lvalues[iterval];
                    }

                    if ( (side != PastixLCoef) &&
                         (rownum > (solvcblk->fcolnum + itercoltab)) )
                    {
#if defined(PRECISION_z) || defined(PRECISION_c)
                        if (bcsc->mtxtype == PastixHermitian)
                            ucoeftab[coefindx] = conj(Uvalues[iterval]);
                        else
#endif
                            ucoeftab[coefindx] = Uvalues[iterval];
                    }
                }
                else {
#if defined(PASTIX_DEBUG_COEFTAB)
                    fprintf(stderr, "cpucblk_zfillin: drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n",
                            (long)solvcblk->fcolnum + itercoltab, (long)itercoltab,
                            (long)rownum, (long)iterval, (long)itercblk,
                            (long)solvcblk->fcolnum, (long)solvcblk->lcolnum );
#endif
                }
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Initialize the coeftab structure from the internal bcsc.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Define which side of the matrix must be initialized.
 *          @arg PastixLCoef if lower part only
 *          @arg PastixUCoef if upper part only
 *          @arg PastixLUCoef if both sides.
 *
 * @param[in] solvmtx
 *          PaStiX structure to store numerical data and flags
 *
 * @param[in] bcsc
 *          The internal bcsc structure that hold the graph with permutation
 *          stored by cblk.
 *
 * @param[in] itercblk
 *          The index of the cblk to fill in both bcsc and solvmtx structures.
 *
 *******************************************************************************/
void
cpucblk_zfillin( pastix_coefside_t    side,
                 const SolverMatrix  *solvmtx,
                 const pastix_bcsc_t *bcsc,
                 pastix_int_t         itercblk )
{
    if ( (solvmtx->cblktab + itercblk)->cblktype & CBLK_COMPRESSED ) {
        cpucblk_zfillin_lr( side, solvmtx, bcsc, itercblk );
    }
    else {
        cpucblk_zfillin_fr( side, solvmtx, bcsc, itercblk );
    }
}
