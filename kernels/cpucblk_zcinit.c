/**
 *
 * @file cpucblk_zcinit.c
 *
 * Mixed-Precision dependent coeficient array initialization routines.
 *
 * @copyright 2015-2024 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.2
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @author Esragul Korkmaz
 * @author Tony Delarue
 * @author Brieuc Nicolas
 * @date 2023-07-21
 *
 * @precisions mixed zc -> ds
 *
 **/
#include "common/common.h"
#include "blend/solver.h"
#include "bcsc/bcsc.h"
#include "pastix_zcores.h"
#include "pastix_zlrcores.h"
#include <lapacke.h>

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#if defined(PRECISION_zc)
#define cpucblk_zccheck_overflow( value, overflow )                                                       \
{                                                                                                         \
    double rvalue = fabs( creal( value ) );                                                               \
    double ivalue = fabs( cimag( value ) );                                                               \
    if ( pastix_unlikely( (rvalue > overflow) || (ivalue > overflow) ) ) {                                \
        pastix_print_warning( "cpucblk_zccheck_overflow, Incorrect value overflow for mixed precision" ); \
        return 1;                                                                                         \
    }                                                                                                     \
}
#else
#define cpucblk_dscheck_overflow( value, overflow )                                                       \
{                                                                                                         \
    double valabs = fabs( value );                                                                        \
    if ( pastix_unlikely( valabs > overflow ) ) {                                                         \
        pastix_print_warning( "cpucblk_dscheck_overflow, Incorrect value overflow for mixed precision" ); \
        return 1;                                                                                         \
    }                                                                                                     \
}
#endif
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

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
 *******************************************************************************
 *
 * @return 0 on success, 1 in case of overflow during initialization.
 *
 *******************************************************************************/
static inline int
cpucblk_zcfillin_fr( pastix_coefside_t    side,
                     const SolverMatrix  *solvmtx,
                     const pastix_bcsc_t *bcsc,
                     pastix_int_t         itercblk )
{
    SolverCblk *solvcblk       = solvmtx->cblktab + itercblk;
    const bcsc_cblk_t *csccblk = bcsc->cscftab + solvcblk->bcscnum;
    SolverBlok *solvblok;
    SolverBlok *lsolvblok = (solvcblk+1)->fblokptr;
    pastix_complex32_t *lcoeftab = solvcblk->lcoeftab;
    pastix_complex32_t *ucoeftab = solvcblk->ucoeftab;
    pastix_complex64_t *Lvalues  = bcsc->Lvalues;
    pastix_complex64_t *Uvalues  = bcsc->Uvalues;
    pastix_int_t ldd = solvcblk->stride;
    pastix_int_t itercoltab, iterval, coefindx;
    int is2d = solvcblk->cblktype & CBLK_LAYOUT_2D;
    double overflow = LAPACKE_slamch( 'o' );

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
                        cpucblk_zccheck_overflow( Lvalues[iterval], overflow )
                        lcoeftab[coefindx] = (pastix_complex32_t) Lvalues[iterval];
                    }

                    if ( (side != PastixLCoef) &&
                         (rownum > (solvcblk->fcolnum + itercoltab)) )
                    {
                        cpucblk_zccheck_overflow( Uvalues[iterval], overflow );
#if defined(PRECISION_zc)
                        if (bcsc->mtxtype == PastixHermitian) {
                            ucoeftab[coefindx] = (pastix_complex32_t) conj(Uvalues[iterval]);
                        }
                        else
#endif
                        {
                            ucoeftab[coefindx] = (pastix_complex32_t) Uvalues[iterval];
                        }
                    }
                }
#if defined(PASTIX_DEBUG_COEFTAB)
                else
                {
                    fprintf( stderr, "cpucblk_zfillin: drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n",
                             (long)solvcblk->fcolnum + itercoltab, (long)itercoltab,
                             (long)rownum, (long)iterval, (long)itercblk,
                             (long)solvcblk->fcolnum, (long)solvcblk->lcolnum );
                }
#endif
            }
        }
    }
    (void)overflow;
    return 0;
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
 *******************************************************************************
 *
 * @return 0 on success, 1 in case of overflow during initialization.
 *
 *******************************************************************************/
static inline int
cpucblk_zcfillin_lr( pastix_coefside_t    side,
                     const SolverMatrix  *solvmtx,
                     const pastix_bcsc_t *bcsc,
                     pastix_int_t         itercblk )
{
    SolverCblk *solvcblk       = solvmtx->cblktab + itercblk;
    const bcsc_cblk_t *csccblk = bcsc->cscftab + solvcblk->bcscnum;
    SolverBlok *solvblok;
    SolverBlok *lsolvblok = (solvcblk+1)->fblokptr;
    pastix_complex32_t *lcoeftab, *ucoeftab;
    pastix_complex64_t *Lvalues  = bcsc->Lvalues;
    pastix_complex64_t *Uvalues  = bcsc->Uvalues;
    pastix_int_t itercoltab, iterval, coefindx, ldd;
    double overflow = LAPACKE_slamch( 'o' );

    assert( solvcblk->cblktype & CBLK_LAYOUT_2D );

    for (itercoltab=0; itercoltab<csccblk->colnbr; itercoltab++)
    {
        pastix_int_t frow = csccblk->coltab[itercoltab];
        pastix_int_t lrow = csccblk->coltab[itercoltab+1];

        solvblok = solvcblk->fblokptr;
        ldd = blok_rownbr( solvblok );
        lcoeftab = (pastix_complex32_t*)(solvblok->LRblock[0]->u);
        ucoeftab = (pastix_complex32_t*)(solvblok->LRblock[1]->u);

        for (iterval=frow; iterval<lrow; iterval++)
        {
            pastix_int_t rownum = bcsc->rowtab[iterval];

            /* If values in the lower part of the matrix */
            if (rownum >= (solvcblk->fcolnum+itercoltab))
            {
                while ((solvblok < lsolvblok) &&
                       ((solvblok->lrownum < rownum) ||
                        (solvblok->frownum > rownum)))
                {
                    solvblok++;
                    ldd = blok_rownbr( solvblok );
                    lcoeftab = (pastix_complex32_t*)(solvblok->LRblock[0]->u);
                    ucoeftab = (pastix_complex32_t*)(solvblok->LRblock[1]->u);
                }

                if ( solvblok < lsolvblok )
                {
                    coefindx  = rownum - solvblok->frownum; /* Row shift    */
                    coefindx += itercoltab * ldd;           /* Column shift */
                    pastix_cblk_lock( solvcblk );
                    solvblok->iluklvl = 0;
                    pastix_cblk_unlock( solvcblk );

                    if ( side != PastixUCoef ) {
                        cpucblk_zccheck_overflow( Lvalues[iterval], overflow );
                        lcoeftab[coefindx] = (pastix_complex32_t) Lvalues[iterval];
                    }

                    if ( (side != PastixLCoef) &&
                         (rownum > (solvcblk->fcolnum + itercoltab)) )
                    {
                        cpucblk_zccheck_overflow( Uvalues[iterval], overflow);
#if defined(PRECISION_zc)
                        if (bcsc->mtxtype == PastixHermitian) {
                            ucoeftab[coefindx] = (pastix_complex32_t) conj(Uvalues[iterval]);
                        }
                        else
#endif
                        {
                            ucoeftab[coefindx] = (pastix_complex32_t) Uvalues[iterval];
                        }
                    }
                }
#if defined(PASTIX_DEBUG_COEFTAB)
                else
                {
                    fprintf( stderr, "cpucblk_zfillin: drop coeff from CSC c=%ld(%ld) l=%ld(%ld) cblk=%ld fcol=%ld lcol=%ld\n",
                             (long)solvcblk->fcolnum + itercoltab, (long)itercoltab,
                             (long)rownum, (long)iterval, (long)itercblk,
                             (long)solvcblk->fcolnum, (long)solvcblk->lcolnum );
                }
#endif
            }
        }
    }
    (void)overflow;
    return 0;
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
 *******************************************************************************
 *
 * @return 0 on success, 1 in case of overflow during initialization.
 *
 *******************************************************************************/
int
cpucblk_zcfillin( pastix_coefside_t    side,
                  const SolverMatrix  *solvmtx,
                  const pastix_bcsc_t *bcsc,
                  pastix_int_t         itercblk )
{
    int rc;
    if ( (solvmtx->cblktab + itercblk)->cblktype & CBLK_COMPRESSED ) {
        rc = cpucblk_zcfillin_lr( side, solvmtx, bcsc, itercblk );
    }
    else {
        rc = cpucblk_zcfillin_fr( side, solvmtx, bcsc, itercblk );
    }
    return rc;
}
