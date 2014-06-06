/************************************************************/
/**                                                        **/
/**   NAME       : solver.h                                **/
/**                                                        **/
/**   AUTHORS    : David GOUDIN                            **/
/**                Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                Pierre RAMET                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the solver matrix.                  **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     28 oct 1998     **/
/**                # Version 1.0  : from : 06 jun 2002     **/
/**                                 to     06 jun 2002     **/
/**                                                        **/
/************************************************************/

#ifndef SOLVER_H
#define SOLVER_H

#include "solver_struct.h"

/*
**  The type and structure definitions.
*/

#define COMP_1D                     0
#define DIAG                        1
#define E1                          2
#define E2                          3
#define DRUNK                       4

pastix_int_t
sizeofsolver(const SolverMatrix *solvptr,
             pastix_int_t *iparm );

/**
 * Indicates whether a column block is in halo.
 *
 * @param datacode SolverMatrix structure.
 * @param column block SolverCblk structure to test.
 *
 * @retval API_YES if the column block is in halo.
 * @retval API_NO  if the column block is not in halo.
 */
static inline
int cblk_ishalo( SolverMatrix * datacode,
                 SolverCblk   * cblk ) {
#ifdef PASTIX_WITH_STARPU
    if ((size_t)cblk >= (size_t)datacode->hcblktab &&
        (size_t)cblk < (size_t)(datacode->hcblktab+datacode->hcblknbr)) {
        return API_YES;
    }
#endif
    return API_NO;
}

/**
 * Indicates whether a column block is a fanin column block.
 *
 * @param datacode SolverMatrix structure.
 * @param column block SolverCblk structure to test.
 *
 * @retval API_YES if the column block is a fanin column block.
 * @retval API_NO  if the column block is not a fanin column block.
 */
static inline
int cblk_isfanin( SolverMatrix * datacode,
                  SolverCblk   * cblk ) {
    pastix_int_t clustnum;
#ifdef PASTIX_WITH_STARPU
    for (clustnum = 0; clustnum < datacode->clustnbr; clustnum++) {
        if ((size_t)cblk >= (size_t)datacode->fcblktab[clustnum] &&
            (size_t)cblk < (size_t)(datacode->fcblktab[clustnum] +
                                    datacode->fcblknbr[clustnum])) {
            return API_YES;
        }
    }
#endif
    return API_NO;
}


/**
 * Indicates whether a column block is a local column block.
 *
 * @param datacode SolverMatrix structure.
 * @param column block SolverCblk structure to test.
 *
 * @retval API_YES if the column block is a local column block.
 * @retval API_NO  if the column block is not a local column block.
 */
static inline
int cblk_islocal( SolverMatrix * datacode,
                  SolverCblk   * cblk ) {
    if ((size_t)cblk >= (size_t)datacode->cblktab &&
        (size_t)cblk < (size_t)(datacode->cblktab+datacode->cblknbr)) {
        return API_YES;
    }
    return API_NO;
}


/**
 * Get the index of a local column block.
 *
 * @param datacode SolverMatrix structure.
 * @param column block SolverCblk structure to test.
 *
 * @returns the index of the column block.
 */
static inline
pastix_int_t cblk_getnum( SolverMatrix * datacode,
                          SolverCblk   * cblk ) {
    assert(cblk_islocal(datacode, cblk) == API_YES);
    return cblk - datacode->cblktab;
}


/**
 * Get the index of a halo column block.
 *
 * @param datacode SolverMatrix structure.
 * @param column block SolverCblk structure to test.
 *
 * @returns the index of the column block.
 */
static inline
pastix_int_t hcblk_getnum( SolverMatrix * datacode,
                           SolverCblk   * cblk ) {
#ifdef PASTIX_WITH_STARPU
    assert(cblk_ishalo(datacode, cblk) == API_YES);
    return cblk - datacode->hcblktab;
#else
    return -1;
#endif
}


/**
 * Get the index of a fanin column block.
 *
 * @param datacode SolverMatrix structure.
 * @param column block SolverCblk structure to test.
 *
 * @returns the index of the column block.
 */
static inline
pastix_int_t fcblk_getnum( SolverMatrix * datacode,
                           SolverCblk   * cblk,
                           pastix_int_t   procnum ) {
#ifdef PASTIX_WITH_STARPU
    assert(cblk_isfanin(datacode, cblk) == API_YES);
    return cblk - datacode->fcblktab[procnum];

#else
    return -1;
#endif
}

static inline
pastix_int_t fcblk_getorigin( SolverMatrix * datacode,
                              SolverMatrix * cblk ) {
#ifdef PASTIX_WITH_STARPU
    pastix_int_t clustnum;
    for (clustnum = 0; clustnum < datacode->clustnbr; clustnum++) {
        if ((size_t)cblk >= (size_t)datacode->fcblktab[clustnum] &&
            (size_t)cblk < (size_t)(datacode->fcblktab[clustnum] +
                                    datacode->fcblknbr[clustnum])) {
            return clustnum;
        }
    }
#endif
    return -1;
}
/**
 * Get the number of columns of a column block.
 *
 * @param column block SolverCblk structure.
 *
 * @returns the number of columns in the column block.
 */
static inline
pastix_int_t cblk_colnbr( SolverCblk * cblk ) {
    return cblk->lcolnum - cblk->fcolnum + 1;
}


static inline
pastix_int_t cblk_save( SolverCblk * cblk, char * name, pastix_float_t * coef) {
#ifdef PASTIX_DUMP_CBLK
    pastix_int_t i,j;
    SolverBlok *b;
    FILE * file = fopen(name, "w");
    for ( b = cblk->fblokptr; b < cblk[1].fblokptr; b++) {
        fprintf(file, "%ld %ld\n", b->frownum, b->lrownum);
    }
    for (j = 0; j < cblk->stride; j++) {
        for (i = 0; i < cblk_colnbr(cblk); i++) {
            fprintf(file, "%10.5lg ", coef[j+i*cblk->stride]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
#endif
    return PASTIX_SUCCESS;
}
/**
 * Indicate if a blok is included inside an other block.
 * i.e. indicate if the row range of the first block is included in the
 * one of the second.
 *
 * @param first block SolverBlok structure to test.
 * @param second block SolverBlok structure to test.
 *
 * @retval true   if the first block is     included in the second one.
 * @retval false  if the first block is not included in the second one.
 */
#  if defined(NAPA_SOPALIN)
static inline int is_block_inside_fblock( SolverBlok *blok, SolverBlok *fblok ) {
    return (((blok->frownum >= fblok->frownum) &&
             (blok->lrownum <= fblok->lrownum)) ||
            ((blok->frownum <= fblok->frownum) &&
             (blok->lrownum >= fblok->lrownum)) ||
            ((blok->frownum <= fblok->frownum) &&
             (blok->lrownum >= fblok->frownum)) ||
            ((blok->frownum <= fblok->lrownum) &&
             (blok->lrownum >= fblok->lrownum)));
}
#  else
static inline int is_block_inside_fblock( SolverBlok *blok, SolverBlok *fblok ) {
    return ((blok->frownum >= fblok->frownum) &&
            (blok->lrownum <= fblok->lrownum));
}
#  endif /* defined(NAPA_SOPALIN) */

#endif /* SOLVER_H */
