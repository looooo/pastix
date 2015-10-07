/**
 *
 * @file panel_trsm_cpu.h
 *
 *  Factoization of a column block on a CPU.
 *
 * @version 6.0.0
 * @author Xavier Lacoste
 * @date 2014-04-09
 *
 *
 **/
#ifndef _PANEL_TRSM_CPU_H
#define _PANEL_TRSM_CPU_H

#include "common.h"
#include "z_solver.h"
#include "compute_trsm.h"

/**
 * Compute the factorization of a column block diagonal block in the genral case.
 *
 * @param[in] cblk
 *      Column block to factorize.
 * @param[in] criteria
 *      Criteria for pivoting.
 *
 * @returns Number of pivoting value detected.
 */
static inline
int
panel_trsm_cpu(const z_SolverCblk * cblk,
               pastix_complex64_t   * buffer) {
    pastix_int_t dima = cblk->lcolnum - cblk->fcolnum + 1;
    pastix_int_t dimb = cblk->stride - dima;
    pastix_complex64_t *lExtraDiag = cblk->coeftab + (cblk->fblokptr+1)->coefind;
#if (defined CHOL_SOPALIN && defined SOPALIN_LU)
    pastix_complex64_t *uExtraDiag = cblk->ucoeftab + (cblk->fblokptr+1)->coefind;
#endif
    kernel_trsm(dimb, dima,
                cblk->coeftab,
#if (defined CHOL_SOPALIN && defined SOPALIN_LU)
                cblk->ucoeftab,
#endif
                cblk->stride,
                lExtraDiag,
#if (defined CHOL_SOPALIN && defined SOPALIN_LU)
                uExtraDiag,
#endif
#ifndef CHOL_SOPALIN
                buffer,
#endif
                cblk->stride);
    return PASTIX_SUCCESS;

}

#endif /* _PANEL_TRSM_CPU_H */
