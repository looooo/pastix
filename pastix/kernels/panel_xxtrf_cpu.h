/**
 *
 * @file panel_xxtrf_cpu.h
 *
 *  Factoization of a column block on a CPU.
 *
 * @version 6.0.0
 * @author Xavier Lacoste
 * @date 2014-04-09
 *
 *
 **/
#ifndef _PANEL_XXTRF_CPU_H
#define _PANEL_XXTRF_CPU_H

#include "common.h"
#include "solver.h"
#include "sopalin_compute.h"
#include "compute_diag.h"

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
panel_getrf_cpu(const SolverCblk * cblk,
                double criteria) {
    pastix_int_t npiv = 0;
    pastix_int_t dima = cblk->lcolnum - cblk->fcolnum + 1;

    /* Add U diagonal updates into L */
    SOPALIN_GEAM("T", "N", dima, dima, 1.0,
                 cblk->ucoeftab, cblk->stride,
                 cblk->coeftab, cblk->stride);
    /* Factorize diagonal block (two terms version with workspace) */
    PASTIX_getrf_block(cblk->coeftab, dima, dima, cblk->stride,
                       &npiv, criteria);
    /* Transpose L_diag in U_diag Matrix */
    DimTrans(cblk->coeftab, cblk->stride, dima, cblk->ucoeftab);
    return npiv;
}


/**
 * Compute the factorization of a column block diagonal block in the SDP case.
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
panel_potrf_cpu(const SolverCblk * cblk,
                double criteria) {
    pastix_int_t npiv = 0;
    pastix_int_t dima = cblk->lcolnum - cblk->fcolnum + 1;

    PASTIX_potrf_block(cblk->coeftab, dima, cblk->stride,
                       &npiv, criteria);
    return npiv;
}


/**
 * Compute the factorization of a column block diagonal block in the LDLt case.
 *
 * @param[in] cblk
 *      Column block to factorize.
 * @param[in] criteria
 *      Criteria for pivoting.
 * @param[out] buffer
 *      Byffer used to store temporary data.
 *
 * @returns Number of pivoting value detected.
 */
static inline
int
panel_sytrf_cpu(const SolverCblk * cblk,
                double criteria,
                pastix_float_t * buffer) {
    pastix_int_t npiv = 0;
    pastix_int_t dima = cblk->lcolnum - cblk->fcolnum + 1;

    PASTIX_sytrf_block(cblk->coeftab, dima, cblk->stride,
                       &npiv, criteria, buffer);
    return npiv;
}


/**
 * Compute the factorization of a column block diagonal block hermitian case.
 *
 * @param[in] cblk
 *      Column block to factorize.
 * @param[in] criteria
 *      Criteria for pivoting.
 * @param[out] buffer
 *      Byffer used to store temporary data.
 *
 * @returns Number of pivoting value detected.
 */
static inline
int
panel_hetrf_cpu(const SolverCblk * cblk,
                double criteria,
                pastix_float_t * buffer) {
    pastix_int_t npiv = 0;
    pastix_int_t dima = cblk->lcolnum - cblk->fcolnum + 1;

    PASTIX_hetrf_block(cblk->coeftab, dima, cblk->stride,
                       &npiv, criteria, buffer);
    return npiv;
}

/**
 * Compute the factorization of a column block diagonal block automatic detection of the case.
 *
 * @param[in] cblk
 *      Column block to factorize.
 * @param[in] criteria
 *      Criteria for pivoting.
 * @param[out] buffer
 *      Byffer used to store temporary data.
 *
 * @returns Number of pivoting value detected.
 */
static inline
int
panel_xxtrf_cpu(const SolverCblk * cblk,
                double criteria,
                pastix_float_t * buffer) {
#ifdef CHOL_SOPALIN
#  ifdef SOPALIN_LU
    return panel_getrf_cpu(cblk, criteria);
#  else
    return panel_potrf_cpu(cblk, criteria);
#  endif
#else
#  ifdef HERMITIAN
    return panel_hetrf_cpu(cblk, criteria, buffer);
#  else
    return panel_sytrf_cpu(cblk, criteria, buffer);
#  endif
#endif
}

#endif /* _PANEL_XXTRF_CPU_H */
