/**
 *
 * @file symbol_cost_flops.c
 *
 *  PaStiX symbol functions to compute the number of flops induced by the chosen
 *  symbolic structure.
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author David Goudin
 * @author Francois Pelegrin
 * @author Mathieu Faverge
 * @author Pascal Henon
 * @author Pierre Ramet
 * @date 2013-06-24
 *
 **/
#include "common.h"
#include "symbol_cost.h"
#include "flops.h"

/**
 * Computations flops of diagonal blocks
 */
static double
flops_zgetrf_diag( pastix_int_t N ) {
    return FLOPS_ZGETRF( N, N );
}

static inline double
flops_dgetrf_diag( pastix_int_t N ) {
    return FLOPS_DGETRF( N, N );
}

static inline double
flops_zpotrf_diag( pastix_int_t N ) {
    return FLOPS_ZPOTRF( N );
}

static inline double
flops_dpotrf_diag( pastix_int_t N ) {
    return FLOPS_DPOTRF( N );
}

static inline double
flops_zsytrf_diag( pastix_int_t N ) {
    return FLOPS_ZSYTRF( N );
}

static inline double
flops_dsytrf_diag( pastix_int_t N ) {
    return FLOPS_DSYTRF( N );
}

/**
 * Computations flops of the solve step
 */
static inline double
flops_zpotrf_trsm( pastix_int_t M, pastix_int_t N ) {
    return FLOPS_ZTRSM( PlasmaRight, M, N );
}

static inline double
flops_dpotrf_trsm( pastix_int_t M, pastix_int_t N ) {
    return FLOPS_DTRSM( PlasmaRight, M, N );
}

static double
flops_zgetrf_trsm( pastix_int_t M, pastix_int_t N ) {
    return 2. * FLOPS_ZTRSM( PlasmaRight, M, N );
}

static inline double
flops_dgetrf_trsm( pastix_int_t M, pastix_int_t N ) {
    return 2. * FLOPS_DTRSM( PlasmaRight, M, N );
}

static inline double
flops_zsytrf_trsm( pastix_int_t M, pastix_int_t N ) {
    return FLOPS_ZTRSM( PlasmaRight, M, N ) + 6. * (double)N * (double)M;
}

static inline double
flops_dsytrf_trsm( pastix_int_t M, pastix_int_t N ) {
    return FLOPS_DTRSM( PlasmaRight, M, N ) + (double)N * (double)M;
}

/**
 * Theroretical computation flops of the update step
 */
static inline double
flops_zpotrf_update( pastix_int_t K, pastix_int_t M ) {
    return FLOPS_ZHERK( K, M );
}

static inline double
flops_dpotrf_update( pastix_int_t K, pastix_int_t M ) {
    return FLOPS_DSYRK( K, M );
}

static double
flops_zgetrf_update( pastix_int_t K, pastix_int_t M ) {
    return FLOPS_ZGEMM( K, M, M );
}

static inline double
flops_dgetrf_update( pastix_int_t K, pastix_int_t M ) {
    return FLOPS_DGEMM( K, M, M );
}

static inline double
flops_zsytrf_update( pastix_int_t K, pastix_int_t M ) {
    return FLOPS_ZSYRK( K, M ) + 6. * (double)M * (double)M;
}

static inline double
flops_dsytrf_update( pastix_int_t K, pastix_int_t M ) {
    return FLOPS_DSYRK( K, M ) + (double)M * (double)M;
}

/**
 * Real computation flops of the update step
 */
static inline double
flops_zpotrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    return FLOPS_ZGEMM( M, N, K ) + 2. * (double)M * (double)N;
}

static inline double
flops_dpotrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    return FLOPS_DGEMM( M, N, K ) + (double)M * (double)N;
}

static double
flops_zgetrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    return FLOPS_ZGEMM( M, N, K ) + FLOPS_ZGEMM( M-N, N, K )
        + 2. * (double)M * (double)N + 2. * (double)(M-N) * (double)(N); /* Add step */
}

static inline double
flops_dgetrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    return FLOPS_DGEMM( M, N, K ) + FLOPS_DGEMM( M-N, N, K )
        + (double)M * (double)N + (double)(M-N) * (double)(N); /* Add step */
}

static inline double
flops_zsytrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    /* If we consider that we stored the D * A^t somewhere */
#if 0
    return FLOPS_ZGEMM( M, N, K )
        + 2. * (double)M * (double)N;
#else
    /* If not, as it is the case in the runtime */
    return FLOPS_ZGEMM( M, N, K )
        + 2. * (double)M * (double)N   /* Add step   */
        + 6. * (double)M * (double)N;  /* Scale step */
#endif
}

static inline double
flops_dsytrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    /* If we consider that we stored the D * A^t somewhere */
#if 0
    return FLOPS_DGEMM( M, N, K )
        + (double)M * (double)N;
#else
    /* If not, as it is the case in the runtime */
    return FLOPS_DGEMM( M, N, K )
        + (double)M * (double)N   /* Add step   */
        + (double)M * (double)N;  /* Scale step */
#endif
}

symbol_function_t flopstable[2][4] = {
    {
        {flops_zgetrf_diag, flops_zgetrf_trsm, flops_zgetrf_update, flops_zgetrf_blkupdate },
        {flops_zpotrf_diag, flops_zpotrf_trsm, flops_zpotrf_update, flops_zpotrf_blkupdate },
        {flops_zsytrf_diag, flops_zsytrf_trsm, flops_zsytrf_update, flops_zsytrf_blkupdate },
        {flops_zsytrf_diag, flops_zsytrf_trsm, flops_zsytrf_update, flops_zsytrf_blkupdate }
    },
    {
        {flops_dgetrf_diag, flops_dgetrf_trsm, flops_dgetrf_update, flops_dgetrf_blkupdate },
        {flops_dpotrf_diag, flops_dpotrf_trsm, flops_dpotrf_update, flops_dpotrf_blkupdate },
        {flops_dsytrf_diag, flops_dsytrf_trsm, flops_dsytrf_update, flops_dsytrf_blkupdate },
        {flops_dsytrf_diag, flops_dsytrf_trsm, flops_dsytrf_update, flops_dsytrf_blkupdate }
    }
};
