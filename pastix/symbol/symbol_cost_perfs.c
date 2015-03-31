/**
 *
 * @file symbol_cost_perfs.c
 *
 *  PaStiX symbol functions to compute the computational time induced by the chosen
 *  symbolic structure with a given performance model.
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
#include "perf.h"

/**
 * Time models of diagonal blocks
 */
static inline double
perfs_dpotrf_diag( pastix_int_t N ) {
    double total = PERF_POTRF(N);
    assert( N > 0 );
    return (total > 0.) ? total : 0.;
}

static inline double
perfs_dsytrf_diag( pastix_int_t N ) {
    double total = PERF_SYTRF(N) + (double)N * PERF_COPY(N);
    assert( N > 0 );
    return (total > 0.) ? total : 0.;
}

static inline double
perfs_dgetrf_diag( pastix_int_t N ) {
    /* Approximate to 2 times potrf */
    return 2. * perfs_dpotrf_diag( N );
}

/**
 * Time models of the solve step
 */
static inline double
perfs_dpotrf_trsm( pastix_int_t M, pastix_int_t N ) {
    double total = PERF_TRSM( N, M );
    assert( (M > 0) && (N > 0) );
    return (total > 0.) ? total : 0.;
}

static inline double
perfs_dsytrf_trsm( pastix_int_t M, pastix_int_t N ) {
    double total = PERF_TRSM( N, M )
        +  (double)N * (PERF_SCAL(M) + PERF_COPY(M));
    assert( (M > 0) && (N > 0) );
    return (total > 0.) ? total : 0.;
}

static inline double
perfs_dgetrf_trsm( pastix_int_t M, pastix_int_t N ) {
    double total = 2. * PERF_TRSM( N, M );
    assert( (M > 0) && (N > 0) );
    return (total > 0.) ? total : 0.;
}

/**
 * Time model of the update step per block
 */
static inline double
perfs_dpotrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    double total = PERF_GEMM( M, N, K ) + PERF_GEAM( M, N );
    assert( (M > 0) && (N > 0) && (K > 0) );
    return (total > 0.) ? total : 0.;
}

static inline double
perfs_dsytrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    double total = PERF_GEMM( M, N, K ) + PERF_GEAM( M, N );
    assert( (M > 0) && (N > 0) && (K > 0) );
    return (total > 0.) ? total : 0.;
}

static inline double
perfs_dgetrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    double total = PERF_GEMM(     M, N, K ) + PERF_GEAM(     M, N )
        +          PERF_GEMM( (M-N), N, K ) + PERF_GEAM( (M-N), N );
    assert( (M > 0) && (N > 0) && (K > 0) );
    return (total > 0.) ? total : 0.;
}

symbol_function_t perfstable[2][4] = {
    {
        {perfs_dgetrf_diag, perfs_dgetrf_trsm, NULL, perfs_dgetrf_blkupdate },
        {perfs_dpotrf_diag, perfs_dpotrf_trsm, NULL, perfs_dpotrf_blkupdate },
        {perfs_dsytrf_diag, perfs_dsytrf_trsm, NULL, perfs_dsytrf_blkupdate },
        {perfs_dsytrf_diag, perfs_dsytrf_trsm, NULL, perfs_dsytrf_blkupdate }
    },
    {
        {perfs_dgetrf_diag, perfs_dgetrf_trsm, NULL, perfs_dgetrf_blkupdate },
        {perfs_dpotrf_diag, perfs_dpotrf_trsm, NULL, perfs_dpotrf_blkupdate },
        {perfs_dsytrf_diag, perfs_dsytrf_trsm, NULL, perfs_dsytrf_blkupdate },
        {perfs_dsytrf_diag, perfs_dsytrf_trsm, NULL, perfs_dsytrf_blkupdate }
    }
};
