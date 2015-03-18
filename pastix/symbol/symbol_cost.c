#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "common.h"
#include "dof.h"
#include "symbol.h"
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
    return FLOPS_DSYRK( K, M ) + 6. * (double)M * (double)M;
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
    return FLOPS_ZGEMM( M, N, K ) + FLOPS_ZGEMM( M-N, N, K );
    + 2. * (double)M * (double)N + 2. * (double)(M-N) * (double)(N); /* Add step */
}

static inline double
flops_dgetrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    return FLOPS_DGEMM( M, N, K ) + FLOPS_DGEMM( M-N, N, K );
    + (double)M * (double)N + (double)(M-N) * (double)(N); /* Add step */
}

static inline double
flops_zsytrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    /* If we consider that we stored the D * A^t somewhere */
    //return FLOPS_ZGEMM( M, N, K ) + 2. * (double)M * (double)N;
    /* If not, as it is the case in the runtime */
    return FLOPS_ZGEMM( M, N, K )
        + 2. * (double)M * (double)N   /* Add step   */
        + 6. * (double)M * (double)N;  /* Scale step */
}

static inline double
flops_dsytrf_blkupdate( pastix_int_t M, pastix_int_t N, pastix_int_t K )
{
    /* If we consider that we stored the D * A^t somewhere */
    //return FLOPS_ZGEMM( M, N, K ) + (double)M * (double)N;
    /* If not, as it is the case in the runtime */
    return FLOPS_DGEMM( M, N, K )
        + (double)M * (double)N   /* Add step   */
        + (double)M * (double)N;  /* Scale step */
}

typedef struct flops_function_s {
    double (*diag     )(pastix_int_t);
    double (*trsm     )(pastix_int_t, pastix_int_t);
    double (*update   )(pastix_int_t, pastix_int_t);
    double (*blkupdate)(pastix_int_t, pastix_int_t, pastix_int_t);
} flops_function_t;


static flops_function_t flopstable[2][4] = {
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

static double
sum1d(const flops_function_t *fptr,
            pastix_int_t      cblknum,
      const SymbolMatrix     *symbmtx,
      const Dof              *dofptr)
{
    double M, N, K;
    pastix_int_t k;
    double nbops = 0.;

    /*
     * Size of the factorization kernel (square)
     */
    N = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);

    /*
     * Height of the TRSM to which apply the TRSM
     */
    M = 0;
    for(k = symbmtx->cblktab[cblknum].bloknum+1;
        k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
        M += (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);
    }

#ifdef DOF_CONSTANT
    N *= (double)dofptr->noddval;
    M *= (double)dofptr->noddval;
#endif

    nbops = fptr->diag( N );
    if( M > 0 ) {
        nbops += fptr->trsm( M, N );
        nbops += fptr->update( N, M );
    }

    return nbops;
}

static double
sum2d(const flops_function_t *fptr,
            pastix_int_t      cblknum,
      const SymbolMatrix     *symbmtx,
      const Dof              *dofptr)
{
    double M, N, K;
    pastix_int_t k;
    double nbops = 0.;

    /*
     * Size of the factorization kernel (square)
     */
    N = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);

    /*
     * Height of the TRSM to which apply the TRSM
     */
    M = 0;
    for(k = symbmtx->cblktab[cblknum].bloknum+1;
        k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
        M += (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);
    }

#ifdef DOF_CONSTANT
    N *= (double)dofptr->noddval;
    M *= (double)dofptr->noddval;
#endif

    nbops = fptr->diag( N );
    if( M > 0 ) {
        nbops += fptr->trsm( M, N );
    }

    /*
     * Compute the cost of each GEMM
     */
    K = N;
    for(k = symbmtx->cblktab[cblknum].bloknum+1;
        k < symbmtx->cblktab[cblknum+1].bloknum-1; k++)
    {
        N = (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k  ].frownum + 1);

#ifdef DOF_CONSTANT
        N *= (double)dofptr->noddval;
#endif
        nbops += fptr->blkupdate( K, M, N );

        M -= N;
    }
    return nbops;
}

static double
recursive_sum(pastix_int_t a, pastix_int_t b,
              double (*fval)(const flops_function_t *, pastix_int_t, const SymbolMatrix *, const Dof *),
              flops_function_t *fptr, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
    if(a != b)
        return recursive_sum(        a, (a+b)/2, fval, fptr, symbmtx, dofptr)
            +  recursive_sum((a+b)/2+1,       b, fval, fptr, symbmtx, dofptr);

    return fval(fptr, a, symbmtx, dofptr);
}

void
symbolCost(const SymbolMatrix *symbmtx, const Dof *dofptr,
           pastix_coeftype_t flttype, pastix_factotype_t factotype,
           pastix_int_t *nnz, double *thflops, double *rlflops )
{
    /* Compute NNZ */
    if ( nnz != NULL ) {
        *nnz = symbolGetNNZ( symbmtx );
    }

    /* Compute theroritical flops */
    if ( thflops != NULL ) {
        *thflops = recursive_sum(0, symbmtx->cblknbr-1, sum1d,
                                 &(flopstable[flttype][factotype]),
                                 symbmtx, dofptr);
    }

    /* Compute real flops */
    if ( rlflops != NULL ) {
        *rlflops = recursive_sum(0, symbmtx->cblknbr-1, sum2d,
                                 &(flopstable[flttype][factotype]),
                                 symbmtx, dofptr);
    }
}
