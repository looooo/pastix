#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "common.h"
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
      const SymbolMatrix     *symbmtx)
{
    double M, N;
    pastix_int_t k;
    double nbops = 0.;
    double dof = (double)(symbmtx->dof);

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

    if ( dof > 0. ) {
        N *= dof;
        M *= dof;
    }

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
      const SymbolMatrix     *symbmtx)
{
    double M, N, K;
    pastix_int_t k;
    double nbops = 0.;
    double dof = (double)(symbmtx->dof);

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

    N *= dof;
    M *= dof;

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
        N *= dof;

        nbops += fptr->blkupdate( K, M, N );

        M -= N;
    }
    return nbops;
}

static double
recursive_sum(pastix_int_t a, pastix_int_t b,
              double (*fval)(const flops_function_t *, pastix_int_t, const SymbolMatrix *),
              flops_function_t *fptr, const SymbolMatrix * symbmtx)
{
    if(a != b)
        return recursive_sum(        a, (a+b)/2, fval, fptr, symbmtx)
            +  recursive_sum((a+b)/2+1,       b, fval, fptr, symbmtx);

    return fval(fptr, a, symbmtx);
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbol
 *
 * symbolGetNNZ - Computes the number of non-zero elements stored in the
 * symbol matrix in order to compute the fill-in. This routines returns the
 * number of non-zero of the strictly lower part of the matrix.
 *
 *******************************************************************************
 *
 * @param[in] symbptr
 *          The symbol structure to study.
 *
 *******************************************************************************
 *
 * @return
 *          The number of non zero elements in the strictly lower part of the
 *          full symbol matrix.
 *
 *******************************************************************************/
pastix_int_t
symbolGetNNZ(const SymbolMatrix *symbptr)
{
    SymbolCblk *cblk;
    SymbolBlok *blok;
    pastix_int_t itercblk;
    pastix_int_t cblknbr;
    pastix_int_t nnz = 0;
    pastix_int_t dof = symbptr->dof;

    cblknbr = symbptr->cblknbr;
    cblk    = symbptr->cblktab;
    blok    = symbptr->bloktab;

    for(itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
    {
        pastix_int_t iterblok = cblk[0].bloknum + 1;
        pastix_int_t lbloknum = cblk[1].bloknum;

        pastix_int_t colnbr = dof * (cblk->lcolnum - cblk->fcolnum + 1);

        /* Diagonal block */
        blok++;
        nnz += ( colnbr * (colnbr+1) ) / 2 - colnbr;

        /* Off-diagonal blocks */
        for( ; iterblok < lbloknum; iterblok++, blok++)
        {
            pastix_int_t rownbr = (blok->lrownum - blok->frownum + 1) * dof;

            nnz += rownbr * colnbr;
        }
    }

    return nnz;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_symbol
 *
 * symbolGetFlops - Computes the number of theoritical and real flops to
 * factorized the given symbolic matrix with the specified type of factorization
 * and floating type.
 *
 *******************************************************************************
 *
 * @param[in] symbptr
 *          The symbol structure to study.
 *
 * @param[in] flttype
 *          The floating type of the elements in the matrix.
 *          PastixPattern, PastixFloat, PastixDouble, PastixComplex32 or
 *          PastixComplex64. In case of PastixPattern, values for PastixDouble
 *          are returned.
 *
 * @param[in] factotype
 *          The factorization algorithm to perform: PastixFactLLT,
 *          PastixFactLDLT, PastixFactLDLH or PastixFactLU.
 *
 * @param[out] thflops
 *          Returns the number of theoretical flops to perform.
 *          NULL if not asked.
 *
 * @param[out] rlflops
 *          Returns the number of real flops to perform, taking into account
 *          copies and scatter operations.
 *          NULL if not asked.
 *
 *******************************************************************************/
void
symbolGetFlops(const SymbolMatrix *symbmtx,
               pastix_coeftype_t flttype, pastix_factotype_t factotype,
               double *thflops, double *rlflops )
{
    int iscomplex = (flttype == PastixComplex32) || (flttype == PastixComplex64);

    /* Compute NNZ */
    if ( nnz != NULL ) {
        *nnz = symbolGetNNZ( symbmtx );
    }

    /* Compute theroritical flops */
    if ( thflops != NULL ) {
        *thflops = recursive_sum(0, symbmtx->cblknbr-1, sum1d,
                                 &(flopstable[iscomplex][factotype]),
                                 symbmtx);
    }

    /* Compute real flops */
    if ( rlflops != NULL ) {
        *rlflops = recursive_sum(0, symbmtx->cblknbr-1, sum2d,
                                 &(flopstable[iscomplex][factotype]),
                                 symbmtx);
    }
}
