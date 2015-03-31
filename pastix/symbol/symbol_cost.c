#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "common.h"
#include "symbol.h"
#include "symbol_cost.h"

static double
sum1d(const symbol_function_t *fptr,
      const SymbolMatrix     *symbmtx,
            pastix_int_t      cblknum)
{
    SymbolCblk *cblk = symbmtx->cblktab + cblknum;
    pastix_int_t M, N, k;
    double nbops = 0.;
    double dof = (double)(symbmtx->dof);

    /*
     * Size of the factorization kernel (square)
     */
    N = (cblk->lcolnum - cblk->fcolnum + 1);

    /*
     * Height of the TRSM to which apply the TRSM
     */
    M = 0;
    for(k = cblk[0].bloknum+1; k < cblk[1].bloknum; k++)
    {
        M += (symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);
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
sum2d(const symbol_function_t *fptr,
      const SymbolMatrix     *symbmtx,
            pastix_int_t      cblknum)
{
    SymbolCblk *cblk = symbmtx->cblktab + cblknum;
    pastix_int_t M, N, K, l;
    double nbops = 0.;
    double dof = (double)(symbmtx->dof);

    /*
     * Size of the factorization kernel (square)
     */
    N = (cblk->lcolnum - cblk->fcolnum + 1);

    /*
     * Height of the TRSM to which apply the TRSM
     */
    M = 0;
    for(l = cblk[0].bloknum+1; l < cblk[1].bloknum; l++)
    {
        M += (symbmtx->bloktab[l].lrownum - symbmtx->bloktab[l].frownum + 1);
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
    for(l = cblk[0].bloknum+1; l < cblk[1].bloknum; l++)
    {
        N = (symbmtx->bloktab[l].lrownum - symbmtx->bloktab[l].frownum + 1);
        N *= dof;

        nbops += fptr->blkupdate( K, M, N );

        M -= N;
    }
    return nbops;
}

static double
sum2dext(const symbol_function_t *fptr,
         const SymbolMatrix     *symbmtx,
               pastix_int_t      cblknum,
               double           *blokcost)
{
    SymbolCblk *cblk = symbmtx->cblktab + cblknum;
    pastix_int_t M, N, K, l;
    double nbops = 0.;
    double dof = (double)(symbmtx->dof);

    /*
     * Size of the factorization kernel (square)
     */
    N = (cblk->lcolnum - cblk->fcolnum + 1);

    /*
     * Height of the TRSM to which apply the TRSM
     */
    M = 0;
    for(l = cblk[0].bloknum+1; l < cblk[1].bloknum; l++)
    {
        M += (symbmtx->bloktab[l].lrownum - symbmtx->bloktab[l].frownum + 1);
    }

    N *= dof;
    M *= dof;

    nbops = fptr->diag( N );
    if( M > 0 ) {
        nbops += fptr->trsm( M, N );
    }
    *blokcost = nbops;
    blokcost++;

    /*
     * Compute the cost of each GEMM
     */
    K = N;
    for(l = cblk[0].bloknum+1; l < cblk[1].bloknum; l++, blokcost++)
    {
        N = (symbmtx->bloktab[l].lrownum - symbmtx->bloktab[l].frownum + 1);
        N *= dof;

        *blokcost = fptr->blkupdate( K, M, N );
        nbops += *blokcost;

        M -= N;
    }
    return nbops;
}

static double
recursive_sum(pastix_int_t a, pastix_int_t b,
              double (*fval)(const symbol_function_t *, const SymbolMatrix *, pastix_int_t),
              symbol_function_t *fptr, const SymbolMatrix * symbmtx)
{
    if(a != b)
        return recursive_sum(        a, (a+b)/2, fval, fptr, symbmtx)
            +  recursive_sum((a+b)/2+1,       b, fval, fptr, symbmtx);

    return fval(fptr, symbmtx, a);
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
symbolGetTimes(const SymbolMatrix *symbmtx,
               pastix_coeftype_t flttype, pastix_factotype_t factotype,
               double *cblkcost, double *blokcost )
{
    symbol_function_t *f;
    double *cblkptr, *blokptr;
    pastix_int_t i;
    int iscomplex = (flttype == PastixComplex32) || (flttype == PastixComplex64);
    f = &(perfstable[iscomplex][factotype]);

    /* Initialize costs */
    cblkptr = cblkcost;
    blokptr = blokcost;

    for(i=0; i<symbmtx->cblknbr; i++, cblkptr++) {
        *cblkptr = sum2dext( f, symbmtx, i, blokptr );

        blokptr += symbmtx->cblktab[i+1].bloknum
            -      symbmtx->cblktab[i  ].bloknum;
    }

    assert( ( cblkptr - cblkcost ) == symbmtx->cblknbr );
    assert( ( blokptr - blokcost ) == symbmtx->bloknbr );
}
