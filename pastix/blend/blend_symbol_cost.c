#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>

#include "common.h"
#include "symbol.h"
#include "dof.h"
#include "perf.h"
#include "symbol_cost.h"

#define PlasmaLeft  141
#define PlasmaRight 142
#include "flops.h"

double flops_zgetrf(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr);
double flops_dgetrf(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr);
double flops_zpotrf(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr);
double flops_dpotrf(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr);

double recursive_sum(pastix_int_t a, pastix_int_t b,
                     double (*fval)(pastix_int_t, const SymbolMatrix *, const Dof *),
                     const SymbolMatrix * symbmtx, const Dof * dofptr)
{
  if(a != b)
    return recursive_sum(        a, (a+b)/2, fval, symbmtx, dofptr)
         + recursive_sum((a+b)/2+1,       b, fval, symbmtx, dofptr);

  return fval(a, symbmtx, dofptr);
}

double crout_2t(pastix_int_t cblknum, const SymbolMatrix *symbmtx, const Dof * dofptr)
{
    pastix_int_t i;
    double gk = 0;
    double lk = 0;

    /* lk is the dimension of the diagonal blok */
    lk = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);
    lk = lk * dofptr->noddval;

    /* gk is the height of off-diag bloks */
    for(i=symbmtx->cblktab[cblknum].bloknum+1; i<symbmtx->cblktab[cblknum+1].bloknum; i++)
        gk += (double)(symbmtx->bloktab[i].lrownum - symbmtx->bloktab[i].frownum +1);

    gk = gk * dofptr->noddval;

    return( (2*lk*lk*lk + (6*gk+3)*lk*lk + (6*gk*gk+6*gk-5)*lk)/6);

}


double crout_3t(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
    pastix_int_t i;
    double gk = 0;
    double lk = 0;

    /* lk is the dimension of the diagonal blok */
    lk = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);
    lk = lk * dofptr->noddval;

    /* gk is the height of off-diag bloks */
    for(i=symbmtx->cblktab[cblknum].bloknum+1;i<symbmtx->cblktab[cblknum+1].bloknum;i++)
        gk += (double)(symbmtx->bloktab[i].lrownum - symbmtx->bloktab[i].frownum +1);

    gk = gk * dofptr->noddval;

    return( (lk*lk*lk + (3*gk+1)*lk*lk + (3*gk*gk-2+2*gk)*lk)/2 );

}


double crout_hyb(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
    pastix_int_t i;
    double gk = 0;
    double lk = 0;

    /* lk is the dimension of the diagonal blok */
    lk =(double)( symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);
    lk = lk * dofptr->noddval;

    /* gk is the height of off-diag bloks */
    for(i=symbmtx->cblktab[cblknum].bloknum+1;i<symbmtx->cblktab[cblknum+1].bloknum;i++)
        gk += (double)(symbmtx->bloktab[i].lrownum - symbmtx->bloktab[i].frownum +1);

    gk = gk * dofptr->noddval;

    return( (lk*lk*lk + 3*(gk+1)*lk*lk + (3*gk*gk -4 + 6*gk)*lk)/3 );
}

double cholesky(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
    pastix_int_t i;
    double gk = 0;
    double lk = 0;

    /* lk is the dimension of the diagonal blok */
    lk = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);
    lk = lk * dofptr->noddval;

    /* gk is the height of off-diag bloks */
    for(i=symbmtx->cblktab[cblknum].bloknum+1;i<symbmtx->cblktab[cblknum+1].bloknum;i++)
        gk += (double)(symbmtx->bloktab[i].lrownum - symbmtx->bloktab[i].frownum +1);

    gk = gk * dofptr->noddval;

    return( (  2*lk*lk*lk +
               (6*gk-3)*lk*lk +
               (6*gk   *gk + 1-6*gk)*lk) /6 );
}

/*******************************************/
/* Number of non zero  extradiagonal terms */
/*******************************************/

double nnz(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
    pastix_int_t i;
    pastix_int_t noddval = ( dofptr == NULL ) ? 1 : dofptr->noddval;
    double gk = 0;
    double lk = 0;

    /* lk is the dimension of the diagonal blok */
    lk = (double)(symbmtx->cblktab[cblknum].lcolnum - symbmtx->cblktab[cblknum].fcolnum + 1);
    lk = lk * noddval;

    /* gk is the height of off-diag bloks */
    for(i=symbmtx->cblktab[cblknum].bloknum+1;i<symbmtx->cblktab[cblknum+1].bloknum;i++)
        gk +=(double)( symbmtx->bloktab[i].lrownum - symbmtx->bloktab[i].frownum +1);

    gk = gk * noddval;

    return( lk*(lk+1)/2 + gk*lk - lk);
}


double crout_blok(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
    double M, N, K;
    pastix_int_t k;
    double nbops = 0;

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

    assert( N > 0 );

#ifdef DOF_CONSTANT
    N *= (double)dofptr->noddval;
    M *= (double)dofptr->noddval;
#endif

    nbops = FLOPS_DSYTRF(N);
    if( M > 0 )
        nbops += FLOPS_DTRSM( PlasmaRight, M, N ) + N * (double)(OPS_SCAL(M));

    /*
     * Compute the cost of each GEMM
     */
    K = N;
    for(k = symbmtx->cblktab[cblknum].bloknum+1;
        k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
        N = (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);

#ifdef DOF_CONSTANT
        N *= (double)dofptr->noddval;
#endif

        nbops += /*K * (double)(OPS_SCAL(M)) + */
            FLOPS_DGEMM( M, N, K ) +
            (double)(OPS_GEAM( M, N ));

        M -= N;
    }

    return nbops;
}

double flops_zgetrf(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
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

    nbops  = FLOPS_ZGETRF( N, N );
    if( M > 0 )
        nbops += 2. * FLOPS_ZTRSM( PlasmaRight, M, N );

    /*
     * Compute the cost of each GEMM
     */
    K = N;
    for(k = symbmtx->cblktab[cblknum].bloknum+1;
        k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
        N = (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);

#ifdef DOF_CONSTANT
        N *= (double)dofptr->noddval;
#endif

        nbops += 2. * FLOPS_ZGEMM( M, N, K );

        M -= N;
    }

    return nbops;
}

double flops_dgetrf(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
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

    nbops  = FLOPS_DGETRF( N, N );
    if( M > 0 )
        nbops += 2. * FLOPS_DTRSM( PlasmaRight, M, N );

    /*
     * Compute the cost of each GEMM
     */
    K = N;
    for(k = symbmtx->cblktab[cblknum].bloknum+1;
        k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
        N = (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);

#ifdef DOF_CONSTANT
        N *= (double)dofptr->noddval;
#endif

        nbops += 2. * FLOPS_DGEMM( M, N, K );

        M -= N;
    }

    return nbops;
}

double flops_zpotrf(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
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

    nbops  = FLOPS_ZPOTRF( N );
    if( M > 0 )
        nbops += FLOPS_ZTRSM( PlasmaRight, M, N );

    /*
     * Compute the cost of each GEMM
     */
    K = N;
    for(k = symbmtx->cblktab[cblknum].bloknum+1;
        k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
        N = (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);

#ifdef DOF_CONSTANT
        N *= (double)dofptr->noddval;
#endif

        nbops += FLOPS_ZGEMM( M, N, K );

        M -= N;
    }

    return nbops;
}

double flops_dpotrf(pastix_int_t cblknum, const SymbolMatrix * symbmtx, const Dof * dofptr)
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

    nbops  = FLOPS_DPOTRF( N );
    if( M > 0 )
        nbops += FLOPS_DTRSM( PlasmaRight, M, N );

    /*
     * Compute the cost of each GEMM
     */
    K = N;
    for(k = symbmtx->cblktab[cblknum].bloknum+1;
        k < symbmtx->cblktab[cblknum+1].bloknum; k++)
    {
        N = (double)(symbmtx->bloktab[k].lrownum - symbmtx->bloktab[k].frownum + 1);

#ifdef DOF_CONSTANT
        N *= (double)dofptr->noddval;
#endif

        nbops += FLOPS_DGEMM( M, N, K );

        M -= N;
    }

    return nbops;
}

void symbCost(pastix_int_t *iparm, double *dparm, const SymbolMatrix * symbmtx, const Dof * dofptr)
{
    double flops = 0.;
    printf("SymbolCost: number of operations Cholesky %g \n",
           recursive_sum(0, symbmtx->cblknbr-1, cholesky, symbmtx, dofptr));
    printf("SymbolCost: number of operations Crout2t  %g \n",
           recursive_sum(0, symbmtx->cblknbr-1, crout_2t, symbmtx, dofptr));
    printf("SymbolCost: number of operations Crout3t  %g \n",
           recursive_sum(0, symbmtx->cblknbr-1, crout_3t, symbmtx, dofptr));
    printf("SymbolCost: number of operations CroutHyb %g \n",
           recursive_sum(0, symbmtx->cblknbr-1, crout_hyb, symbmtx, dofptr));
    printf("SymbolCost: number of operations CroutHyb blok %g \n",
           recursive_sum(0, symbmtx->cblknbr-1, crout_blok, symbmtx, dofptr));
    printf("SymbolCost: number of non-zero   %g \n",
           recursive_sum(0, symbmtx->cblknbr-1, nnz, symbmtx, dofptr));

    set_iparm(iparm, IPARM_NNZEROS,   (pastix_int_t)recursive_sum(0, symbmtx->cblknbr-1, nnz,        symbmtx, dofptr));

    if ( iparm[IPARM_FACTORIZATION] == API_FACT_LU ) {
        if ( (iparm[IPARM_FLOAT] == API_COMPLEXDOUBLE) ||
             (iparm[IPARM_FLOAT] == API_COMPLEXSINGLE) ) {
            flops = recursive_sum(0, symbmtx->cblknbr-1, flops_zgetrf, symbmtx, dofptr);
        }
        else {
            flops = recursive_sum(0, symbmtx->cblknbr-1, flops_dgetrf, symbmtx, dofptr);
        }
    } else {
        if ( (iparm[IPARM_FLOAT] == API_COMPLEXDOUBLE) ||
             (iparm[IPARM_FLOAT] == API_COMPLEXSINGLE) ) {
            flops = recursive_sum(0, symbmtx->cblknbr-1, flops_zpotrf, symbmtx, dofptr);
        }
        else {
            flops = recursive_sum(0, symbmtx->cblknbr-1, flops_dpotrf, symbmtx, dofptr);
        }
    }
    printf("SymbolCost: number of operations with new formula %g \n", flops );
    set_dparm(dparm, DPARM_FACT_FLOPS, flops);
}
