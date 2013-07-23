
/************************************************************/
/**                                                        **/
/**   NAME       : amalgamate.c                            **/
/**                                                        **/
/**   AUTHOR     : Pascal HENON                            **/
/**                                                        **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 15/08/2006      **/
/**                                                        **/
/**                                                        **/
/************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "common.h"
#include "queue.h"
#include "perf.h"

/*#define BLAS_GAIN*/    /** Amalgamation use the best ratio time/nnzadd
                         to merge cblk **/
/*#define BLAS_SOLVE*/   /** Amalgamation seeks to optimize the triangular
                         solve **/

#ifdef BLAS_SOLVE
#define CBLKTIME  cblk_time_fact
#else
#define CBLKTIME  cblk_time_solve
#endif

/** 2 percents **/
#define RAT_CBLK 0.02

#define INFINI 10e6

#define print_one(fmt, ...)    if( procnum == 0) fprintf(stdout, fmt, ##__VA_ARGS__)

double cblk_time_fact (pastix_int_t n, pastix_int_t *ja, pastix_int_t colnbr);
double cblk_time_solve(pastix_int_t n, pastix_int_t *ja, pastix_int_t colnbr);
void   get_son(pastix_int_t node, pastix_int_t *sonindex, pastix_int_t *sontab, pastix_int_t *colweight, pastix_int_t *ns, pastix_int_t *list);


void get_son( pastix_int_t  node,
              pastix_int_t *sonindex,
              pastix_int_t *sontab,
              pastix_int_t *colweight,
              pastix_int_t *ns,
              pastix_int_t *list )
{
    pastix_int_t i, s;
    pastix_int_t nss;
    pastix_int_t ind;
    ind = 0;
    for(i=sonindex[node];i<sonindex[node+1];i++)
    {
        s = sontab[i];
        if(colweight[s] <= 0)
        {
            get_son(s, sonindex, sontab, colweight, &nss, list+ind);
            ind += nss;
        }
        else
            list[ind++] = s;
    }
    *ns = ind;
}

double cblk_time_fact(pastix_int_t n, pastix_int_t *ja, pastix_int_t colnbr)
{
  /*******************************************/
  /* Compute the time to compute a cblk      */
  /* according to the BLAS modelization      */
  /*******************************************/
  double cost;
  pastix_int_t i;
  pastix_int_t L, G, H;

  /** The formula are based on the costfunc.c in blend **/
  /** @@@ OIMBE: il faudra faire les DOF_CONSTANT ***/

  /** Diagonal factorization and TRSM **/
  L = colnbr;
  G = n-L;
#define CHOLESKY
#ifndef CHOLESKY
  cost = (double)(L*PERF_COPY(L)+ PERF_PPF(L) + PERF_TRSM(L, G) + L*PERF_SCAL(G)
                  + L*PERF_COPY(G));
#else
  cost = (double)(PERF_POF(L) + PERF_TRSM(L, G)) ;
#endif

  /** Contributions **/
  i = colnbr;
  while(i<n)
    {
      H = 1;
      i++;
      while(i<n && ja[i] == ja[i-1]+1)
        {
          i++;
          H++;
        }

      cost += (double)(PERF_GEMM(G, H, L));
      G -= H;

    }

  return cost;
}


double cblk_time_solve(pastix_int_t n, pastix_int_t *ja, pastix_int_t colnbr)
{
  double cost;
  pastix_int_t L;
  (void)ja;

  L = colnbr;

  cost = (double)PERF_TRSV(L) + (double) PERF_GEMV(L, n-L);
  return cost;
}
