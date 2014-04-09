#ifndef _COMPUTE_TRSM_H
#define _COMPUTE_TRSM_H

#include "common.h"
#include "sopalin_compute.h"

/*
 * Definition of the 3 kernels for the 3 different factorizations
 *
 *    m   - The number of rows of L (and U)
 *    n   - The number of columns of L (and U)
 *    dL  - Diagonal block of matrix L
 *    dU  - Diagonal block of matrix U
 *    ldd - Leading dimension of matrix dL (and dU)
 *    L   - extra diagonal block of L
 *    U   - extra diagonal block of U
 *    ldl - Leading dimension of matrix L (and U)
 */
#define kernel_trsm API_CALL(kernel_trsm)
#ifdef CHOL_SOPALIN
#  ifdef SOPALIN_LU

static inline
void kernel_trsm(pastix_int_t m, pastix_int_t n,
               pastix_float_t *dL, pastix_float_t *dU, pastix_int_t ldd,
               pastix_float_t *L,  pastix_float_t *U,  pastix_int_t ldl )
{
  SOPALIN_TRSM("R", "U", "N", "N", m, n,
               fun, dL, ldd, L, ldl);

  SOPALIN_TRSM("R", "U", "N", "U", m, n,
               fun, dU, ldd, U, ldl);
}

#  else

static inline
void kernel_trsm(pastix_int_t m, pastix_int_t n,
               pastix_float_t *dL, pastix_int_t ldd,
               pastix_float_t *L,  pastix_int_t ldl )
{
  SOPALIN_TRSM("R", "L",
               "T",
               "N", m, n,
         fun, dL, ldd, L, ldl);
}
#  endif /* SOPALIN_LU */

#else

static inline
void kernel_trsm(pastix_int_t m, pastix_int_t n,
               pastix_float_t *dL, pastix_int_t ldd,
               pastix_float_t *L,  pastix_float_t *L2, pastix_int_t ldl )
{
  pastix_int_t k;
#ifdef HERMITIAN
  SOPALIN_TRSM("R", "L",
               "C",
               "U", m, n,
               fun, dL, ldd, L, ldl);
#else
  SOPALIN_TRSM("R", "L",
               "T",
               "U", m, n,
               fun, dL, ldd, L, ldl);
#endif

  for (k=0; k<n; k++)
  {
    pastix_float_t alpha;
# ifdef COMPUTE
    ASSERTDBG(dL[k+k*ldd] != 0., MOD_SOPALIN);
    alpha = fun / dL[k+k*ldd];
# endif
    SOPALIN_COPY(m, &(L[ k*ldl]), iun,
            &(L2[k*m  ]), iun);
    SOPALIN_SCAL(m, alpha, &(L[k*ldl]), iun);
  }
}

#endif /* CHOL_SOPALIN */
#endif _COMPUTE_TRSM_H
