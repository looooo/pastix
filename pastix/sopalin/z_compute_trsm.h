/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#ifndef Z_COMPUTE_TRSM_H
#define Z_COMPUTE_TRSM_H

#include "common.h"
#include "z_sopalin_compute.h"

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
#define z_kernel_trsm API_CALL(z_kernel_trsm)
#ifdef CHOL_SOPALIN
#  ifdef SOPALIN_LU

static inline
void z_kernel_trsm(pastix_int_t m, pastix_int_t n,
               pastix_complex64_t *dL, pastix_complex64_t *dU, pastix_int_t ldd,
               pastix_complex64_t *L,  pastix_complex64_t *U,  pastix_int_t ldl )
{
  SOPALIN_TRSM("R", "U", "N", "N", m, n,
               fun, dL, ldd, L, ldl);

  SOPALIN_TRSM("R", "U", "N", "U", m, n,
               fun, dU, ldd, U, ldl);
}

#  else

static inline
void z_kernel_trsm(pastix_int_t m, pastix_int_t n,
               pastix_complex64_t *dL, pastix_int_t ldd,
               pastix_complex64_t *L,  pastix_int_t ldl )
{
  SOPALIN_TRSM("R", "L",
               "T",
               "N", m, n,
         fun, dL, ldd, L, ldl);
}
#  endif /* SOPALIN_LU */

#else

static inline
void z_kernel_trsm(pastix_int_t m, pastix_int_t n,
               pastix_complex64_t *dL, pastix_int_t ldd,
               pastix_complex64_t *L,  pastix_complex64_t *L2, pastix_int_t ldl )
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
    pastix_complex64_t alpha;
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
#endif /*Z_COMPUTE_TRSM_H */
