/*
  File: compute_trsm.c

  Computation functions.

  Pierre Ramet    : fev 2003
  Mathieu Faverge
  Xavier Lacoste

*/

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
#ifdef CHOL_SOPALIN
#  ifdef SOPALIN_LU

static inline
void API_CALL(kernel_trsm)(pastix_int_t m, pastix_int_t n,
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
void API_CALL(kernel_trsm)(pastix_int_t m, pastix_int_t n,
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
void API_CALL(kernel_trsm)(pastix_int_t m, pastix_int_t n,
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

/*********************************************************************
 *
 * Trsm executed in 1D version, or in 1D+esp:
 *   all the blocks are grouped together to do only one TRSM
 *
 *   sopalin_data - pointer on data associated to the factorization
 *   me           - Thread id
 *   c            - Block column id
 *
 */
void API_CALL(factor_trsm1d)(Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t c)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  pastix_float_t *dL, *L;
#ifdef SOPALIN_LU
  pastix_float_t *dU, *U;
#endif
#ifndef CHOL_SOPALIN
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  pastix_float_t *L2;
#endif
  pastix_int_t dima, dimb, stride, offsetD, offsetED;
  (void)me;

  offsetD  = SOLV_COEFIND(SYMB_BLOKNUM(c));
  offsetED = SOLV_COEFIND(SYMB_BLOKNUM(c)+1);

  /* diagonal column block address */
  dL = &(SOLV_COEFTAB(c)[offsetD]);

  /* first extra-diagonal bloc in column block address */
  L  = &(SOLV_COEFTAB(c)[offsetED]);

  stride = SOLV_STRIDE(c);

  /* horizontal dimension */
  dima = SYMB_LCOLNUM(c) - SYMB_FCOLNUM(c) + 1;
  /* vertical dimension */
  dimb = stride - dima;

#ifdef CHOL_SOPALIN
#  ifdef SOPALIN_LU
  dU = &(SOLV_UCOEFTAB(c)[offsetD ]);
  U  = &(SOLV_UCOEFTAB(c)[offsetED]);
  API_CALL(kernel_trsm)(dimb, dima, dL, dU, stride, L,  U, stride);
#  else
  API_CALL(kernel_trsm)(dimb, dima, dL,     stride, L,     stride);
#  endif /* SOPALIN_LU */
#else
  L2 = thread_data->maxbloktab1;
  ASSERTDBG(SOLV_COEFMAX >= dimb*dima, MOD_SOPALIN);
  API_CALL(kernel_trsm)(dimb, dima, dL,     stride, L, L2, stride);
#endif /* CHOL_SOPALIN */
}
