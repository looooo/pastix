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


/*********************************************************************
 *
 * Trsm executed in 2D version:
 *   all the blocks are grouped together to do only one TRSM
 *
 *   sopalin_data - pointer on data associated to the factorization
 *   me           - Thread id
 *   task         - Task id
 *
 */
void API_CALL(factor_trsm2d)(Sopalin_Data_t *sopalin_data, pastix_int_t me, pastix_int_t task)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  pastix_float_t *dL, *L, *L2;
#ifdef SOPALIN_LU
  pastix_float_t *dU, *U, *U2;
#endif
  pastix_int_t c, b, dima, dimb, stride, offsetED, size;
  (void)me;

  c = TASK_CBLKNUM(task);
  b = TASK_BLOKNUM(task);

  offsetED = SOLV_COEFIND(b);

  /* Address of diagonal block */
  /* Even if the data is local dL = RTASK_COEFTAB */
  dL = (pastix_float_t *)RTASK_COEFTAB(task);

  /* Adress of extra-diagonal block */
  L = &(SOLV_COEFTAB(c)[offsetED]);

  /* horizontal dimension / also leading dimension of dL since dL is square */
  dima = SYMB_LCOLNUM(c) - SYMB_FCOLNUM(c) + 1;
  /* vertical dimension */
  dimb = SYMB_LROWNUM(b) - SYMB_FROWNUM(b) + 1;

  /* Leading dimension of b */
  stride = SOLV_STRIDE(c);

  /* Size of L2/U2 */
#ifdef SOPALIN_LU
  size = dima * dimb * 2;
#else
  size = dima * dimb;
#endif
  MALLOC_INTERN(L2, size, pastix_float_t);

  print_debug(DBG_SOPALIN_ALLOC, "alloc block coeff %x\n",
        (unsigned int)(intptr_t)L2);

  STATS_ADD(size);

  /*
   * Resolution of L_kk * A_jk^t (ga, gb, stride, dimb, dima);
   */
#ifdef CHOL_SOPALIN
#  ifdef SOPALIN_LU
  dU = dL + dima*dima;
  U2 = L2 + dima*dimb;
  U  = &(SOLV_UCOEFTAB(c)[offsetED]);

  API_CALL(kernel_trsm)(dimb, dima, dL, dU, dima, L,  U, stride);

  SOPALIN_LACPY(dimb, dima, L, stride, L2, dimb);
  SOPALIN_LACPY(dimb, dima, U, stride, U2, dimb);
#  else
  API_CALL(kernel_trsm)(dimb, dima, dL,     dima, L,     stride);
  SOPALIN_LACPY(dimb, dima, L, stride, L2, dimb);
#  endif /* SOPALIN_LU */
#else
  API_CALL(kernel_trsm)(dimb, dima, dL,     dima, L, L2, stride);
#endif /* CHOL_SOPALIN */

  /* Save pointer on the copy for future update and send */
  MUTEX_LOCK(&(sopalin_data->mutex_task[task]));
  STASK_COEFTAB(task) = L2;
  MUTEX_UNLOCK(&(sopalin_data->mutex_task[task]));

  /* Free E2 tasks waiting on previous pointer */
  pthread_cond_broadcast(&(sopalin_data->cond_task[task]));
}
