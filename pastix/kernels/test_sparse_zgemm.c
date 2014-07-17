/**
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <inttypes.h>
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas.h>
#ifdef WITH_MAGMABLAS
#include <magmablas.h>
#endif
#include <math.h>
#include <float.h>

#include "common.h"
#include <cblas.h>
#include "pastix_zcores.h"
#include "timing.h"
/* #include "sparse_zgemm.h" */
#ifdef PRECISIONS_z
#  define TYPE_COMPLEX
#  define PREC_DOUBLE
#elif PRECISIONS_d
#  ifdef TYPE_COMPLEX
#    undef TYPE_COMPLEX
#  endif
#  define PREC_DOUBLE
#elif PRECISIONS_c
#  define TYPE_COMPLEX
#  ifdef PREC_DOUBLE
#    undef PREC_DOUBLE
#  endif
#elif PRECISIONS_s
#  ifdef TYPE_COMPLEX
#    undef TYPE_COMPLEX
#  endif
#  ifdef PREC_DOUBLE
#    undef PREC_DOUBLE
#  endif
#endif
#include "sopalin_compute.h"
#include "pastix_cuda_helper.h"

#if (CUDA_SM_VERSION >= 20)
#include "sparse_zgemm_fermi.h"
#endif

#define TO_STR(x) STR(x)
#define STR(x) #x
#if (defined PRECISION_z || defined PRECISION_d)
#  define FLOAT_EPSILON DBL_EPSILON
#else
#  define FLOAT_EPSILON FLT_EPSILON
#endif

#define SAVE_MIN_TIME(time)                                             \
  do {                                                                  \
    if (run_idx == 0)                                                   \
      min_ ## time [run2_idx] = time;                                   \
    else                                                                \
      min_ ## time [run2_idx] = MIN(time, min_ ## time [run2_idx]);     \
  } while (0)

void usage(char * name)
{
  fprintf(stdout, "usage: %s <nrowsA> <ncolsA> <nrowsA11> [nruns_by_dist n_dist]\n", name);
}

#define READ_INT(m, i) do {                     \
    char *       endptr;                        \
    m = strtol(argv[i], &endptr, 0);            \
    if (argv[i] == endptr)                      \
      {                                         \
        usage(argv[0]);                         \
        return 1;                               \
      }                                         \
  } while(0)

#define FILL(A,n) do {                                                 \
        int fill_i;                                                     \
        for (fill_i = 0; fill_i < n; fill_i++) {                        \
            A[fill_i] = (pastix_complex64_t)(((double)rand())/          \
                                             ((double)RAND_MAX));       \
            /* fprintf(stdout, "%s[%d] = %g\n", #A, fill_i, A[fill_i]); */    \
        }                                                               \
    } while(0)


#define CUDA_CALL(x) do {                               \
    cudaError_t CUAD_CALL_err;                          \
    if (cudaSuccess != (CUAD_CALL_err = x))             \
      {                                                 \
        errorPrint("%s %s (%s,%d)\n",                   \
                   cudaGetErrorString(CUAD_CALL_err),   \
                   #x, __FILE__,__LINE__);              \
      }                                                 \
  } while(0)


#define CUDA_SPARSE_GEMM(TRANSA, TRANSB,                                \
                         dimi, dimj, dima,                              \
                         alpha,                                         \
                         A,  stride_A,                                  \
                         B,  stride_B,                                  \
                         beta,                                          \
                         C, stride_C,                                   \
                         blocknbr, blocktab,                            \
                         fblocknbr, fblocktab)                          \
  do {                                                                  \
      magmablas_sparse_zgemm_kernel_N_T_64_16_4_16_4((int)dimi,         \
                                                     (int)dimj,         \
                                                     (int)dima,         \
                                                     (float)alpha,      \
                                                     A,                 \
                                                     (int)stride_A,     \
                                                      B,               \
                                                     (int)stride_B,    \
                                                     (float)beta,      \
                                                     C,                \
                                                     (int)stride_C,    \
                                                     blocknbr,         \
                                                     blocktab,         \
                                                     fblocknbr,        \
                                                     fblocktab);       \
  } while(0)

#define COMPARE_RES(B1,B2) do {                                         \
    int cmpres_i;                                                       \
    int cmpres_p = 0;                                                   \
    double cmpres_norm =  0.0, cmpres_sum = 0.0;                        \
    double cmpres_maxdiff = 0.0;                                        \
    for (cmpres_i = 0; cmpres_i < ldb*n; cmpres_i++)                    \
      {                                                                 \
        double cmpres_diff = (double)((B1[cmpres_i]-B2[cmpres_i])*      \
                                      conj(B1[cmpres_i]-B2[cmpres_i])); \
        cmpres_norm += cmpres_diff;                                     \
        cmpres_maxdiff = MAX(cmpres_maxdiff, cmpres_diff);              \
        cmpres_sum  += B1[cmpres_i]*conj(B1[cmpres_i]);                 \
        if (cmpres_p < 10 && sqrt(cmpres_diff) > 0.01)                  \
          {                                                             \
            fprintf(stdout,                                             \
                    "%s:%d %s[%d] %.10g %.10g !="                       \
                    " %s[%d] %.10g %.10g (%.10g)\n",                    \
                    __FILE__, __LINE__,                                 \
                    #B1, cmpres_i,                                      \
                    creal(B1[cmpres_i]), cimag(B1[cmpres_i]),           \
                    #B2, cmpres_i,                                      \
                    creal(B2[cmpres_i]), cimag(B2[cmpres_i]),           \
                    cmpres_diff);                                       \
            cmpres_p++;                                                 \
          }                                                             \
      }                                                                 \
    fprintf(stdout, "%d.%d: norm2 = %e\n",                              \
            run2_idx, run_idx, sqrt(cmpres_norm/cmpres_sum));           \
    fprintf(stdout, "%d.%d: normM = %e\n",                              \
            run2_idx, run_idx, sqrt(cmpres_maxdiff));                   \
    fprintf(stdout, "%d.%d: normM/(sum(normA,At,C)*MAX(m,n)*EPS) = %e\n", \
            run2_idx, run_idx, sqrt(cmpres_maxdiff)/                    \
            ((norm_A+norm_At+norm_B)*n*FLOAT_EPSILON));                 \
  } while(0)

#define PRINT_TIME(str, time, ops) do {                         \
    fprintf(stdout,  "%d.%d: " str " : %.2g s, %.2g GFLOPS\n",  \
            run2_idx, run_idx, time, ops/(time*(1<<30)));       \
  } while (0)
#define COMPARE_TIME(str, time, ops, time_ref) do {                     \
    fprintf(stdout,                                                     \
            "%d.%d: " str " : %.2g s, %.2f GFLOPS, Acceleration %.2f\n", \
            run2_idx, run_idx, time, ops/(time*(1<<30)), time_ref/time); \
  } while (0)

static pastix_complex64_t zzero =  0.;

int sparse_zgemm_cpu( int transa, int transb,
                      int m, int n, int k,
                      pastix_complex64_t alpha,
                      const pastix_complex64_t * a, int lda,
                      const pastix_complex64_t * b, int ldb,
                      pastix_complex64_t beta,
                      pastix_complex64_t       * c, unsigned int ldc,
                      int blocknbr,  const int * blocktab,
                      int fblocknbr, const int * fblocktab,
                      pastix_complex64_t *work, int worksize)
{
  int col;
  int   C_index = 0;
  int   W_index = 0;
  const int * fblock;
  const int * lfblock = &(fblocktab[2*(fblocknbr-1)]);
  const int * block;
  const int * lblock  = &(blocktab[2*(blocknbr-1)]);
  (void)worksize;

  cblas_zgemm( CblasColMajor, transa, transb,
               m, n, k,
               CBLAS_SADDR(beta),  a,  lda,
               b,  ldb,
               CBLAS_SADDR(zzero), work, m  );

  for (block = blocktab, fblock = fblocktab;
       block <= lblock;
       W_index += block[1]-block[0]+1, block+=2) {

      int size;
      /* search for facing block */
      while ( ( !( block[0] >= fblock[0] &&
                   block[1] <= fblock[1])) &&
              lfblock >= fblock) {
          C_index += fblock[1]-fblock[0]+1;
          fblock+=2;
      }

      /* Block not found */
      if (lfblock < fblock) {
          fprintf(stderr, "block [%d, %d] not found in facing column.",
                  block[0], block[1]);
          return BADPARAMETER_ERR;
      }

      size = block[1]-block[0]+1;
      /* fprintf(stdout, "%d %d %d %g %d %d %g\n", C_index, block[0], fblock[0], work[W_index], size, C_index + block[0]-fblock[0], alpha); */
      /* fprintf(stdout, "%s:%d sparse_zgemm_cpu %g %g\n", __FILE__, __LINE__, work[W_index], c[(C_index + block[0]-fblock[0])]); */
      /* fprintf(stdout, "%s:%d transa %d %d\n", __FILE__, __LINE__, transa, sizeof(transa)); */
      /* fprintf(stdout, "%s:%d size %d %d\n",   __FILE__, __LINE__, size, sizeof(size)); */
      /* fprintf(stdout, "%s:%d n %d %d\n",      __FILE__, __LINE__, n, sizeof(n)); */
      /* fprintf(stdout, "%s:%d alpha %g %d\n",  __FILE__, __LINE__, alpha, sizeof(alpha)); */
      /* fprintf(stdout, "%s:%d work %p %d\n",  __FILE__, __LINE__, work+W_index, sizeof(work)); */
      /* fprintf(stdout, "%s:%d c %p %d\n",  __FILE__, __LINE__, c+(C_index + block[0]-fblock[0]), sizeof(c)); */
      /* fprintf(stdout, "%s:%d m %d %d\n",      __FILE__, __LINE__, m, sizeof(m)); */
      /* fprintf(stdout, "%s:%d ldc %d %d\n",      __FILE__, __LINE__, ldc, sizeof(ldc)); */
      core_zgeadd( transa, size, n, alpha,
                   &(work[W_index]), m,
                   &(c[C_index + block[0]-fblock[0]]), ldc);
      /* for (col = 0; col < n; col++) */
      /*   SOPALIN_SCAL(size, beta, &(c[C_index + block[0]-fblock[0] + col*ldc]), 1); */
      /* SOPALIN_GEAM("N", "N", size, n, 1.0, */
      /*              &(work[W_index]), m, */
      /*              &(c[C_index + block[0]-fblock[0]]),  ldc); */
      /* fprintf(stdout, "%g\n", c[C_index + block[0]-fblock[0]]); */
    }

  return NO_ERR;
}


int sparse_zgemdm_cpu( int transa, int transb,
                       int m, int n, int k,
                       pastix_complex64_t alpha,
                       const pastix_complex64_t * a, int lda,
                       const pastix_complex64_t * d, int ldd,
                       const pastix_complex64_t * b, int ldb,
                       pastix_complex64_t beta,
                       pastix_complex64_t       * c, unsigned int ldc,
                       int blocknbr,  const int * blocktab,
                       int fblocknbr, const int * fblocktab,
                       pastix_complex64_t *work, int worksize,
                       pastix_complex64_t *work2, int worksize2)
{
  int col;
  int   C_index = 0;
  int   W_index = 0;
  const int * fblock;
  const int * lfblock = &(fblocktab[2*(fblocknbr-1)]);
  const int * block;
  const int * lblock  = &(blocktab[2*(blocknbr-1)]);
  (void)worksize;

  core_zgemdm( transa, transb,
               m, n, k,
               beta,
               a,  lda,
               b,  ldb,
               zzero,
               work, m,
               d,    ldd,
               work2, worksize2);

  for (block = blocktab, fblock = fblocktab;
       block <= lblock;
       W_index += block[1]-block[0]+1, block+=2) {

      int size;
      /* search for facing block */
      while ( ( !( block[0] >= fblock[0] &&
                   block[1] <= fblock[1])) &&
              lfblock >= fblock) {
          C_index += fblock[1]-fblock[0]+1;
          fblock+=2;
      }

      /* Block not found */
      if (lfblock < fblock) {
          fprintf(stderr, "block [%d, %d] not found in facing column.",
                  block[0], block[1]);
          return BADPARAMETER_ERR;
      }

      size = block[1]-block[0]+1;
      /* fprintf(stdout, "%d %d %d %g %d %d %g\n", C_index, block[0], fblock[0], work[W_index], size, C_index + block[0]-fblock[0], alpha); */
      /* fprintf(stdout, "%s:%d sparse_zgemm_cpu %g %g\n", __FILE__, __LINE__, work[W_index], c[(C_index + block[0]-fblock[0])]); */
      /* fprintf(stdout, "%s:%d transa %d %d\n", __FILE__, __LINE__, transa, sizeof(transa)); */
      /* fprintf(stdout, "%s:%d size %d %d\n",   __FILE__, __LINE__, size, sizeof(size)); */
      /* fprintf(stdout, "%s:%d n %d %d\n",      __FILE__, __LINE__, n, sizeof(n)); */
      /* fprintf(stdout, "%s:%d alpha %g %d\n",  __FILE__, __LINE__, alpha, sizeof(alpha)); */
      /* fprintf(stdout, "%s:%d work %p %d\n",  __FILE__, __LINE__, work+W_index, sizeof(work)); */
      /* fprintf(stdout, "%s:%d c %p %d\n",  __FILE__, __LINE__, c+(C_index + block[0]-fblock[0]), sizeof(c)); */
      /* fprintf(stdout, "%s:%d m %d %d\n",      __FILE__, __LINE__, m, sizeof(m)); */
      /* fprintf(stdout, "%s:%d ldc %d %d\n",      __FILE__, __LINE__, ldc, sizeof(ldc)); */
      core_zgeadd( transa, size, n, alpha,
                   &(work[W_index]), m,
                   &(c[C_index + block[0]-fblock[0]]), ldc);
      /* for (col = 0; col < n; col++) */
      /*   SOPALIN_SCAL(size, beta, &(c[C_index + block[0]-fblock[0] + col*ldc]), 1); */
      /* SOPALIN_GEAM("N", "N", size, n, 1.0, */
      /*              &(work[W_index]), m, */
      /*              &(c[C_index + block[0]-fblock[0]]),  ldc); */
      /* fprintf(stdout, "%g\n", c[C_index + block[0]-fblock[0]]); */
    }

  return NO_ERR;
}


double compute_norme(pastix_complex64_t * A1, int lda, int ncol, int nrow)
{
  int i,j;
  double norm = 0.0;
  for (i = 0; i <ncol; i++)
    {
      for(j = 0; j < nrow; j++)
        {
          norm = MAX(norm, sqrt(A1[i*lda + j]*conj(A1[i*lda + j])));
        }
    }
  return norm;
}

#ifdef PASTIX_WITH_MAGMABLAS
#  define CUBLAS_GEMM magmablas_zgemm
#else
#  define CUBLAS_GEMM cublasZgemm
#endif
int
main(int argc, char ** argv) {
  unsigned int  iseed = 0;//(unsigned int)time(NULL);
  long          m,n,k,nrowsA11;
  int  *        blocks;
  int  *        fblocks;
  int  *        d_blocks;
  int  *        d_fblocks;
  int           before_size;
  int           after_size;
  int           nb_blocks_in;
  int           b;
  int           block_size;
  pastix_complex64_t        *A1, *D;
  int           lda;
  pastix_complex64_t        *B0, *B_ref, *B_ref_gemdm, *B_res;
  pastix_complex64_t        *work, *work2;
  int           worksize, worksize2;
  int           ldb;
  int           ldd;
  double        clk[2];
  double        clk_wt[2];
#define myClockInit(clock)  clock[0] = clockGet()
#define myClockStart(clock)  clock[0] = clockGet()
#define myClockStop(clock)  clock[1] = clockGet()
#define myClockVal(clock)   (clock[1] - clock[0])

  pastix_complex64_t alpha = -1.0;
  pastix_complex64_t beta  = 1.0;

  CU_FLOAT cu_alpha;
  CU_FLOAT cu_beta;
  CU_FLOAT *d_A, *d_B, *d_D;

  double time_sparse_CPU,          *min_time_sparse_CPU;
  double time_sparse_GPU,          *min_time_sparse_GPU;
  double time_sparse_GPU_wt,       *min_time_sparse_GPU_wt;
#if (CUDA_SM_VERSION >= 20)
  double time_sparse_GPU_FERMI,    *min_time_sparse_GPU_FERMI;
  double time_sparse_GPU_FERMI_wt, *min_time_sparse_GPU_FERMI_wt;
  double time_sparse_GPU_FERMI_GEMDM,    *min_time_sparse_GPU_FERMI_GEMDM;
  double time_sparse_GPU_FERMI_GEMDM_wt, *min_time_sparse_GPU_FERMI_GEMDM_wt;
#endif
  double time_dense_CPU,           *min_time_dense_CPU;
  double time_dense_GPU,           *min_time_dense_GPU;
  double time_dense_GPU_wt,        *min_time_dense_GPU_wt;
  double ops;
  double norm_A;
  double norm_B;
  double norm_At;
  int    nruns = 1, nruns2 = 1, run_idx, run2_idx;

  cu_alpha = CU_FLOAT_INIT(creal(alpha), cimag(alpha));
  cu_beta  = CU_FLOAT_INIT(creal(beta),  cimag(beta));
  srand (iseed);

  if (argc != 4 && argc != 6)
    {
      usage(argv[0]);
      return 1;
    }

  READ_INT(m, 1);
  READ_INT(k, 2);
  READ_INT(nrowsA11, 3);
  assert(nrowsA11 <= m);
  if (argc == 6)
    {
      READ_INT(nruns, 4);
      READ_INT(nruns2, 5);
    }
  MALLOC_INTERN(min_time_sparse_CPU,    nruns2, double);
  MALLOC_INTERN(min_time_sparse_GPU,    nruns2, double);
  MALLOC_INTERN(min_time_sparse_GPU_wt, nruns2, double);
#if (CUDA_SM_VERSION >= 20)
  MALLOC_INTERN(min_time_sparse_GPU_FERMI,    nruns2, double);
  MALLOC_INTERN(min_time_sparse_GPU_FERMI_wt, nruns2, double);
  MALLOC_INTERN(min_time_sparse_GPU_FERMI_GEMDM,    nruns2, double);
  MALLOC_INTERN(min_time_sparse_GPU_FERMI_GEMDM_wt, nruns2, double);
#endif
  MALLOC_INTERN(min_time_dense_CPU,    nruns2, double);
  MALLOC_INTERN(min_time_dense_GPU,    nruns2, double);
  MALLOC_INTERN(min_time_dense_GPU_wt, nruns2, double);

  for (run2_idx = 0; run2_idx < nruns2; run2_idx++)
    {
      int           nb_blocks  = 0;
      int           size       = 0;
      int           last       = 0;
      int           nb_fblocks = 0;
      /* Build a sparse block column */
      /* The total size of all blocks must be m,
       so there are m blocks max. */
      MALLOC_INTERN(blocks, 2*m, int);
      do {
        block_size = (int)(((double)m*(double)rand())/
                           ((double)RAND_MAX)/10)+1;
        if (nb_blocks == 0)
          block_size =nrowsA11;

        if (size + block_size > m)
          block_size = m - block_size;

        blocks[2*nb_blocks]   = last + (int)(((double)m*(double)rand())/
                                             ((double)RAND_MAX)) +1;
        blocks[2*nb_blocks+1] = blocks[2*nb_blocks] + block_size - 1;
        if ( size + blocks[2*nb_blocks+1] - blocks[2*nb_blocks] + 1 > m)
          blocks[2*nb_blocks+1] = blocks[2*nb_blocks] + m - size - 1;
        fprintf(stdout, "block [%d, %d] ([%d, %d])\n",
                blocks[2*nb_blocks], blocks[2*nb_blocks+1],
                size, size + block_size - 1);
        size += blocks[2*nb_blocks+1] - blocks[2*nb_blocks] + 1;
        last = blocks[2*nb_blocks+1];
        nb_blocks++;
      } while(size < m);

      MALLOC_INTERN(fblocks, 2*m, int);
      b = 0;

      /* Build a facing sparse block column */
      size = 0;
      do {
        nb_blocks_in = (int)((nb_blocks*(double)rand())/
                             ((double)RAND_MAX));

        if (nb_fblocks == 0)
          before_size = blocks[0];
        else
          before_size = blocks[2*b]  - fblocks[2*nb_fblocks-1] - 1;
        before_size = (int)rint(((double)before_size*(double)rand())/
                                ((double)RAND_MAX));
        if (nb_blocks_in == 0 || b + nb_blocks_in >= nb_blocks)
          {
            nb_blocks_in = nb_blocks - b;
            after_size = 100;
          }
        else
          after_size = blocks[2*(b+nb_blocks_in)] -
            blocks[2*(b+nb_blocks_in-1)+1]-1;
        after_size = (int)rint(((double)after_size*(double)rand())/
                               ((double)RAND_MAX));

        fblocks[2*nb_fblocks]   = blocks[2*b] - before_size;
        block_size = blocks[2*(b+nb_blocks_in-1)+1] - blocks[2*b] + 1;
        fblocks[2*nb_fblocks+1] = fblocks[2*nb_fblocks] +
          before_size + block_size + after_size - 1;
        fprintf(stdout, "fblock [%d, %d] ([%d, %d]) %d %d %d %d\n",
                fblocks[2*nb_fblocks], fblocks[2*nb_fblocks+1],
                size , size + before_size + block_size + after_size - 1,
                nb_blocks_in, before_size, block_size, after_size);
	size += before_size + block_size + after_size; 
        b += nb_blocks_in;
        nb_fblocks++;
      } while(b != nb_blocks);

      /* check */
      {
        int fb = 0;
        for (b = 0; b < nb_blocks; b++)
          {
            while (!(blocks[2*b] >= fblocks[2*fb] &&
                     blocks[2*b+1] <= fblocks[2*fb+1]))
              {
                fb++;
                assert(fb < nb_fblocks);
              }
          }
      }

      /* allocate and fill the matrices */
      lda = m;
      ldd = lda;
      ldb = fblocks[2*(nb_fblocks-1)+1]-fblocks[0]+1;
      n = blocks[1] - blocks[0] +1;
      MALLOC_INTERN(A1, lda*k, pastix_complex64_t);
      MALLOC_INTERN(D,  k*ldd,     pastix_complex64_t);
      MALLOC_INTERN(B0, ldb*n, pastix_complex64_t);
      memset(B0, 0, ldb*n*sizeof(pastix_complex64_t));
      MALLOC_INTERN(B_ref, ldb*n, pastix_complex64_t);
      MALLOC_INTERN(B_ref_gemdm, ldb*n, pastix_complex64_t);
      MALLOC_INTERN(B_res, ldb*n, pastix_complex64_t);
      MALLOC_INTERN(work, m*n, pastix_complex64_t);
      worksize = m*n;
      worksize2 = (MAX(m,n)+1)*k;
      MALLOC_INTERN(work2, worksize2, pastix_complex64_t);
      FILL(A1, lda*k);
      FILL(B0, ldb*n);
      norm_A  = compute_norme(A1, lda, k, lda);
      norm_At = compute_norme(A1, lda, k, blocks[1]-blocks[0]+1);
      norm_B  = compute_norme(B0, ldb, n, ldb);
      {
        int i,j;
        for (i = 0; i < ldd; i++)
	  for (j = 0; j < k; j++)
	    D[j*ldd + i] = 2*(i==j);
      }
      
      for (run_idx = 0; run_idx < nruns; run_idx++)
        {
          memcpy(B_ref, B0, ldb*n*sizeof(pastix_complex64_t));
          myClockInit(clk);
          myClockStart(clk);
          sparse_zgemm_cpu(CblasNoTrans, CblasConjTrans,
                           m, n, k,
                           alpha,
                           A1, lda,
                           A1, lda,
                           beta,
                           B_ref, ldb,
                           nb_blocks,  blocks,
                           nb_fblocks, fblocks,
                           work, worksize);
          myClockStop(clk);
          time_sparse_CPU = myClockVal(clk);
          ops = m*n*k*2;
          SAVE_MIN_TIME(time_sparse_CPU);
          PRINT_TIME("sparse ZGEMM on CPU", time_sparse_CPU, ops);



          memcpy(B_ref_gemdm, B0, ldb*n*sizeof(pastix_complex64_t));
          myClockInit(clk);
          myClockStart(clk);
          sparse_zgemdm_cpu(CblasNoTrans, CblasConjTrans,
                            m, n, k,
                            alpha,
                            A1, lda,
                            D, ldd+1,
                            A1, lda,
                            beta,
                            B_ref_gemdm, ldb,
                            nb_blocks,  blocks,
                            nb_fblocks, fblocks,
                            work, worksize,
                            work2, worksize2);
          myClockStop(clk);
          time_sparse_CPU = myClockVal(clk);
          ops = m*n*k*2;
          SAVE_MIN_TIME(time_sparse_CPU);
          PRINT_TIME("sparse ZGEMDM on CPU", time_sparse_CPU, ops);


          memcpy(B_res, B0, ldb*n*sizeof(pastix_complex64_t));
          if (CUDA_SUCCESS != cuInit(0))
            {
              errorPrint("cuInit()");
              assert(0);
            }
          CUDA_CALL(cudaSetDevice(0));

          CUDA_CALL(cudaMalloc((void*)&(d_blocks),
                               2*nb_blocks*sizeof(int)));
          CUDA_CALL(cudaMemcpy((void*)d_blocks, blocks,
                               2*nb_blocks*sizeof(int),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaMalloc((void*)&(d_fblocks),
                               2*nb_fblocks*sizeof(int)));
          CUDA_CALL(cudaMemcpy((void*)d_fblocks, fblocks,
                               2*nb_fblocks*sizeof(int),
                               cudaMemcpyHostToDevice));
/* #if (CUDA_SM_VERSION >= 20 || !(defined PREC_DOUBLE && defined TYPE_COMPLEX)) */
/*           if (nruns <  10) */
/*             fprintf(stdout, ">>> %s <<<\n", */
/*                     TO_STR(PASTIX_PREFIX_F(magmablas_sparse_zgemm_kernel_N_T_64_16_4_16_4))); */

/*           myClockInit(clk_wt); */
/*           myClockStart(clk_wt); */
/*           CUDA_CALL(cudaMalloc((void*)&(d_A), */
/*                                lda*k*sizeof(pastix_complex64_t))); */
/*           CUDA_CALL(cudaMemcpy((void*)d_A, A1, */
/*                                lda*k*sizeof(pastix_complex64_t), */
/*                                cudaMemcpyHostToDevice)); */
/*           CUDA_CALL(cudaMalloc((void*)&(d_B), */
/*                                ldb*n*sizeof(pastix_complex64_t))); */
/*           CUDA_CALL(cudaMemcpy((void*)d_B, B0, */
/*                                ldb*n*sizeof(pastix_complex64_t), */
/*                                cudaMemcpyHostToDevice)); */

/*           CUDA_CALL(cudaThreadSynchronize()); */
/*           myClockInit(clk); */
/*           myClockStart(clk); */
/*           CUDA_SPARSE_GEMM("N", "C", */
/*                            m, n, k, */
/*                            alpha, */
/*                            (pastix_complex64_t*)d_A, lda, */
/*                            (pastix_complex64_t*)d_A, lda, */
/*                            beta, */
/*                            (pastix_complex64_t*)d_B, ldb, */
/*                            nb_blocks,  d_blocks, */
/*                            nb_fblocks, d_fblocks); */
/*           CUDA_CALL(cudaThreadSynchronize()); */

/*           myClockStop(clk); */

/*           CUDA_CALL(cudaMemcpy((void*)B_res, d_B, */
/*                                ldb*n*sizeof(pastix_complex64_t), */
/*                                cudaMemcpyDeviceToHost)); */
/*           CUDA_CALL(cudaFree(d_A)); */
/*           CUDA_CALL(cudaFree(d_B)); */
/*           myClockStop(clk_wt); */

/*           time_sparse_GPU = myClockVal(clk); */
/*           SAVE_MIN_TIME(time_sparse_GPU); */
/*           time_sparse_GPU_wt = myClockVal(clk_wt); */
/*           SAVE_MIN_TIME(time_sparse_GPU_wt); */
/*           COMPARE_TIME("sparse ZGEMM on GPU (SM < 20)", */
/*                        time_sparse_GPU, ops, time_sparse_CPU); */
/*           COMPARE_TIME("sparse ZGEMM on GPU (SM < 20) with transfer", */
/*                        time_sparse_GPU_wt, ops, time_sparse_CPU); */
/*           COMPARE_RES(B_ref, B_res); */

/* #endif */




#if (CUDA_SM_VERSION >= 20)
          if (nruns <  10)
            fprintf(stdout, ">>> %s <<< (FERMI)\n",
                    TO_STR(GENERATE_SM_VERSION_NAME(gemm)));

          /* FERMI kernel */
          myClockInit(clk_wt);
          myClockStart(clk_wt);
          CUDA_CALL(cudaMalloc((void*)&(d_A),
                               lda*k*sizeof(pastix_complex64_t)));
          CUDA_CALL(cudaMemcpy((void*)d_A, A1,
                               lda*k*sizeof(pastix_complex64_t),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaMalloc((void*)&(d_B),
                               ldb*n*sizeof(pastix_complex64_t)));
          CUDA_CALL(cudaMemcpy((void*)d_B, B0,
                               ldb*n*sizeof(pastix_complex64_t),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaThreadSynchronize());
          myClockInit(clk);
          myClockStart(clk);
          GENERATE_SM_VERSION_NAME(gemm)('N',
                                         'C',
                                         m, n, k,
                                         cu_alpha,
                                         d_A, lda,
                                         d_A, lda,
                                         cu_beta,
                                         d_B, ldb,
                                         nb_blocks,  d_blocks,
                                         nb_fblocks, d_fblocks,
                                         0 );
          CUDA_CALL(cudaThreadSynchronize());
          myClockStop(clk);
          CUDA_CALL(cudaMemcpy((void*)B_res, d_B,
                               ldb*n*sizeof(pastix_complex64_t),
                               cudaMemcpyDeviceToHost));
          CUDA_CALL(cudaFree(d_A));
          CUDA_CALL(cudaFree(d_B));
          myClockStop(clk_wt);
          time_sparse_GPU_FERMI = myClockVal(clk);
          SAVE_MIN_TIME(time_sparse_GPU_FERMI);
          time_sparse_GPU_FERMI_wt = myClockVal(clk_wt);
          SAVE_MIN_TIME(time_sparse_GPU_FERMI_wt);
          COMPARE_TIME("sparse ZGEMM on GPU (Fermi)",
                       time_sparse_GPU_FERMI, ops, time_sparse_CPU);
          COMPARE_TIME("sparse ZGEMM on GPU (Fermi) with transfert",
                       time_sparse_GPU_FERMI_wt, ops, time_sparse_CPU);
          COMPARE_RES(B_ref,B_res);


          if (nruns <  10)
            fprintf(stdout, ">>> %s <<< (FERMI)\n",
                    TO_STR(GENERATE_SM_VERSION_NAME(gemdm)));



          /* FERMI kernel */
          myClockInit(clk_wt);
          myClockStart(clk_wt);
          CUDA_CALL(cudaMalloc((void*)&(d_A),
                               lda*k*sizeof(pastix_complex64_t)));
          CUDA_CALL(cudaMemcpy((void*)d_A, A1,
                               lda*k*sizeof(pastix_complex64_t),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaMalloc((void*)&(d_D),
                               k*ldd*sizeof(pastix_complex64_t)));
          CUDA_CALL(cudaMemcpy(d_D, D,
                               k*ldd*sizeof(pastix_complex64_t),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaMalloc((void*)&(d_B),
                               ldb*n*sizeof(pastix_complex64_t)));
          CUDA_CALL(cudaMemcpy((void*)d_B, B0,
                               ldb*n*sizeof(pastix_complex64_t),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaThreadSynchronize());
          myClockInit(clk);
          myClockStart(clk);
          GENERATE_SM_VERSION_NAME(gemdm)('N',
                                          'C',
                                          m, n, k,
                                          cu_alpha,
                                          d_A, lda,
                                          d_D, ldd,
                                          d_A, lda,
                                          cu_beta,
                                          d_B, ldb,
                                          nb_blocks,  d_blocks,
                                          nb_fblocks, d_fblocks,
                                          0 );
	  CUDA_CALL(cudaThreadSynchronize());
          myClockStop(clk);
          CUDA_CALL(cudaMemcpy((void*)B_res, d_B,
                               ldb*n*sizeof(pastix_complex64_t),
                               cudaMemcpyDeviceToHost));
          CUDA_CALL(cudaFree(d_A));
          CUDA_CALL(cudaFree(d_B));
          myClockStop(clk_wt);
          time_sparse_GPU_FERMI_GEMDM = myClockVal(clk);
          SAVE_MIN_TIME(time_sparse_GPU_FERMI_GEMDM);
          time_sparse_GPU_FERMI_GEMDM_wt = myClockVal(clk_wt);
          SAVE_MIN_TIME(time_sparse_GPU_FERMI_GEMDM_wt);
          COMPARE_TIME("sparse ZGEMDM on GPU (Fermi)",
                       time_sparse_GPU_FERMI_GEMDM, ops, time_sparse_CPU);
          COMPARE_TIME("sparse ZGEMDM on GPU (Fermi) with transfert",
                       time_sparse_GPU_FERMI_GEMDM_wt, ops, time_sparse_CPU);
          COMPARE_RES(B_ref_gemdm,B_res);
#endif




          CUDA_CALL(cudaFree(d_blocks));
          CUDA_CALL(cudaFree(d_fblocks));
          /* Just for timing, perform dense GEMM */

          myClockInit(clk);
          myClockStart(clk);
          cblas_zgemm( CblasColMajor, CblasNoTrans, CblasConjTrans,
                        m, n, k,
                        CBLAS_SADDR(alpha),
                        A1, lda,
                        A1, lda,
                        CBLAS_SADDR(beta),
                        B_ref, ldb);
          myClockStop(clk);
          time_dense_CPU = myClockVal(clk);
          SAVE_MIN_TIME(time_dense_CPU);
          PRINT_TIME("dense ZGEMM on CPU", time_dense_CPU, ops);

          myClockInit(clk_wt);
          myClockStart(clk_wt);
          CUDA_CALL(cudaMalloc((void*)&(d_A),
                               lda*k*sizeof(pastix_complex64_t)));
          CUDA_CALL(cudaMemcpy((void*)d_A, A1,
                               lda*k*sizeof(pastix_complex64_t),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaMalloc((void*)&(d_B),
                               ldb*n*sizeof(pastix_complex64_t)));
          CUDA_CALL(cudaMemcpy((void*)d_B, B0,
                               ldb*n*sizeof(pastix_complex64_t),
                               cudaMemcpyHostToDevice));
          CUDA_CALL(cudaThreadSynchronize());

          myClockInit(clk);
          myClockStart(clk);
          CUDA_CALL(cudaThreadSynchronize());
          CUBLAS_GEMM('N', 'C', m, n, k, cu_alpha,
                      (CU_FLOAT*)d_A, lda, (CU_FLOAT*)d_A, lda,
                      cu_beta, (CU_FLOAT*)d_B, ldb);
          CUDA_CALL(cudaThreadSynchronize());
          myClockStop(clk);


          CUDA_CALL(cudaMemcpy((void*)B_res, d_B,
                               ldb*n*sizeof(pastix_complex64_t),
                               cudaMemcpyDeviceToHost));
          CUDA_CALL(cudaFree(d_A));
          CUDA_CALL(cudaFree(d_B));
          myClockStop(clk_wt);

          time_dense_GPU = myClockVal(clk);
          SAVE_MIN_TIME(time_dense_GPU);
          time_dense_GPU_wt = myClockVal(clk_wt);
          SAVE_MIN_TIME(time_dense_GPU_wt);
#ifdef WITH_MAGMABLAS
          COMPARE_TIME("dense magZGEMM on GPU", time_dense_GPU, ops, time_dense_CPU);
          COMPARE_TIME("dense magZGEMM on GPU with transfert", time_dense_GPU_wt, ops, time_dense_CPU);
#else
          COMPARE_TIME("dense cuZGEMM on GPU", time_dense_GPU, ops, time_dense_CPU);
          COMPARE_TIME("dense cuZGEMM on GPU with transfert", time_dense_GPU_wt, ops, time_dense_CPU);
#endif


        }
      PRINT_TIME("(min) sparse ZGEMM on CPU", min_time_sparse_CPU[run2_idx], ops);
      /* COMPARE_TIME("(min) sparse ZGEMM in GPU", min_time_sparse_GPU[run2_idx], ops, min_time_sparse_CPU[run2_idx]); */
      /* COMPARE_TIME("(min) sparse ZGEMM in GPU with transfert", min_time_sparse_GPU_wt[run2_idx], ops, min_time_sparse_CPU[run2_idx]); */
#if (CUDA_SM_VERSION >= 20)
      COMPARE_TIME("(min) sparse ZGEMM in GPU (FERMI)", min_time_sparse_GPU_FERMI[run2_idx], ops, min_time_sparse_CPU[run2_idx]);
      COMPARE_TIME("(min) sparse ZGEMM in GPU (FERMI) with transfert", min_time_sparse_GPU_FERMI_wt[run2_idx], ops, min_time_sparse_CPU[run2_idx]);
      COMPARE_TIME("(min) sparse ZGEMDM in GPU (FERMI)", min_time_sparse_GPU_FERMI_GEMDM[run2_idx], ops, min_time_sparse_CPU[run2_idx]);
      COMPARE_TIME("(min) sparse ZGEMDM in GPU (FERMI) with transfert", min_time_sparse_GPU_FERMI_GEMDM_wt[run2_idx], ops, min_time_sparse_CPU[run2_idx]);
#endif
      PRINT_TIME("(min) dense ZGEMM on CPU", min_time_dense_CPU[run2_idx], ops);
      COMPARE_TIME("(min) dense ZGEMM on GPU", min_time_dense_GPU[run2_idx], ops, min_time_dense_CPU[run2_idx]);
      COMPARE_TIME("(min) dense ZGEMM on GPU with transfert", min_time_dense_GPU_wt[run2_idx], ops, min_time_dense_CPU[run2_idx]);


      memFree_null(A1);
      memFree_null(B0);
      memFree_null(B_res);
      memFree_null(B_ref);
      memFree_null(B_ref_gemdm);
      memFree_null(blocks);
      memFree_null(fblocks);
    }
  memFree_null(min_time_sparse_CPU);
  /* memFree_null(min_time_sparse_GPU); */
  /* memFree_null(min_time_sparse_GPU_wt); */
#if (CUDA_SM_VERSION >= 20)
  memFree_null(min_time_sparse_GPU_FERMI);
  memFree_null(min_time_sparse_GPU_FERMI_wt);
#endif
  memFree_null(min_time_dense_CPU);
  memFree_null(min_time_dense_GPU);
  memFree_null(min_time_dense_GPU_wt);
  
  return EXIT_SUCCESS;
}
