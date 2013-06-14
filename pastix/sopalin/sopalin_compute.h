/*
 * File: sopalin_compute.h
 *
 * Blas interface definition for LINUX/GNU operating systems.
 *
 * Sensible to those defines :
 *   PREC_DOUBLE    - Determine if coefficients are doubles or not.
 *   TYPE_COMPLEX   - Determine if coefficients are complex or not.
 *   pastix_int_t            - Defined in common_pastix.h to determine integer type.
 *   CBLAS          - Determine if we want to use cblas or not.
 *   DEBUG_NAN      - If we want to check for NaN values in BLAS operartions,
 *                    must be used with DBG_SOPALIN_NAN set to 1.
 *   BLAS_USE_COPY  - To use memcpy instead of BLAS copy function.
 *   BLAS_USE_SCAL  - To use "by hand" scal macro instead of BLAS one.
 *
 * Authors:
 *   Ramet  Pierre  - ramet@labri.fr
 *   Xavier Lacoste - lacoste@labri.fr
 */
#ifndef SOPALIN_COMPUTE_H
#define SOPALIN_COMPUTE_H
#ifndef X_INCLUDE_ESSL
/*#define CBLAS*/
#endif /* X_INCLUDE_ESSL */

/*
 *  Define: BLAS_INT
 *  integer to use in blas calls.
 */
#define BLAS_INT pastix_int_t

/*
 Define: BLAS macros

 Defines to set the correct blas depending on the
 coefficients type and the set of cblas usage or not.

 */
#ifdef X_INCLUDE_ESSL
#  ifdef PREC_DOUBLE
#    ifdef TYPE_COMPLEX
#      ifdef CBLAS
#        define BLAS_CALL(func) cblas_z ## func
#      else
#        define BLAS_CALL(func)       z ## func
#      endif
#    else
#      ifdef CBLAS
#        define BLAS_CALL(func) cblas_d ## func
#      else
#        define BLAS_CALL(func)       d ## func
#      endif
#    endif
#  else
#    ifdef TYPE_COMPLEX
#      ifdef CBLAS
#        define BLAS_CALL(func) cblas_c ## func
#      else
#        define BLAS_CALL(func)       c ## func
#      endif
#    else
#      ifdef CBLAS
#        define BLAS_CALL(func) cblas_s ## func
#      else
#        define BLAS_CALL(func)       s ## func
#      endif
#    endif
#  endif
#else /* not X_INCLUDE_ESSL */
#  ifdef PREC_DOUBLE
#    ifdef TYPE_COMPLEX
#      ifdef CBLAS
#        define BLAS_CALL(func) cblas_z ## func
#      else
#        define BLAS_CALL(func)       z ## func ## _
#      endif
#    else
#      ifdef CBLAS
#        define BLAS_CALL(func) cblas_d ## func
#      else
#        define BLAS_CALL(func)       d ## func ## _
#      endif
#    endif
#  else
#    ifdef TYPE_COMPLEX
#      ifdef CBLAS
#        define BLAS_CALL(func) cblas_c ## func
#      else
#        define BLAS_CALL(func)       c ## func ## _
#      endif
#    else
#      ifdef CBLAS
#        define BLAS_CALL(func) cblas_s ## func
#      else
#        define BLAS_CALL(func)       s ## func ## _
#      endif
#    endif
#  endif
#endif /* not X_INCLUDE_ESSL */

#ifdef X_INCLUDE_ESSL
#define     GESUB BLAS_CALL(gesub)
#define     GEADD BLAS_CALL(geadd)
#  define   AXPY(n,al,x,ix,y,iy)                \
  BLAS_CALL(axpy)((*n),(*al),x,(*ix),y,(*iy))
#  define   GEMM(ta,tb,m,n,k,al,a,lda,b,ldb,be,c,ldc) \
  BLAS_CALL(gemm)(ta,tb,(*m),(*n),(*k),(*al),a,(*lda),b,(*ldb),(*be),c,(*ldc))
#  define   GEMV(t,m,n,al,a,lda,x,ix,be,y,iy)   \
  BLAS_CALL(gemv)(t,(*m),(*n),(*al),a,(*lda),x,(*ix),(*be),y,(*iy))
#  define   TRSM(s,u,t,d,m,n,al,a,lda,b,ldb)    \
  BLAS_CALL(trsm)(s,u,t,d,(*m),(*n),(*al),a,(*lda),b,(*ldb))
#  define   TRSV(u,t,d,n,a,lda,x,ix)            \
  BLAS_CALL(trsv)(u,t,d,(*n),a,(*lda),x,(*ix))
#  define   POF(c,p,i1,i2)                      \
  BLAS_CALL(pof)(c,p,(*i1),(*i2))
#  define   PPF(f,i,j)                          \
  BLAS_CALL(ppf)(f,(*i),(*j)
#  ifdef TYPE_COMPLEX
#    define GER(m,n,al,x,ix,y,iy,a,lda)         \
  BLAS_CALL(geru)((*m),(*n),(*al),x,(*ix),y,(*iy),a,(*lda))
#    define SYR(u,n,al,x,ix,a,lda)              \
  BLAS_CALL(her)(u,(*n),(*al),x,(*ix),a,(*lda))
#    define SYRK(u,t,n,k,al,a,lda,be,c,ldc)     \
  BLAS_CALL(herk)(u,t,(*n),(*k),(*al),a,(*lda),(*be),c,(*ldc))
#  else  /* not TYPE_COMPLEX */
#    define GER(m,n,al,x,ix,y,iy,a,lda)         \
  BLAS_CALL(ger)((*m),(*n),(*al),x,(*ix),y,(*iy),a,(*lda))
#    define SYR(u,n,al,x,ix,a,lda)              \
  BLAS_CALL(syr)(u,(*n),(*al),x,(*ix),a,(*lda))
#    define SYRK(u,t,n,k,al,a,lda,be,c,ldc)     \
  BLAS_CALL(syrk)(u,t,(*n),(*k),(*al),a,(*lda),(*be),c,(*ldc))
#  endif /* not TYPE_COMPLEX */
#  define   COPY(n,x,ix,y,iy)                   \
  BLAS_CALL(copy)((*n),x,(*ix),y,(*iy))
#  define   SCAL(n,al,x,ix)                     \
  BLAS_CALL(scal)((*n),(*al),x,(*ix))
#  define DEADD BLAS_CALL(geadd)
#else /* not X_INCLUDE_ESSL */
#  define   AXPY  BLAS_CALL(axpy)
#  define   GEMM  BLAS_CALL(gemm)
#  define   GEMV  BLAS_CALL(gemv)
#  define   TRSM  BLAS_CALL(trsm)
#  define   TRSV  BLAS_CALL(trsv)
#  define   POF   BLAS_CALL(pof)
#  define   PPF   BLAS_CALL(ppf)
#  ifdef TYPE_COMPLEX
#    define GER   BLAS_CALL(geru)
#    define SYR   BLAS_CALL(her)
#    define SYRK  BLAS_CALL(herk)
#  else  /* not TYPE_COMPLEX */
#    define GER   BLAS_CALL(ger)
#    define SYR   BLAS_CALL(syr)
#    define SYRK  BLAS_CALL(syrk)
#  endif /* not TYPE_COMPLEX */
#  define   COPY  BLAS_CALL(copy)
#  define   SCAL  BLAS_CALL(scal)
extern void SYR(char * i, BLAS_INT * n , BLAS_REAL * x, BLAS_FLOAT * a,
                BLAS_INT * u, BLAS_FLOAT * b, BLAS_INT * v);
extern void SYRK (char * i, char * t,  BLAS_INT * n, BLAS_INT * k, BLAS_REAL * x,
                  BLAS_FLOAT * a, BLAS_INT * u, BLAS_REAL * y, BLAS_FLOAT * b, BLAS_INT *v);
extern void COPY(BLAS_INT *n, BLAS_FLOAT *x, BLAS_INT * u, BLAS_FLOAT * y, BLAS_INT * v);
extern void SCAL(BLAS_INT * n , BLAS_FLOAT * a, BLAS_FLOAT * x, BLAS_INT * u);
extern void TRSM (char * s, char * p, char * t, char * d, BLAS_INT * m,
                  BLAS_INT * n , BLAS_FLOAT * x, BLAS_FLOAT * a, BLAS_INT * u,
                  BLAS_FLOAT * b, BLAS_INT * v);
extern void TRSV(char *, char *, char *, BLAS_INT *, BLAS_FLOAT *, BLAS_INT *,
                 BLAS_FLOAT *, BLAS_INT *);
extern void GEMM(char * i, char * j, BLAS_INT * m, BLAS_INT * n,
                 BLAS_INT * k, BLAS_FLOAT * x, BLAS_FLOAT * a, BLAS_INT * u,
                 BLAS_FLOAT * b, BLAS_INT * v, BLAS_FLOAT * y, BLAS_FLOAT * c, BLAS_INT * w);
extern void GEMV(char *, BLAS_INT *, BLAS_INT *, BLAS_FLOAT*, BLAS_FLOAT*, BLAS_INT*,
                 BLAS_FLOAT*,BLAS_INT*,BLAS_FLOAT*,BLAS_FLOAT*,BLAS_INT*);
extern void AXPY(BLAS_INT * n , BLAS_FLOAT * a, BLAS_FLOAT * x, BLAS_INT *  ix,
                 BLAS_FLOAT * y, BLAS_INT * iy);
extern void GER(BLAS_INT * m, BLAS_INT * n, BLAS_FLOAT * x, BLAS_FLOAT * a,
                BLAS_INT * u, BLAS_FLOAT * b, BLAS_INT * v, BLAS_FLOAT * c,
                BLAS_INT * w );
#endif /* not X_INCLUDE_ESSL */

void dim_dgeam(char *transa,char *transb,
               pastix_int_t m,pastix_int_t n,pastix_float_t alpha,
               pastix_float_t *a, pastix_int_t lda,
               pastix_float_t *b,pastix_int_t ldb);



/*
 * Macro: SOPALIN_GEAM
 *
 * Computes $ b = x \times a + b $.
 *
 * Parameters:
 *    i - indicates if a needs to be transposed.
 *    j - indicates if b needs to be transposed.
 *    m - number of row in a and b.
 *    n - number of colonnes in a and b.
 *    x - scalar.
 *    a - Matrix a.
 *    u - Stride between 2 columns of a.
 *    b - Matrix b.
 *    v - Stride between 2 columns of b.
 */
#ifdef X_INCLUDE_ESSL
#  define SOPALIN_GEAM(i,j,m,n,x,a,u,b,v) do {  \
    TAB_CHECK_NAN(((pastix_float_t*)a), n*m);            \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*m);            \
    GEADD((BLAS_FLOAT *) a,u,i,                 \
          (BLAS_FLOAT *) b,v,j,                 \
          (BLAS_FLOAT *) b,v,m,n);              \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*m);            \
  } while (0)

#else
/* Pb pas de truc d'addition C=A+B */
#  define SOPALIN_GEAM(i,j,m,n,x,a,u,b,v) do {  \
    TAB_CHECK_NAN(((pastix_float_t*)a), n*m);            \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*m);            \
    dim_dgeam(i,j,m,n,x,a,u,b,v);               \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*m);            \
  } while (0)
#endif

/* Macro: SOPALIN_GESM
 *
 * Computes $ b = x \times a - b $.
 *
 * Parameters:
 *    i - indicates if a needs to be transposed.
 *    j - indicates if b needs to be transposed.
 *    m - number of row in a and b.
 *    n - number of colonnes in a and b.
 *    x - scalar.
 *    a - Matrix a.
 *    u - Stride between 2 columns of a.
 *    b - Matrix b.
 *    v - Stride between 2 columns of b.
 */
#ifdef X_INCLUDE_ESSL
#  define SOPALIN_GESM(i,j,m,n,x,a,u,b,v) do {  \
    TAB_CHECK_NAN(((pastix_float_t*)a), m*n);            \
    TAB_CHECK_NAN(((pastix_float_t*)b), m*n);            \
    GESUB((BLAS_FLOAT *)b,v,j,                  \
          (BLAS_FLOAT *)a,u,i,                  \
          (BLAS_FLOAT *)b,v,m,n);               \
    TAB_CHECK_NAN(((pastix_float_t*)b), m*n);            \
  } while (0)
#else
/* Pb pas de truc de soustraction C=A-B */
#  define SOPALIN_GESM(i,j,m,n,x,a,u,b,v) do {  \
    TAB_CHECK_NAN(((pastix_float_t*)a), m*n);            \
    TAB_CHECK_NAN(((pastix_float_t*)b), m*n);            \
    dim_dgeam(i,j,m,n,-(x),a,u,b,v);            \
    TAB_CHECK_NAN(((pastix_float_t*)b), m*n);            \
  } while (0)
#endif
/*
 *  Macro: SOPALIN_COPY
 *
 *  Copy *x* vector to *y* vector.
 *
 *  Care: if BLAS_USE_COPY is defined, only works
 *    with *u* and *v* equal to 1.
 *
 *  Parameters:
 *    n - Size of the vectors.
 *    x - Source vector
 *    u - *x* vector stride between to cefficients.
 *    y - Destination vector.
 *    v - *y* vector stride between to cefficients.
 */
#ifdef BLAS_USE_COPY
#  define SOPALIN_COPY(n,x,u,y,v) {             \
    memcpy((void *)(y), (void *)(x),            \
           ((pastix_int_t)(n))*sizeof(pastix_float_t));           \
  }
#else /* BLAS_USE_COPY */

#  define SOPALIN_COPY(n,x,u,y,v) {             \
    BLAS_INT varin = (BLAS_INT)(n);             \
    BLAS_INT variu = (BLAS_INT)(u);             \
    BLAS_INT variv = (BLAS_INT)(v);             \
    TAB_CHECK_NAN(((pastix_float_t*)x), n);              \
    TAB_CHECK_NAN(((pastix_float_t*)y), n);              \
    COPY(&varin,                                \
         (BLAS_FLOAT*)(x), &variu,              \
         (BLAS_FLOAT*)(y), &variv);             \
    TAB_CHECK_NAN(((pastix_float_t*)y), n);              \
  }
#endif /* BLAS_USE_COPY */

/*
 * Macro: SOPALIN_LACPY
 *
 * Equivalent of lacpy form LAPACK with uplo = 'General'
 *
 * Care: if BLAS_USE_COPY is defined, only works
 *   with *lda* and *ldb* equal to 1.
 *
 * Parameters:
 *   m   - Number of rows to copy
 *   n   - Number of columns to copy
 *   A   - lda-by-n Matrix A
 *   lda - Leading dimension of A (lda >= m)
 *   B   - Destination ldb-by-n matrix B
 *   ldb - Leading dimension of B (ldb >= m)
 */
#define SOPALIN_LACPY(m, n, A, lda, B, ldb)                     \
  {                                                             \
    pastix_int_t k;                                                      \
    for (k=0; k<(n); k++)                                       \
    {                                                           \
      SOPALIN_COPY(m, &(A[k*(lda)]), iun, &(B[k*(ldb)]), iun);  \
    }                                                           \
  }


/*
 * Macro: SOPALIN_SCAL
 *
 * Multiply a vector by a scalar.
 *
 * Care: if BLAS_USE_SCAL is defined, only works
 *   with *u* equal to 1.
 *
 * Parameters:
 *   n - Size of the vector
 *   a - Scalar to multiply the vector by.
 *   x - Vector.
 *   u - stride between to element of the vector.
 *
 */
#ifdef BLAS_USE_SCAL
#  define SOPALIN_SCAL(n,a,x,u) {               \
    pastix_int_t i; pastix_float_t *pt,*p=(x);                    \
    for((pt)=(p+n);(p)<(pt);)                   \
      *((p)++)*= (a);                           \
  }
#else /* BLAS_USE_SCAL */
#  define SOPALIN_SCAL(n,a,x,u) {                                   \
    BLAS_INT varin     = (BLAS_INT)(n);                             \
    BLAS_INT variu     = (BLAS_INT)(u);                             \
    pastix_float_t    float_tmp = a;                                         \
    BLAS_FLOAT varia   = *((BLAS_FLOAT*) ((void*)(&(float_tmp))));  \
    TAB_CHECK_NAN(((pastix_float_t*)x), n);                                  \
    SCAL(&varin, &varia, (BLAS_FLOAT*)(x), &variu);                 \
    TAB_CHECK_NAN(((pastix_float_t*)x), n);                                  \
  }
#endif

/*
 * Macro: SOPALIN_TRSM
 *
 * Solves $ A  \times  X = x  \times  B $ or $ X  \times  A = x  \times  B $ depending on
 * s value.
 * A is a triangular matrix.
 *
 * Parameters:
 *   s - Position of *a* 'L' or 'l' for left, 'R' or 'r' for right.
 *   p - Indicate if the matrice is upper ('U') or lower ('L')
 *       triangular matrix.
 *   t - Indicate if user wants to transpose *a*.
 *   d - Indicate if diagonal is united ('U') or not ('N').
 *   m - Number of lines in *b*.
 *   n - Number of columns in *b*.
 *   x - Coefficient to multiply right-hand-side by.
 *   a - *a* matrix (MxM or NxN depending on s).
 *   u - Stride for *a*
 *   b - *b* MxN matrix.
 *   v - Stride for *b*.
 */
#define SOPALIN_TRSM(s,p,t,d,m,n,x,a,u,b,v) {                     \
    BLAS_INT   varim     = (BLAS_INT)(m);                         \
    BLAS_INT   varin     = (BLAS_INT)(n);                         \
    BLAS_INT   variu     = (BLAS_INT)(u);                         \
    BLAS_INT   variv     = (BLAS_INT)(v);                         \
    pastix_float_t      float_tmp = x;                                     \
    BLAS_FLOAT varix     = *((BLAS_FLOAT*)((void*)&(float_tmp))); \
    TAB_CHECK_NAN(((pastix_float_t*)b), m*n);                              \
    if ( *s == 'L' || *s == 'l' ) {                               \
      TAB_CHECK_NAN(((pastix_float_t*)a), m*m);                            \
    }                                                             \
    else {                                                        \
      TAB_CHECK_NAN(((pastix_float_t*)a), n*n);                            \
    }                                                             \
    TRSM((s), (p), (t), (d), &varim, &varin,                      \
         &varix, (BLAS_FLOAT*)(a), &variu,                        \
         (BLAS_FLOAT*)(b), &variv);                               \
                                                                  \
    TAB_CHECK_NAN(((pastix_float_t*)b), m*n);                              \
  }

/*
 * Macro: SOPALIN_TRSV
 *
 * Solves a system of linear equation with
 * a triangular matrix.
 *
 * Parameters:
 *   p - Indicate if the matrix is upper ('U') or lower ('L') triangular.
 *   t - Indicate if user wants to transpose ('T') the matrix or not ('N').
 *   d - Indicate if the matrix has a united diagonal ('U') ot not ('N').
 *   n - Size of the matrix and vector.
 *   a - Triangular matrix NxN.
 *   u - Stride for *a*.
 *   b - Right-hand-side member.
 *   v - Stride for *b*.
 */
#define SOPALIN_TRSV(p,t,d,n,a,u,b,v) {                     \
    BLAS_INT varin = (BLAS_INT)(n);                         \
    BLAS_INT variu = (BLAS_INT)(u);                         \
    BLAS_INT variv = (BLAS_INT)(v);                         \
    TAB_CHECK_NAN(((pastix_float_t*)a), n*n);                        \
    TAB_CHECK_NAN(((pastix_float_t*)b), n);                          \
                                                            \
    TRSV((p), (t), (d), &varin,                             \
         (BLAS_FLOAT*)(a), &variu,                          \
         (BLAS_FLOAT*)(b), &variv);                         \
                                                            \
    TAB_CHECK_NAN(((pastix_float_t*)b), n);                          \
  }

/*
 * Macro: SOPALIN_GEMM
 *
 * Computes $ c = x \times a \times b + y \times c $.
 *
 * Parameters:
 *   i - indicate if user wants to transpose *a*.
 *   j - indicate if user wants to transpose *b*.
 *   m - Number of lines in *c* and *a*.
 *   n - Number of columns in *c* and *b*.
 *   k - Number of lines in *b* and columns in *a*.
 *   x - coefficient for $ a \times b$
 *   a - $m \times k$ Matrix.
 *   u - Stride for *a*.
 *   b - $k \times n$ Matrix.
 *   v - Stride for *b*.
 *   y - Coefficient for *c*
 *   c - $m \times n$ matrix.
 *   w - Stride for *c*
 */
#define SOPALIN_GEMM(i,j,m,n,k,x,a,u,b,v,y,c,w) {                     \
    BLAS_INT   varim = (BLAS_INT)(m);                                 \
    BLAS_INT   varin = (BLAS_INT)(n);                                 \
    BLAS_INT   varik = (BLAS_INT)(k);                                 \
    BLAS_INT   variu = (BLAS_INT)(u);                                 \
    BLAS_INT   variv = (BLAS_INT)(v);                                 \
    BLAS_INT   variw = (BLAS_INT)(w);                                 \
    pastix_float_t      float_tmp;                                             \
    BLAS_FLOAT varix;                                                 \
    BLAS_FLOAT variy;                                                 \
    float_tmp = x;                                                    \
    varix     = *((BLAS_FLOAT*) ((void*)&(float_tmp)));               \
    float_tmp = y;                                                    \
    variy     = *((BLAS_FLOAT*) ((void*)&(float_tmp)));               \
    TAB_CHECK_NAN(((pastix_float_t*)a), k*m);                                  \
    TAB_CHECK_NAN(((pastix_float_t*)b), k*n);                                  \
    TAB_CHECK_NAN(((pastix_float_t*)c), m*n);                                  \
                                                                      \
    GEMM((i), (j), &varim, &varin, &varik, &varix,                    \
         (BLAS_FLOAT*)(a), &variu,                                    \
         (BLAS_FLOAT*)(b), &variv, &variy,                            \
         (BLAS_FLOAT*)(c), &variw);                                   \
                                                                      \
    TAB_CHECK_NAN(((pastix_float_t*)c), m*n);                                  \
  }

#define SOPALIN_GEMV(i,m,n,x,a,u,b,v,y,c,w) {           \
    BLAS_INT   varim = (BLAS_INT)(m);                   \
    BLAS_INT   varin = (BLAS_INT)(n);                   \
    BLAS_INT   variu = (BLAS_INT)(u);                   \
    BLAS_INT   variv = (BLAS_INT)(v);                   \
    BLAS_INT   variw = (BLAS_INT)(w);                   \
    pastix_float_t      float_tmp;                               \
    BLAS_FLOAT varix;                                   \
    BLAS_FLOAT variy;                                   \
    float_tmp = x;                                      \
    varix     = *((BLAS_FLOAT*) ((void*)&(float_tmp))); \
    float_tmp = y;                                      \
    variy     = *((BLAS_FLOAT*) ((void*)&(float_tmp))); \
                                                        \
    TAB_CHECK_NAN(((pastix_float_t*)a), m*n);                    \
    if (*i == 'N' || *i == 'n')   {                     \
      TAB_CHECK_NAN(((pastix_float_t*)b), n);                    \
      TAB_CHECK_NAN(((pastix_float_t*)c), m);                    \
    }                                                   \
    else {                                              \
      TAB_CHECK_NAN(((pastix_float_t*)b), m);                    \
      TAB_CHECK_NAN(((pastix_float_t*)c), n);                    \
    }                                                   \
    GEMV((i), &varim, &varin, &varix,                   \
         (BLAS_FLOAT*)(a), &variu,                      \
         (BLAS_FLOAT*)(b), &variv,                      \
         &variy, (BLAS_FLOAT*)(c), &variw);             \
    if (*i == 'N' || *i == 'n')   {                     \
      TAB_CHECK_NAN(((pastix_float_t*)c), m);                    \
    }                                                   \
    else {                                              \
      TAB_CHECK_NAN(((pastix_float_t*)c), n);                    \
    }                                                   \
  }


#  define SOPALIN_AXPY(n,a,x,ix,y,iy) {                            \
    BLAS_INT   varin     = (BLAS_INT)(n);                          \
    BLAS_INT   varix     = (BLAS_INT)(ix);                         \
    BLAS_INT   variy     = (BLAS_INT)(iy);                         \
    pastix_float_t      float_tmp = a;                                      \
    BLAS_FLOAT varia     = *((BLAS_FLOAT*) ((void*)&(float_tmp))); \
    TAB_CHECK_NAN(((pastix_float_t*)x), n);                                 \
    TAB_CHECK_NAN(((pastix_float_t*)y), n);                                 \
    AXPY(&varin, &varia,                                           \
         (BLAS_FLOAT*)(x), &varix,                                 \
         (BLAS_FLOAT*)(y), &variy);                                \
    TAB_CHECK_NAN(((pastix_float_t*)y), n);                                 \
  }

#ifdef CPLX
#  define SOPALIN_SYR(i,n,x,a,u,b,v) {                                  \
    BLAS_INT   varin     = (BLAS_INT)(n);                               \
    BLAS_INT   variu     = (BLAS_INT)(u);                               \
    BLAS_INT   variv     = (BLAS_INT)(v);                               \
    pastix_float_t      float_tmp = x;                                           \
    BLAS_FLOAT varix     = *((BLAS_FLOAT*) ((void*)&(float_tmp)));      \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*n);                                    \
    TAB_CHECK_NAN(((pastix_float_t*)a), n);                                      \
    GER(&varin, &varin, &varix,                                         \
        (BLAS_FLOAT*)(a), &variu,                                       \
        (BLAS_FLOAT*)(a), &variu,                                       \
        (BLAS_FLOAT*)(b), &variv);                                      \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*n);                                    \
  }
#  define SOPALIN_HER(i,n,x,a,u,b,v) {                               \
    BLAS_INT   varin = (BLAS_INT)(n);                                \
    BLAS_INT   variu = (BLAS_INT)(u);                                \
    BLAS_INT   variv = (BLAS_INT)(v);                                \
    BLAS_REAL  varix = (BLAS_REAL)(x);                               \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*n);                                 \
    TAB_CHECK_NAN(((pastix_float_t*)a), n);                                   \
    SYR((i), &varin, &varix,                                         \
        (BLAS_FLOAT*)(a), &variu,                                    \
        (BLAS_FLOAT*)(b), &variv);                                   \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*n);                                 \
  }
#else  /* CPLX */
#  define SOPALIN_SYR(i,n,x,a,u,b,v) {                            \
    BLAS_INT varin     = (BLAS_INT)(n);                           \
    BLAS_INT variu     = (BLAS_INT)(u);                           \
    BLAS_INT variv     = (BLAS_INT)(v);                           \
    pastix_float_t    float_tmp = x;                                       \
    BLAS_FLOAT varix   = *((BLAS_FLOAT*) ((void*)&(float_tmp)));  \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*n);                              \
    TAB_CHECK_NAN(((pastix_float_t*)a), n);                                \
    SYR((i), &varin, &varix,                                      \
        (BLAS_FLOAT*)(a), &variu,                                 \
        (BLAS_FLOAT*)(b), &variv);                                \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*n);                              \
  }
#  define SOPALIN_HER(i,n,x,a,u,b,v) {                             \
    BLAS_INT   varin     = (BLAS_INT)(n);                          \
    BLAS_INT   variu     = (BLAS_INT)(u);                          \
    BLAS_INT   variv     = (BLAS_INT)(v);                          \
    pastix_float_t      float_tmp = x;                                      \
    BLAS_FLOAT varix     = *((BLAS_FLOAT*) ((void*)&(float_tmp))); \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*n);                               \
    TAB_CHECK_NAN(((pastix_float_t*)a), n);                                 \
    SYR((i), &varin, &varix,                                       \
        (BLAS_FLOAT*)(a), &variu,                                  \
        (BLAS_FLOAT*)(b), &variv);                                 \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*n);                               \
  }
#endif
#define SOPALIN_SYRK(i,t,n,k,x,a,u,y,b,v) {            \
    BLAS_INT   varin = (BLAS_INT)(n);                  \
    BLAS_INT   varik = (BLAS_INT)(k);                  \
    BLAS_INT   variu = (BLAS_INT)(u);                  \
    BLAS_INT   variv = (BLAS_INT)(v);                  \
    BLAS_REAL  varix = (BLAS_REAL) (x);                \
    BLAS_REAL  variy = (BLAS_REAL) (y);                \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*n);                   \
    if (*t == 'N' || *t == 'N') {                      \
      TAB_CHECK_NAN(((pastix_float_t*)a), k*n);                 \
    }                                                  \
    else {                                             \
      TAB_CHECK_NAN(((pastix_float_t*)a), n*k);                 \
    }                                                  \
    SYRK((i), (t), &varin, &varik, &varix,             \
         (BLAS_FLOAT*)(a), &variu,                     \
         &variy, (BLAS_FLOAT*)(b), &variv);            \
    TAB_CHECK_NAN(((pastix_float_t*)b), n*n);                   \
  }

#define SOPALIN_GER(m,n,x,a,u,b,v,c,w) {                                \
    BLAS_INT varim     = (BLAS_INT)(m);                                 \
    BLAS_INT varin     = (BLAS_INT)(n);                                 \
    BLAS_INT variu     = (BLAS_INT)(u);                                 \
    BLAS_INT variv     = (BLAS_INT)(v);                                 \
    BLAS_INT variw     = (BLAS_INT)(w);                                 \
    pastix_float_t    float_tmp = x;                                             \
    BLAS_FLOAT varix   = *((BLAS_FLOAT*) ((void*)&(float_tmp)));        \
    TAB_CHECK_NAN(((pastix_float_t*)c), n*m);                                    \
    TAB_CHECK_NAN(((pastix_float_t*)a), n);                                      \
    TAB_CHECK_NAN(((pastix_float_t*)b), n);                                      \
    GER(&varim, &varin, &varix, (BLAS_FLOAT*)(a), &variu,               \
        (BLAS_FLOAT*)(b), &variv,                                       \
        (BLAS_FLOAT*)(c), &variw);                                      \
    TAB_CHECK_NAN(((pastix_float_t*)c), n*m);                                    \
  }

#define SOPALIN_POF(i,a,u,n)  POF((i), (a), (u), (n))
#define SOPALIN_PPF(a,n,f)    PPF((a), (n), (f))
#endif /* SOPALIN_COMPUTE_H */
