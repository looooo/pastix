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
#define CS_VER 1    /* CSparse Version 1.2.0 */
#define CS_SUBVER 2
#define CS_SUBSUB 0
#define CS_DATE "Mar 6, 2006"    /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006"

typedef struct z_cs_sparse    /* matrix in compressed-column or triplet form */
{
  pastix_int_t nzmax ;    /* maximum number of entries */
  pastix_int_t m ;    /* number of rows */
  pastix_int_t n ;    /* number of columns */
  pastix_int_t *p ;    /* column poINTers (size n+1) or col indices (size nzmax) */
  pastix_int_t *i ;    /* row indices, size nzmax */
  pastix_complex64_t *x ;    /* numerical values, size nzmax */
  pastix_int_t nz ;    /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

/* keep all large entries */
pastix_int_t z_cs_droptol (cs *A, pastix_complex64_t tol) ;
/* keep all nonzero entries */
pastix_int_t z_cs_dropzeros (cs *A) ;
/* C = alpha*A + beta*B */
cs *z_cs_add (const cs *A, const cs *B, pastix_complex64_t alpha, pastix_complex64_t beta) ;
/* removes duplicate entries from A */
pastix_int_t z_cs_dupl (cs *A) ;
/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
pastix_int_t z_cs_entry (cs *T, pastix_int_t i, pastix_int_t j, pastix_complex64_t x) ;
/* drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
pastix_int_t z_cs_fkeep (cs *A, pastix_int_t (*fkeep) (pastix_int_t, pastix_int_t, pastix_complex64_t, void *), void *other) ;
/* y = A*x+y */
pastix_int_t z_cs_gaxpy (const cs *A, const pastix_complex64_t *x, pastix_complex64_t *y) ;
/* C = A*B */
cs *z_cs_multiply (const cs *A, const cs *B) ;
/* 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum */
pastix_complex64_t z_cs_norm (const cs *A) ;
/* C = A(P,Q) where P and Q are permutations of 0..m-1 and 0..n-1. */
cs *z_cs_permute (const cs *A, const pastix_int_t *P, const pastix_int_t *Q, pastix_int_t values) ;
/* Pinv = P', or P = Pinv' */
pastix_int_t *z_cs_pinv (const pastix_int_t *P, pastix_int_t n) ;
/* C = A' */
cs *z_cs_transpose (const cs *A, pastix_int_t values) ;
/* C = compressed-column form of a triplet matrix T */
cs *z_cs_triplet (const cs *T) ;
/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
pastix_int_t z_cs_scatter (const cs *A, pastix_int_t j, pastix_complex64_t beta, pastix_int_t *w, pastix_complex64_t *x, pastix_int_t mark,
                cs *C, pastix_int_t nz) ;
/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
pastix_int_t z_cs_cumsum (pastix_int_t *p, pastix_int_t *c, pastix_int_t n) ;

/* utilities */
/* wrapper for malloc */
void *z_cs_malloc (pastix_int_t n, size_t size) ;
/* wrapper for calloc */
void *z_cs_calloc (pastix_int_t n, size_t size) ;
/* wrapper for free */
void *z_cs_free (void *p) ;
/* wrapper for realloc */
void *z_cs_realloc (void *p, pastix_int_t n, size_t size, pastix_int_t *ok) ;
/* allocate a sparse matrix (triplet form or compressed-column form) */
cs *z_cs_spalloc (pastix_int_t m, pastix_int_t n, pastix_int_t nzmax, pastix_int_t values, pastix_int_t triplet) ;
/* change the max # of entries sparse matrix */
pastix_int_t z_cs_sprealloc (cs *A, pastix_int_t nzmax) ;
/* free a sparse matrix */
cs *z_cs_spfree (cs *A) ;
/* free workspace and return a sparse matrix result */
cs *z_cs_done (cs *C, void *w, void *x, pastix_int_t ok) ;

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(Ap,j) (Ap [j] < 0)
#define CS_MARK(Ap,j) { Ap [j] = CS_FLIP (Ap [j]) ; }
#define CS_OVERFLOW(n,size) (n > INT_MAX / (pastix_int_t) size)
