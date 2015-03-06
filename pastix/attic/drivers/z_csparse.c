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

#include "z_pastix.h"

static pastix_int_t z_cs_tol (pastix_int_t i, pastix_int_t j, pastix_complex64_t aij, void *tol)
{
  return (fabs (aij) > *((pastix_complex64_t *) tol)) ;
}

pastix_int_t z_cs_droptol (cs *A, pastix_complex64_t tol)
{
  return (z_cs_fkeep (A, &z_cs_tol, &tol)) ;    /* keep all large entries */
}

static pastix_int_t z_cs_nonzero (pastix_int_t i, pastix_int_t j, pastix_complex64_t aij, void *other)
{
  return (aij != 0) ;
}

pastix_int_t z_cs_dropzeros (cs *A)
{
  return (z_cs_fkeep (A, &z_cs_nonzero, NULL)) ;/* keep all nonzero entries */
}

/* C = alpha*A + beta*B */
cs *z_cs_add ( const cs *A, const cs *B, pastix_complex64_t alpha, pastix_complex64_t beta )
{
  pastix_int_t p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values ;
  pastix_complex64_t *x, *Bx, *Cx ;
  cs *C ;
  if (!A || !B) return (NULL) ;/* check inputs */
  m = A->m ; anz = A->p [A->n] ;
  n = B->n ; Bp = B->p ; Bx = B->x ; bnz = Bp [n] ;
  w = z_cs_calloc (m, sizeof (pastix_int_t)) ;
  values = (A->x != NULL) && (Bx != NULL) ;
  x = values ? z_cs_malloc (m, sizeof (pastix_complex64_t)) : NULL ;
  C = z_cs_spalloc (m, n, anz + bnz, values, 0) ;
  if (!C || !w || (values && !x)) return (z_cs_done (C, w, x, 0)) ;
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (j = 0 ; j < n ; j++)
    {
      Cp [j] = nz ;/* column j of C starts here */
      nz = z_cs_scatter (A, j, alpha, w, x, j+1, C, nz) ;   /* alpha*A(:,j)*/
      nz = z_cs_scatter (B, j, beta, w, x, j+1, C, nz) ;    /* beta*B(:,j) */
      if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
  Cp [n] = nz ;/* finalize the last column of C */
  z_cs_sprealloc (C, 0) ;/* remove extra space from C */
  return (z_cs_done (C, w, x, 1)) ;/* success; free workspace, return C */
}

/* removes duplicate entries from A */
pastix_int_t z_cs_dupl (cs *A)
{
  pastix_int_t i, j, p, q, nz = 0, n, m, *Ap, *Ai, *w ;
  pastix_complex64_t *Ax ;
  if (!A) return (0) ;/* check inputs */
  m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  w = z_cs_malloc (m, sizeof (pastix_int_t)) ;/* get workspace */
  if (!w) return (0) ;/* out of memory */
  for (i = 0 ; i < m ; i++) w [i] = -1 ;/* row i not yet seen */
  for (j = 0 ; j < n ; j++)
    {
      q = nz ;/* column j will start at q */
      for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	  i = Ai [p] ;/* A(i,j) is nonzero */
	  if (w [i] >= q)
	    {
	      Ax [w [i]] += Ax [p] ;/* A(i,j) is a duplicate */
	    }
	  else
	    {
	      w [i] = nz ;/* record where row i occurs */
	      Ai [nz] = i ;/* keep A(i,j) */
	      Ax [nz++] = Ax [p] ;
	    }
	}
      Ap [j] = q ;/* record start of column j */
    }
  Ap [n] = nz ;/* finalize A */
  z_cs_free (w) ;/* free workspace */
  return (z_cs_sprealloc (A, 0)) ;/* remove extra space from A */
}

/* add an entry to a triplet matrix; return 1 if ok, 0 otherwise */
pastix_int_t z_cs_entry (cs *T, pastix_int_t i, pastix_int_t j, pastix_complex64_t x)
{
  if (!T || (T->nz >= T->nzmax && !z_cs_sprealloc (T, 2*(T->nzmax)))) return(0);
  if (T->x) T->x [T->nz] = x ;
  T->i [T->nz] = i ;
  T->p [T->nz++] = j ;
  T->m = CS_MAX (T->m, i+1) ;
  T->n = CS_MAX (T->n, j+1) ;
  return (1) ;
}

/* drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 */
pastix_int_t z_cs_fkeep (cs *A, pastix_int_t (*fkeep) (pastix_int_t, pastix_int_t, pastix_complex64_t, void *), void *other)
{
  pastix_int_t baseval = 1 ;
  pastix_int_t j, p, nz = 0, n, *Ap, *Ai ;
  pastix_complex64_t *Ax ;
  if (!A || !fkeep) return (-1) ;    /* check inputs */
  n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  for (j = 0 ; j < n ; j++)
    {
      p = Ap [j] - baseval ;    /* get current location of col j */
      Ap [j] = nz + baseval ;    /* record new location of col j */
      for ( ; p < Ap [j+1] - baseval ; p++)
	{
	  if (fkeep (Ai [p], j, Ax ? Ax [p] : 1, other))
	    {
	      if (Ax) Ax [nz] = Ax [p] ;  /* keep A(i,j) */
	      Ai [nz++] = Ai [p] ;
	    }
	  else printf("drop %ld,%ld\n",(long)j,(long)p);
	}
    }
  /* finalize A and return nnz(A) */
  return (Ap [n] = nz + baseval) ;
}

/* y = A*x+y */
pastix_int_t z_cs_gaxpy (const cs *A, const pastix_complex64_t *x, pastix_complex64_t *y)
{
  pastix_int_t p, j, n, *Ap, *Ai ;
  pastix_complex64_t *Ax ;
  if (!A || !x || !y) return (0) ;    /* check inputs */
  n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  for (j = 0 ; j < n ; j++)
    {
      for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	  y [Ai [p]] += Ax [p] * x [j] ;
	}
    }
  return (1) ;
}


/* C = A*B */
cs *z_cs_multiply (const cs *A, const cs *B)
{
  pastix_int_t p, j, nz = 0, anz, *Cp, *Ci, *Bp, m, n, bnz, *w, values, *Bi ;
  pastix_complex64_t *x, *Bx, *Cx ;
  cs *C ;
  if (!A || !B) return (NULL) ;/* check inputs */
  m = A->m ; anz = A->p [A->n] ;
  n = B->n ; Bp = B->p ; Bi = B->i ; Bx = B->x ; bnz = Bp [n] ;
  w = z_cs_calloc (m, sizeof (pastix_int_t)) ;
  values = (A->x != NULL) && (Bx != NULL) ;
  x = values ? z_cs_malloc (m, sizeof (pastix_complex64_t)) : NULL ;
  C = z_cs_spalloc (m, n, anz + bnz, values, 0) ;
  if (!C || !w || (values && !x)) return (z_cs_done (C, w, x, 0)) ;
  Cp = C->p ;
  for (j = 0 ; j < n ; j++)
    {
      if (nz + m > C->nzmax && !z_cs_sprealloc (C, 2*(C->nzmax)+m))
	{
	  return (z_cs_done (C, w, x, 0)) ;/* out of memory */
	} 
      Ci = C->i ; Cx = C->x ;/* C may have been reallocated */
      Cp [j] = nz ;/* column j of C starts here */
      for (p = Bp [j] ; p < Bp [j+1] ; p++)
	{
	  nz = z_cs_scatter (A, Bi [p], Bx ? Bx [p] : 1, w, x, j+1, C, nz) ;
	}
      if (values) for (p = Cp [j] ; p < nz ; p++) Cx [p] = x [Ci [p]] ;
    }
  Cp [n] = nz ;/* finalize the last column of C */
  z_cs_sprealloc (C, 0) ;/* remove extra space from C */
  return (z_cs_done (C, w, x, 1)) ;/* success; free workspace, return C */
}

/* 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum */
pastix_complex64_t z_cs_norm (const cs *A)
{
  pastix_int_t p, j, n, *Ap ;
  pastix_complex64_t *Ax,  norm = 0, s ;
  if (!A || !A->x) return (-1) ;/* check inputs */
  n = A->n ; Ap = A->p ; Ax = A->x ;
  for (j = 0 ; j < n ; j++)
    {
      for (s = 0, p = Ap [j] ; p < Ap [j+1] ; p++) s += fabs (Ax [p]) ;
      norm = CS_MAX (norm, s) ;
    }
  return (norm) ;
}

/* C = A(P,Q) where P and Q are permutations of 0..m-1 and 0..n-1. */
cs *z_cs_permute (const cs *A, const pastix_int_t *Pinv, const pastix_int_t *Q, pastix_int_t values)
{
  pastix_int_t p, j, k, nz = 0, m, n, *Ap, *Ai, *Cp, *Ci ;
  pastix_complex64_t *Cx, *Ax ;
  cs *C ;
  if (!A) return (NULL) ;/* check inputs */
  m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  C = z_cs_spalloc (m, n, Ap [n], values && Ax != NULL, 0) ;
  if (!C) return (z_cs_done (C, NULL, NULL, 0)) ;   /* out of memory */
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (k = 0 ; k < n ; k++)
    {
      Cp [k] = nz ;/* column k of C is column Q[k] of A */
      j = Q ? (Q [k]) : k ;
      for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	  if (Cx) Cx [nz] = Ax [p] ;/* row i of A is row Pinv[i] of C */
	  Ci [nz++] = Pinv ? (Pinv [Ai [p]]) : Ai [p] ;
	}
    }
  Cp [n] = nz ;/* finalize the last column of C */
  return (z_cs_done (C, NULL, NULL, 1)) ;
}

/* Pinv = P', or P = Pinv' */
pastix_int_t *z_cs_pinv (pastix_int_t const *P, pastix_int_t n)
{
  pastix_int_t k, *Pinv ;
  if (!P) return (NULL) ;/* P = NULL denotes identity */
  Pinv = z_cs_malloc (n, sizeof (pastix_int_t)) ;/* allocate resuult */
  if (!Pinv) return (NULL) ;/* out of memory */
  for (k = 0 ; k < n ; k++) Pinv [P [k]] = k ;/* invert the permutation */
  return (Pinv) ;/* return result */
}

/* C = A' */
cs *z_cs_transpose (const cs *A, pastix_int_t values)
{
  pastix_int_t p, q, j, *Cp, *Ci, n, m, *Ap, *Ai, *w ;
  pastix_complex64_t *Cx, *Ax ;
  cs *C ;
  if (!A) return (NULL) ;
  m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
  C = z_cs_spalloc (n, m, Ap [n], values && Ax, 0) ;   /* allocate result */
  w = z_cs_calloc (m, sizeof (pastix_int_t)) ;
  if (!C || !w) return (z_cs_done (C, w, NULL, 0)) ;   /* out of memory */
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (p = 0 ; p < Ap [n] ; p++) w [Ai [p]]++ ;   /* row counts */
  z_cs_cumsum (Cp, w, m) ;   /* row poINTers */
  for (j = 0 ; j < n ; j++)
    {
      for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	  Ci [q = w [Ai [p]]++] = j ;/* place A(i,j) as entry C(j,i) */
	  if (Cx) Cx [q] = Ax [p] ;
	}
    }
  return (z_cs_done (C, w, NULL, 1)) ;/* success; free w and return C */
}

/* C = compressed-column form of a triplet matrix T */
cs *z_cs_triplet (const cs *T)
{
  pastix_int_t m, n, nz, p, k, *Cp, *Ci, *w, *Ti, *Tj ;
  pastix_complex64_t *Cx, *Tx ;
  cs *C ;
  if (!T) return (NULL) ;/* check inputs */
  m = T->m ; n = T->n ; Ti = T->i ; Tj = T->p ; Tx = T->x ; nz = T->nz ;
  C = z_cs_spalloc (m, n, nz, Tx != NULL, 0) ;/* allocate result */
  w = z_cs_calloc (n, sizeof (pastix_int_t)) ;/* get workspace */
  if (!C || !w) return (z_cs_done (C, w, NULL, 0)) ;/* out of memory */
  Cp = C->p ; Ci = C->i ; Cx = C->x ;
  for (k = 0 ; k < nz ; k++) w [Tj [k]]++ ;/* column counts */
  z_cs_cumsum (Cp, w, n) ;/* column poINTers */
  for (k = 0 ; k < nz ; k++)
    {
      Ci [p = w [Tj [k]]++] = Ti [k] ;    /* A(i,j) is the pth entry in C */
      if (Cx) Cx [p] = Tx [k] ;
    }
  return (z_cs_done (C, w, NULL, 1)) ;    /* success; free w and return C */
}

/* x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse */
pastix_int_t z_cs_scatter (const cs *A, pastix_int_t j, pastix_complex64_t beta, pastix_int_t *w, pastix_complex64_t *x, pastix_int_t mark,
		cs *C, pastix_int_t nz)
{
  pastix_int_t i, p, *Ap, *Ai, *Ci ;
  pastix_complex64_t *Ax ;
  if (!A || !w || !C) return (-1) ;/* ensure inputs are valid */
  Ap = A->p ; Ai = A->i ; Ax = A->x ; Ci = C->i ;
  for (p = Ap [j] ; p < Ap [j+1] ; p++)
    {
      i = Ai [p] ;/* A(i,j) is nonzero */
      if (w [i] < mark)
	{
	  w [i] = mark ;/* i is new entry in column j */
	  Ci [nz++] = i ;/* add i to pattern of C(:,j) */
	  if (x) x [i] = beta * Ax [p] ;/* x(i) = beta*A(i,j) */
	}
      else if (x) x [i] += beta * Ax [p] ;/* i exists in C(:,j) already */
    }
  return (nz) ;
}

/* p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c */
pastix_int_t z_cs_cumsum (pastix_int_t *p, pastix_int_t *c, pastix_int_t n)
{
  pastix_int_t i, nz = 0 ;
  if (!p || !c) return (-1) ;    /* check inputs */
  for (i = 0 ; i < n ; i++)
    {
      p [i] = nz ;
      nz += c [i] ;
      c [i] = p [i] ;
    }
  p [n] = nz ;
  return (nz) ;    /* return sum (c [0..n-1]) */
}

/* wrapper for malloc */
void *z_cs_malloc (pastix_int_t n, size_t size)
{
  return (CS_OVERFLOW (n,size) ? NULL : malloc (CS_MAX (n,1) * size)) ;
}

/* wrapper for calloc */
void *z_cs_calloc (pastix_int_t n, size_t size)
{
  return (CS_OVERFLOW (n,size) ? NULL : calloc (CS_MAX (n,1), size)) ;
}

/* wrapper for free */
void *z_cs_free (void *p)
{
  if (p) free (p) ;    /* free p if it is not already NULL */
  return (NULL) ;    /* return NULL to simplify the use of z_cs_free */
}

/* wrapper for realloc */
void *z_cs_realloc (void *p, pastix_int_t n, size_t size, pastix_int_t *ok)
{
  void *p2 ;
  *ok = !CS_OVERFLOW (n,size) ;    /* guard against pastix_int_t overflow */
  if (!(*ok)) return (p) ;    /* p unchanged if n too large */
  p2 = realloc (p, CS_MAX (n,1) * size) ; /* realloc the block */
  *ok = (p2 != NULL) ;
  return ((*ok) ? p2 : p) ;    /* return original p if failure */
}

/* allocate a sparse matrix (triplet form or compressed-column form) */
cs *z_cs_spalloc (pastix_int_t m, pastix_int_t n, pastix_int_t nzmax, pastix_int_t values, pastix_int_t triplet)
{
  cs *A = z_cs_calloc (1, sizeof (cs)) ;    /* allocate the cs struct */
  if (!A) return (NULL) ;    /* out of memory */
  A->m = m ;    /* define dimensions and nzmax */
  A->n = n ;
  A->nzmax = nzmax = CS_MAX (nzmax, 1) ;
  A->nz = triplet ? 0 : -1 ;    /* allocate triplet or comp.col */
  A->p = z_cs_malloc (triplet ? nzmax : n+1, sizeof (pastix_int_t)) ;
  A->i = z_cs_malloc (nzmax, sizeof (pastix_int_t)) ;
  A->x = values ? z_cs_malloc (nzmax, sizeof (pastix_complex64_t)) : NULL ;
  return ((!A->p || !A->i || (values && !A->x)) ? z_cs_spfree (A) : A) ;
}

/* change the max # of entries sparse matrix */
pastix_int_t z_cs_sprealloc (cs *A, pastix_int_t nzmax)
{
  pastix_int_t ok, oki, okj = 1, okx = 1 ;
  if (!A) return (0) ;
  nzmax = (nzmax <= 0) ? (A->p [A->n]) : nzmax ;
  A->i = z_cs_realloc (A->i, nzmax, sizeof (pastix_int_t), &oki) ;
  if (A->nz >= 0) A->p = z_cs_realloc (A->p, nzmax, sizeof (pastix_int_t), &okj) ;
  if (A->x) A->x = z_cs_realloc (A->x, nzmax, sizeof (pastix_complex64_t), &okx) ;
  ok = (oki && okj && okx) ;
  if (ok) A->nzmax = nzmax ;
  return (ok) ;
}

/* free a sparse matrix */
cs *z_cs_spfree (cs *A)
{
  if (!A) return (NULL) ;/* do nothing if A already NULL */
  z_cs_free (A->p) ;
  z_cs_free (A->i) ;
  z_cs_free (A->x) ;
  return (z_cs_free (A)) ;/* free the cs struct and return NULL */
}

/* free workspace and return a sparse matrix result */
cs *z_cs_done (cs *C, void *w, void *x, pastix_int_t ok)
{
  z_cs_free (w) ;/* free workspace */
  z_cs_free (x) ;
  return (ok ? C : z_cs_spfree (C)) ;/* return result if OK, else free it */
}
