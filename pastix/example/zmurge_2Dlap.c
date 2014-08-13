/**
 *  File: zmurge_2Dlap.c
 *
 *  Example that generate A 2D Laplacian with multiple
 *  degrees of freedom and solves it.
 *
 *  That example only works with the -DDISTRIBUTED option.
 *
 *  This exameple is based on PetSC example src/ksp/ksp/examples/tutorials/ex3.c
 *
 *  Usage:
 *  > ./murge <size> <DofNbr> <mode>
 *
 *  Authors:
 *    Xavier LACOSTE - lacoste@labri.fr
 *
 *  @version 1.0.0
 *  @author Mathieu Faverge
 *  @author Pierre Ramet
 *  @author Xavier Lacoste
 *  @date 2011-11-11
 *  @precisions normal z -> c d s
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#ifdef TYPE_COMPLEX
#include <complex.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef FORCE_NOMPI
#  include <mpi.h>
#endif

#include "z_pastix.h"
#include "zmurge.h"
#include "zmurge_pastix.h"
#include "mem_trace.h"
extern char *optarg;

#define MEMORY_WRITE(mem) ( ((mem) < 1<<10) ?                          \
                            ( (double)(mem) ) :                         \
                            ( ( (mem) < 1<<20 ) ?                       \
                              ( (double)(mem)/(double)(1<<10) ) :       \
                              ( ((mem) < 1<<30 ) ?                      \
                                ( (double)(mem)/(double)(1<<20) ) :     \
                                ( (double)(mem)/(double)(1<<30) ))))
#define MEMORY_UNIT_WRITE(mem) (((mem) < 1<<10) ?                       \
                                "o" :                                   \
                                ( ( (mem) < 1<<20 ) ?                   \
                                  "Ko" :                                \
                                  ( ( (mem) < 1<<30 ) ?                 \
                                    "Mo" :                              \
                                    "Go" )))

#define MEMORY_PRINT(mem) MEMORY_WRITE(mem), MEMORY_UNIT_WRITE(mem)

/* enable some timing outputs */
#ifndef MURGE_TIME
#  define MURGE_TIME
#endif
#ifdef MURGE_TIME
#  define START_TIMER(t1) do {                                          \
    MPI_Barrier(MPI_COMM_WORLD);                                        \
    gettimeofday(&tv, NULL);                                            \
    t1 = ((double) tv.tv_sec + (double) tv.tv_usec * (double)1.0e-6L);  \
  } while(0)

#  define STOP_TIMER(str, t1) do {                                        \
    MPI_Barrier(MPI_COMM_WORLD);                                        \
    gettimeofday(&tv, NULL);                                            \
    t2 = ((double) tv.tv_sec + (double) tv.tv_usec * (double)1.0e-6L);  \
    if (rank == 0)                                                      \
      fprintf(stdout, " $$ time for '%-40s' %.2e s $$\n", str, t2-t1);  \
    mem = z_pastix_getMemoryUsage();                                    \
    maxmem = z_pastix_getMaxMemoryUsage();                              \
    if (rank == 0)                                                      \
      fprintf(stdout, " @@ memory usage after '%-40s'"                  \
              " %.3g %s, max %.3g %s @@\n", str,                        \
              MEMORY_PRINT(mem), MEMORY_PRINT(maxmem));                 \
  } while(0)
#else /* not MURGE_TIME */
#  define START_TIMER(t1)
#  define STOP_TIMER(str, t1)
#endif /* MURGE_TIME */

#define SUFFIX(name) name ## _2Dlap
#define MURGE_UserData_t SUFFIX(MURGE_UserData_t)
#define MURGE_UserData_  SUFFIX(MURGE_UserData_)
typedef struct MURGE_UserData_ {
  int m;
} MURGE_UserData_t;
#undef MURGE_user_data_t
#undef MURGE_user_data_
#define VERT_PER_ELEMENT(d) 4
#define GET_VERTICES(i, idx, d)                                       \
  do {                                                                \
    idx[0] = (d->m+1)*(i/d->m) + (i % d->m);                          \
    idx[1] = idx[0]+1;                                                \
    idx[2] = idx[1] + d->m + 1;                                       \
    idx[3] = idx[2] - 1;                                              \
  } while (0)
#include "ZMURGE_GetLocalElementList.c"
#define ZMURGE_GetLocalElementNbr  SUFFIX(ZMURGE_GetLocalElementNbr)
#define ZMURGE_GetLocalElementList SUFFIX(ZMURGE_GetLocalElementList)
#define MURGE_UserData_t          SUFFIX(MURGE_UserData_t)
#define MURGE_UserData_           SUFFIX(MURGE_UserData_)
#undef SUFIX
#undef VERT_PER_ELEMENT
#undef GET_VERTICES

#define CHKERRQ(ierr)                                      \
  do {                                                     \
    if (ierr != MURGE_SUCCESS)                             \
      {                                                    \
        fprintf(stderr, "%s:%d ERROR %ld in MURGE call\n", \
                __FILE__, __LINE__, (long)ierr);           \
        abort();                                           \
      }                                                    \
  } while(0)

#ifdef FORCE_NOMPI
#  define MPICHKERRQ(ierr)
#else
#  define MPICHKERRQ(ierr)                                 \
  do {                                                     \
    if (ierr != MPI_SUCCESS)                               \
      {                                                    \
        fprintf(stderr, "%s:%d ERROR %ld in MURGE call\n", \
                __FILE__, __LINE__, (long)ierr);           \
        abort();                                           \
      }                                                    \
  } while(0)
#endif


extern int trace_task;

/* element stiffness for Laplacian */
static inline
int FormElementStiffness(REAL H,COEF *Ke, int dof)
{
  int k;
  memset(Ke, 0, dof*dof*16*sizeof(COEF));
  for(k = 0; k < dof; k++)
    {
      int offset = k*dof+k;
      Ke += offset;
      Ke[0*dof*dof]  = H/6.0;
      Ke[4*dof*dof]  = -H/8.0;
      Ke[8*dof*dof]  = H/12.0;
      Ke[12*dof*dof] = -H/8;

      Ke[1*dof*dof]  = -H/8.0;
      Ke[5*dof*dof]  = H/6.0;
      Ke[9*dof*dof]  = -H/8.0;
      Ke[13*dof*dof]  = H/12.0;

      Ke[2*dof*dof]  = H/12.0;
      Ke[6*dof*dof]  = -H/8.0;
      Ke[10*dof*dof] = H/6.0;
      Ke[14*dof*dof] = -H/8.0;

      Ke[3*dof*dof]  = -H/8.0;
      Ke[7*dof*dof] = H/12.0;
      Ke[11*dof*dof] = -H/8.0;
      Ke[15*dof*dof] = H/6.0;
      Ke -= offset;
    }
  return 0;
}

/* --------------------------------------------------------------------- */
static inline
int FormElementRhs(REAL x,REAL y,REAL H,COEF *r, int dof)
{
   int k;
   for (k =0; k < dof; k++) {
     r[0*dof+k] = 0.;
     r[1*dof+k] = 0.;
     r[2*dof+k] = 0.;
     r[3*dof+k] = 0.;
   }
   return 0;
}

int main(int argc,char **argv)
{
  COEF    *u,*b,*ustar; /* approx solution, RHS, exact solution */
  INTS     N;           /* dimension of system (global) */
  INTS     M;           /* number of elements (global) */
  int      rank;        /* processor rank */
  int      size;        /* size of communicator */
  int      required, provided;
  COEF    *Ke;      /* element matrix */
  COEF     r[4];        /* element vector */
  REAL     h;           /* mesh width */
  REAL     x,y;
  COEF     val;
  INTS     ierr;
  INTS     idx[4],count,*rows,i, ie, m = 5,dof=1, start,end,j, k, l;
  INTS     id;
  INTL     nnz;
  INTS     localElementNbr;
  INTS    *localElements  = NULL;
  REAL     prec;
  INTS     nb_threads;
  INTS     iter_blocks = 0;
  INTS     n_blocks = 200;
  int      solv, n_solv = 1;
  int      fact, n_fact = 1;
  INTS    *indexes  = NULL;
  COEF    *values   = NULL;
  MURGE_UserData_t d;
  int      mode = 0;
  int      opt;
  char     out_string[40];
#ifdef MURGE_TIME
  struct timeval tv;
  double t0, t1, t2;
  unsigned long mem, maxmem;
#endif

#ifdef FORCE_NOMPI
  rank = 0;
  size = 1;
#  define MPI_COMM_WORLD 0
#else
  required=MPI_THREAD_MULTIPLE;
  MPI_Init_thread(&argc, &argv, required, &provided);
  START_TIMER(t0);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPICHKERRQ(ierr);
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&size);MPICHKERRQ(ierr);
#endif /* FORCE_NOMPI */
  memtrace_start();
  m = 100;
  while ((opt = getopt(argc, argv, "n:d:m:b:f:s:h")) != -1) {
    switch (opt) {
    case 'n':
      m = atoi(optarg);
      if (m < 1) goto usage;
      break;
    case 'd':
      dof = atoi(optarg);
      if (dof < 1) goto usage;
      break;
    case 'm':
      mode = atoi(optarg);
      if (mode > 3 || mode < 0) goto usage;
      break;
    case 'b':
      n_blocks = atoi(optarg);
      if (n_blocks < 1) goto usage;
      break;
    case 's':
      n_solv = atoi(optarg);
      if (n_solv < 1) goto usage;
      break;
    case 'f':
      n_fact = atoi(optarg);
      if (n_fact < 1) goto usage;
      break;
    case 'h':
    default:
    usage:
      printf("unprocessed option -%c %s\n", opt, optarg);
      printf("usage: %s -n <size> -d <DofNumber> -m <mode> -b <nblocks> -s <solve number> -f <factorization number> ]\n", argv[0]);
      return 1;
    }
  }

  if (mode == 3) {
    indexes = malloc(n_blocks*4*sizeof(INTS));
    values  = malloc(n_blocks*4*4*dof*dof*sizeof(COEF));
  }

  N = (m+1)*(m+1);
  M = m*m;
  h = 1.0/m;
  /* Starting MURGE*/
  ierr = ZMURGE_Initialize((INTS)1);
  if (ierr != MURGE_SUCCESS) {
    fprintf(stderr, "Error %ld in ZMURGE_Initialize\n", (long)ierr);
    return 1;
  }
  id = 0;

  /* Set Options */
  prec = 1e-7;

  ZMURGE_SetDefaultOptions(id, 0);
  ZMURGE_SetOptionINT(id, IPARM_VERBOSE, API_VERBOSE_NO);
  ZMURGE_SetOptionINT(id, IPARM_MATRIX_VERIFICATION, API_YES);

  nb_threads = 1;
#ifdef _OPENMP
#pragma omp parallel shared(nb_threads)
  {
    nb_threads = omp_get_num_threads();
  }
#endif /* _OPENMP */

  if (rank == 0) {
    fprintf(stdout, "Running on %ld threads and %d MPI Tasks\n",
            (long)nb_threads, size);
  }
  ZMURGE_SetOptionINT(id, IPARM_THREAD_NBR, nb_threads);
  ZMURGE_SetOptionINT(id, MURGE_IPARAM_DOF, dof);
  ZMURGE_SetOptionINT(id, MURGE_IPARAM_SYM, MURGE_BOOLEAN_FALSE);
  ZMURGE_SetOptionINT(id, MURGE_IPARAM_BASEVAL, 0);

  ZMURGE_SetOptionREAL(id, MURGE_RPARAM_EPSILON_ERROR, prec);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Au = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create stiffness matrix
  */
  ierr = ZMURGE_SetCommunicator(id, MPI_COMM_WORLD); CHKERRQ(ierr);


  /*
     Mark rows for for Dirichlet boundary conditions
  */
  rows = malloc(N*sizeof(INTS));
  memset(rows, 0, N*sizeof(INTS));

  for (i=0; i<m+1; i++) {
    rows[i] = 1; /* bottom */
    rows[m*(m+1) + i] = 1; /* top */
  }
  for (i=m+1; i<m*(m+1); i+= m+1) {
    rows[i] = 1;
  }
  for (i=2*m+1; i<m*(m+1); i+= m+1) {
    rows[i] = 1;
  }

  ZMURGE_SetDropRows(id, N, rows);
  fprintf(stdout, "mode %d\n", mode);
  if (mode == 2) {
    d.m = m;
    trace_task++;
    START_TIMER(t1);
    ierr = ZMURGE_GetLocalElementNbr(id,
                                    N,
                                    M,
                                    &localElementNbr,
                                    MURGE_DISTRIBUTE_ELEMENTS,
                                    &d); CHKERRQ(ierr);
    localElements = malloc(localElementNbr*sizeof(INTS));
    ierr = ZMURGE_GetLocalElementList(id, localElements); CHKERRQ(ierr);
    STOP_TIMER("getting element list", t1);
    trace_task++;

    start = 0;
    end   = localElementNbr;
  }
  else {
    start = rank*(M/size) + ((M%size) < rank ? (M%size) : rank);
    end   = start + M/size + ((M%size) > rank);
  }
  for (fact = 0; fact < n_fact; fact ++) {
  /*
     Assemble matrix
  */
  nnz = 4*4*(end-start)*dof*dof;

  Ke = malloc(16*dof*dof*sizeof(COEF));
  ierr = FormElementStiffness(h*h,Ke,dof);
  trace_task++;
  MPI_Barrier(MPI_COMM_WORLD);
  START_TIMER(t1);
  ierr = ZMURGE_AssemblyBegin(id, N, nnz,
                             MURGE_ASSEMBLY_ADD,  /* What to do if one entry appears twice */
                             MURGE_ASSEMBLY_ADD,  /* What to do if an entry appears on 2 processors */
                             MURGE_ASSEMBLY_FOOL, /* Do we respect the solver distribution */
                             MURGE_BOOLEAN_FALSE); CHKERRQ(ierr);
  for (ie=start; ie<end; ie++) {
    if (mode == 2)
      i = localElements[ie];
    else
      i = ie;
    /* location of lower left corner of element */
    x = h*(i % m); y = h*(i/m);
    /* node numbers for the four corners of element */
    idx[0] = (m+1)*(i/m) + (i % m);
    idx[1] = idx[0]+1;
    idx[2] = idx[1] + m + 1;
    idx[3] = idx[2] - 1;
    if (mode == 0) {
      for (j = 0; j < 4; j++)
        for (k =0; k < 4; k ++)
          ierr = ZMURGE_AssemblySetNodeValues(id, idx[j], idx[k], &(Ke[(j+k*4)*dof*dof])); CHKERRQ(ierr);
    }
    else {
      if (mode == 3) {
        memcpy(&(indexes[4*iter_blocks]),           idx, 4*sizeof(INTS));
        memcpy(&(values[4*4*dof*dof*iter_blocks]),  Ke,  4*4*dof*dof*sizeof(COEF));
        iter_blocks++;
        if (iter_blocks == n_blocks) {
          ierr = ZMURGE_AssemblySetListOfBlockValues(id, n_blocks,
                                                    4, indexes,
                                                    4, indexes,
                                                    values); CHKERRQ(ierr);
          iter_blocks = 0;
        }
      } else {
      ierr = ZMURGE_AssemblySetBlockValues(id, 4, idx, 4, idx, Ke); CHKERRQ(ierr);
    }
  }
  }
  if (mode == 3 && iter_blocks != 0) {
    ierr = ZMURGE_AssemblySetListOfBlockValues(id, iter_blocks,
                                              4, indexes,
                                              4, indexes,
                                              values); CHKERRQ(ierr);
  }
    sprintf(out_string, "First assembly loop (fact %d)", fact); 
    STOP_TIMER(out_string, t1);
    ierr = ZMURGE_AssemblyEnd(id); CHKERRQ(ierr);
    sprintf(out_string, "First assembly end (fact %d)", fact); 
    STOP_TIMER(out_string, t1);
  trace_task++;
  free(Ke);
  /*
     Create right-hand-side and solution vectors
  */
  /*
     Assemble right-hand-side vector
  */
  START_TIMER(t1);
  b = malloc(N*sizeof(COEF)*dof);
  memset(b, 0, N*sizeof(COEF)*dof);
  for (ie=start; ie<end; ie++) {
    if (mode == 2)
      i = localElements[ie];
    else
      i = ie;
    /* location of lower left corner of element */
    x = h*(i % m);
    y = h*(i/m);
    /* node numbers for the four corners of element */
    idx[0] = (m+1)*(i/m) + (i % m);
    idx[1] = idx[0]+1;
    idx[2] = idx[1] + m + 1;
    idx[3] = idx[2] - 1;
    ierr = FormElementRhs(x,y,h*h,r,1);
    for (k = 0; k < 4; k++)
      {
        for (l = 0; l  < dof; l++)
          b[idx[k]*dof+l] += r[k];
      }
  }
#ifndef FORCE_NOMPI
  {
    COEF * b_recv = malloc(N*sizeof(COEF)*dof);
    memset(b_recv, 0, N*sizeof(COEF)*dof);
    MPI_Allreduce(b, b_recv, N*dof, MURGE_MPI_COEF, MPI_SUM, MPI_COMM_WORLD);
    free(b);
    b = b_recv;
  }
#endif
  rows = malloc(4*m*sizeof(INTS));
  for (i=0; i<m+1; i++) {
    rows[i] = i; /* bottom */
    rows[3*m - 1 +i] = m*(m+1) + i; /* top */
  }
  count = m+1; /* left side */
  for (i=m+1; i<m*(m+1); i+= m+1) {
    rows[count++] = i;
  }
  count = 2*m; /* left side */
  for (i=2*m+1; i<m*(m+1); i+= m+1) {
    rows[count++] = i;
  }


  for (i=0; i<4*m; i++) {
    x = h*(rows[i] % (m+1));
    y = h*(rows[i]/(m+1));
    val = y;
    for (k = 0; k < dof; k++) {
      b[rows[i]*dof+k] = val;
    }
  }

    sprintf(out_string, "RHS assembly (fact %d)", fact); 
    STOP_TIMER(out_string, t1);
  START_TIMER(t1);
  trace_task++;
  ierr = ZMURGE_AssemblyBegin(id, N, 4*m*dof*dof,
                             MURGE_ASSEMBLY_OVW,  /* What to do if one entry appears twice */
                             MURGE_ASSEMBLY_OVW,  /* What to do if an entry appears on 2 processors */
			       MURGE_ASSEMBLY_DROPNONLOCAL, /* Do we respect the solver distribution */
                             MURGE_BOOLEAN_FALSE); CHKERRQ(ierr);
  {
    COEF * vals = malloc(dof*dof*sizeof(COEF));
    memset(vals, 0, dof*dof*sizeof(COEF));
    for (k = 0; k < dof;k++)
      vals[k+k*dof] = 1;
    for (i=0; i<4*m; i++) {
      ZMURGE_AssemblySetNodeValues(id, rows[i], rows[i], vals);
    }
    free(vals);
  }
    sprintf(out_string, "Second assembly loop (fact %d)", fact); 
    STOP_TIMER(out_string, t1);
    ierr = ZMURGE_AssemblyEnd(id); CHKERRQ(ierr);
    sprintf(out_string, "Second assembly end (fact %d)", fact); 
    STOP_TIMER(out_string, t1);
  trace_task++;

    for (solv = 0; solv < n_solv; solv++) {
  START_TIMER(t1);
  trace_task++;
  ZMURGE_SetGlobalRHS(id, b, -1, MURGE_ASSEMBLY_OVW);
  trace_task++;
      sprintf(out_string, "Setting RHS (fact %d, solv %d)", fact, solv); 
      STOP_TIMER(out_string, t1);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 -                     Solve the linear system                     -
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  START_TIMER(t1);
  u = (COEF*)malloc(N*dof*sizeof(COEF));
  trace_task++;
  {
    INTS root = -1; /* All processors */
    ZMURGE_GetGlobalSolution(id, u, root);
  }

      sprintf(out_string, "Getting Solution (fact %d, solv %d)", fact, solv); 
      STOP_TIMER(out_string, t1);
  trace_task++;
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 -                Check solution and clean up                      -
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  START_TIMER(t1);
  /* Check error */
  ustar = (COEF*)malloc(N*dof*sizeof(COEF));
  for (i=0; i<N; i++) {
    x = h*(i % (m+1)); y = h*(i/(m+1));
    val = y;
    for (k = 0; k < dof; k++) {
      ustar[i*dof+k] = val;
    }
  }

  {
    FILE * f = NULL;
    double sum1 = 0;
    double sum2 = 0;
    if (rank == 0) {
      f = fopen("plot.dat", "w");
    }
    for (k = 0; k < N*dof; k++) {
      if (N*dof < 20)
        fprintf(stdout, "u[%ld] = %g, ustar[%ld] = %g\n",
                (long)k, u[k], (long)k, ustar[k]);
      if (rank == 0) {
        fprintf(f, "%ld %ld %lg\n", (long)k%(m+1), (long)k/(m+1), u[k]);
      }
      sum1 += ustar[k]*ustar[k];
      sum2 += (ustar[k] - u[k])*(ustar[k]-u[k]);
    }
    if (rank == 0) {
      fprintf(stdout, "||u* - u||/||u*||  : %.15g\n", sqrt(sum2/sum1));
      fclose(f);
    }
  }
      sprintf(out_string, "Checking Solution (fact %d, solv %d)", fact, solv); 
      STOP_TIMER(out_string, t1);
      ZMURGE_RHSReset(id);
    }
    ZMURGE_MatrixReset(id);
  }
  /*
     Free work space.
  */
  if (localElements) free(localElements);
  free(rows);
  free(ustar);
  free(u);
  free(b);
  START_TIMER(t1);
  ZMURGE_Clean(id);
  ZMURGE_Finalize();
  STOP_TIMER("Cleaning", t1);
  STOP_TIMER("Total run", t0);
#ifndef FORCE_NOMPI
  MPI_Finalize();
#endif
  memtrace_stop();
  return 0;
}
