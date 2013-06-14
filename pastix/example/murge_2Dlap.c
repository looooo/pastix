/**
 *  File: murge_2Dlap.c
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
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#ifdef TYPE_COMPLEX
#include <complex.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef FORCE_NOMPI
#  include <mpi.h>
#endif

#include "pastix.h"
#include "murge.h"
#include "mem_trace.h"

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
#  define START_TIMER do {                                              \
    MPI_Barrier(MPI_COMM_WORLD);                                        \
    gettimeofday(&tv, NULL);                                            \
    t1 = ((double) tv.tv_sec + (double) tv.tv_usec * (double)1.0e-6L);  \
  } while(0)

#  define STOP_TIMER(str) do {                                          \
    MPI_Barrier(MPI_COMM_WORLD);                                        \
    gettimeofday(&tv, NULL);                                            \
    t2 = ((double) tv.tv_sec + (double) tv.tv_usec * (double)1.0e-6L);  \
    if (rank == 0)                                                      \
      fprintf(stdout, " $$ time for '%-40s' %.2e s $$\n", str, t2-t1);  \
    mem = pastix_getMemoryUsage();                                      \
    maxmem = pastix_getMaxMemoryUsage();                                \
    if (rank == 0)                                                      \
      fprintf(stdout, " @@ memory usage after '%-40s'"                  \
              " %.3g %s, max %.3g %s @@\n", str,                        \
              MEMORY_PRINT(mem), MEMORY_PRINT(maxmem));                 \
  } while(0)
#else /* not MURGE_TIME */
#  define START_TIMER
#  define STOP_TIMER(str)
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
#include "MURGE_GetLocalElementList.c"
#define MURGE_GetLocalElementNbr  SUFFIX(MURGE_GetLocalElementNbr)
#define MURGE_GetLocalElementList SUFFIX(MURGE_GetLocalElementList)
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
int FormElementStiffness(REAL H,COEF *Ke, int dof)
{
  int k;
  memset(Ke, 0, dof*dof*16*sizeof(COEF));
  for(k = 0; k < dof; k++)
    {
      int offset = k*dof+k;
      Ke += offset;
      Ke[0*dof*dof]  = H/6.0;
      Ke[4*dof*dof]  = -.125*H;
      Ke[8*dof*dof]  = H/12.0;
      Ke[12*dof*dof]  = -.125*H;

      Ke[1*dof*dof]  = -.125*H;
      Ke[5*dof*dof]  = H/6.0;
      Ke[9*dof*dof]  = -.125*H;
      Ke[13*dof*dof]  = H/12.0;

      Ke[2*dof*dof]  = H/12.0;
      Ke[6*dof*dof]  = -.125*H;
      Ke[10*dof*dof] = H/6.0;
      Ke[14*dof*dof] = -.125*H;

      Ke[3*dof*dof] = -.125*H;
      Ke[7*dof*dof] = H/12.0;
      Ke[11*dof*dof] = -.125*H;
      Ke[15*dof*dof] = H/6.0;
      Ke -= offset;
    }
  return 0;
}

/* --------------------------------------------------------------------- */
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
  COEF    *Ke, *Ke_cpy;      /* element matrix */
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
  MURGE_UserData_t d;
  int      mode = 0;
#ifdef MURGE_TIME
  struct timeval tv;
  double t1, t2;
  unsigned long mem, maxmem;
#endif

#ifdef FORCE_NOMPI
  rank = 0;
  size = 1;
#  define MPI_COMM_WORLD 0
#else
  required=MPI_THREAD_MULTIPLE;
  MPI_Init_thread(&argc, &argv, required, &provided);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPICHKERRQ(ierr);
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&size);MPICHKERRQ(ierr);
#endif /* FORCE_NOMPI */
  memtrace_start();

  if (argc > 1) {
    m = atoi(argv[1]);
    if (argc > 2) {
      dof = atoi(argv[2]);
      if (dof < 1)
        goto usage;
    }
    if (argc == 4) {
      mode = atoi(argv[3]);
      if (mode > 2 || mode < 0)
        goto usage;
    }
  }
  else {
    usage:
    if (rank == 0) {
      fprintf(stderr, "Usage: %s <size> [<DofNumber> <mode>]\n", argv[0]);
      return 1;
    }
  }

  N = (m+1)*(m+1);
  M = m*m;
  h = 1.0/m;
  /* Starting MURGE*/
  ierr = MURGE_Initialize((INTS)1);
  if (ierr != MURGE_SUCCESS) {
    fprintf(stderr, "Error %ld in MURGE_Initialize\n", (long)ierr);
    return 1;
  }
  id = 0;

  /* Set Options */
  prec = 1e-7;

  MURGE_SetDefaultOptions(id, 0);
  MURGE_SetOptionINT(id, IPARM_VERBOSE, API_VERBOSE_NO);
  MURGE_SetOptionINT(id, IPARM_MATRIX_VERIFICATION, API_YES);

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
  MURGE_SetOptionINT(id, IPARM_THREAD_NBR, nb_threads);
  MURGE_SetOptionINT(id, MURGE_IPARAM_DOF, dof);
  MURGE_SetOptionINT(id, MURGE_IPARAM_SYM, MURGE_BOOLEAN_FALSE);
  MURGE_SetOptionINT(id, MURGE_IPARAM_BASEVAL, 0);

  MURGE_SetOptionREAL(id, MURGE_RPARAM_EPSILON_ERROR, prec);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Compute the matrix and right-hand-side vector that define
         the linear system, Au = b.
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
     Create stiffness matrix
  */
  ierr = MURGE_SetCommunicator(id, MPI_COMM_WORLD); CHKERRQ(ierr);


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

  MURGE_SetDropNodes(id, N, rows);

  if (mode == 2) {
    d.m = m;
    trace_task++;
    START_TIMER;
    ierr = MURGE_GetLocalElementNbr(id,
                                    N,
                                    M,
                                    &localElementNbr,
                                    MURGE_DISTRIBUTE_ELEMENTS,
                                    &d); CHKERRQ(ierr);
    localElements = malloc(localElementNbr*sizeof(INTS));
    ierr = MURGE_GetLocalElementList(id, localElements); CHKERRQ(ierr);
    STOP_TIMER("getting element list");
    trace_task++;

    start = 0;
    end   = localElementNbr;
  }
  else {
    start = rank*(M/size) + ((M%size) < rank ? (M%size) : rank);
    end   = start + M/size + ((M%size) > rank);
  }
  /*
     Assemble matrix
  */
  nnz = 4*4*(end-start)*dof*dof;

  Ke = malloc(16*dof*dof*sizeof(COEF));
  Ke_cpy = malloc(16*dof*dof*sizeof(COEF));
  ierr = FormElementStiffness(h*h,Ke,dof);
  trace_task++;
  MPI_Barrier(MPI_COMM_WORLD);
  START_TIMER;
  ierr = MURGE_AssemblyBegin(id, N, nnz,
                             MURGE_ASSEMBLY_ADD,  /* What to do if one entry appears twice */
                             MURGE_ASSEMBLY_ADD,  /* What to do if an entry appears on 2 processors */
                             MURGE_ASSEMBLY_FOOL, /* Do we respect the solver distribution */
                             MURGE_BOOLEAN_FALSE); CHKERRQ(ierr);
  for (ie=start; ie<end; ie++) {
    if (mode == 2)
      i = localElements[ie];
    else
      i = ie;
    memcpy(Ke_cpy, Ke, 16*dof*dof*sizeof(COEF));
    /* location of lower left corner of element */
     x = h*(i % m); y = h*(i/m);
     /* node numbers for the four corners of element */
     idx[0] = (m+1)*(i/m) + (i % m);
     for (k =0; k < dof; k++)
       {
         for (l = 0; l < dof; l++)
           {
             Ke_cpy += l*dof+k;
             if (rows[idx[0]] > 0)
               {
                 Ke_cpy[0] = 0.;
                 Ke_cpy[4*dof*dof] = 0.;
                 Ke_cpy[8*dof*dof] = 0.;
                 Ke_cpy[12*dof*dof] = 0.;
               }
             idx[1] = idx[0]+1;
             if (rows[idx[1]] > 0)
               {
                 Ke_cpy[1*dof*dof] = 0.;
                 Ke_cpy[5*dof*dof] = 0.;
                 Ke_cpy[9*dof*dof] = 0.;
                 Ke_cpy[13*dof*dof] = 0.;
               }
             idx[2] = idx[1] + m + 1;
             if (rows[idx[2]] > 0)
               {
                 Ke_cpy[2*dof*dof] = 0.;
                 Ke_cpy[6*dof*dof] = 0.;
                 Ke_cpy[10*dof*dof] = 0.;
                 Ke_cpy[14*dof*dof] = 0.;
               }
             idx[3] = idx[2] - 1;
             if (rows[idx[3]] > 0)
               {
                 Ke_cpy[3*dof*dof] = 0.;
                 Ke_cpy[7*dof*dof] = 0.;
                 Ke_cpy[11*dof*dof] = 0.;
                 Ke_cpy[15*dof*dof] = 0.;
               }
             Ke_cpy -= l*dof+k;
           }
       }
     if (mode == 0) {
       for (j = 0; j < 4; j++)
         for (k =0; k < 4; k ++)
           ierr = MURGE_AssemblySetNodeValues(id, idx[j], idx[k], &(Ke_cpy[(j+k*4)*dof*dof])); CHKERRQ(ierr);
     }
     else {
       ierr = MURGE_AssemblySetBlockValues(id, 4, idx, 4, idx, Ke_cpy); CHKERRQ(ierr);
     }
  }
  STOP_TIMER("First assembly loop");
  ierr = MURGE_AssemblyEnd(id); CHKERRQ(ierr);
  STOP_TIMER("First assembly end");
  trace_task++;
  free(rows);
  free(Ke);
  free(Ke_cpy);
  /*
     Create right-hand-side and solution vectors
  */
  /*
     Assemble right-hand-side vector
  */
  START_TIMER;
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

  STOP_TIMER("RHS Assembly");
  START_TIMER;
  trace_task++;
  ierr = MURGE_AssemblyBegin(id, N, 4*m*dof*dof,
                             MURGE_ASSEMBLY_OVW,  /* What to do if one entry appears twice */
                             MURGE_ASSEMBLY_OVW,  /* What to do if an entry appears on 2 processors */
                             MURGE_ASSEMBLY_FOOL, /* Do we respect the solver distribution */
                             MURGE_BOOLEAN_FALSE); CHKERRQ(ierr);
  {
    COEF * vals = malloc(dof*dof*sizeof(COEF));
    memset(vals, 0, dof*dof*sizeof(COEF));
    for (k = 0; k < dof;k++)
      vals[k+k*dof] = 1;
    for (i=0; i<4*m; i++) {
      MURGE_AssemblySetNodeValues(id, rows[i], rows[i], vals);
    }
    free(vals);
  }
  STOP_TIMER("Second assembly loop");
  ierr = MURGE_AssemblyEnd(id); CHKERRQ(ierr);
  STOP_TIMER("Second assembly end");
  trace_task++;

  START_TIMER;
  trace_task++;
  MURGE_SetGlobalRHS(id, b, -1, MURGE_ASSEMBLY_OVW);
  trace_task++;
  STOP_TIMER("Setting RHS");
  free(rows);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the linear system
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  START_TIMER;
  u = (COEF*)malloc(N*dof*sizeof(COEF));
  trace_task++;
  {
    INTS root = -1; /* All processors */
    MURGE_GetGlobalSolution(id, u, root);
  }
  STOP_TIMER("Getting Solution");
  trace_task++;
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Check solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  START_TIMER;
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
    FILE * f;
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
  STOP_TIMER("CHecking solution");
  /*
     Free work space.
  */
  if (localElements) free(localElements);
  free(ustar);
  free(u);
  free(b);
  START_TIMER;
  MURGE_Clean(id);
  MURGE_Finalize();
  STOP_TIMER("Cleaning");
#ifndef FORCE_NOMPI
  MPI_Finalize();
#endif
  memtrace_stop();
  return 0;
}
