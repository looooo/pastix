/**
 * sopalin3d.h -- 
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/

#ifndef SOPALIN3D_H
#define SOPALIN3D_H

/************************************************/
/*  Structure pour le AllReduce Funneled        */
/************************************************/
/*
  Struct: Pastix_Allreduce_

  Structure used for MPI_Allreduce operations in Thread Funneled mode.
*/
typedef struct Pastix_Allreduce_ {
  void *           sendbuf;        /*+ Sending buffer                  +*/
  void *           recvbuf;        /*+ receiving buffer                +*/
  int              count;          /*+ Number of elements to reduce    +*/
  MPI_Datatype     datatype;       /*+ MPI Datatype to reduce          +*/
  MPI_Op           op;             /*+ MPI operation                   +*/
} Pastix_Allreduce_t;

/*
  Enum: SOPALIN_TASK

  tasks wich will be computed in current PaStiX call

  constants:

  SOPALIN_ONLY       - Only runs factorisation.
  SOPALIN_UPDO       - Runs factorisation and up down.
  SOPALIN_UPDO_GMRES - Runs factorisation, up down and GMRES.
  SOPALIN_UPDO_GRAD  - Runs factorisation, up down and conjugate gradient.
  SOPALIN_UPDO_PIVOT - Runs factorisation, up down and simple iterative refinement.
  UPDO_ONLY          - Only up down.
  REFINE_GMRES       - Only GMRES.
  REFINE_GRAD        - Only conjugate gradient.
  REFINE_PIVOT       - Only simple iterative refinement.
  SOPALIN_NBTASKS    - Number of existing tasks.
*/
enum SOPALIN_TASK {
  SOPALIN_ONLY = 0,
  SOPALIN_UPDO,
  SOPALIN_UPDO_GMRES,
  SOPALIN_UPDO_GRAD,
  SOPALIN_UPDO_PIVOT,
  UPDO_ONLY,
  REFINE_GMRES,
  REFINE_GRAD,
  REFINE_PIVOT,
  SOPALIN_NBTASKS
};

typedef struct starpu_task_stats_ {
  double        delay_sum;
  double        length_sum;
  unsigned int  cnt;
  unsigned long ops;
} starpu_task_stats_t;


#endif /* SOPALIN3D_H*/
