/**
 *  File: murge_3Dlap.c
 *
 *  Example that generate A 3D Laplacian with multiple
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
 *    Algiane Froehly - lacoste@labri.fr
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
#ifdef PRECISION_z
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

#define MEMORY_WRITE(mem) ( ((mem) < 1<<10) ?                           \
                            ( (double)(mem) ) :                         \
                            ( ( (mem) < 1<<20 ) ?                       \
                              ( (double)(mem)/(double)(1<<10) ) :       \
                              ( ((mem) < 1<<30 ) ?                      \
                                ( (double)(mem)/(double)(1<<20) ) :     \
                                ( (double)(mem)/(double)(1<<30) ))))
#define MEMORY_UNIT_WRITE(mem) (((mem) < 1<<10) ?       \
                                "o" :                   \
                                ( ( (mem) < 1<<20 ) ?   \
                                  "Ko" :                \
                                  ( ( (mem) < 1<<30 ) ? \
                                    "Mo" :              \
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

#  define STOP_TIMER(str, t1) do {                                      \
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
#define VERT_PER_ELEMENT(d) 8
#define GET_VERTICES(i, idx, d)                         \
  do {                                                  \
    int plan_idx = (i) / ((d)->m*(d)->m);               \
    int row_idx  = (i % ((d)->m*(d)->m)) / (d)->m;      \
    int col_idx  = (i) % (d)->m;                        \
    (idx)[0] = plan_idx*((d)->m+1)*((d)->m+1) +         \
      row_idx * ((d)->m+1) + col_idx;                   \
    (idx)[1] = (idx)[0] + 1;                            \
    (idx)[2] = (idx)[1] + (d)->m + 1;                   \
    (idx)[3] = (idx)[2] - 1;                            \
    (idx)[4] = (idx)[0] + ((d)->m + 1)*((d)->m + 1);    \
    (idx)[5] = (idx)[1] + ((d)->m + 1)*((d)->m + 1);    \
    (idx)[6] = (idx)[2] + ((d)->m + 1)*((d)->m + 1);    \
    (idx)[7] = (idx)[3] + ((d)->m + 1)*((d)->m + 1);    \
  } while (0)
#include "ZMURGE_GetLocalElementList.c"
#define ZMURGE_GetLocalElementNbr  SUFFIX(ZMURGE_GetLocalElementNbr)
#define ZMURGE_GetLocalElementList SUFFIX(ZMURGE_GetLocalElementList)
#define MURGE_UserData_t          SUFFIX(MURGE_UserData_t)
#define MURGE_UserData_           SUFFIX(MURGE_UserData_)
#undef SUFIX
//#undef VERT_PER_ELEMENT
//#undef GET_VERTICES

#define CHKERRQ(ierr)                                           \
  do {                                                          \
    if (ierr != MURGE_SUCCESS)                                  \
      {                                                         \
        fprintf(stderr, "%s:%d ERROR %ld in MURGE call\n",      \
                __FILE__, __LINE__, (long)ierr);                \
        abort();                                                \
      }                                                         \
  } while(0)

#ifdef FORCE_NOMPI
#  define MPICHKERRQ(ierr)
#else
#  define MPICHKERRQ(ierr)                                      \
  do {                                                          \
    if (ierr != MPI_SUCCESS)                                    \
      {                                                         \
        fprintf(stderr, "%s:%d ERROR %ld in MURGE call\n",      \
                __FILE__, __LINE__, (long)ierr);                \
        abort();                                                \
      }                                                         \
  } while(0)
#endif


extern int trace_task;

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 ** Pour mon pti chou: les fonctions qui te manquent pour le calcul du laplacien 3d/2d
 **^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/


/*  Pour approcher l integrale de f sur l element de reference, on la remplace par : sum_{i=1,nquad} weight_i*f(xi_i) avec weight_i les poids et xi_i les points de quadratures */

/** alloue xi et weight, des tableaux de double */

/**
 * \brief Valeur des poids et points de quadrature pour approcher les integrales
 * \param[out] xi Coordonnees des points de quadrature
 * \param[out] weight Valeur des poids
 * \param[out] nPtQuad nombre de points de quadrature
 * \param[in] dim Dimension de travail
 **/
int quadrature(double** xi, double** weight, int* nPtQuad, int dim) {

  *nPtQuad = pow(3,dim);

  (*xi) = (double*)malloc((*nPtQuad)*dim*sizeof(double));
  if ( !(*xi) ) {
    perror("  ## Memory problem: malloc");
    exit(EXIT_FAILURE);
  }
  (*weight) = (double*)malloc((*nPtQuad)*sizeof(double));
  if ( !(*xi) ) {
    perror("  ## Memory problem: malloc");
    exit(EXIT_FAILURE);
  }

  if ( dim == 1 ) {
    /* points de quadrature de Gauss-Legendre pour le segment [0;1]
     ** (3 points, methode d ordre 5) */

    (*xi)[0] = 0.5;
    (*xi)[1] = 0.5-sqrt(3./20.);
    (*xi)[2] = 0.5+sqrt(3./20.);

    (*weight)[1] = 5./18.;
    (*weight)[2] = 5./18.;
    (*weight)[0] = 1.-(*weight)[1]-(*weight)[2];

    return 1;
  }
  else if ( dim == 2 ) {
    /* points de quadraturede Gauss-Legendre pour le rectangle [0;1]x[0;1]
     ** (9 points, methode d ordre 5)) */

    (*xi)[0]  = 0.5;              (*xi)[1]  = 0.5;
    (*xi)[2]  = 0.5-sqrt(3./20.); (*xi)[3]  = 0.5-sqrt(3./20.);
    (*xi)[4]  = 0.5-sqrt(3./20.); (*xi)[5]  = 0.5+sqrt(3./20.);
    (*xi)[6]  = 0.5+sqrt(3./20.); (*xi)[7]  = 0.5-sqrt(3./20.);
    (*xi)[8]  = 0.5+sqrt(3./20.); (*xi)[9]  = 0.5+sqrt(3./20.);
    (*xi)[10] = 0.5-sqrt(3./20.); (*xi)[11] = 0.5;
    (*xi)[12] = 0.5+sqrt(3./20.); (*xi)[13] = 0.5;
    (*xi)[14] = 0.5;              (*xi)[15] = 0.5-sqrt(3./20.);
    (*xi)[16] = 0.5;              (*xi)[17] = 0.5+sqrt(3./20.);

    (*weight)[0] = 16./81.;
    (*weight)[1] = 25./324.;
    (*weight)[2] = 25./324.;
    (*weight)[3] = 25./324.;
    (*weight)[4] = 25./324.;
    (*weight)[5] = 10./81.;
    (*weight)[6] = 10./81.;
    (*weight)[7] = 10./81.;
    (*weight)[8] = 10./81.;

    return 1;
  }
  else if ( dim == 3 ) {
    /* points de quadrature de Gauss-Legendre pour le parallellepipede [0;1]x[0;1]x[0;1]
     ** (27 points, meth. d ordre 5)*/
    (*xi) = (double*)malloc(27*dim*sizeof(double));
    (*weight) = (double*)malloc(27*sizeof(double));

    (*xi)[0]  = 0.5;              (*xi)[1]  = 0.5;              (*xi)[2]  = 0.5;
    (*xi)[3]  = 0.5-sqrt(3./20.); (*xi)[4]  = 0.5-sqrt(3./20.); (*xi)[5]  = 0.5;
    (*xi)[6]  = 0.5-sqrt(3./20.); (*xi)[7]  = 0.5+sqrt(3./20.); (*xi)[8]  = 0.5;
    (*xi)[9]  = 0.5+sqrt(3./20.); (*xi)[10] = 0.5-sqrt(3./20.); (*xi)[11] = 0.5;
    (*xi)[12] = 0.5+sqrt(3./20.); (*xi)[13] = 0.5+sqrt(3./20.); (*xi)[14] = 0.5;
    (*xi)[15] = 0.5-sqrt(3./20.); (*xi)[16] = 0.5;              (*xi)[17] = 0.5;
    (*xi)[18] = 0.5+sqrt(3./20.); (*xi)[19] = 0.5;              (*xi)[20] = 0.5;
    (*xi)[21] = 0.5;              (*xi)[22] = 0.5-sqrt(3./20.); (*xi)[23] = 0.5;
    (*xi)[24] = 0.5;              (*xi)[25] = 0.5+sqrt(3./20.); (*xi)[26] = 0.5;

    (*xi)[27] = 0.5;              (*xi)[28] = 0.5;              (*xi)[29] = 0.5-sqrt(3./20.);
    (*xi)[30] = 0.5-sqrt(3./20.); (*xi)[31] = 0.5-sqrt(3./20.); (*xi)[32] = 0.5-sqrt(3./20.);
    (*xi)[33] = 0.5-sqrt(3./20.); (*xi)[34] = 0.5+sqrt(3./20.); (*xi)[35] = 0.5-sqrt(3./20.);
    (*xi)[36] = 0.5+sqrt(3./20.); (*xi)[37] = 0.5-sqrt(3./20.); (*xi)[38] = 0.5-sqrt(3./20.);
    (*xi)[39] = 0.5+sqrt(3./20.); (*xi)[40] = 0.5+sqrt(3./20.); (*xi)[41] = 0.5-sqrt(3./20.);
    (*xi)[42] = 0.5-sqrt(3./20.); (*xi)[43] = 0.5;              (*xi)[44] = 0.5-sqrt(3./20.);
    (*xi)[45] = 0.5+sqrt(3./20.); (*xi)[46] = 0.5;              (*xi)[47] = 0.5-sqrt(3./20.);
    (*xi)[48] = 0.5;              (*xi)[49] = 0.5-sqrt(3./20.); (*xi)[50] = 0.5-sqrt(3./20.);
    (*xi)[51] = 0.5;              (*xi)[52] = 0.5+sqrt(3./20.); (*xi)[53] = 0.5-sqrt(3./20.);

    (*xi)[54] = 0.5;              (*xi)[55] = 0.5;              (*xi)[56] = 0.5+sqrt(3./20.);
    (*xi)[57] = 0.5-sqrt(3./20.); (*xi)[58] = 0.5-sqrt(3./20.); (*xi)[59] = 0.5+sqrt(3./20.);
    (*xi)[60] = 0.5-sqrt(3./20.); (*xi)[61] = 0.5+sqrt(3./20.); (*xi)[62] = 0.5+sqrt(3./20.);
    (*xi)[63] = 0.5+sqrt(3./20.); (*xi)[64] = 0.5-sqrt(3./20.); (*xi)[65] = 0.5+sqrt(3./20.);
    (*xi)[66] = 0.5+sqrt(3./20.); (*xi)[67] = 0.5+sqrt(3./20.); (*xi)[68] = 0.5+sqrt(3./20.);
    (*xi)[69] = 0.5-sqrt(3./20.); (*xi)[70] = 0.5;              (*xi)[71] = 0.5+sqrt(3./20.);
    (*xi)[72] = 0.5+sqrt(3./20.); (*xi)[73] = 0.5;              (*xi)[74] = 0.5+sqrt(3./20.);
    (*xi)[75] = 0.5;              (*xi)[76] = 0.5-sqrt(3./20.); (*xi)[77] = 0.5+sqrt(3./20.);
    (*xi)[78] = 0.5;              (*xi)[79] = 0.5+sqrt(3./20.); (*xi)[80] = 0.5+sqrt(3./20.);


    (*weight)[0]  =  64./729.;
    (*weight)[1]  =  25./729.;
    (*weight)[2]  =  25./729.;
    (*weight)[3]  =  25./729.;
    (*weight)[4]  =  25./729.;
    (*weight)[5]  =  40./729.;
    (*weight)[6]  =  40./729.;
    (*weight)[7]  =  40./729.;
    (*weight)[8]  =  40./729.;
    (*weight)[9]  =  40./729.;
    (*weight)[10] = 125./5832.;
    (*weight)[11] = 125./5832.;
    (*weight)[12] = 125./5832.;
    (*weight)[13] = 125./5832.;
    (*weight)[14] =  25./729.;
    (*weight)[15] =  25./729.;
    (*weight)[16] =  25./729.;
    (*weight)[17] =  25./729.;
    (*weight)[18] =  40./729.;
    (*weight)[19] = 125./5832.;
    (*weight)[20] = 125./5832.;
    (*weight)[21] = 125./5832.;
    (*weight)[22] = 125./5832.;
    (*weight)[23] =  25./729.;
    (*weight)[24] =  25./729.;
    (*weight)[25] =  25./729.;
    (*weight)[26] =  25./729.;

    return 1;
  }
  printf("wrong value for the dimension %d\n",dim);
  return 0;
}

/**
 * \brief Basis functions
 * \param[in] xi Coordinates where the function is evaluated
 * \param[in] dim Working dimension
 * \param[out] basis Values of the basis functions
 **/
int basis(double *xi, int dim, double *basis) {

  if ( dim == 1 ) {
    basis[0] = 1-xi[0];
    basis[1] = xi[0];
  }
  else if ( dim == 2 ) {
    basis[0] = (1-xi[0])*(1-xi[1]);
    basis[1] =    xi[0] *(1-xi[1]);
    basis[2] =    xi[0] *   xi[1];
    basis[3] = (1-xi[0])*   xi[1];
  }
  else if ( dim == 3 ) {
    basis[3] = (1-xi[0])*(1-xi[1]) * (1-xi[2]);
    basis[2] =    xi[0] *(1-xi[1]) * (1-xi[2]);
    basis[6] =    xi[0] *   xi[1]  * (1-xi[2]);
    basis[7] = (1-xi[0])*   xi[1]  * (1-xi[2]);
    basis[0] = (1-xi[0])*(1-xi[1]) *    xi[2];
    basis[1] =    xi[0] *(1-xi[1]) *    xi[2];
    basis[5] =    xi[0] *   xi[1]  *    xi[2];
    basis[4] = (1-xi[0])*   xi[1]  *    xi[2];
  }

  return 1;
}

/**
 * \brief Gradient des fonctions de base
 * \param[in] xi Coordonnees du point d'evaluation du gradient
 * \param[in] dim Dimension de travail
 * \param[out] gradBasis Valeurs du gradient des fonctions de base (gradient = (d_x,d_y))
 **/
int gradient_basis(double *xi, int dim, double **gradBasis) {

  if ( dim == 1 ) {
    /* phi_1 = 1-x
     ** phi_2 = x   */
    gradBasis[0][0] = -1;
    gradBasis[0][1] =  1;
  }
  else if ( dim == 2 ) {
    /* phi_1 = (1-x)*(1-y)
     ** phi_2 =    x *(1-y)
     ** phi_3 =    x *   y
     ** phi_4 = (1-x)*   y  */
    gradBasis[0][0] = xi[1]-1.; gradBasis[1][0] = xi[0]-1.;
    gradBasis[0][1] = 1.-xi[1]; gradBasis[1][1] = -xi[0];
    gradBasis[0][2] =    xi[1]; gradBasis[1][2] =  xi[0];
    gradBasis[0][3] =   -xi[1]; gradBasis[1][3] = 1.-xi[0];

  }
  else if ( dim == 3 ) {
    /* phi_3 = (1-x)*(1-y) * (1-z)
     ** phi_2 =    x *(1-y) * (1-z)
     ** phi_6 =    x *   y  * (1-z)
     ** phi_7 = (1-x)*   y  * (1-z)
     ** phi_0 = (1-x)*(1-y) *    z
     ** phi_1 =    x *(1-y) *    z
     ** phi_5 =    x *   y  *    z
     ** phi_4 = (1-x)*   y  *    z
     */
    gradBasis[0][3] = (xi[1]-1.)*(1.-xi[2]);
    gradBasis[1][3] = (xi[0]-1.)*(1.-xi[2]);
    gradBasis[2][3] = (xi[0]-1.)*(1.-xi[1]);

    gradBasis[0][2] = (1.-xi[1])*(1.-xi[2]);
    gradBasis[1][2] = -xi[0]    *(1.-xi[2]);
    gradBasis[2][2] = -xi[0]    *(1.-xi[1]);

    gradBasis[0][6] =  xi[1]*(1.-xi[2]);
    gradBasis[1][6] =  xi[0]*(1.-xi[2]);
    gradBasis[2][6] =  -xi[0]*xi[1];

    gradBasis[0][7] =    -xi[1] *(1.-xi[2]);
    gradBasis[1][7] = (1.-xi[0])*(1.-xi[2]);;
    gradBasis[2][7] = -(1.-xi[0])*xi[1];

    gradBasis[0][0] = (xi[1]-1.)*xi[2];
    gradBasis[1][0] = (xi[0]-1.)*xi[2];
    gradBasis[2][0] = (1.-xi[0])*(1.-xi[1]);

    gradBasis[0][1] = (1.-xi[1])*xi[2];
    gradBasis[1][1] = -xi[0]    *xi[2];
    gradBasis[2][1] =  xi[0]    *(1.-xi[1]);

    gradBasis[0][5] =  xi[1]*xi[2];
    gradBasis[1][5] =  xi[0]*xi[2];
    gradBasis[2][5] =  xi[0]*xi[1];

    gradBasis[0][4] = -xi[1]*xi[2];
    gradBasis[1][4] = (1.-xi[0])*xi[2];
    gradBasis[2][4] = (1.-xi[0])*xi[1];
  }

  return 1;
}

/**
 * \brief Determinant de la matrice jacobienne de la transformation element_reference->element_reel
 * \param[in] coorEl Coordonnees de l'element reel
 * \param[in] gradBasis Valeurs du gradient des fonctions de base
 * \param[in] nodesPerEl Nombre de noeuds par element
 * \param[in] dim Dimension de travail
 * \result Valeur du determinant
 **/
double det_jacobian(double **coorEl,double **gradBasis,int nodesPerEl,int dim){
  double det, jaco[dim][dim];
  int i,j,k;

  /* init */
  det = 0;
  for ( i=0; i< dim; i++ ) {
    for ( j=0; j< dim; j++ ) {
      jaco[i][j] = 0.;
    }
  }

  /* jacobian */
  for ( i=0; i< dim; i++ ) {
    for ( j=0; j< dim; j++ ) {
      for ( k=0; k< nodesPerEl; k++ ) {
        jaco[i][j] = jaco[i][j] + coorEl[i][k]*gradBasis[j][k];
      }
    }
  }

  /* determinant */
  if ( dim > 1 ) {
    det =  jaco[0][0]*jaco[1][1]-jaco[1][0]*jaco[0][1];
    if ( dim==3 ) {
      det =  det*jaco[2][2]-jaco[2][1]*( jaco[0][0]*jaco[1][2]-jaco[1][0]*jaco[0][2])
        +jaco[2][0]*( jaco[0][1]*jaco[1][2]-jaco[1][1]*jaco[0][2]);
    }
  }

  return det;
}


/**
 * \brief Le laplacien sur segment, quadrangle 2d ou pave 3d suivant dim
 * (resp. matrice 2x2, 4x4, 8x8)
 * \param[in] coorEl Coordonnees des sommets de l'element (reel)
 * \param[in] dim Dimension de travail
 * \param[in] nodesPerEl Nombre de noeuds par elements
 * \param[out] matLapla matrice laplacienne
 **/
int laplacian(double** coorEl,int dim,int nodesPerEl,COEF* matLapla, int dofnbr) {

  /* for quadrature */
  double *xi=NULL,*weight=NULL;
  int nPtQuad;

  /* to store the gradient of basis functions */
  double **gradBasis=NULL;


  /* determinant of the jacobian of the transformation */
  double detJac;

  int i,j,d,idof, iQuad, addr,adr_quad;

  if ( nodesPerEl != pow(2,dim) ) {
    printf("ERROR: prevu seulement pour des elements P1");
    exit(EXIT_FAILURE);
  }


  /* Allocation of the matrix of the basis functions gradients */
  gradBasis = (double**)malloc(dim*sizeof(double*));
  if ( !gradBasis ) {
    perror("  ## Memory problem: malloc");
    exit(EXIT_FAILURE);
  }
  for ( i=0; i < dim; i++ ) {
    gradBasis[i] = (double*)malloc(nodesPerEl*sizeof(double));
    if ( !gradBasis[i] ) {
      perror("  ## Memory problem: malloc");
      exit(EXIT_FAILURE);
    }
  }


  if ( !quadrature(&xi,&weight,&nPtQuad,dim) ) exit(EXIT_FAILURE);

  /* init */
  addr = 0;
  for( i=0; i<nodesPerEl; i++ ) {
    for( j=0; j<nodesPerEl; j++ ) {
      matLapla[addr] = 0.0;
      addr++;
    }
  }

  /* laplacian computation :
   **  laplacien_ij = sum_{k=1,nPtQuad}  (  weight_k *
   **      (d_x basisFunction_i(xi_k) *  d_x basisFunction_j(xi_k) +
   **       d_y basisFunction_i(xi_k) *  d_y basisFunction_j(xi_k) )
   **      * det(JacobianTransfo(xi_k)   )  */
  for ( iQuad=0; iQuad < nPtQuad; iQuad++ ) {
    addr=0;
    adr_quad = iQuad*dim;

    if ( !gradient_basis(&(xi[adr_quad]),dim,gradBasis) )
      exit(EXIT_FAILURE);

    detJac = det_jacobian(coorEl,gradBasis,nodesPerEl,dim);

    for( i=0; i<nodesPerEl; i++ ) {
      for( j=0; j<nodesPerEl; j++ ) {
        for( d=0; d<dim; d++) {
          for (idof=0; idof < dofnbr; idof++) {
            matLapla[addr*dofnbr*dofnbr+idof*(dofnbr+1)] =
              matLapla[addr*dofnbr*dofnbr+idof*(dofnbr+1)] + detJac * weight[iQuad] * gradBasis[d][i] * gradBasis[d][j];
          }
        }
        addr++;
      }
    }
  }


  /* Memory free */
  free(xi);
  xi = NULL;
  free(weight);
  weight = NULL;

  for ( i=0; i<dim; i++ ) {
    free(gradBasis[i]);
    gradBasis[i] = NULL;
  }
  free(gradBasis);
  gradBasis = NULL;

  return(1);
}

/**
 * \brief Call finite element laplacian for pastix
 * \param[in] h Grid step
 * \param[in,out] Ke Elementary laplacian matrix
 * \param[in] dof Number of degrees of freedom per node
 * \param[in] i Index of the element
 * \param[m] m Number of elements per dimension
 **/
int FormElementStiffness(double h, COEF *Ke, int dof,int i,int m){
  /* h = pas de la grille
   ** Ke = matrice elementaire du laplacien
   ** dof = nombre de dofs par inconnue
   ** i = indice element
   ** m = nb element par dimension */
  int dim = 3, nodesPerEl = 8;

  double ** coorEl;

  /** note : h is unused */
  /* if ( dof !=1 ) { */
  /*   printf( " dof !=1: cas non prevu\n" ); */
  /*   return 0; */
  /* } */

  coorEl = (double**)malloc(dim*sizeof(double*));
  if ( !coorEl ) {
    perror("  ## Memory problem: malloc");
    exit(EXIT_FAILURE);
  }
  for ( i=0; i < dim; i++ ) {
    coorEl[i] = (double*)malloc(nodesPerEl*sizeof(double));
    if ( !coorEl[i] ) {
      perror("  ## Memory problem: malloc");
      exit(EXIT_FAILURE);
    }
  }


  /*
   **            4.-----------.5
   **            /|           /|
   **           / |          / |
   **         0.------------.1 |
   **          |  |         |  |
   **          |  .7--------|--.6
   **          | /          | /
   **          |/           |/
   **         3.-----------.2
   */

  /** Coordonnees des points de l element */
  /* Pt 0 */
  coorEl[0][0] = ((double)(i%m)) * h;
  coorEl[1][0] = ((double)(i/(m*m))) * h;
  coorEl[2][0] = ((double)(m-(i%(m*m))/m)) * h;

  /* Pt 1 */
  coorEl[0][1] = coorEl[0][0] + h;
  coorEl[1][1] = coorEl[1][0];
  coorEl[2][1] = coorEl[2][0];

  /* Pt 2 */
  coorEl[0][2] = coorEl[0][0] + h;
  coorEl[1][2] = coorEl[1][0];
  coorEl[2][2] = coorEl[2][0] - h;

  /* Pt 3 */
  coorEl[0][3] = coorEl[0][0];
  coorEl[1][3] = coorEl[1][0];
  coorEl[2][3] = coorEl[2][0] - h;

  /* Pt 4 */
  coorEl[0][4] = coorEl[0][0];
  coorEl[1][4] = coorEl[1][0] + h;
  coorEl[2][4] = coorEl[2][0];

  /* Pt 5 */
  coorEl[0][5] = coorEl[0][0] + h;
  coorEl[1][5] = coorEl[1][0] + h;
  coorEl[2][5] = coorEl[2][0];

  /* Pt 6 */
  coorEl[0][6] = coorEl[0][0] + h;
  coorEl[1][6] = coorEl[1][0] + h;
  coorEl[2][6] = coorEl[2][0] - h;

  /* Pt 7 */
  coorEl[0][7] = coorEl[0][0];
  coorEl[1][7] = coorEl[1][0] + h;
  coorEl[2][7] = coorEl[2][0] - h;



  if ( !laplacian(coorEl,dim,nodesPerEl,Ke, dof) ) return 0;

  for (i=0;i<dim;i++){
    free(coorEl[i]);
    coorEl[i] = NULL;
  }
  free(coorEl);
  coorEl=NULL;
  return 1;
}

/* --------------------------------------------------------------------- */
static inline
int FormElementRhs(REAL x,REAL y,REAL H,COEF *r, int dof, int i_elt, int m)
{
  int k;
  int plan_idx = (i_elt) / (m*m);
  int row_idx  = (i_elt % (m*m)) / m;
  int col_idx  = (i_elt) % m;
  for (k =0; k < dof; k++) {
    r[0*dof+k] = 0.;
    r[1*dof+k] = 0.;
    r[2*dof+k] = 0.;
    r[3*dof+k] = 0.;
    r[4*dof+k] = 0.;
    r[5*dof+k] = 0.;
    r[6*dof+k] = 0.;
    r[7*dof+k] = 0.;
  }

  if (plan_idx == 0) {
    r[0] += 1;
    r[1] += 1;
    r[2] += 1;
    r[3] += 1;
  }
  if (plan_idx == m-1) {
    r[4+0] += 1;
    r[4+1] += 1;
    r[4+2] += 1;
    r[4+3] += 1;
  }
  if (col_idx == 0) {
    int i;
    int idx[4] = {0, 3, 4, 7};
    for (i = 0; i < 4; i++) {
      r[idx[i]] += 1;
    }
  }
  if (col_idx == m-1) {
    int i;
    int idx[4] = {1, 2, 5, 6};
    for (i = 0; i < 4; i++) {
      r[idx[i]] += 1;
    }
  }
  if (row_idx == 0) {
    int i;
    int idx[4] = {0, 1, 4, 5};
    for (i = 0; i < 4; i++) {
      r[idx[i]] += 1;
    }
  }
  if (row_idx == m-1) {
    int i;
    int idx[4] = {3, 2, 7, 6};
    for (i = 0; i < 4; i++) {
      r[idx[i]] += 1;
    }
  }

  return 0;
}

/* --------------------------------------------------------------------- */
/** WARNING : DON'T WORK FOR DIMENSION 1 OR 2
 * \brief Compute elementary right hand side for the finite element method
 * \param[in] x ??
 * \param[in] y ??
 * \param[in] H Grid Step
 * \param[out] *r Right hand side
 * \param[in] dof Number of degrees of freedom per node
 * \param[in] i_elt Index of the element
 * \param[in] m Number of element per dimensions
 **/
static inline
int FormFiniteElementRhs(REAL x,REAL y,REAL H,COEF *r, int dof, int i_elt, int m)
{
  int k,l;
  int plan_idx = (i_elt) / (m*m);
  int row_idx  = (i_elt % (m*m)) / m;
  int col_idx  = (i_elt) % m;
  double *weight=NULL,*xi=NULL,bas[4], **gradBasis=NULL, coorEl[2][4];
  double sumGradient[2], normale[2];
  int nPtQuad,nodesPerEl=4,dim=3;

  switch (dim) {
  case 1:
  case 2:
    printf("Case not available (dimension 1 or 2)\n");
    exit(EXIT_FAILURE);
    break;
  case 3:
    break;
  default:
    printf("wrong dimension\n");
    exit(EXIT_FAILURE);
  }

  /* Allocation of the matrix of the basis functions gradients */
  gradBasis = (double**)malloc((dim-1)*sizeof(double*));
  if ( !gradBasis ) {
    perror("  ## Memory problem: malloc");
    exit(EXIT_FAILURE);
  }
  for ( k=0; k < dim-1; k++ ) {
    gradBasis[k] = (double*)malloc(nodesPerEl*sizeof(double));
    if ( !gradBasis[k] ) {
      perror("  ## Memory problem: malloc");
      exit(EXIT_FAILURE);
    }
  }


  /* We integrate on the boundaries so we work in dim-1 */
  if ( !quadrature(&xi,&weight,&nPtQuad,dim-1) ) exit( EXIT_FAILURE );

  for (k =0; k < dof; k++) {
    r[0*dof+k] = 0.;
    r[1*dof+k] = 0.;
    r[2*dof+k] = 0.;
    r[3*dof+k] = 0.;
    r[4*dof+k] = 0.;
    r[5*dof+k] = 0.;
    r[6*dof+k] = 0.;
    r[7*dof+k] = 0.;
  }

  // On the boundaries :
  // r[i] = sum_k [ w_k basis[i](xi_k) *
  //                [ (sum_l Y_l * gradient_basis[1][l](xi_k))*(sum_l gradient_basis[0][l])
  //                 -(sum_l X_l * gradient_basis[0][l](xi_k))*(sum_l gradient_basis[1][l]) ]

  if (plan_idx == 0) {  // y=0,
    // to work in two dimension, we store the z coordinate in coorEl[1][.]

    /*
     **         0.-----------.1
     **          |3         2|
     **          |           |
     **          |           |
     **          |0         1|
     **         3.-----------.2
     */

    // Local coordinates on the surface
    coorEl[0][3] = ((double)(i_elt%m)) * H;           // Pt 0 in 3d numerotation
    coorEl[1][3] = ((double)(m-(i_elt%(m*m))/m)) * H; // Pt 0 in 3d numerotation
    coorEl[0][2] = coorEl[0][3] + H; // Pt 1 in 3d numerotation
    coorEl[1][2] = coorEl[1][3];     // Pt 1 in 3d numerotation
    coorEl[0][1] = coorEl[0][3] + H; // Pt 2 in 3d numerotation
    coorEl[1][1] = coorEl[1][3] - H; // Pt 2 in 3d numerotation
    coorEl[0][0] = coorEl[0][3];     // Pt 3 in 3d numerotation
    coorEl[1][0] = coorEl[1][3] - H; // Pt 3 in 3d numerotation

    sumGradient[0] = 0; sumGradient[1] = 0;
    normale[0] = 0; normale[1] = 0;

    for ( k=0; k < nPtQuad; k++ ) {

      int addr = k*(dim-1);

      if ( !basis(&(xi[addr]),dim-1,bas) ) exit( EXIT_FAILURE );
      if ( !gradient_basis(&(xi[addr]),dim-1,gradBasis) ) exit( EXIT_FAILURE );

      for ( l = 0; l<nodesPerEl; l++ ) {
        sumGradient[0] += gradBasis[0][l];
        sumGradient[1] += gradBasis[1][l];

        normale[0] += coorEl[1][l] * gradBasis[1][l];
        normale[1] -= coorEl[0][l] * gradBasis[0][l];

      }
      r[0] *= weight[k]*bas[3]*(normale[0]*sumGradient[0]+
                                normale[1]*sumGradient[1]);
      r[1] *= weight[k]*bas[2]*(normale[0]*sumGradient[0]+
                                normale[1]*sumGradient[1]);
      r[2] *= weight[k]*bas[1]*(normale[0]*sumGradient[0]+
                                normale[1]*sumGradient[1]);
      r[3] *= weight[k]*bas[0]*(normale[0]*sumGradient[0]+
                                normale[1]*sumGradient[1]);
    }
  }

  if (plan_idx == m-1) { // y = H
    // to work in two dimension, we store the z coordinate in coorEl[1][.]
    /*
     **         5.-----------.4
     **          |3         2|
     **          |           |
     **          |           |
     **          |0         1|
     **         6.-----------.7
     */

    // Local coordinates on surface
    coorEl[0][2] = ((double)(i_elt%m)) * H;           // Pt 4 in 3d numerotation
    coorEl[1][2] = ((double)(m-(i_elt%(m*m))/m)) * H; // Pt 4 in 3d numerotation
    coorEl[0][3] = coorEl[0][2] + H; // Pt 5 in 3d numerotation
    coorEl[1][3] = coorEl[1][2];     // Pt 5 in 3d numerotation
    coorEl[0][0] = coorEl[0][2] + H; // Pt 6 in 3d numerotation
    coorEl[1][0] = coorEl[1][2] - H; // Pt 6 in 3d numerotation
    coorEl[0][1] = coorEl[0][2];     // Pt 7 in 3d numerotation
    coorEl[1][1] = coorEl[1][2] - H; // Pt 7 in 3d numerotation

    sumGradient[0] = 0; sumGradient[1] = 0;
    normale[0] = 0; normale[1] = 0;

    for ( k=0; k < nPtQuad; k++ ) {

      int addr = k*(dim-1);

      if ( !basis(&(xi[addr]),dim-1,bas) ) exit( EXIT_FAILURE );
      if ( !gradient_basis(&(xi[addr]),dim-1,gradBasis) ) exit( EXIT_FAILURE );

      for ( l = 0; l<nodesPerEl; l++ ) {
        sumGradient[0] += gradBasis[0][l];
        sumGradient[1] += gradBasis[1][l];

        normale[0] += coorEl[1][l] * gradBasis[1][l];
        normale[1] -= coorEl[0][l] * gradBasis[0][l];

      }
      r[4+0] *= weight[k]*bas[2]*(normale[0]*sumGradient[0]+
                                  normale[1]*sumGradient[1]);
      r[4+1] *= weight[k]*bas[3]*(normale[0]*sumGradient[0]+
                                  normale[1]*sumGradient[1]);
      r[4+2] *= weight[k]*bas[0]*(normale[0]*sumGradient[0]+
                                  normale[1]*sumGradient[1]);
      r[4+3] *= weight[k]*bas[1]*(normale[0]*sumGradient[0]+
                                  normale[1]*sumGradient[1]);
    }
  }

  if (col_idx == 0) {  // x=0
    int idx[4] = {0, 3, 4, 7};
    // to work in two dimension, we store the (y,z) coordinate in (coorEl[0][.],coorEl[1][.])
    /*
     **         4.-----------.0
     **          |3         2|
     **          |           |
     **          |           |
     **          |0         1|
     **         7.-----------.3
     */

    // Local coordinates on surface
    coorEl[0][2] = ((double)(i_elt/(m*m))) * H;       // Pt 0 in 3d numerotation
    coorEl[1][2] = ((double)(m-(i_elt%(m*m))/m)) * H; // Pt 0 in 3d numerotation
    coorEl[0][1] = coorEl[0][2];     // Pt 3 in 3d numerotation
    coorEl[1][1] = coorEl[1][2] - H; // Pt 3 in 3d numerotation
    coorEl[0][0] = coorEl[0][2] + H; // Pt 7 in 3d numerotation
    coorEl[1][0] = coorEl[1][2] - H; // Pt 7 in 3d numerotation
    coorEl[0][3] = coorEl[0][2] + H; // Pt 4 in 3d numerotation
    coorEl[1][3] = coorEl[1][2];     // Pt 4 in 3d numerotation

    sumGradient[0] = 0; sumGradient[1] = 0;
    normale[0] = 0; normale[1] = 0;

    for ( k=0; k < nPtQuad; k++ ) {

      int addr = k*(dim-1);

      if ( !basis(&(xi[addr]),dim-1,bas) ) exit( EXIT_FAILURE );
      if ( !gradient_basis(&(xi[addr]),dim-1,gradBasis) ) exit( EXIT_FAILURE );

      for ( l = 0; l<nodesPerEl; l++ ) {
        sumGradient[0] += gradBasis[0][l];
        sumGradient[1] += gradBasis[1][l];

        normale[0] += coorEl[1][l] * gradBasis[1][l];
        normale[1] -= coorEl[0][l] * gradBasis[0][l];

      }
      r[idx[0]] *= weight[k]*bas[2]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
      r[idx[1]] *= weight[k]*bas[1]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
      r[idx[2]] *= weight[k]*bas[3]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
      r[idx[3]] *= weight[k]*bas[0]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
    }
  }

  if (col_idx == m-1) { //x=H
    // to work in two dimension, we store the (y,z) coordinate in (coorEl[0][.],coorEl[1][.])
    int idx[4] = {1, 2, 5, 6};

    /*
     **         1.-----------.5
     **          |3         2|
     **          |           |
     **          |           |
     **          |0         1|
     **         2.-----------.6
     */


    // Local coordinates on surface
    coorEl[0][3] = ((double)(i_elt/(m*m))) * H;       // Pt 1 in 3d numerotation
    coorEl[1][3] = ((double)(m-(i_elt%(m*m))/m)) * H; // Pt 1 in 3d numerotation
    coorEl[0][2] = coorEl[0][3] + H; // Pt 5 in 3d numerotation
    coorEl[1][2] = coorEl[1][3];     // Pt 5 in 3d numerotation
    coorEl[0][1] = coorEl[0][3] + H; // Pt 6 in 3d numerotation
    coorEl[1][1] = coorEl[1][3] - H; // Pt 6 in 3d numerotation
    coorEl[0][0] = coorEl[0][3];     // Pt 2 in 3d numerotation
    coorEl[1][0] = coorEl[1][3] - H; // Pt 2 in 3d numerotation

    sumGradient[0] = 0; sumGradient[1] = 0;
    normale[0] = 0; normale[1] = 0;

    for ( k=0; k < nPtQuad; k++ ) {

      int addr = k*(dim-1);

      if ( !basis(&(xi[addr]),dim-1,bas) ) exit( EXIT_FAILURE );
      if ( !gradient_basis(&(xi[addr]),dim-1,gradBasis) ) exit( EXIT_FAILURE );

      for ( l = 0; l<nodesPerEl; l++ ) {
        sumGradient[0] += gradBasis[0][l];
        sumGradient[1] += gradBasis[1][l];

        normale[0] += coorEl[1][l] * gradBasis[1][l];
        normale[1] -= coorEl[0][l] * gradBasis[0][l];

      }
      r[idx[0]] *= weight[k]*bas[3]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
      r[idx[1]] *= weight[k]*bas[0]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
      r[idx[2]] *= weight[k]*bas[2]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
      r[idx[3]] *= weight[k]*bas[1]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
    }
  }

  if (row_idx == 0) {  // z=0
    // to work in two dimension, we store the (x,y) coordinate in (coorEl[0][.],coorEl[1][.])
    int idx[4] = {0, 1, 4, 5};

    /*
     **         4.-----------.5
     **          |3         2|
     **          |           |
     **          |           |
     **          |0         1|
     **         0.-----------.1
     */


    // Local coordinates on surface
    coorEl[0][0] = ((double)(i_elt%m)) * H;     // Pt 0 in 3d numerotation
    coorEl[1][0] = ((double)(i_elt/(m*m))) * H; // Pt 0 in 3d numerotation
    coorEl[0][1] = coorEl[0][0] + H; // Pt 1 in 3d numerotation
    coorEl[1][1] = coorEl[1][0];     // Pt 1 in 3d numerotation
    coorEl[0][2] = coorEl[0][0] + H; // Pt 5 in 3d numerotation
    coorEl[1][2] = coorEl[1][0] + H; // Pt 5 in 3d numerotation
    coorEl[0][3] = coorEl[0][0];     // Pt 4 in 3d numerotation
    coorEl[1][3] = coorEl[1][0] + H; // Pt 4 in 3d numerotation

    sumGradient[0] = 0; sumGradient[1] = 0;
    normale[0] = 0; normale[1] = 0;

    for ( k=0; k < nPtQuad; k++ ) {

      int addr = k*(dim-1);

      if ( !basis(&(xi[addr]),dim-1,bas) ) exit( EXIT_FAILURE );
      if ( !gradient_basis(&(xi[addr]),dim-1,gradBasis) ) exit( EXIT_FAILURE );

      for ( l = 0; l<nodesPerEl; l++ ) {
        sumGradient[0] += gradBasis[0][l];
        sumGradient[1] += gradBasis[1][l];

        normale[0] += coorEl[1][l] * gradBasis[1][l];
        normale[1] -= coorEl[0][l] * gradBasis[0][l];

      }
      r[idx[0]] *= weight[k]*bas[0]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
      r[idx[1]] *= weight[k]*bas[1]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
      r[idx[2]] *= weight[k]*bas[3]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
      r[idx[3]] *= weight[k]*bas[2]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
    }
  }

  if (row_idx == m-1) { // z=H
    // to work in two dimension, we store the (x,y) coordinate in (coorEl[0][.],coorEl[1][.])

    int idx[4] = {3, 2, 7, 6};

    /*
     **         3.-----------.2
     **          |3         2|
     **          |           |
     **          |           |
     **          |0         1|
     **         7.-----------.6
     */


    // Local coordinates on surface
    coorEl[0][3] =  ((double)(i_elt%m)) * H;    ; // Pt 3 in 3d numerotation
    coorEl[1][3] =  ((double)(i_elt/(m*m))) * H;; // Pt 3 in 3d numerotation
    coorEl[0][2] =  coorEl[0][3] + H; // Pt 2 in 3d numerotation
    coorEl[1][2] =  coorEl[1][3];     // Pt 2 in 3d numerotation
    coorEl[0][1] =  coorEl[0][3] + H; // Pt 6 in 3d numerotation
    coorEl[1][1] =  coorEl[1][3] + H; // Pt 6 in 3d numerotation
    coorEl[0][0] =  coorEl[0][3];     // Pt 7 in 3d numerotation
    coorEl[1][0] =  coorEl[1][3] + H; // Pt 7 in 3d numerotation

    sumGradient[0] = 0; sumGradient[1] = 0;
    normale[0] = 0; normale[1] = 0;

    for ( k=0; k < nPtQuad; k++ ) {

      int addr = k*(dim-1);

      if ( !basis(&(xi[addr]),dim-1,bas) ) exit( EXIT_FAILURE );
      if ( !gradient_basis(&(xi[addr]),dim-1,gradBasis) ) exit( EXIT_FAILURE );

      for ( l = 0; l<nodesPerEl; l++ ) {
        sumGradient[0] += gradBasis[0][l];
        sumGradient[1] += gradBasis[1][l];

        normale[0] += coorEl[1][l] * gradBasis[1][l];
        normale[1] -= coorEl[0][l] * gradBasis[0][l];

      }
      r[idx[0]] *= weight[k]*bas[3]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
      r[idx[1]] *= weight[k]*bas[2]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
      r[idx[2]] *= weight[k]*bas[0]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
      r[idx[3]] *= weight[k]*bas[1]*(normale[0]*sumGradient[0]+
                                     normale[1]*sumGradient[1]);
    }
  }

  /* Memory free */
  free(xi);
  xi = NULL;
  free(weight);
  weight = NULL;

  for ( k=0; k<dim-1; k++ ) {
    free(gradBasis[k]);
    gradBasis[k] = NULL;
  }
  free(gradBasis);
  gradBasis = NULL;

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
  COEF     r[8];        /* element vector */
  REAL     h;           /* mesh width */
  REAL     x,y;
  COEF     val;
  INTS     ierr;
  INTS     idx[8],/* count, */ *rows, i, ie, m = 5,dof=1, start,end,j, k, l;
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
  if (rank == 0) {
    switch (provided) {
    case MPI_THREAD_SINGLE:
      printf("MPI_Init_thread level = MPI_THREAD_SINGLE\n");
      break;
    case MPI_THREAD_FUNNELED:
      printf("MPI_Init_thread level = MPI_THREAD_FUNNELED\n");
      break;
    case MPI_THREAD_SERIALIZED:
      printf("MPI_Init_thread level = MPI_THREAD_SERIALIZED\n");
      break;
    case MPI_THREAD_MULTIPLE:
      printf("MPI_Init_thread level = MPI_THREAD_MULTIPLE\n");
      break;
    default:
      printf("MPI_Init_thread level = ???\n");
    }
  }
#endif /* FORCE_NOMPI */
  memtrace_start();
  m = 100;
  while ((opt = getopt(argc, argv, "n:d:m:b:f:s:t:h")) != -1) {
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
    case 't':
      nb_threads = atoi(optarg);
      if (nb_threads < 1) goto usage;
      break;
    case 'h':
    default:
    usage:
      printf("unprocessed option -%c %s\n", opt, optarg);
      printf("usage: %s -n <size> -d <DofNumber> -m <mode> -b <nblocks> -s <solve number> -f <factorization number> -t <number of threads>]\n", argv[0]);
      return 1;
    }
  }

  if (mode == 3) {
    indexes = malloc(n_blocks*8*sizeof(INTS));
    values  = malloc(n_blocks*8*8*dof*dof*sizeof(COEF));
  }

  N = (m+1)*(m+1)*(m+1);
  M = m*m*m;
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
  ZMURGE_SetOptionINT(id, IPARM_VERBOSE, API_VERBOSE_YES);
  ZMURGE_SetOptionINT(id, IPARM_MATRIX_VERIFICATION, API_NO);

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
    /* CSCd Required for product in verification */
  ZMURGE_SetOptionINT(id, IPARM_FREE_CSCUSER, API_CSC_PRESERVE);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Compute the matrix and right-hand-side vector that define
   the linear system, Au = b.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  /*
   Create stiffness matrix
   */
#ifndef FORCE_NOMPI
  ierr = ZMURGE_SetCommunicator(id, MPI_COMM_WORLD); CHKERRQ(ierr);
#endif

  /*
     Mark rows for for Dirichlet boundary conditions
  */
  rows = malloc(N*sizeof(INTS));
  memset(rows, 0, N*sizeof(INTS));

  for (i=0; i<m+1; i++) {
    for (j=0; j < m+1; j++) {
      rows[j*(m+1) + i] = 1;                 /* front */
      rows[m*(m+1)*(m+1) + j*(m+1) + i] = 1; /* back */
      rows[j*(m+1)*(m+1) + i] = 1;           /* top */
      rows[j*(m+1)*(m+1) + m*(m+1) + i] = 1; /* bottom */
      /* fprintf(stdout, "Drop %d %d %d %d\n", j*(m+1) + i, */
      /*         m*(m+1)*(m+1) + j*(m+1) + i, */
      /*         j*(m+1)*(m+1) + i, */
      /*         j*(m+1)*(m+1) + m*(m+1) + i); */
    }
  }
  ZMURGE_SetDropRows(id, N, rows);
  fprintf(stdout, "mode %d\n", mode);
  d.m = m;
  if (mode == 2) {
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
    nnz = 8*8*(end-start)*dof*dof;

    Ke = malloc(8*8*dof*dof*sizeof(COEF));
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
      ierr = FormElementStiffness(h,Ke,dof,i,m);
      /* { */
      /*   int i,j; */
      /*   for (i = 0; i < 8; i++) { */
      /*     for (j = 0; j < 8; j++) */
      /*       fprintf(stdout, "%12.5g ", Ke[j+i*8]); */
      /*     fprintf(stdout, "\n"); */
      /*   } */
      /* } */

      /* location of lower left corner of element */
      x = h*(i % m); y = h*(i/m);
      /* node numbers for the four corners of element */
      GET_VERTICES(i, idx, &d);
      if (mode == 0) {
        for (j = 0; j < 8; j++)
          for (k =0; k < 8; k ++) {
            ierr = ZMURGE_AssemblySetNodeValues(id, idx[j], idx[k], &(Ke[(j+k*8)*dof*dof])); CHKERRQ(ierr);
          }
      }
      else {
        if (mode == 3) {
          memcpy(&(indexes[8*iter_blocks]),           idx, 8*sizeof(INTS));
          memcpy(&(values[8*8*dof*dof*iter_blocks]),  Ke,  8*8*dof*dof*sizeof(COEF));
          iter_blocks++;
          if (iter_blocks == n_blocks) {
            ierr = ZMURGE_AssemblySetListOfBlockValues(id, n_blocks,
                                                      8, indexes,
                                                      8, indexes,
                                                      values); CHKERRQ(ierr);
            iter_blocks = 0;
          }
        } else {
          ierr = ZMURGE_AssemblySetBlockValues(id, 8, idx, 8, idx, Ke); CHKERRQ(ierr);
        }
      }
    }
    if (mode == 3 && iter_blocks != 0) {
      ierr = ZMURGE_AssemblySetListOfBlockValues(id, iter_blocks,
                                                8, indexes,
                                                8, indexes,
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
      GET_VERTICES(i, idx, &d);
      ierr = FormFiniteElementRhs(x,y,h*h,r,1,i,m);
      /* { */
      /*   int i,j; */
      /*   for (i = 0; i < 1; i++) { */
      /*     for (j = 0; j < 8; j++) */
      /*       fprintf(stdout, "r[%d] %12.5g \n", j+i*8, r[j+i*8]); */
      /*     fprintf(stdout, "\n"); */
      /*   } */
      /* } */

      for (k = 0; k < 8; k++)
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
    /* rows = malloc(8*m*sizeof(INTS)); */
    /* for (i=0; i<m+1; i++) { */
    /*   rows[i] = i; /\* bottom *\/ */
    /*   rows[3*m - 1 +i] = m*(m+1) + i; /\* top *\/ */
    /* } */
    /* count = m+1; /\* left side *\/ */
    /* for (i=m+1; i<m*(m+1); i+= m+1) { */
    /*   rows[count++] = i; */
    /* } */
    /* count = 2*m; /\* left side *\/ */
    /* for (i=2*m+1; i<m*(m+1); i+= m+1) { */
    /*   rows[count++] = i; */
    /* } */

    ierr = ZMURGE_AssemblyBegin(id, N, (2*(m+1)*(m+1) + 2*(m-1)*(m+1))*dof*dof,
                               MURGE_ASSEMBLY_OVW,  /* What to do if one entry appears twice */
                               MURGE_ASSEMBLY_OVW,  /* What to do if an entry appears on 2 processors */
                               MURGE_ASSEMBLY_FOOL, /* Do we respect the solver distribution */
                               MURGE_BOOLEAN_FALSE); CHKERRQ(ierr);

    {
      COEF * vals = malloc(dof*dof*sizeof(COEF));
      memset(vals, 0, dof*dof*sizeof(COEF));
      for (k = 0; k < dof;k++)
	vals[k+k*dof] = 1;
      for (i=0; i<N; i++) {
        if (rows[i]) {
          x = h*(i % (m+1));
          y = h*(i/(m+1));
          val = 1;
          for (k = 0; k < dof; k++) {
            b[i*dof+k] = val;
          }
          ZMURGE_AssemblySetNodeValues(id, i, i, vals);
          /* fprintf(stdout, "A[%d, %d] = %g\n", i, i, *vals); */
        }
      }
      free(vals);
    }

    sprintf(out_string, "RHS assembly (fact %d)", fact);
    STOP_TIMER(out_string, t1);
    START_TIMER(t1);
    trace_task++;
    /* { */
    /*   COEF * vals = malloc(dof*dof*sizeof(COEF)); */
    /*   memset(vals, 0, dof*dof*sizeof(COEF)); */
    /*   for (k = 0; k < dof;k++) */
    /*     vals[k+k*dof] = 1; */
    /*   for (i=0; i<8*m; i++) { */
    /*     ZMURGE_AssemblySetNodeValues(id, rows[i], rows[i], vals); */
    /*   } */
    /*   free(vals); */
    /* } */
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
        val = 1;
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
  /* free(rows); */
  fprintf(stdout, "> Check solution without refinement <\n");

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
