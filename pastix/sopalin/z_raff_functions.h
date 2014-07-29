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
/*
 * File: z_raff_functions.h
 *
 * Functions computing operations for reffinement methods
 *
 */

#ifndef Z_RAFF_FUNCTIONS_H
#define Z_RAFF_FUNCTIONS_H

typedef PASTIX_INT z_RAFF_INT;
typedef pastix_complex64_t z_RAFF_FLOAT;

#define MONO_BEGIN(arg) if(solveur.me(arg)==0){
#define MONO_END(arg)   }

#define SYNC_COEF(var) do {                      \
    MONOTHREAD_BEGIN;                            \
    sopalin_data->common_flt[0] = var;           \
    MONOTHREAD_END;                              \
    SYNCHRO_THREAD;                              \
    var = sopalin_data->common_flt[0];           \
    /* To be sure noone uses common_flt[0] */    \
    SYNCHRO_THREAD;                              \
  } while(0)
#define SYNC_REAL(var) do {                      \
    MONOTHREAD_BEGIN;                            \
    sopalin_data->common_dbl[0] = var;           \
    MONOTHREAD_END;                              \
    SYNCHRO_THREAD;                              \
    var = sopalin_data->common_dbl[0];           \
    /* To be sure noone uses common_dbl[0] */    \
    SYNCHRO_THREAD;                              \
  } while(0)
#ifdef SMP_RAFF
#  define MULTITHREAD_BEGIN
#  define MULTITHREAD_END(sync)
#  define NOSMP_SYNC_COEF(var) do {} while(0)
#  define NOSMP_SYNC_REAL(var) do {} while(0)
#else /* SMP_RAFF */
#  define MULTITHREAD_BEGIN if (me == 0) {
#  define MULTITHREAD_END(sync) } if (sync) {SYNCHRO_THREAD;}
#  define NOSMP_SYNC_COEF(var) SYNC_COEF(var)
#  define NOSMP_SYNC_REAL(var) SYNC_REAL(var)
#endif /* SMP_RAFF */

#define SYNCHRO(arg)                                                    \
  do {                                                                  \
    z_Sopalin_Data_t * sopalin_data;                                      \
    sopalin_data = (z_Sopalin_Data_t*)((sopthread_data_t *)arg)->data;    \
    z_SolverMatrix     *datacode     = sopalin_data->datacode;            \
    SYNCHRO_THREAD;                                                     \
  } while(0)

/*** ALLOCATIONS ET SYNCHRONISATIONS ***/

/* Synchronise le vecteur x dans la nb-ieme variable de la structure */
#define z_Pastix_Synchro_Vect API_CALL(z_Pastix_Synchro_Vect)
pastix_complex64_t *z_Pastix_Synchro_Vect(void *, void *, int);

/* Alloue un vecteur de taille size octets */
#define z_Pastix_Malloc API_CALL(z_Pastix_Malloc)
void *z_Pastix_Malloc(void *, size_t );

/* Libere un vecteur */
#define z_Pastix_Free API_CALL(z_Pastix_Free)
void z_Pastix_Free(void *, void *);


/*** GESTION DE L'INTERFACE ***/

/* Affichage à chaque itération et communication de certaines informations à la structure */
#define z_Pastix_Verbose API_CALL(z_Pastix_Verbose)
void z_Pastix_Verbose(void *, double, double, double, PASTIX_INT);

/* Affichage final */
#define z_Pastix_End API_CALL(z_Pastix_End)
void z_Pastix_End(void*, pastix_complex64_t, PASTIX_INT, double, pastix_complex64_t *);

/* Vecteur solution X */
#define Pastix_X API_CALL(Pastix_X)
void Pastix_X(void *, pastix_complex64_t *);

/* Taille d'un vecteur */
#define z_Pastix_n API_CALL(z_Pastix_n)
PASTIX_INT z_Pastix_n(void *);

/* Nombre de second membres */
#define z_Pastix_m API_CALL(z_Pastix_m)
PASTIX_INT z_Pastix_m(void *);

/* Second membre */
#define z_Pastix_B API_CALL(z_Pastix_B)
void z_Pastix_B(void *, pastix_complex64_t *);

/* Epsilon */
#define z_Pastix_Eps API_CALL(z_Pastix_Eps)
pastix_complex64_t z_Pastix_Eps(void *);

/* Itermax */
#define z_Pastix_Itermax API_CALL(z_Pastix_Itermax)
PASTIX_INT z_Pastix_Itermax(void *);


/* Itermax */
#define z_Pastix_Krylov_Space API_CALL(z_Pastix_Krylov_Space)
PASTIX_INT z_Pastix_Krylov_Space(void *);

/*** OPERATIONS DE BASE ***/
/* Multiplication pour plusieurs second membres */
#define z_Pastix_Mult API_CALL(z_Pastix_Mult)
void z_Pastix_Mult(void *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *, int);

/* Division pour plusieurs second membres */
#define z_Pastix_Div API_CALL(z_Pastix_Div)
void z_Pastix_Div(void *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *, int);

/* Calcul de la norme de frobenius */
#define z_Pastix_Norm2 API_CALL(z_Pastix_Norm2)
pastix_complex64_t z_Pastix_Norm2(void *, pastix_complex64_t *);

/* Copie d'un vecteur */
#define z_Pastix_Copy API_CALL(z_Pastix_Copy)
void z_Pastix_Copy(void *, pastix_complex64_t *, pastix_complex64_t *, int);

/* Application du préconditionneur */
#define z_Pastix_Precond API_CALL(z_Pastix_Precond)
void z_Pastix_Precond(void *, pastix_complex64_t *, pastix_complex64_t *, int);

/* Calcul de alpha * x */
#define z_Pastix_Scal API_CALL(z_Pastix_Scal)
void z_Pastix_Scal(void *, pastix_complex64_t, pastix_complex64_t *, int);

/* Calcul du produit scalaire */
#define z_Pastix_Dotc API_CALL(z_Pastix_Dotc)
void z_Pastix_Dotc(void *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *, int);

#define z_Pastix_Dotc_Gmres API_CALL(z_Pastix_Dotc_Gmres)
void z_Pastix_Dotc_Gmres(void *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *, int);

/* Produit matrice vecteur */
#define z_Pastix_Ax API_CALL(z_Pastix_Ax)
void z_Pastix_Ax(void *, pastix_complex64_t *, pastix_complex64_t *);


/*** A MODIFIER! ***/
#define z_Pastix_bMAx API_CALL(z_Pastix_bMAx)
void z_Pastix_bMAx(void *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *);

#define z_Pastix_BYPX API_CALL(z_Pastix_BYPX)
void z_Pastix_BYPX(void *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *, int);

#define z_Pastix_AXPY API_CALL(z_Pastix_AXPY)
void z_Pastix_AXPY(void *, double, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *, int);

#define z_Pastix_me API_CALL(z_Pastix_me)
PASTIX_INT z_Pastix_me(void *);

struct z_solver
{
  /*** ALLOCATIONS ET SYNCHRONISATIONS ***/
  pastix_complex64_t* (* Synchro)(void *, void *, int);
  void* (* Malloc)(void*, size_t);
  void (* Free)(void*, void*);

  /*** GESTION DE L'INTERFACE ***/
  void (* Verbose)(void *, double, double, double, PASTIX_INT);
  void (* End)(void* , pastix_complex64_t, PASTIX_INT, double, pastix_complex64_t*);
  void (* X)(void *, pastix_complex64_t*);
  PASTIX_INT (* N)(void *);
  void (* B)(void *, pastix_complex64_t*);
  pastix_complex64_t (* Eps)(void *);
  PASTIX_INT (* Itermax)(void *);
  PASTIX_INT (* Krylov_Space)(void *);
  PASTIX_INT (* me)(void *);


  /*** OPERATIONS DE BASE ***/
  void (* Mult)(void *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *, int);
  void (* Div)(void *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *, int);
  void (* Dotc_Gmres)(void *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *, int);

  pastix_complex64_t (* Norm)(void* , pastix_complex64_t *);
  void (* Copy)(void *, pastix_complex64_t *, pastix_complex64_t *, int);
  void (* Precond)(void *, pastix_complex64_t *, pastix_complex64_t *, int);

  void (* Scal)(void *, pastix_complex64_t, pastix_complex64_t *, int);
  void (* Dotc)(void *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *, int);
  void (* Ax)(void *, pastix_complex64_t *, pastix_complex64_t *);

  void (* bMAx)(void *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *);
  void (* BYPX)(void *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *, int);
  void (* AXPY)(void *, double, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *, int);
};

#define z_Pastix_Solveur API_CALL(z_Pastix_Solveur)
void z_Pastix_Solveur(struct z_solver *);

/*
 ** Section: Function creating threads
 */
/*
 Function: method)

 Launch sopaparam->nbthrdcomm threads which will compute
 <method_smp)>.

 Parameters:
 datacode  - PaStiX <z_SolverMatrix> structure.
 sopaparam - <z_SopalinParam> parameters structure.
 */
#define z_raff_thread API_CALL(z_raff_thread)
void z_raff_thread(z_SolverMatrix *, z_SopalinParam *, void*(*)(void *));

#endif
