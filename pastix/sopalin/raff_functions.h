/*
 * File: raff_functions.h
 *
 * Functions computing operations for reffinement methods
 *
 */

#ifndef RAFF_FUNCTIONS_H
#define RAFF_FUNCTIONS_H

typedef PASTIX_INT RAFF_INT;
typedef pastix_float_t RAFF_FLOAT;

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
    Sopalin_Data_t * sopalin_data;                                      \
    sopalin_data = (Sopalin_Data_t*)((sopthread_data_t *)arg)->data;    \
    SolverMatrix     *datacode     = sopalin_data->datacode;            \
    SYNCHRO_THREAD;                                                     \
  } while(0)

/*** ALLOCATIONS ET SYNCHRONISATIONS ***/

/* Synchronise le vecteur x dans la nb-ieme variable de la structure */
#define Pastix_Synchro_Vect API_CALL(Pastix_Synchro_Vect)
pastix_float_t *Pastix_Synchro_Vect(void *, void *, int);

/* Alloue un vecteur de taille size octets */
#define Pastix_Malloc API_CALL(Pastix_Malloc)
void *Pastix_Malloc(void *, size_t );

/* Libere un vecteur */
#define Pastix_Free API_CALL(Pastix_Free)
void Pastix_Free(void *, void *);


/*** GESTION DE L'INTERFACE ***/

/* Affichage à chaque itération et communication de certaines informations à la structure */
#define Pastix_Verbose API_CALL(Pastix_Verbose)
void Pastix_Verbose(void *, double, double, double, PASTIX_INT);

/* Affichage final */
#define Pastix_End API_CALL(Pastix_End)
void Pastix_End(void*, pastix_float_t, PASTIX_INT, double, pastix_float_t *);

/* Vecteur solution X */
#define Pastix_X API_CALL(Pastix_X)
void Pastix_X(void *, pastix_float_t *);

/* Taille d'un vecteur */
#define Pastix_n API_CALL(Pastix_n)
PASTIX_INT Pastix_n(void *);

/* Nombre de second membres */
#define Pastix_m API_CALL(Pastix_m)
PASTIX_INT Pastix_m(void *);

/* Second membre */
#define Pastix_B API_CALL(Pastix_B)
void Pastix_B(void *, pastix_float_t *);

/* Epsilon */
#define Pastix_Eps API_CALL(Pastix_Eps)
pastix_float_t Pastix_Eps(void *);

/* Itermax */
#define Pastix_Itermax API_CALL(Pastix_Itermax)
PASTIX_INT Pastix_Itermax(void *);


/* Itermax */
#define Pastix_Krylov_Space API_CALL(Pastix_Krylov_Space)
PASTIX_INT Pastix_Krylov_Space(void *);

/*** OPERATIONS DE BASE ***/
/* Multiplication pour plusieurs second membres */
#define Pastix_Mult API_CALL(Pastix_Mult)
void Pastix_Mult(void *, pastix_float_t *, pastix_float_t *, pastix_float_t *, int);

/* Division pour plusieurs second membres */
#define Pastix_Div API_CALL(Pastix_Div)
void Pastix_Div(void *, pastix_float_t *, pastix_float_t *, pastix_float_t *, int);

/* Calcul de la norme de frobenius */
#define Pastix_Norm2 API_CALL(Pastix_Norm2)
pastix_float_t Pastix_Norm2(void *, pastix_float_t *);

/* Copie d'un vecteur */
#define Pastix_Copy API_CALL(Pastix_Copy)
void Pastix_Copy(void *, pastix_float_t *, pastix_float_t *, int);

/* Application du préconditionneur */
#define Pastix_Precond API_CALL(Pastix_Precond)
void Pastix_Precond(void *, pastix_float_t *, pastix_float_t *, int);

/* Calcul de alpha * x */
#define Pastix_Scal API_CALL(Pastix_Scal)
void Pastix_Scal(void *, pastix_float_t, pastix_float_t *, int);

/* Calcul du produit scalaire */
#define Pastix_Dotc API_CALL(Pastix_Dotc)
void Pastix_Dotc(void *, pastix_float_t *, pastix_float_t *, pastix_float_t *, int);

#define Pastix_Dotc_Gmres API_CALL(Pastix_Dotc_Gmres)
void Pastix_Dotc_Gmres(void *, pastix_float_t *, pastix_float_t *, pastix_float_t *, int);

/* Produit matrice vecteur */
#define Pastix_Ax API_CALL(Pastix_Ax)
void Pastix_Ax(void *, pastix_float_t *, pastix_float_t *);


/*** A MODIFIER! ***/
#define Pastix_bMAx API_CALL(Pastix_bMAx)
void Pastix_bMAx(void *, pastix_float_t *, pastix_float_t *, pastix_float_t *);

#define Pastix_BYPX API_CALL(Pastix_BYPX)
void Pastix_BYPX(void *, pastix_float_t *, pastix_float_t *, pastix_float_t *, int);

#define Pastix_AXPY API_CALL(Pastix_AXPY)
void Pastix_AXPY(void *, double, pastix_float_t *, pastix_float_t *, pastix_float_t *, int);

#define Pastix_me API_CALL(Pastix_me)
PASTIX_INT Pastix_me(void *);

struct solver
{
  /*** ALLOCATIONS ET SYNCHRONISATIONS ***/
  pastix_float_t* (* Synchro)(void *, void *, int);
  void* (* Malloc)(void*, size_t);
  void (* Free)(void*, void*);

  /*** GESTION DE L'INTERFACE ***/
  void (* Verbose)(void *, double, double, double, PASTIX_INT);
  void (* End)(void* , pastix_float_t, PASTIX_INT, double, pastix_float_t*);
  void (* X)(void *, pastix_float_t*);
  PASTIX_INT (* N)(void *);
  void (* B)(void *, pastix_float_t*);
  pastix_float_t (* Eps)(void *);
  PASTIX_INT (* Itermax)(void *);
  PASTIX_INT (* Krylov_Space)(void *);
  PASTIX_INT (* me)(void *);


  /*** OPERATIONS DE BASE ***/
  void (* Mult)(void *, pastix_float_t *, pastix_float_t *, pastix_float_t *, int);
  void (* Div)(void *, pastix_float_t *, pastix_float_t *, pastix_float_t *, int);
  void (* Dotc_Gmres)(void *, pastix_float_t *, pastix_float_t *, pastix_float_t *, int);

  pastix_float_t (* Norm)(void* , pastix_float_t *);
  void (* Copy)(void *, pastix_float_t *, pastix_float_t *, int);
  void (* Precond)(void *, pastix_float_t *, pastix_float_t *, int);

  void (* Scal)(void *, pastix_float_t, pastix_float_t *, int);
  void (* Dotc)(void *, pastix_float_t *, pastix_float_t *, pastix_float_t *, int);
  void (* Ax)(void *, pastix_float_t *, pastix_float_t *);

  void (* bMAx)(void *, pastix_float_t *, pastix_float_t *, pastix_float_t *);
  void (* BYPX)(void *, pastix_float_t *, pastix_float_t *, pastix_float_t *, int);
  void (* AXPY)(void *, double, pastix_float_t *, pastix_float_t *, pastix_float_t *, int);
};

#define Pastix_Solveur API_CALL(Pastix_Solveur)
void Pastix_Solveur(struct solver *);

/*
 ** Section: Function creating threads
 */
/*
 Function: method)

 Launch sopaparam->nbthrdcomm threads which will compute
 <method_smp)>.

 Parameters:
 datacode  - PaStiX <SolverMatrix> structure.
 sopaparam - <SopalinParam> parameters structure.
 */
#define raff_thread API_CALL(raff_thread)
void raff_thread(SolverMatrix *, SopalinParam *, void*(*)(void *));

#endif
