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
 * File: zmurge.c
 *
 * This file implements <Murge> interface.
 *
 * About: Authors
 *   Mathieu Faverge - faverge@labri.fr
 *   Xavier Lacoste  - xavier.lacoste@inria.fr
 */
#ifdef _WIN32
#  define MURGE_DLL_EXPORT
/* Required for CreateSymbolicLink() call */
#  define _WIN32_WINNT 0x0600
#endif

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#ifndef FORCE_NOSMP
#  include <pthread.h>
#endif
#ifdef FORCE_NOMPI
#  include "nompi.h"
#else /* not FORCE_NOMPI */
#  include <mpi.h>
#endif /* not FORCE_NOMPI */
#ifdef WITH_SEM_BARRIER
#  include <semaphore.h>
#endif
#include "common.h"
#include "z_tools.h"
#include "sopalin_define.h"

#ifdef WITH_SCOTCH
#  ifdef    DISTRIBUTED
#    include "ptscotch.h"
#  else
#    include "scotch.h"
#  endif /* DISTRIBUTED */
#endif /* WITH_SCOTCH */

#include "z_ftgt.h"
#include "symbol.h"
#include "z_csc.h"
#include "z_updown.h"
#include "queue.h"
#include "bulles.h"
#include "z_solver.h"
#include "order.h"
#include "sopalin_thread.h"
#include "stack.h"
#include "z_sopalin3d.h"
#include "z_pastixdata.h"
#include "z_pastix.h"
#include "z_pastix_internal.h"
#include "z_cscd_utils.h"
#include "z_cscd_utils_intern.h"
#include "z_csc_intern_compute.h"
#include "z_csc_intern_updown.h"
#include "z_sopalin_init.h"
#include "z_sopalin_compute.h"

#if (defined PRECISION_z || defined PRECISION_d)
#  define PREC_DOUBLE
#endif
#if (defined PRECISION_z || defined PRECISION_c)
#  define TYPE_COMPLEX
#endif
#define OVERWRITE_MURGE_CALLS
#include "zmurge.h"
#include "zmurge_pastix.h"
#include "zmurge_defines.h"

// ALWAYS_INLINE ---------------------------------------------//
// Macro to use in place of 'inline' to force a function to be inline
#define ALWAYS_INLINE
#if !defined(ALWAYS_INLINE)
#  if defined(_MSC_VER)
#    define ALWAYS_INLINE __forceinline
#  elif defined(__GNUC__) && __GNUC__ > 3
// Clang also defines __GNUC__ (as 4)
#    define ALWAYS_INLINE __attribute__ ((__always_inline__))
#  else
#    define ALWAYS_INLINE
#  endif
#endif
/* #define ALWAYS_INLINE __attribute__((always_inline)) */

/* internal function declarations */
static inline
INTS ZMURGE_GraphSetEdge_ (INTS id, INTS ROW, INTS COL) ALWAYS_INLINE;

static inline
INTS ZMURGE_AssemblySetValue_    (INTS id, INTS ROW, INTS COL, COEF value) ALWAYS_INLINE;
static inline
INTS ZMURGE_AssemblySetNodeValues_ (INTS id, INTS ROW, INTS COL, COEF *values) ALWAYS_INLINE;
static inline
INTS ZMURGE_AssemblySetBlockValues_(INTS id, INTS nROW, INTS *ROWlist,
                                   INTS nCOL, INTS *COLlist, COEF *values) ALWAYS_INLINE;

static inline
INTS ZMURGE_GraphSetEdge_ (INTS id, INTS ROW, INTS COL) ALWAYS_INLINE;


/******************************************************************************/
/***                           Section: Structures                          ***/
/******************************************************************************/

/*
 * Structure: ijv_
 *
 * Structure to represente coefficients.
 *
 * Contains:
 *   i     - row
 *   j     - column
 *   v     - pointer to the value array (can be several degree of freedom)
 *   owner - process which own the coefficient.
 */
struct ijv_ {
  pastix_int_t    i;
  pastix_int_t    j;
  pastix_complex64_t* v;
#ifdef MURGE_FOLLOW_INPUT_ORDER
  pastix_int_t    idx;
#endif
  int    owner;
};

/*
 * Typedef: ijv_t
 *
 * Alias to structure <ijv_>.
 */
typedef struct ijv_ ijv_t;

/*
 * struct: z_murge_seq_t
 *
 * Structure used to store assembly sequence.
 *
 * Contains:
 *   indexes              - Order sequences of indexes to be set.
 *   coefnbr              - Number of entries in the sequence.
 *   recv_nbr             - Number of entries to receive from each processor.
 *   recv_indexes         - Indexes of entries which will be received.
 *   mode                 - Local entries or communicating mode
 *   op_local_entries     - Operation to perform when a coefficient appear twice
 *   op_dist_entries      - Operation to perform when a coefficient appear
 *                            twice, given by two processors.
 *   ijv_size             - size of the required array to store
 *                            not local entries.
 *   nodes                - 0 entries are entered value by value,
 *                          1 entries are entries node by node.
 *   next                 - Next entry in the list of sequences.
 */
typedef struct z_murge_seq_ z_murge_seq_t;
struct z_murge_seq_ {
  INTL        * indexes;
  INTL          coefnbr;
  INTL        * recv_nbr;
  INTL       ** recv_indexes;
  INTS          mode;
  INTS          op_local_entries;
  INTS          op_dist_entries;
  INTL          ijv_size;
  INTS          nodes;
  z_murge_seq_t * next;
  pastix_int_t           ID;
};

typedef struct variable_csc_s * variable_csc_t;


/*
 * Typedef: z_murge_data_t
 *
 * alias to structure <murge_data_>.
 */
typedef struct z_murge_data_ z_murge_data_t;

/*
 * struct: z_murge_product_data_
 *
 * thread_id - thread identification number
 * solver    - ponter to <z_murge_data_t> structure.
 * t_prod    - buffer for product computation (shared).
 * ret       - return value.
 */
struct z_murge_product_data_
{
  INTS           thread_id;
  z_murge_data_t * solver;
  COEF         * t_prod;
  INTS           ret;
  INTS           counter;
  INTS           first;
  INTS           last;
};

/*
 * typedef: z_murge_product_data_t
 *
 * alias to struct <z_murge_product_data_>
 */
typedef struct z_murge_product_data_ z_murge_product_data_t;

struct z_murge_thread_data_ {
    z_murge_product_data_t * pdata;
    z_murge_data_t         * solver;
};

/*
 * typedef: z_murge_thread_data_t
 *
 * alias to struct <z_murge_thread_data_>
 */
typedef struct z_murge_thread_data_ z_murge_thread_data_t;

typedef enum murge_thread_state_ {
  MURGE_THREAD_INIT,
  MURGE_THREAD_WAIT,
  MURGE_THREAD_PRODUCT,
  MURGE_THREAD_END
} murge_thread_state_t;


/*
 * struct: z_murge_data_t
 *
 * Structure used to store murge data
 *
 * Contains:
 *   pastix_data    - Pointer to the <z_pastix_data_t> associated to
 *                    the solver instance
 *   n              - Number of local column indicated by murge user
 *   N              - Number of global column indicated by murge user
 *   colptr         - Colptr in murge's user CSCd
 *   rows           - Rows in murge's user CSCd
 *   values         - Values in murge's user CSCd
 *   l2g            - Local to global column number in murge's user CSCd
 *   g2l            - Global to local column number in murge's user CSCd
 *   perm           - Permtab for murge's user
 *   b              - Right-hand-side member(s) given by murge's user
 *   nrhs           - Number of right-hand-side member(s) given by murge's user
 *   cnt            - Iterator for number of entered edges
 *   edgenbr        - Number of edges
 *   state          - State of the solver
 *   mode           - Local entries or communicating mode
 *   op             - Operation to perform when a coefficient appear twice
 *   op2            - Operation to perform when a coefficient appear twice,
 *                    given by two processors.
 *   sym            - Indicate if we have to check that the matrix is symmetric
 *   sequences      - List of saved sequences.
 *   malloc_size    - Current allocated memory.
 *   malloc_maxsize - Maximum allocated memory.
 *   threadnbr      - Number of thread launched.
 *   threads        - Array of pthread_t.
 *   threads_data   - data associated to each thread.
 *   thread_state   - flag to control the threads.
 *   barrier        - barrier used for the threads.
 */
struct z_murge_data_ {
  z_pastix_data_t          *pastix_data;
  pastix_int_t              n;
  pastix_int_t              N;
  pastix_int_t             *colptr;
  pastix_int_t             *rows;
  pastix_complex64_t           *values;
  pastix_int_t             *l2g;
  pastix_int_t             *g2l;
  pastix_int_t             *perm;
#ifdef CENTRALISED
  pastix_int_t             *invp;
#endif
  pastix_complex64_t           *b;
  pastix_int_t              nrhs;
  variable_csc_t         vcsc;
#ifdef MURGE_THREADSAFE
  pthread_mutex_t         mutex_tmpmatrix;
#endif
  pastix_int_t              cnt;
  pastix_int_t              cnt_zero;
  pastix_int_t              cnt_node;
  pastix_int_t              edgenbr;
  pastix_int_t              coefnbr;
  pastix_int_t              nodenbr;
  int                     state;
  int                     mode;
  int                     op;
  int                     op2;
  int                     sym;
  int                     dynamic;
  z_murge_seq_t            *sequences;
  pastix_int_t              seq_ID;
  int                     ndump;
  char                   *dropmask;
  char                   *droprows;
  char                   *dropcols;
  int64_t                 malloc_size;
  int64_t                 malloc_maxsize;
  int                     threadnbr;
  pthread_t              *threads;
  z_murge_thread_data_t  *threads_data;
  murge_thread_state_t    threads_state;
  sopthread_barrier_t     barrier;
  pthread_mutex_t         mutex_state;
  pthread_cond_t          cond_state;
};

/******************************************************************************/
/***                           Section: Global variables                    ***/
/******************************************************************************/


/*
 * Variables: Global variables
 *
 * idnbr   - Number of solvers instances.
 * solvers - Murge solver instances array (<z_murge_data_t>).
 */
INTS           idnbr   = 0;
z_murge_data_t **zmurge_solvers = NULL;
#include "variable_csc.c"


/******************************************************************************/
/***                           Section: Functions                           ***/
/******************************************************************************/


/*******************************************************************************
 * Group: Auxilary functions
 */

#ifdef CENTRALISED
#  define ALLOC_INVP                                                      \
  if (NULL == zmurge_solvers[id]->invp)                                        \
    MURGE_MEMALLOC(zmurge_solvers[id]->invp, zmurge_solvers[id]->n, pastix_int_t);
#  define INVP zmurge_solvers[id]->invp
#else
#  define ALLOC_INVP {};
#  define INVP NULL
#endif


#define z_product_thread PASTIX_PREFIX_F(z_product_thread)
void* z_product_thread(void * data) {
  INTS first, last;
  INTS gfirst, glast;
  INTS nnz, nnz_per_thread;
  INTS itercol;
  INTS iterrows;
  INTS dof, baseval;
  INTS row;
  COEF * mat = NULL;
  z_murge_product_data_t * pdata = (z_murge_product_data_t *)data;
  INTS id = &(pdata->solver) - zmurge_solvers;
  INTS threadnbr = pdata->solver->pastix_data->iparm[IPARM_THREAD_NBR];

  baseval = pdata->solver->pastix_data->iparm[IPARM_BASEVAL];
  dof = pdata->solver->pastix_data->iparm[IPARM_DOF_NBR];
  pdata->ret = MURGE_SUCCESS;
  if (pdata->t_prod == NULL)
      MURGE_MEMALLOC_RET(pdata->solver,
                         pdata->t_prod, (pdata->last-pdata->first)*dof, COEF,
                         pdata->ret);
  for (iterrows = 0; iterrows < pdata->last-pdata->first; iterrows++) {
    INTS i;
    for (i = 0; i < dof; i++) {
      pdata->t_prod[(iterrows)*dof+i] = 0;
  }
  }
  SYNCHRO_X_THREAD(threadnbr, pdata->solver->barrier);

  for (itercol = 0; itercol < pdata->solver->n; itercol++) {
    for (iterrows = pdata->solver->colptr[itercol]-baseval;
         iterrows < pdata->solver->colptr[itercol+1]-baseval;
           iterrows++) {
        row = pdata->solver->rows[iterrows]-baseval;
      if (row < pdata->first) continue;
      if (row >= pdata->last) break;
        mat = &(pdata->solver->values[iterrows*dof*dof]);
        SOPALIN_GEMV("N", dof, dof, 1.0, mat, dof,
                     &(pdata->solver->b[itercol*dof]), 1, 1.0,
                   &(pdata->t_prod[(row-pdata->first)*dof]), 1);
    }
  }

  SYNCHRO_X_THREAD(threadnbr, pdata->solver->barrier);
  if (pdata->thread_id == 0) {
    pthread_mutex_lock(&(pdata->solver->mutex_state));
    pdata->solver->threads_state = MURGE_THREAD_WAIT;
    pthread_mutex_unlock(&(pdata->solver->mutex_state));
    pthread_cond_broadcast(&(pdata->solver->cond_state));
  }
  SYNCHRO_X_THREAD(threadnbr, pdata->solver->barrier);
  return 0;
}


void* z_control_thread(void * data) {
  z_murge_thread_data_t * tdata= (z_murge_thread_data_t*)data;
  INTS threadnbr = tdata->solver->pastix_data->iparm[IPARM_THREAD_NBR];
  if (tdata->pdata->thread_id == 0) {
    pthread_mutex_lock(&(tdata->solver->mutex_state));
    tdata->solver->threads_state = MURGE_THREAD_WAIT;
    pthread_mutex_unlock(&(tdata->solver->mutex_state));
  }
  SYNCHRO_X_THREAD(threadnbr, tdata->solver->barrier);
  if (tdata->pdata->thread_id == 0) {
    pthread_cond_broadcast(&(tdata->solver->cond_state));
  }
  pthread_mutex_lock(&(tdata->solver->mutex_state));
  while (tdata->solver->threads_state != MURGE_THREAD_END) {
    switch(tdata->solver->threads_state) {
    case MURGE_THREAD_END:
      pthread_mutex_unlock(&(tdata->solver->mutex_state));
      return NULL;
      break;
    case MURGE_THREAD_WAIT:
      pthread_cond_wait(&(tdata->solver->cond_state),
                        &(tdata->solver->mutex_state));
      break;
    case MURGE_THREAD_PRODUCT:
      pthread_mutex_unlock(&(tdata->solver->mutex_state));
      z_product_thread(tdata->pdata);
      pthread_mutex_lock(&(tdata->solver->mutex_state));
      break;
    default:
      errorPrint("Undefined state : %d", tdata->solver->threads_state);
      break;
    }
  }
  pthread_mutex_unlock(&(tdata->solver->mutex_state));
  return NULL;
}

INTS z_start_threads(INTS id) {
  INTS iter;
  int ret;
  COEF ** all_prods;
  zmurge_solvers[id]->threadnbr = zmurge_solvers[id]->pastix_data->iparm[IPARM_THREAD_NBR];
  pthread_mutex_init(&(zmurge_solvers[id]->mutex_state), NULL);
  pthread_cond_init(&(zmurge_solvers[id]->cond_state), NULL);
  MURGE_MEMALLOC(zmurge_solvers[id]->threads_data,
                 zmurge_solvers[id]->threadnbr,
                 z_murge_thread_data_t);
  MURGE_MEMALLOC(zmurge_solvers[id]->threads,
                 zmurge_solvers[id]->threadnbr,
                 pthread_t);
  for (iter = 0; iter < zmurge_solvers[id]->threadnbr; iter++) {
    MURGE_MEMALLOC(zmurge_solvers[id]->threads_data[iter].pdata,
                   zmurge_solvers[id]->threadnbr,
                   z_murge_product_data_t);
    zmurge_solvers[id]->threads_data[iter].solver = zmurge_solvers[id];
    {
      /* Initialize pdata */
      z_murge_product_data_t * pdata       = zmurge_solvers[id]->threads_data[iter].pdata;
      INTS                   N           = zmurge_solvers[id]->N;
      INTS                   N_perthread = N/zmurge_solvers[id]->threadnbr;
      INTS                   size        = N_perthread;

      pdata->thread_id = iter;
      pdata->solver = zmurge_solvers[id];
      pdata->t_prod = NULL;
      pdata->counter = 0;
      pdata->first = pdata->thread_id*N_perthread;
      if (pdata->thread_id == zmurge_solvers[id]->threadnbr-1)
        size = N-pdata->first;
      pdata->last=pdata->first+size;
    }
  }
  zmurge_solvers[id]->threads_state = MURGE_THREAD_INIT;
  for (iter = 0; iter < zmurge_solvers[id]->threadnbr; iter++) {
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    ret = pthread_create(&(zmurge_solvers[id]->threads[iter]), &attr,
                         z_control_thread, (void *)&(zmurge_solvers[id]->threads_data[iter]));
    if (ret) {errorPrint("thread create."); EXIT(MOD_SOPALIN,THREAD_ERR);}
  }
  pthread_mutex_lock(&(zmurge_solvers[id]->mutex_state));
  while (zmurge_solvers[id]->threads_state != MURGE_THREAD_WAIT) {
    pthread_cond_wait(&(zmurge_solvers[id]->cond_state),
                      &(zmurge_solvers[id]->mutex_state));
  }
  pthread_mutex_unlock(&(zmurge_solvers[id]->mutex_state));
  return MURGE_SUCCESS;
}

INTS z_stop_threads(INTS id) {
  int ret;
  INTS iter;
  pthread_mutex_lock(&(zmurge_solvers[id]->mutex_state));
  zmurge_solvers[id]->threads_state = MURGE_THREAD_END;
  pthread_mutex_unlock(&(zmurge_solvers[id]->mutex_state));
  pthread_cond_broadcast(&(zmurge_solvers[id]->cond_state));
  for (iter = 0; iter < zmurge_solvers[id]->threadnbr; iter++) {
    ret = pthread_join(zmurge_solvers[id]->threads[iter],(void**)NULL);
    if (ret) {errorPrint("thread join."); EXIT(MOD_SOPALIN,THREAD_ERR);}
    if (zmurge_solvers[id]->threads_data[iter].pdata &&
        zmurge_solvers[id]->threads_data[iter].pdata->t_prod)
    MURGE_FREE(zmurge_solvers[id]->threads_data[iter].pdata->t_prod);
    MURGE_FREE(zmurge_solvers[id]->threads_data[iter].pdata);
  }
  MURGE_FREE(zmurge_solvers[id]->threads);
  MURGE_FREE(zmurge_solvers[id]->threads_data);
  pthread_mutex_destroy(&(zmurge_solvers[id]->mutex_state));
  pthread_cond_destroy(&(zmurge_solvers[id]->cond_state));
  zmurge_solvers[id]->threadnbr = 0;
  return MURGE_SUCCESS;
}

static inline
int murge_redistribute_data(INTS id) {
  z_pastix_data_t *pastix_data = zmurge_solvers[id]->pastix_data;
  PASTIX_INT    *iparm       = pastix_data->iparm;
  PASTIX_INT     newn;
  PASTIX_INT    *newl2g      = NULL;
  PASTIX_INT    *newcolptr   = NULL;
  PASTIX_INT    *newrows     = NULL;
  COEF          *newvals     = NULL, **newvalsptr = NULL;
  COEF          *newb        = NULL, **newbptr    = NULL;
#ifdef MURGE_TIME
  Clock clock;
#endif
  CLOCK_INIT;
  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_VALUES_OK)) {
    newvalsptr = &newvals;
  }
  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_RHS_OK)) {
    newbptr = &newb;
  }

  newn = z_pastix_getLocalNodeNbr(&(pastix_data));
  MURGE_FREE(zmurge_solvers[id]->perm);
  MURGE_MEMALLOC(zmurge_solvers[id]->perm, newn, PASTIX_INT);
  memset(zmurge_solvers[id]->perm, 0, newn*sizeof(PASTIX_INT));
  MURGE_MEMALLOC(newl2g,            newn, PASTIX_INT);
  z_pastix_getLocalNodeLst(&(zmurge_solvers[id]->pastix_data),
                           newl2g);
  CLOCK_PRINT("pastix_getLocalNodeLst");

  z_cscd_redispatch_int(zmurge_solvers[id]->n,
                      zmurge_solvers[id]->colptr,
                      zmurge_solvers[id]->rows,
                      zmurge_solvers[id]->values,
                      zmurge_solvers[id]->b, zmurge_solvers[id]->nrhs,  zmurge_solvers[id]->l2g,
                      newn,
                      &newcolptr,
                      &newrows,
                      newvalsptr,
                      newbptr, newl2g, API_YES,
                      zmurge_solvers[id]->pastix_data->pastix_comm,
                      iparm[IPARM_DOF_NBR]);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(newcolptr), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(newrows), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(newvals), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(newb), char);
  CLOCK_PRINT("cscd_redispatch_int");

  MURGE_FREE(zmurge_solvers[id]->l2g);
  if (zmurge_solvers[id]->g2l)
    MURGE_FREE(zmurge_solvers[id]->g2l);
  MURGE_FREE(zmurge_solvers[id]->colptr);
  MURGE_FREE(zmurge_solvers[id]->rows);
  if (zmurge_solvers[id]->values)
    MURGE_FREE(zmurge_solvers[id]->values);
  if (zmurge_solvers[id]->b)
    MURGE_FREE(zmurge_solvers[id]->b);
  zmurge_solvers[id]->n      = newn;
  zmurge_solvers[id]->l2g    = newl2g;
  zmurge_solvers[id]->colptr = newcolptr;
  zmurge_solvers[id]->rows   = newrows;
  zmurge_solvers[id]->values = newvals;
  zmurge_solvers[id]->b      = newb;
  return MURGE_SUCCESS;
}


/*
 * Function: check_preprocessing
 *
 * Checks if preprocessing (blend) has been called.
 *
 * If it hasn't, it will allocate permutation tabular
 * and call preprocessing step.
 *
 * After calling preprocessing, it will set local number of column
 * and local to global column number tabular to their new values.
 *
 * Colptr and rows will be destroyed because it is obsolete,
 * and state will be set to indicate that preprocessing has been performed.
 *
 * Parameters:
 *   id - Solver instance ID we want to check
 *
 * Returns:
 *   MURGE_ERR_ALLOCATE - If any allocation error occurs.
 */
static inline
int check_preprocessing(int id) {
  z_pastix_data_t   *pastix_data = zmurge_solvers[id]->pastix_data;
  pastix_int_t             *iparm       = pastix_data->iparm;
#ifdef MURGE_TIME
  Clock            clock;
#endif

  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK) &&
      !(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_BLEND_OK))) {
    /* Si il n'est pas fait on effectue le pretraitement */
    if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_SYMB_OK)))
      iparm[IPARM_START_TASK]  = API_TASK_ORDERING;
    else
      iparm[IPARM_START_TASK]  = API_TASK_BLEND;
    iparm[IPARM_END_TASK]    = API_TASK_BLEND;

    CLOCK_INIT;

    if (iparm[IPARM_MATRIX_VERIFICATION] == API_YES) {
      pastix_int_t err;
#ifdef MEMORY_USAGE
      size_t new_size, old_size = PTR_MEMSIZE(zmurge_solvers[id]->rows);
#endif
      if (NO_ERR !=
          ((err =
            z_pastix_checkMatrix_int(pastix_data->pastix_comm,
                                     iparm[IPARM_VERBOSE],
                                     iparm[IPARM_SYM],
                                     API_YES,
                                     zmurge_solvers[id]->n,
                                     &(zmurge_solvers[id]->colptr),
                                     &(zmurge_solvers[id]->rows),
                                     NULL,
                                     &(zmurge_solvers[id]->l2g),
                                     iparm[IPARM_DOF_NBR],
                                     API_YES)))) {
        errorPrint("z_pastix_checkMatrix : err %ld\n", (long)err);
        return MURGE_ERR_PARAMETER;
      }
#ifdef MEMORY_USAGE
      new_size = PTR_MEMSIZE(zmurge_solvers[id]->rows);
      MURGE_TRACE_MALLOC(new_size-old_size, char);
#endif
    }

    z_pastix_welcome_print(zmurge_solvers[id]->pastix_data,
                           zmurge_solvers[id]->colptr,
                           zmurge_solvers[id]->n);

    if (zmurge_solvers[id]->perm == NULL) {
      MURGE_MEMALLOC(zmurge_solvers[id]->perm, zmurge_solvers[id]->n, pastix_int_t);
      memset(zmurge_solvers[id]->perm, 0, zmurge_solvers[id]->n*sizeof(pastix_int_t));
    }

    ALLOC_INVP;

    DPASTIX(&(pastix_data),
            pastix_data->pastix_comm,
            zmurge_solvers[id]->n,
            zmurge_solvers[id]->colptr,
            zmurge_solvers[id]->rows,
            zmurge_solvers[id]->values,
            zmurge_solvers[id]->l2g,
            zmurge_solvers[id]->perm,
            INVP,
            zmurge_solvers[id]->b,
            zmurge_solvers[id]->nrhs,
            iparm,
            pastix_data->dparm);
    CLOCK_PRINT("Preprocessing in PaStiX");
#ifdef MURGE_INSERT_DIRECTLY
    /* We build the new CSC which will receive the coefficients */
    murge_redistribute_data(id);
#else /* not MURGE_INSERT_DIRECTLY */
    /* On corrige n et l2g */

    zmurge_solvers[id]->n = z_pastix_getLocalNodeNbr(&(pastix_data));
    MURGE_FREE(zmurge_solvers[id]->l2g);
    MURGE_FREE(zmurge_solvers[id]->g2l);
    MURGE_FREE(zmurge_solvers[id]->perm);
    MURGE_FREE(zmurge_solvers[id]->colptr);
    MURGE_FREE(zmurge_solvers[id]->rows);
    MURGE_MEMALLOC(zmurge_solvers[id]->perm, zmurge_solvers[id]->n, pastix_int_t);
    memset(zmurge_solvers[id]->perm, 0, zmurge_solvers[id]->n*sizeof(pastix_int_t));
    MURGE_MEMALLOC(zmurge_solvers[id]->l2g, zmurge_solvers[id]->n, pastix_int_t);
    z_pastix_getLocalNodeLst(&(zmurge_solvers[id]->pastix_data),
                           zmurge_solvers[id]->l2g);
    zmurge_solvers[id]->N=-1;

#endif /* not MURGE_INSERT_DIRECTLY */
    /*
     * Building global to local column number array
     */
    z_cscd_build_g2l(zmurge_solvers[id]->n,
                     zmurge_solvers[id]->l2g,
                     zmurge_solvers[id]->pastix_data->pastix_comm,
                     &zmurge_solvers[id]->N,
                     &zmurge_solvers[id]->g2l);
    MURGE_TRACE_MALLOC(PTR_MEMSIZE(zmurge_solvers[id]->g2l), char);
    MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_BLEND_OK);
    MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_SYMB_OK);
    CLOCK_PRINT("z_cscd_build_g2l");
  }
  return MURGE_SUCCESS;
}

static inline
int check_fact(INTS id) {
  z_pastix_data_t   *pastix_data = zmurge_solvers[id]->pastix_data;
  pastix_int_t      *iparm       = pastix_data->iparm;

  /* Check that values have been filled */
  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_VALUES_OK))) {
    errorPrint("Need to set values before.");
    return MURGE_ERR_ORDER;
  }
  /* Check That RHS has been filled */
  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_RHS_OK))) {
    errorPrint("Need to set right-hand-side member before.");
    return MURGE_ERR_ORDER;
  }
  /* Do we have to perform factorization ¨*/
  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_FACTO_OK))) {
    if (iparm[IPARM_MATRIX_VERIFICATION] == API_YES) {
      pastix_int_t err;
#ifdef MEMORY_USAGE
      size_t new_size;
      size_t old_size = PTR_MEMSIZE(zmurge_solvers[id]->rows) +
	PTR_MEMSIZE(zmurge_solvers[id]->values);
#endif
      if (NO_ERR !=
	  ((err =
	    z_pastix_checkMatrix_int(pastix_data->pastix_comm,
				   iparm[IPARM_VERBOSE],
				   iparm[IPARM_SYM],
				   API_YES,
				   zmurge_solvers[id]->n,
				   &(zmurge_solvers[id]->colptr),
				   &(zmurge_solvers[id]->rows),
				   &(zmurge_solvers[id]->values),
				   &(zmurge_solvers[id]->l2g),
				   iparm[IPARM_DOF_NBR],
				   API_YES)))) {
	errorPrint("z_pastix_checkMatrix : err %ld\n", (long)err);
	return MURGE_ERR_PARAMETER;
      }
#ifdef MEMORY_USAGE
      new_size = PTR_MEMSIZE(zmurge_solvers[id]->rows) +
	PTR_MEMSIZE(zmurge_solvers[id]->values);
      MURGE_TRACE_MALLOC(new_size-old_size, char);
#endif
    }
    CLOCK_PRINT("check_fact -- pastix_checkMatrix");
    /* Do we have to perform preprocessing */
    if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_BLEND_OK))) {
      /* Do we have to perform symbolic factorization */
      if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state,
                               MURGE_SYMB_OK)))
        iparm[IPARM_START_TASK]  = API_TASK_ORDERING;
      else
        iparm[IPARM_START_TASK]  = API_TASK_BLEND;
      if (zmurge_solvers[id]->perm == NULL) {
        MURGE_MEMALLOC(zmurge_solvers[id]->perm, zmurge_solvers[id]->n, pastix_int_t);
        memset(zmurge_solvers[id]->perm, 0, zmurge_solvers[id]->n*sizeof(pastix_int_t));
      }
      z_pastix_welcome_print(zmurge_solvers[id]->pastix_data,
                           zmurge_solvers[id]->colptr,
                           zmurge_solvers[id]->n);

      ALLOC_INVP;
      /* Call preprocessing separatly if required so that we can delete CSC
       * before Factorization call
       */
      iparm[IPARM_END_TASK]    = API_TASK_BLEND;
      DPASTIX(&pastix_data,
	      pastix_data->pastix_comm,
	      zmurge_solvers[id]->n,
	      zmurge_solvers[id]->colptr,
	      zmurge_solvers[id]->rows,
	      zmurge_solvers[id]->values,
	      zmurge_solvers[id]->l2g,
	      zmurge_solvers[id]->perm,
	      NULL,
	      zmurge_solvers[id]->b,
	      zmurge_solvers[id]->nrhs,
	      pastix_data->iparm,
	      pastix_data->dparm);
      murge_redistribute_data(id);
    }
    CLOCK_PRINT("check_fact -- preprocessing");

    /* If not performed yet we
     * - fill the internal CSC
     * - delete murge CSC and
     * - perform factorization
     */
    iparm[IPARM_CSCD_CORRECT] = API_YES;
    PASTIX_FILLIN_CSC(zmurge_solvers[id]->pastix_data,
		      zmurge_solvers[id]->pastix_data->pastix_comm,
		      zmurge_solvers[id]->n,
		      zmurge_solvers[id]->colptr,
		      zmurge_solvers[id]->rows,
		      zmurge_solvers[id]->values,
		      zmurge_solvers[id]->b,
		      zmurge_solvers[id]->nrhs,
		      zmurge_solvers[id]->l2g);
    CLOCK_PRINT("check_fact -- PASTIX_FILLIN_CSC");
    zmurge_solvers[id]->pastix_data->cscInternFilled = API_YES;
    if (zmurge_solvers[id]->pastix_data->iparm[IPARM_FREE_CSCUSER] == API_CSC_FREE) {
      MURGE_FREE(zmurge_solvers[id]->colptr);
      MURGE_FREE(zmurge_solvers[id]->rows);
      MURGE_FREE(zmurge_solvers[id]->values);
    }
    zmurge_solvers[id]->pastix_data->iparm[IPARM_START_TASK]  =
      API_TASK_NUMFACT;
  } else {
    /* Else, solve is enough */
    if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_SOLVE_DONE))) {
      zmurge_solvers[id]->pastix_data->iparm[IPARM_START_TASK]  =
        API_TASK_SOLVE;
    } else {
      zmurge_solvers[id]->pastix_data->iparm[IPARM_START_TASK]  =
        API_TASK_REFINE;
    }
  }

  if (iparm[IPARM_ONLY_RAFF] ==  API_YES) {
    /* When only raff we need to call solve to set RHS */
    iparm[IPARM_END_TASK] = API_TASK_SOLVE;
    DPASTIX(&pastix_data,
            pastix_data->pastix_comm,
            zmurge_solvers[id]->n,
            zmurge_solvers[id]->colptr,
            zmurge_solvers[id]->rows,
            zmurge_solvers[id]->values,
            zmurge_solvers[id]->l2g,
            zmurge_solvers[id]->perm,
            NULL,
            zmurge_solvers[id]->b,
            zmurge_solvers[id]->nrhs,
            pastix_data->iparm,
            pastix_data->dparm);
  }
  if (iparm[IPARM_MURGE_REFINEMENT] == API_YES) {
    iparm[IPARM_END_TASK] = API_TASK_REFINE;
    if (zmurge_solvers[id]->pastix_data->cscInternFilled == API_NO) {
      errorPrint("Trying to refine without internal CSC\n"
                 "\t You need to keep PaStiX internal CSC using"
                 " IPARM_MURGE_MAY_REFINE = API_YES.");
      return MURGE_ERR_PARAMETER;
    }
  } else {
    /* No need to IPARM_FREE_CSCUSER as it is done when facto is performed,
     * after fillin internal CSC, not possible to free it inside PaStiX as it
     * was allocated with memAlloc
     */
    if (iparm[IPARM_MURGE_MAY_REFINE] == API_NO)
      iparm[IPARM_FREE_CSCPASTIX] = API_CSC_FREE;
    iparm[IPARM_END_TASK] = API_TASK_SOLVE;
  }
  DPASTIX(&pastix_data,
          pastix_data->pastix_comm,
          zmurge_solvers[id]->n,
          zmurge_solvers[id]->colptr,
          zmurge_solvers[id]->rows,
          zmurge_solvers[id]->values,
          zmurge_solvers[id]->l2g,
          zmurge_solvers[id]->perm,
          NULL,
          zmurge_solvers[id]->b,
          zmurge_solvers[id]->nrhs,
          pastix_data->iparm,
          pastix_data->dparm);
  CLOCK_PRINT("check_fact -- DPASTIX");

  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_SYMB_OK);
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_SOLVE_DONE);
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_REFINE_DONE);
  if (iparm[IPARM_MURGE_REFINEMENT] == API_NO &&
      iparm[IPARM_MURGE_MAY_REFINE] == API_NO) {
    zmurge_solvers[id]->pastix_data->cscInternFilled = API_NO;
  }
  CLOCK_PRINT("check_fact -- end");

  return MURGE_SUCCESS;
}

static inline
void zmurge_pastix_fillin_csc(INTS            id,
                             z_pastix_data_t * data,
                             MPI_Comm        comm,
                             pastix_int_t      n,
                             pastix_int_t    * colptr,
                             pastix_int_t    * rows,
                             pastix_complex64_t  * values,
                             pastix_complex64_t  * b,
                             pastix_int_t      nrhs,
                             pastix_int_t    * l2g) {
  pastix_int_t tmpN;
  pastix_int_t   *tmpcolptr = NULL, *tmprows = NULL;
  pastix_int_t   *tmpperm = NULL, *tmpinvp = NULL;
  pastix_complex64_t *tmpvalues = NULL, *tmprhs = NULL;
  INTS has_values = API_NO;
  INTS has_rhs = API_NO;

  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_VALUES_OK))
    has_values = API_YES;
  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_RHS_OK))
    has_rhs = API_NO; /* We always enter NULL ptr here... */


  if (data->procnum == 0 && data->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    errorPrintW("To get an optimal MURGE installation"
                " please build PaStiX with -DDISTRIBUTED\n");

  z_cscd2csc_int(n, colptr, rows, values,
                 b, NULL, NULL,
                 &tmpN, &tmpcolptr, &tmprows,
                 (has_values==API_NO)?NULL:&tmpvalues,
                 (has_rhs==API_NO)?NULL:&tmprhs, NULL, NULL,
                 l2g, comm, data->iparm[IPARM_DOF_NBR], API_YES);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpcolptr), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmprows), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpvalues), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpperm), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpinvp), char);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmprhs), char);
  z_pastix_fillin_csc(data,
                    comm,
                    tmpN,
                    tmpcolptr,
                    tmprows,
                    tmpvalues,
                    tmprhs,
                    nrhs,
                    NULL);
  if (NULL != tmpcolptr) MURGE_FREE(tmpcolptr);
  if (NULL != tmprows) MURGE_FREE(tmprows);
  if (NULL != values) MURGE_FREE(tmpvalues);
  if (NULL != tmpperm) MURGE_FREE(tmpperm);
  if (NULL != tmpinvp) MURGE_FREE(tmpinvp);
  if (NULL != b) MURGE_FREE(tmprhs);
}

static inline
int zmurge_dpastix(INTS             id,
                   z_pastix_data_t ** data,
                   MPI_Comm         comm,
                   pastix_int_t       n,
                   pastix_int_t     * colptr,
                   pastix_int_t     * rows,
                   pastix_complex64_t   * values,
                   pastix_int_t     * l2g,
                   pastix_int_t     * perm,
                   pastix_int_t     * invp,
                   pastix_complex64_t   * b,
                   pastix_int_t       nrhs,
                   pastix_int_t     * iparm,
                   double         * dparm) {
  pastix_int_t tmpN;
  pastix_int_t   *tmpcolptr = NULL, *tmprows = NULL;
  pastix_int_t   *tmpperm = NULL, *tmpinvp = NULL;
  pastix_complex64_t *tmpvalues = NULL, *tmprhs = NULL;
  pastix_int_t save_veri = iparm[IPARM_MATRIX_VERIFICATION];
  pastix_int_t save_free = iparm[IPARM_FREE_CSCUSER];
  INTS i;
  INTS has_values = API_NO;
  INTS has_rhs = API_NO;
  INTS has_perm = API_NO;

  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_VALUES_OK))
    has_values = API_YES;
  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_RHS_OK))
    has_rhs = API_YES;
  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK))
    has_perm = API_YES;

  if ((*data)->procnum == 0 && (*data)->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    errorPrintW("To get an optimal MURGE installation"
                " please build PaStiX with -DDISTRIBUTED\n");
  /* After factorization, no more CSC if API_CSC_FREE */
  if ( ( zmurge_solvers[id]->pastix_data->iparm[IPARM_FREE_CSCUSER] == API_CSC_FREE &&
         iparm[IPARM_END_TASK] > API_TASK_NUMFACT )) {
      tmpN = zmurge_solvers[id]->N;
      tmpcolptr = NULL;
      tmprows   = NULL;
      tmpvalues = NULL;
      tmpperm   = NULL;
      tmpinvp   = NULL;
      if (has_rhs) {
          COEF * tmprhs2;
          MURGE_MEMALLOC(tmprhs, zmurge_solvers[id]->N,  pastix_complex64_t);
          MURGE_MEMALLOC(tmprhs2, zmurge_solvers[id]->N, pastix_complex64_t);
          memset(tmprhs, 0, zmurge_solvers[id]->N*sizeof(COEF));
          for (i = 0; i < n; i++) {
              tmprhs2[l2g[i]-1] = b[i];
          }
          MPI_Allreduce(tmprhs2, tmprhs, zmurge_solvers[id]->N, MPI_DOUBLE_COMPLEX, MPI_SUM,
                        zmurge_solvers[id]->pastix_data->pastix_comm);
          MURGE_FREE(tmprhs2);
      }
  } else {
      z_cscd2csc_int(n, colptr, rows, values,
                     b, perm, NULL,
                     &tmpN, &tmpcolptr, &tmprows,
                     (has_values==API_NO)?NULL:&tmpvalues,
                     (has_rhs==API_NO)?NULL:&tmprhs,
                     (has_perm==API_NO)?NULL:&tmpperm, &tmpinvp,
                     l2g, comm, (*data)->iparm[IPARM_DOF_NBR], API_YES);
      MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpcolptr), char);
      MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmprows), char);
      MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpvalues), char);
      MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpperm), char);
      MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpinvp), char);
      MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmprhs), char);
  }

  iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
  iparm[IPARM_FREE_CSCUSER] = API_CSC_PRESERVE;
  z_pastix(data, comm, tmpN, tmpcolptr, tmprows, tmpvalues,
           tmpperm, tmpinvp, tmprhs, nrhs, iparm, dparm);
  iparm[IPARM_MATRIX_VERIFICATION] = save_veri;
  iparm[IPARM_FREE_CSCUSER] = save_free;

  if (has_perm == API_YES)
    for (i = 0; i < n;i++)
      perm[i] = tmpperm[l2g[i]-1];
  if (has_rhs == API_YES)
    for (i = 0; i < n;i++)
      b[i] = tmprhs[l2g[i]-1];
  if (NULL != tmpcolptr) MURGE_FREE(tmpcolptr);
  if (NULL != tmprows) MURGE_FREE(tmprows);
  if (NULL != tmpvalues) MURGE_FREE(tmpvalues);
  if (NULL != tmpperm) MURGE_FREE(tmpperm);
  if (NULL != tmpinvp) MURGE_FREE(tmpinvp);
  if (NULL != tmprhs) MURGE_FREE(tmprhs);
  return MURGE_SUCCESS;
}


/*******************************************************************************
 * Group: Solver setup functions
 */

/*
 * Function: ZMURGE_GetSolver
 *
 * returns ZMURGE_SOLVER_PASTIX
 */
INTS ZMURGE_GetSolver(INTS * solver_id) {
  *solver_id = MURGE_SOLVER_PASTIX;
  return MURGE_SUCCESS;
}

static inline
INTS _ZMURGE_InitId(INTS murge_id) {
  
  zmurge_solvers[murge_id] = (z_murge_data_t*)malloc(sizeof(z_murge_data_t));
  zmurge_solvers[murge_id]->n           = 0;
  zmurge_solvers[murge_id]->N           = 0;
  zmurge_solvers[murge_id]->colptr      = NULL;
  zmurge_solvers[murge_id]->rows        = NULL;
  zmurge_solvers[murge_id]->values      = NULL;
  zmurge_solvers[murge_id]->l2g         = NULL;
  zmurge_solvers[murge_id]->g2l         = NULL;
  zmurge_solvers[murge_id]->perm        = NULL;
#ifdef CENTRALISED
  zmurge_solvers[murge_id]->invp        = NULL;
#endif
  zmurge_solvers[murge_id]->b           = NULL;
  zmurge_solvers[murge_id]->nrhs        = 1;
  zmurge_solvers[murge_id]->state       = MURGE_INIT_OK;
  zmurge_solvers[murge_id]->pastix_data = NULL;
  zmurge_solvers[murge_id]->sequences   = NULL;
  zmurge_solvers[murge_id]->seq_ID      = 0;
  zmurge_solvers[murge_id]->ndump       = 0;
  zmurge_solvers[murge_id]->dropmask    = NULL;
  zmurge_solvers[murge_id]->droprows    = NULL;
  zmurge_solvers[murge_id]->dropcols    = NULL;
  zmurge_solvers[murge_id]->malloc_size = 0;
  zmurge_solvers[murge_id]->malloc_maxsize = 0;
  zmurge_solvers[murge_id]->threadnbr   = 0;
  zmurge_solvers[murge_id]->threads_state = 0;
#ifndef FORCE_NOSMP
  zmurge_solvers[murge_id]->barrier.instance         = 0;
  zmurge_solvers[murge_id]->barrier.blocked_threads  = 0;
  pthread_mutex_init(&(zmurge_solvers[murge_id]->barrier.sync_lock), NULL);
  pthread_cond_init(&(zmurge_solvers[murge_id]->barrier.sync_cond), NULL);
#endif
  z_pastix_task_init(&(zmurge_solvers[murge_id]->pastix_data), MPI_COMM_WORLD, NULL, NULL);
  return MURGE_SUCCESS;
}


/*
 * Function: ZMURGE_Initialize
 *
 * Allocate the instance arrays which will keeps intern data for all
 * solver instances.
 *
 * If user is creating several threads calling the solver, this function
 * has to be called before creating threads to insure solver is thread safe.
 *
 * Parameters:
 *   idnbr - Maximum number of solver instances that will be
 *           launched.
 *
 * Returns:
 *   MURGE_SUCCESS      - If function runned successfully.
 *   MURGE_ERR_ALLOCATE - If for some reason, allocation was not
 *                        successfull.
 */
INTS ZMURGE_Initialize(INTS id_nbr) {
  INTS i;

  print_debug(DBG_MURGE, ">> Zmurge_Initialize\n");

  if (sizeof(COEF) != sizeof(pastix_complex64_t)) {
    errorPrint("Incompatible coefficient type\n");
    return MURGE_ERR_PARAMETER;
  }

  if ( (zmurge_solvers != NULL) ) {
    errorPrint("ZMURGE_Initialize has been already called");
    return MURGE_ERR_ORDER;
  }

  idnbr = id_nbr;

  zmurge_solvers = (z_murge_data_t**)malloc(idnbr*sizeof(z_murge_data_t*));

  for (i=0; i< idnbr; i++) {
    zmurge_solvers[i] = NULL;
    zmurge_solvers[i] = (z_murge_data_t*)malloc(sizeof(z_murge_data_t));

    zmurge_solvers[i]->n           = 0;
    zmurge_solvers[i]->N           = 0;
    zmurge_solvers[i]->colptr      = NULL;
    zmurge_solvers[i]->rows        = NULL;
    zmurge_solvers[i]->values      = NULL;
    zmurge_solvers[i]->l2g         = NULL;
    zmurge_solvers[i]->g2l         = NULL;
    zmurge_solvers[i]->perm        = NULL;
#ifdef CENTRALISED
    zmurge_solvers[i]->invp        = NULL;
#endif
    zmurge_solvers[i]->b           = NULL;
    zmurge_solvers[i]->nrhs        = 1;
    zmurge_solvers[i]->state       = MURGE_INIT_OK;
    zmurge_solvers[i]->pastix_data = NULL;
    zmurge_solvers[i]->sequences   = NULL;
    zmurge_solvers[i]->seq_ID      = 0;
    zmurge_solvers[i]->ndump       = 0;
    zmurge_solvers[i]->dropmask    = NULL;
    zmurge_solvers[i]->droprows    = NULL;
    zmurge_solvers[i]->dropcols    = NULL;
    zmurge_solvers[i]->malloc_size = 0;
    zmurge_solvers[i]->malloc_maxsize = 0;
    zmurge_solvers[i]->threadnbr   = 0;
    zmurge_solvers[i]->threads_state = 0;
#ifndef FORCE_NOSMP
    zmurge_solvers[i]->barrier.instance         = 0;
    zmurge_solvers[i]->barrier.blocked_threads  = 0;
    pthread_mutex_init(&(zmurge_solvers[i]->barrier.sync_lock), NULL);
    pthread_cond_init(&(zmurge_solvers[i]->barrier.sync_cond), NULL);
#endif
    z_pastix_task_init(&(zmurge_solvers[i]->pastix_data), MPI_COMM_WORLD, NULL, NULL);
  }

  print_debug(DBG_MURGE, "<< Zmurge_Initialize\n");
  return MURGE_SUCCESS;
}

/*
 * Function: ZMURGE_SetDefaultOptions
 *
 * Create a solver instance if not created yet.
 *
 * Sets default options, for solver instance number *id*.
 *
 * The default option set correspond to *stratnum* strategy ID,
 * depending on the solver.
 *
 * Needs <ZMURGE_Initialize> to be called before
 * to allocate solver instances array.
 *
 * Parameters:
 *   id       - Solver instance identification number.
 *   stratnum - Strategy for the default option Set.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <MURGE_Initialize> was not called before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *stratnum* is not valid.
 *   MURGE_ERR_ALLOCATE  - If couldn't create solver instance.
 */
INTS ZMURGE_SetDefaultOptions(INTS id, INTS stratnum) {
  print_debug(DBG_MURGE, ">> ZMURGE_SetDefaultOptions\n");
  CHECK_SOLVER_ID(id);

  if (zmurge_solvers[id]->pastix_data->iparm == NULL) {
    MURGE_MEMALLOC_EXT(zmurge_solvers[id]->pastix_data->iparm,
		       IPARM_SIZE, pastix_int_t);
    MURGE_MEMALLOC_EXT(zmurge_solvers[id]->pastix_data->dparm,
		       DPARM_SIZE, double);
#ifdef MURGE_THREADSAFE
    pthread_mutex_init(&zmurge_solvers[id]->mutex_tmpmatrix, NULL);
#endif
  }

  z_pastix_initParam(zmurge_solvers[id]->pastix_data->iparm,
                   zmurge_solvers[id]->pastix_data->dparm);

  zmurge_solvers[id]->pastix_data->iparm[IPARM_PID] =
    zmurge_solvers[id]->pastix_data->pastix_id;
  zmurge_solvers[id]->pastix_data->iparm[IPARM_RHS_MAKING] = API_RHS_B;
  zmurge_solvers[id]->pastix_data->iparm[IPARM_FREE_CSCUSER] = API_YES;
  return MURGE_SUCCESS;
}


/*
 * Function: ZMURGE_SetOptionINT
 *
 * Sets integer option, indicated by *number*, to *value* for the
 * solver instance number *id*.
 *
 * Needs <ZMURGE_SetDefaultOption> to be called before to initiate
 * solver instance data.
 *
 * Parameters:
 *   id     - Solver instance identification number.
 *   number - Identification of the integer parameter.
 *   value  - value to set the parameter to.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <ZMURGE_SetDefaultOption> was not
 *                         called before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
 *                         *number* or *value* are not valid.
 *
 */
INTS ZMURGE_SetOptionINT (INTS id, INTS number, INTS value) {
  pastix_int_t murge_param[64];
  pastix_int_t * iparm = NULL;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  iparm = zmurge_solvers[id]->pastix_data->iparm;

  murge_param[MURGE_IPARAM_BASEVAL       - 1024] =  IPARM_BASEVAL;
  murge_param[MURGE_IPARAM_DOF           - 1024] =  IPARM_DOF_NBR;
  murge_param[MURGE_IPARAM_SYM           - 1024] =  IPARM_SYM;

  if (number >= 1024) {
    if (number == MURGE_IPARAM_SYM) {
      if (value == MURGE_BOOLEAN_TRUE) {
        value = API_SYM_YES;
        if (iparm[IPARM_FACTORIZATION] == API_FACT_LU) {
          iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
        }
      } else {
        if (value == MURGE_BOOLEAN_FALSE) {
          value = API_SYM_NO;
          if (iparm[IPARM_FACTORIZATION] != API_FACT_LU) {
            iparm[IPARM_FACTORIZATION] = API_FACT_LU;
          }
        } else {
          errorPrint("Invalid value");
          return MURGE_ERR_PARAMETER;
        }
      }
    }
    number = murge_param[number-1024];
  }

  if (!( number < IPARM_SIZE )) {
    errorPrint("option number '%d' is too big", number);
    return MURGE_ERR_PARAMETER;
  }

  if (number < 0) {
    errorPrint("option number '%d' is negative", number);
    return MURGE_ERR_PARAMETER;
  }

  /* TODO : Est-ce qu'on ajoute des tests sur les valeurs rentrées ???? */
  iparm[number] = (pastix_int_t)value;
  return MURGE_SUCCESS;
}

/*
 Function: ZMURGE_SetOptionREAL

 Sets real option, indicated by *number*, to *value* for the
 solver instance number *id*.

 Needs <ZMURGE_SetDefaultOption> to be called before to initiate
 solver instance data.

 Parameters:
 id     - Solver instance identification number.
 number - Identification of the integer parameter.
 value  - value to set the parameter to.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_ORDER     - If <ZMURGE_SetDefaultOption> was not
 called before.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
 *number* or *value* are not valid.

 */
INTS ZMURGE_SetOptionREAL(INTS id, INTS number, REAL value) {
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  pastix_int_t murge_param[64];

  murge_param[MURGE_RPARAM_EPSILON_ERROR- 1024] =  DPARM_EPSILON_REFINEMENT;

  if (number >= 1024) {
    number = murge_param[number-1024];
  }

  if (!( number < DPARM_SIZE )) {
    errorPrint("number is too big");
    return MURGE_ERR_PARAMETER;
  }

  if (number < 0) {
    errorPrint("number is negative");
    return MURGE_ERR_PARAMETER;
  }

  /* TODO : Est-ce qu'on ajoute des tests sur les valeurs rentrées ???? */
  zmurge_solvers[id]->pastix_data->dparm[number] = (double)value;
  return MURGE_SUCCESS;
}

/*
 * Function: ZMURGE_SetCommunicator
 *
 * Sets MPI communicator for the given solver instance.
 *
 * Needs <ZMURGE_SetDefaultOption> to be called before to initiate
 * solver instance data.
 *
 * Musn't be called before <MURGE_SAVE>, <MURGE_LOAD>,
 * <ZMURGE_GetLocalNodeNbr> nor <ZMURGE_GetLocalUnknownNbr>
 * because the solver as to be runned with the same MPI
 * communicator all along.
 *
 * If this function is not called, MPI communicator will be
 * *MPI_COMM_WORLD*.
 *
 * This function may not exist if the solver
 * has been compiled without MPI.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   mpicomm - MPI communicator to affect the solver to.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <ZMURGE_SetDefaultOption> was not
 *                         called before or if it is called after
 *                         the solver starts its computing tasks.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range or
 *                         *number* or *value* are not valid.
 */
INTS ZMURGE_SetCommunicator(INTS id, MPI_Comm mpicomm) {
  CHECK_SOLVER_ID(id);
  zmurge_solvers[id]->pastix_data->pastix_comm = mpicomm;
  zmurge_solvers[id]->pastix_data->inter_node_comm = mpicomm;
  zmurge_solvers[id]->pastix_data->intra_node_comm = MPI_COMM_SELF;
  MPI_Comm_size((zmurge_solvers[id]->pastix_data)->inter_node_comm,
                &((zmurge_solvers[id]->pastix_data)->inter_node_procnbr));
  MPI_Comm_rank((zmurge_solvers[id]->pastix_data)->inter_node_comm,
                &((zmurge_solvers[id]->pastix_data)->inter_node_procnum));
  MPI_Comm_size((zmurge_solvers[id]->pastix_data)->intra_node_comm,
                &((zmurge_solvers[id]->pastix_data)->intra_node_procnbr));
  MPI_Comm_rank((zmurge_solvers[id]->pastix_data)->intra_node_comm,
                &((zmurge_solvers[id]->pastix_data)->intra_node_procnum));

  MPI_Comm_size(mpicomm, &(zmurge_solvers[id]->pastix_data->procnbr));
  MPI_Comm_rank(mpicomm, &(zmurge_solvers[id]->pastix_data->procnum));
  return MURGE_SUCCESS;
}


/*******************************************************************************
 * Group: I/O functions
 */

/*
 * Function: ZMURGE_Save
 *
 * Runs preprocessing step, if not done yet, and save the result to disk,
 * into *directory*, so that it can be resume using <ZMURGE_Load>.
 *
 * Needs <ZMURGE_SetDefaultOption> to be called before to initiate
 * solver instance data.
 *
 * Parameters:
 *   id        - Solver instance identification number.
 *   directory - Path to the directory where to save the solver step.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <ZMURGE_SetDefaultOption> was not
 *                         called before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 *   MURGE_ERR_IO        - If file(s) couldn't be writen.
 */
INTS ZMURGE_Save(INTS id, char* directory) {
  char * dest    = NULL;
  char * src     = NULL;
  int    procnum;
  FILE * stream  = NULL;

  CHECK_SOLVER_ID(id);
  procnum = (int)zmurge_solvers[id]->pastix_data->procnum;

  if (zmurge_solvers[id]->pastix_data->iparm == NULL) {
    errorPrint("You need to call ZMURGE_SetDefaultOptions before");
    return MURGE_ERR_ORDER;
  }
  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("You need to set graph before");
    return MURGE_ERR_ORDER;
  }

  zmurge_solvers[id]->pastix_data->iparm[IPARM_IO_STRATEGY] = API_IO_SAVE;
  zmurge_solvers[id]->pastix_data->iparm[IPARM_START_TASK]  = API_TASK_ORDERING;
  zmurge_solvers[id]->pastix_data->iparm[IPARM_END_TASK]    = API_TASK_SYMBFACT;

  if (NULL == zmurge_solvers[id]->perm) {
    MURGE_MEMALLOC(zmurge_solvers[id]->perm, zmurge_solvers[id]->n, pastix_int_t);
    memset(zmurge_solvers[id]->perm, 0, zmurge_solvers[id]->n*sizeof(pastix_int_t));
  }
  z_pastix_welcome_print(zmurge_solvers[id]->pastix_data,
                         zmurge_solvers[id]->colptr,
                         zmurge_solvers[id]->n);
  DPASTIX(&(zmurge_solvers[id]->pastix_data),
          zmurge_solvers[id]->pastix_data->pastix_comm,
          zmurge_solvers[id]->n,
          zmurge_solvers[id]->colptr,
          zmurge_solvers[id]->rows,
          zmurge_solvers[id]->values,
          zmurge_solvers[id]->l2g,
          zmurge_solvers[id]->perm,
          NULL,
          zmurge_solvers[id]->b,
          zmurge_solvers[id]->nrhs,
          zmurge_solvers[id]->pastix_data->iparm,
          zmurge_solvers[id]->pastix_data->dparm);

  MURGE_MEMALLOC(dest, (strlen(directory)+20), char);
  MURGE_MEMALLOC(src, 20, char);

  if (procnum == 0) {
    sprintf(dest,"%s/ordergen", directory);
    RENAME("ordergen", dest);
    sprintf(dest,"%s/symbgen", directory);
    RENAME("symbgen", dest);

    sprintf(dest,"%s/murge.save", directory);
    PASTIX_FOPEN(stream, dest, "w");
    fclose(stream);

  }

  MURGE_FREE(src);
  MURGE_FREE(dest);

  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_BLEND_OK);
  return MURGE_SUCCESS;
}

/*
 * Function: ZMURGE_Load
 *
 * Loads preprocessing result from disk, into *directory*,
 * where it had been saved by <ZMURGE_Save>.
 *
 * If preprocessing data was already computed or loaded, it will
 * be overwriten.
 *
 * Needs <ZMURGE_SetDefaultOption> to be called before to initiate
 * solver instance data.
 *
 * Parameters:
 *   id        - Solver instance identification number.
 *   directory - Path to the directory where to load the solver
 *               preprocessing data.
 *
 * In Fortran, *STR_LEN* is the length of the string directory.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If <ZMURGE_SetDefaultOption> was not
 *                         called before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 *   MURGE_ERR_IO        - If file(s) couldn't be read.
 */
INTS ZMURGE_Load(INTS id, char* directory) {
  char * src     = NULL;
  int    procnum;

  CHECK_SOLVER_ID(id);
  procnum = (int)zmurge_solvers[id]->pastix_data->procnum;

  if (zmurge_solvers[id]->pastix_data->iparm == NULL) {
    errorPrint("You need to call ZMURGE_SetDefaultOptions before");
    return MURGE_ERR_ORDER;
  }

  zmurge_solvers[id]->pastix_data->iparm[IPARM_IO_STRATEGY] = API_IO_LOAD;
  zmurge_solvers[id]->pastix_data->iparm[IPARM_START_TASK]  = API_TASK_ORDERING;
  zmurge_solvers[id]->pastix_data->iparm[IPARM_END_TASK]    = API_TASK_BLEND;

  z_pastix_welcome_print(zmurge_solvers[id]->pastix_data,
                         zmurge_solvers[id]->colptr,
                         zmurge_solvers[id]->n);

  if (NULL == zmurge_solvers[id]->perm) {
    MURGE_MEMALLOC(zmurge_solvers[id]->perm, zmurge_solvers[id]->n, pastix_int_t);
    memset(zmurge_solvers[id]->perm, 0, zmurge_solvers[id]->n*sizeof(pastix_int_t));
  }
  MURGE_MEMALLOC(src, (strlen(directory)+20), char);
  if (procnum == 0) {
    sprintf(src,"%s/ordergen", directory);
    LINK(src, "ordername");
    sprintf(src,"%s/symbgen", directory);
    LINK(src, "symbname");
  }
  MPI_Barrier(zmurge_solvers[id]->pastix_data->pastix_comm);
  MURGE_FREE(src);

  DPASTIX(&(zmurge_solvers[id]->pastix_data),
          zmurge_solvers[id]->pastix_data->pastix_comm,
          zmurge_solvers[id]->n,
          zmurge_solvers[id]->colptr,
          zmurge_solvers[id]->rows,
          zmurge_solvers[id]->values,
          zmurge_solvers[id]->l2g,
          zmurge_solvers[id]->perm,
          NULL,
          zmurge_solvers[id]->b,
          zmurge_solvers[id]->nrhs,
          zmurge_solvers[id]->pastix_data->iparm,
          zmurge_solvers[id]->pastix_data->dparm);
  zmurge_solvers[id]->N = zmurge_solvers[id]->pastix_data->ordemesh->rangtab[zmurge_solvers[id]->pastix_data->ordemesh->cblknbr];
  if (procnum == 0) {
    UNLINK("ordername");
    UNLINK("symbname");
  }


  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_SYMB_OK);
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK);
  return MURGE_SUCCESS;
}


/*******************************************************************************
 * Group: Getting solver's distribution
 */

/*
 * Function: ZMURGE_GetLocalNodeNbr
 *
 * Computes preprocessing step, if not done, and the number of
 * Nodes in the new ditribution of the matrix.
 *
 * Parameters:
 *   id        - Solver instance identification number.
 *   nodenbr   - *INTS* where to store number of nodes.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *nodenbr* is *NULL* (can occur in C).
 */
INTS ZMURGE_GetLocalNodeNbr    (INTS id, INTS *nodenbr){
#ifdef MURGE_TIME
  Clock clock;
#endif
  CLOCK_INIT;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);

  *nodenbr = (INTS)z_pastix_getLocalNodeNbr(&(zmurge_solvers[id]->pastix_data));
  CLOCK_PRINT("CALL z_pastix_getLocalNodeNbr");

  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_NODENBR_OK);
  return MURGE_SUCCESS;
}

/*
 * Function: ZMURGE_GetLocalNodeList
 *
 * Computes the local node list, corresponding to
 * the new distribution, after preprocessing.
 *
 * *nodelist* array has to be allocated before calling
 * this function.
 *
 * As it's result determines the size of *nodelist*
 * array, <ZMURGE_GetLocalNodeNbr> should be run before it.
 *
 * Parameters:
 *   id        - Solver instance identification number.
 *   nodelist  - Array where to store the list of local nodes.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - if <ZMURGE_GetLocalNodeNbr> has not been called
 *                         before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *nodelist* is *NULL* (can occur in C).
 */
INTS ZMURGE_GetLocalNodeList   (INTS id, INTS *nodelist) {
  int ret;
  pastix_int_t nodenbr = 0;
  pastix_int_t i;
  pastix_int_t * intern_nodelist;
#ifdef MURGE_TIME
  Clock clock;
#endif
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_NODENBR_OK))) {
    errorPrint("You need to call ZMURGE_GetLocalNodeNbr before");
    return MURGE_ERR_ORDER;
  }
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_NODELST_OK);

  nodenbr = z_pastix_getLocalNodeNbr(&(zmurge_solvers[id]->pastix_data));
  if (sizeof(pastix_int_t) != sizeof(INTS)) {
    MURGE_MEMALLOC(intern_nodelist, nodenbr, pastix_int_t);
  }
  else
    {
      intern_nodelist = (pastix_int_t*)nodelist;
    }

    CLOCK_INIT;
    if (EXIT_SUCCESS != ( ret =
                          z_pastix_getLocalNodeLst(&(zmurge_solvers[id]->pastix_data),
                                                   intern_nodelist)))
    return MURGE_ERR_SOLVER;
    CLOCK_PRINT("CALL z_pastix_getLocalNodeLst");

  if (sizeof(pastix_int_t) != sizeof(INTS)) {
    for (i = 0; i < nodenbr; i++) {
      nodelist[i] = (INTS) intern_nodelist[i];
    }
    MURGE_FREE(intern_nodelist);
  }
  if (zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
    for (i = 0; i < nodenbr; i++)
      nodelist[i] -= 1;
  }

  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_NODELST_OK);
  return MURGE_SUCCESS;
}

/*
 * Function: ZMURGE_GetLocalUnkownNbr
 *
 * Computes preprocessing step, if not done, and the number of
 * Unkowns in the new ditribution of the matrix.
 *
 * Parameters:
 *   id            - Solver instance identification number.
 *   unkownnbr     - *INTS* where to store number of unkowns.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *unkownnbr* is *NULL* (can occur in C).
 */
INTS ZMURGE_GetLocalUnknownNbr (INTS id, INTS *unkownnbr) {
#ifdef MURGE_TIME
  Clock clock;
#endif
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);
  CLOCK_INIT;
  *unkownnbr = (INTS)z_pastix_getLocalUnknownNbr(&(zmurge_solvers[id]->pastix_data));
  CLOCK_PRINT("CALL z_pastix_getLocalUnknownNbr");
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_NODENBR_OK);
  return MURGE_SUCCESS;
}

/*
 * Function: ZMURGE_GetLocalUnkownList
 *
 * Computes the local unkown list, corresponding to
 * the new distribution, after preprocessing.
 *
 * *unkownlist* array has to be allocated before calling
 * this function.
 *
 * As it's result determines the size of *unkownlist*
 * array, <ZMURGE_GetLocalUnkownNbr> should be run before it.
 *
 * Parameters:
 *   id          - Solver instance identification number.
 *   unkownlist  - Array where to store the list of local unkowns.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - if <ZMURGE_GetLocalUnkownNbr> has not been called
 *                         before.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range
 *                         or *unkownlist* is *NULL* (can occur in C).
 */
INTS ZMURGE_GetLocalUnknownList(INTS id, INTS *unkownlist){
  int ret;
  pastix_int_t nodenbr = 0;
  pastix_int_t i;
  pastix_int_t * intern_nodelist = NULL;
#ifdef MURGE_TIME
  Clock clock;
#endif
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_NODENBR_OK))) {
    errorPrint("You need to call ZMURGE_GetLocalNodeNbr before");
    return MURGE_ERR_ORDER;
  }

  nodenbr = z_pastix_getLocalUnknownNbr(&(zmurge_solvers[id]->pastix_data));
  if (sizeof(pastix_int_t) != sizeof(INTS)) {
    MURGE_MEMALLOC(intern_nodelist, nodenbr, pastix_int_t);
  }
  else
    {
      intern_nodelist = (pastix_int_t*)unkownlist;
    }
  CLOCK_INIT;
  if (EXIT_SUCCESS != ( ret =
                        z_pastix_getLocalUnknownLst(&(zmurge_solvers[id]->pastix_data),
                                                    intern_nodelist)))
    return MURGE_ERR_SOLVER;
  CLOCK_PRINT("CALL z_pastix_getLocalUnknownLst");

  if (sizeof(pastix_int_t) != sizeof(INTS)) {
    for (i = 0; i < nodenbr; i++) {
      unkownlist[i] = (INTS) intern_nodelist[i];
    }
    MURGE_FREE(intern_nodelist);
  }
  if (zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
    for (i = 0; i < nodenbr; i++)
      unkownlist[i] -= 1;
  }
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_NODELST_OK);
  return MURGE_SUCCESS;
}



/*******************************************************************************
 * Group: Graph setup functions
 */

/*
 * Function: ZMURGE_GraphBegin
 *
 * - Allocate temporary structure which will contain graph entries.
 * - Set number of unkowns in the graph.
 * - Set the number of entries that are expected in this building session.
 * - Reset the number of entries for this build session.
 * - Set all states except MURGE_GRAPH_BUILD to FALSE (graph, values, blend,
 * nodelst, nodenbr, facto)
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   N       - Number of unkowns.
 *   edgenbr - Number of edges in this building session.
 *             If edgenbr is negative, PaStiX will perform dynamic
 *             reallocation of the array, with the first allocation of
 *             size -edgenbr.
 *
 * Returns:
 *   MURGE_ERR_ORDER     - ZMURGE_GraphBegin has already been called, or if
 *                         *zmurge_solvers* or *zmurge_solvers[id]* are not allocated,
 *                         or if *iparm* or *dparm* are not allocated.
 *   MURGE_ERR_PARAMETER - If *id* is not in correct range.
 *   MURGE_SUCCESS       - Otherwise.
 */
INTS ZMURGE_GraphBegin(INTS id, INTS N, INTL edgenbr) {
  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_BUILD)) {
    errorPrint("ZMURGE_GraphBegin has been called before");
    return MURGE_ERR_ORDER;
  }

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif
  if (edgenbr < 0) {
    edgenbr = -edgenbr;
    zmurge_solvers[id]->dynamic = API_YES;
  }
  else {
    zmurge_solvers[id]->dynamic = API_NO;
  }
#ifdef HASH_MTX
  MURGE_MEMALLOC(zmurge_solvers[id]->hashgraphtab, edgenbr, hash_graph_entry_t);
  zmurge_solvers[id]->hashgraphtab_size = edgenbr;
#endif
  vcsc_init(&zmurge_solvers[id]->vcsc, N, edgenbr, 0, id);

  zmurge_solvers[id]->N        = N;
  zmurge_solvers[id]->edgenbr  = edgenbr;
  zmurge_solvers[id]->cnt      = 0;
  zmurge_solvers[id]->cnt_zero = 0;

#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_NODELST_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_FACTO_OK);

  MURGE_STATE_TRUE(zmurge_solvers[id]->state,  MURGE_GRAPH_BUILD);
  return MURGE_SUCCESS;
}

/*
 * Function: ZMURGE_GraphSetEdge
 *
 * - Check that the number of entries has not been reach for
 * this session.
 * - Increments ROW and COL if baseval is set to 0.
 * - Checks that ROW and COL ranges are corrects.
 * - Adds an entry to the temporary ijv structure.
 *
 * Parameters:
 *   id  - Solver instance identification number.
 *   ROW - Row of the entry.
 *   COL - Column of the entry.
 *
 * Return:
 *   MURGE_ERR_ORDER     - if we are not in a graph building session, or if
 *                         two many edges have been entered, or if
 *                         *zmurge_solvers* or *zmurge_solvers[id]* are not allocated,
 *                         or if *iparm* or *dparm* are not allocated.
 *   MURGE_ERR_PARAMETER - *ROW* or *COL* are out of range or if *id* is not
 *                         in correct range.
 *   MURGE_SUCCESS       - Otherwise
 */
INTS ZMURGE_GraphSetEdge (INTS id, INTS ROW, INTS COL) {
  return ZMURGE_GraphSetEdge_ (id, ROW, COL);
}
INTS ZMURGE_GraphEdge (INTS id, INTS ROW, INTS COL) {
  return ZMURGE_GraphSetEdge_ (id, ROW, COL);
}
static inline
INTS ZMURGE_GraphSetEdge_ (INTS id, INTS ROW, INTS COL) {

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_BUILD))) {
    errorPrint("Need to call ZMURGE_GraphBegin first");
    return MURGE_ERR_ORDER;
  }
  if (zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
    COL += 1;
    ROW += 1;
  }

  if (zmurge_solvers[id]->dropmask != NULL && ROW != COL)
    if (zmurge_solvers[id]->dropmask[(COL-1)] && zmurge_solvers[id]->dropmask[(ROW-1)]) {
      zmurge_solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
    }

  /* If (ROW, COL) has to be droped */
  if ( ( zmurge_solvers[id]->droprows != NULL && ROW != COL &&
         zmurge_solvers[id]->droprows[(ROW-1)] ) ||
       ( zmurge_solvers[id]->dropcols != NULL && ROW != COL &&
         zmurge_solvers[id]->dropcols[(COL-1)] ) ) {
    /* If (COL, ROW) has to be droped */
    if ( ( zmurge_solvers[id]->droprows != NULL && ROW != COL &&
           zmurge_solvers[id]->droprows[(COL-1)] ) ||
         ( zmurge_solvers[id]->dropcols != NULL && ROW != COL &&
           zmurge_solvers[id]->dropcols[(ROW-1)] ) ) {
      zmurge_solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
    } else {
      /* we don't want to break symmetry */
    }
  }

  if (ROW < 1 || COL < 1 || ROW > zmurge_solvers[id]->N || COL > zmurge_solvers[id]->N) {
    errorPrint("ROW %ld or COL %ld is out of range [%ld-%ld]",
               (long)ROW, (long)COL, (long)1, (long)zmurge_solvers[id]->N);
    return MURGE_ERR_PARAMETER;
  }

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif
  if (zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero + 1 > zmurge_solvers[id]->edgenbr) {
    if (zmurge_solvers[id]->dynamic == API_NO) {
      errorPrint("Too many edges in graph description (%ld > %ld)",
                 zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero + 1, zmurge_solvers[id]->edgenbr);
#ifdef MURGE_THREADSAFE
      pthread_mutex_unlock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif

      return MURGE_ERR_ORDER;
    }
  }
  COEF * VAL = NULL;
  vcsc_add_node(zmurge_solvers[id]->vcsc, COL, ROW, VAL, MURGE_ASSEMBLY_ADD, id);
  zmurge_solvers[id]->cnt++;
#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif
  return MURGE_SUCCESS;
}

INTS ZMURGE_GraphSetBlockEdges(INTS id, INTS nROW, INTS *ROWlist,
                              INTS nCOL, INTS *COLlist) {
  INTS r, c, ret;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  for (r = 0; r < nROW; r++)
    for (c = 0; c < nCOL; c++)
      if (MURGE_SUCCESS !=
          (ret = ZMURGE_GraphSetEdge_(id, ROWlist[r], COLlist[c])))
        return ret;
  return MURGE_SUCCESS;
}


/*
 * Function: ZMURGE_GraphEnd
 *
 * - Sort temporary IJV structure with cols as key.
 * - Distribute columns onto processors.
 * (first column on first proc and so on...)
 * - Build a distributed CSC that will be given to PaStiX.
 *
 * TODO:
 * - In the case of a triangular matrix, count each extra-diagonal twice.
 * - Use initial distribution to compute column distribution,
 * in order to reduce communications.
 *
 * Parameters :
 *   id  - Solver instance identification number.
 *
 * Returns:
 *   MURGE_ERR_ORDER     - if we are not in a graph building session, or if
 *                         all edges have not been entered, or if
 *                         *zmurge_solvers* or *zmurge_solvers[id]* are not allocated,
 *                         or if *iparm* or *dparm* are not allocated.
 *   MURGE_ERR_PARAMETER - *ROW* or *COL* are out of range or if *id* is not
 *                         in correct range.
 *   MURGE_SUCCESS       - Otherwise
 */
INTS ZMURGE_GraphEnd  (INTS id) {
  z_pastix_data_t   *pastix_data = zmurge_solvers[id]->pastix_data;
#ifdef MURGE_TIME
  Clock            clock;
#endif
  CLOCK_INIT;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

#ifdef CENTRALISED
  pastix_data->iparm[IPARM_GRAPHDIST] = API_NO;
#endif

  /*
   * Checking that the function is called at the right time
   */
  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_BUILD))) {
    errorPrint("Need to call ZMURGE_GraphBegin first");
    return MURGE_ERR_ORDER;
  }

  if (zmurge_solvers[id]->cnt_zero != 0) {
    if (zmurge_solvers[id]->pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) {
      fprintf(stdout,
              "%ld (%.2g %%) entries were skipped on proc %ld\n",
              (long)zmurge_solvers[id]->cnt_zero,
              (double)((double)100.0*((double)zmurge_solvers[id]->cnt_zero)/
                       ((double)(zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero))),
              (long)zmurge_solvers[id]->pastix_data->procnum);
    }
    else
      {
        if (zmurge_solvers[id]->pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NO) {
          pastix_int_t nz_glob;
          pastix_int_t zeros_glob;
          pastix_int_t nz = zmurge_solvers[id]->colptr[zmurge_solvers[id]->n]-1;
          MPI_Reduce( &(zmurge_solvers[id]->cnt_zero), &zeros_glob,
                      1, PASTIX_MPI_INT,
                      MPI_SUM, 0, zmurge_solvers[id]->pastix_data->pastix_comm);
          MPI_Reduce( &nz, &nz_glob,
                      1, PASTIX_MPI_INT,
                      MPI_SUM, 0, zmurge_solvers[id]->pastix_data->pastix_comm);
          if (zmurge_solvers[id]->pastix_data->procnum == 0) {
            fprintf(stdout,
                    "%ld entries were skipped"
                    " (from %ld (%.3g%%))\n",
                    (long int)zeros_glob,
                    (long int)nz_glob,
                    100.0*((double)(zeros_glob)/
                           ((double)(nz_glob))));
          }

        }
      }
  }

  if (zmurge_solvers[id]->dynamic == API_NO &&
      zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero < zmurge_solvers[id]->edgenbr) {
    errorPrint("Missing edges entry, expected %ld, entered %ld",
               (long)zmurge_solvers[id]->edgenbr,
               (long)(zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero));
    return MURGE_ERR_ORDER;
  }

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif
  if (zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero > zmurge_solvers[id]->edgenbr) {
    errorPrint("Too many edges entry, expected %ld, entered %ld (%ld dropped)",
               (long)zmurge_solvers[id]->edgenbr,
               (long)(zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero),
               (long)zmurge_solvers[id]->cnt_zero);
#ifdef MURGE_THREADSAFE
    pthread_mutex_unlock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif
    return MURGE_ERR_ORDER;
  }
  vcsc_to_cscd(zmurge_solvers[id]->vcsc, pastix_data->pastix_comm,
               &(zmurge_solvers[id]->n), &(zmurge_solvers[id]->colptr),
               &(zmurge_solvers[id]->rows), NULL,
               &(zmurge_solvers[id]->l2g), &(zmurge_solvers[id]->g2l),
               MURGE_ASSEMBLY_ADD, 0, id);
  vcsc_destroy(zmurge_solvers[id]->vcsc, id);

  MURGE_DUMP_GRAPH;

  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_GRAPH_BUILD);

#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif
  CLOCK_PRINT("ZMURGE_GraphEnd");

  return MURGE_SUCCESS;
}

/*
 * Function: ZMURGE_GraphGlobalCSR
 *
 * Enter the adjency graph in a Column Sparse Row form.
 *
 *
 * If the matrix is symmetric, calls <ZMURGE_GraphGlobalCSC>
 * else uses <ZMURGE_GraphBegin>, <ZMURGE_GraphEdge>,
 * <ZMURGE_GraphEnd> sequence.
 *
 *
 * Parameters:
 *   id     - Solver instance identification number.
 *   N      - Number of rows in the CSR.
 *   rowptr - Indexes of each row in COLS array.
 *   root   - Rank of the processor owning the CSR (-1 for all processors)
 */
INTS ZMURGE_GraphGlobalCSR(INTS id, INTS N, INTL *rowptr, INTS *COLS, INTS root) {
  pastix_int_t  iter;
  pastix_int_t  iter2;
  INTS ret;
  int  baseval;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if ( MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_BUILD) ) {
    errorPrint("Do not call ZMURGE_GraphBegin before");
    return MURGE_ERR_ORDER;
  }

  baseval = zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  /* Si on a un graph symetrique autant faire un ZMURGE_GraphGlobalCSC */
  if (zmurge_solvers[id]->pastix_data->iparm[IPARM_SYM] == API_SYM_YES)
    return ZMURGE_GraphGlobalCSC(id, N, rowptr, COLS, root);


  if (zmurge_solvers[id]->pastix_data->procnum == root || root == -1) {
    if (MURGE_SUCCESS != (ret = ZMURGE_GraphBegin(id, N, rowptr[N]- baseval)))
      return ret;

    for (iter = 0; iter < N; iter++) {
      for (iter2 = rowptr[iter]; iter2 < rowptr[iter+1]; iter2++) {
        ret =
          ZMURGE_GraphEdge(id,
                          iter + baseval,
                          COLS[iter2 - baseval]);
        if (MURGE_SUCCESS != ret)
          return ret;
      }
    }
    if (MURGE_SUCCESS != (ret = ZMURGE_GraphEnd(id)))
      return ret;
  }
  else
    {
      if (MURGE_SUCCESS != (ret = ZMURGE_GraphBegin(id, N, 0)))
        return ret;
      if (MURGE_SUCCESS != (ret = ZMURGE_GraphEnd(id)))
        return ret;
    }

  MURGE_STATE_TRUE (zmurge_solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_NODELST_OK);

  return MURGE_SUCCESS;
}

/*
 * Function: ZMURGE_GraphGlobalCSC
 *
 * Distribute the CSC on the processors and use it for PaStiX calls.
 *
 * Parameters:
 *   id     - Solver instance identification number.
 *   N      - Number of columns in the CSR.
 *   colptr - Indexes of each columns in ROWS array.
 *   root   - Rank of the processor owning the CSR (-1 for all processors)
 */
INTS ZMURGE_GraphGlobalCSC(INTS id, INTS N, INTL *colptr, INTS *ROWS, INTS root) {
  pastix_int_t          globaledgenbr;
  pastix_int_t          averageedgenbr;
  pastix_int_t          localedgenbr;
  pastix_int_t          firstcol;
  pastix_int_t          lastcol         = 0;
  pastix_int_t          iter;
  pastix_int_t          procnum;
  pastix_int_t          firstlast[2];
  INTS        *tmpj            = NULL;
  MPI_Request *requests_fl     = NULL;
  MPI_Request *requests_colptr = NULL;
  MPI_Request *requests_rows   = NULL;
  MPI_Status   status;
  int          baseval;             /* User baseval               */

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if ( MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_BUILD) ) {
    errorPrint("Do not call ZMURGE_GraphBegin before");
    return MURGE_ERR_ORDER;
  }

  zmurge_solvers[id]->N = N;
  baseval        = zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL];

  if (root == zmurge_solvers[id]->pastix_data->procnum || root == -1) {
    globaledgenbr   = colptr[N] - baseval;
    averageedgenbr  = globaledgenbr / zmurge_solvers[id]->pastix_data->procnbr;
    firstcol        = 0;
    if (root != -1) {
      MURGE_MEMALLOC(requests_fl,
                     zmurge_solvers[id]->pastix_data->procnbr,
                     MPI_Request);
      MURGE_MEMALLOC(requests_colptr,
                     zmurge_solvers[id]->pastix_data->procnbr,
                     MPI_Request);
      MURGE_MEMALLOC(requests_rows,
                     zmurge_solvers[id]->pastix_data->procnbr,
                     MPI_Request);
    }
    /* Pour chaque prrocesseur,
     on attribue au processeur un certain nombre de colonnes et donc
     d'arrêtes

     Si le processeur est local on construit le loc2glob
     On construit le colptr
     On copie rows

     Sinon, on envoi les numéros de première et dernière colonnes et
     les morceaux de colptr et rows correspondants.

     */
    for (procnum = 0; procnum <  zmurge_solvers[id]->pastix_data->procnbr; procnum++) {
      while (lastcol < N - 1&&
             colptr[lastcol+1] - colptr[firstcol]< averageedgenbr)
        lastcol++;
      localedgenbr =  colptr[lastcol+1] - colptr[firstcol];

      if (procnum == zmurge_solvers[id]->pastix_data->procnum) {
        zmurge_solvers[id]->n = lastcol-firstcol+1;

        MURGE_MEMALLOC(zmurge_solvers[id]->l2g, zmurge_solvers[id]->n, pastix_int_t);

        for (iter = 0; iter < zmurge_solvers[id]->n; iter++) {
          zmurge_solvers[id]->l2g[iter] = firstcol+iter+1;
        }

        MURGE_MEMALLOC(zmurge_solvers[id]->colptr, zmurge_solvers[id]->n+1, pastix_int_t);

        for (iter = 0; iter < zmurge_solvers[id]->n+1; iter++) {
          zmurge_solvers[id]->colptr[iter] =
            colptr[firstcol+iter] - colptr[firstcol]+1;
        }

        MURGE_MEMALLOC(zmurge_solvers[id]->rows, localedgenbr, pastix_int_t);

        for (iter = 0; iter < localedgenbr; iter++)
          zmurge_solvers[id]->rows[iter] = ROWS[colptr[firstcol]+iter-1];

      }
      else
        {
          if (root != -1) {
            firstlast[0] = firstcol;
            firstlast[1] = lastcol;


            MPI_Isend(firstlast,
                      2, PASTIX_MPI_INT, procnum,
                      TAG_FL, zmurge_solvers[id]->pastix_data->pastix_comm,
                      &requests_fl[procnum]);

            MPI_Isend(&colptr[firstcol],
                      lastcol-firstcol+2,
                      PASTIX_MPI_INT, procnum,
                      TAG_COL, zmurge_solvers[id]->pastix_data->pastix_comm,
                      &requests_colptr[procnum]);

            MPI_Isend(&ROWS[colptr[firstcol]],
                      localedgenbr*sizeof(INTS),
                      MPI_BYTE, procnum,
                      TAG_ROW, zmurge_solvers[id]->pastix_data->pastix_comm,
                      &requests_rows[procnum]);

          }
        }
      firstcol = lastcol + 1;
    }
    if (root != -1) {
      for (procnum = 0;
           procnum < zmurge_solvers[id]->pastix_data->procnbr;
           procnum++) {
        if (procnum != zmurge_solvers[id]->pastix_data->procnum) {
          MPI_Wait(&requests_fl[procnum], &status);
          MPI_Wait(&requests_colptr[procnum], &status);
          MPI_Wait(&requests_rows[procnum], &status);
        }
      }
    }
    MURGE_FREE(requests_rows);
    MURGE_FREE(requests_colptr);
    MURGE_FREE(requests_fl);
  }
  else
    {
      /* Si on est pas le processeur racine

       On recoit les numeros de première et dernière colonnes
       On en déduit le loca2glob
       On recoit les parties locales de colptr et rows
       On construit le colptr local et rows local.
       */
      MPI_Recv(firstlast,
               2, PASTIX_MPI_INT, root,
               TAG_FL, zmurge_solvers[id]->pastix_data->pastix_comm,
               &status);
      firstcol = firstlast[0];
      lastcol  = firstlast[1];

      zmurge_solvers[id]->n = lastcol-firstcol+1;

      MURGE_MEMALLOC(zmurge_solvers[id]->l2g, lastcol-firstcol+1, pastix_int_t);

      for (iter = 0; iter < lastcol-firstcol+1; iter++) {
        zmurge_solvers[id]->l2g[iter] = firstcol+iter;
      }

      MURGE_MEMALLOC(zmurge_solvers[id]->colptr, lastcol-firstcol+2, pastix_int_t);

      MPI_Recv(zmurge_solvers[id]->colptr,
               lastcol-firstcol+2, PASTIX_MPI_INT, root,
               TAG_COL, zmurge_solvers[id]->pastix_data->pastix_comm,
               &status);


      for (iter = 0; iter < lastcol-firstcol+2; iter++) {
        zmurge_solvers[id]->colptr[lastcol-firstcol+1 - iter] -=
          zmurge_solvers[id]->colptr[0];
      }

      localedgenbr = zmurge_solvers[id]->colptr[lastcol-firstcol+1]-1;

      MURGE_MEMALLOC(tmpj, localedgenbr, INTS);

      MPI_Recv(tmpj,
               localedgenbr*sizeof(INTS), MPI_BYTE, root,
               TAG_ROW, zmurge_solvers[id]->pastix_data->pastix_comm,
               &status);

      if (sizeof(INTS) == sizeof(pastix_int_t)) {
        zmurge_solvers[id]->rows = (pastix_int_t*)tmpj;
      }
      else
        {
          MURGE_MEMALLOC(zmurge_solvers[id]->rows, localedgenbr, pastix_int_t);

          for (iter = 0; iter < localedgenbr; iter++)
            zmurge_solvers[id]->rows[iter] = (pastix_int_t)tmpj[iter];
          MURGE_FREE(tmpj);
        }
    }

  MURGE_DUMP_GRAPH;

  MURGE_STATE_TRUE (zmurge_solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_NODELST_OK);

  return MURGE_SUCCESS;
}
/*
 * Function: ZMURGE_GraphGlobalIJV
 *
 * Distribute the graph on the processors, compress the columns
 * array and use the built CSCd to call PaStiX.
 *
 *
 * Parameters:
 *   id     - Solver instance identification number.
 *   N      - Number of columns in the CSR.
 *   NNZ    - Number of non-zeros in the matrix.
 *   ROWS   - Rows array.
 *   COLS   - Columns array.
 *   root   - Rank of the processor owning the CSR (-1 for all processors)
 */
INTS ZMURGE_GraphGlobalIJV(INTS id, INTS N, INTL NNZ, INTS *ROWS,
                          INTS *COLS, INTS root) {
  pastix_int_t       *localn     = NULL;   /* Number of local column on each proc  */
  pastix_int_t       *localedges = NULL;   /* Number of local edges on each proc   */
  pastix_int_t        lnnz;                /* Local number of edges                */
  pastix_int_t        sizes[2];            /* Array to send n and nnz to each proc */
  INTS      *tmprows  = NULL;     /* Temporary local rows tabular         */
  INTS      *tmpcols  = NULL;     /* Temporary local columns tabular      */
  pastix_int_t       *sizecols = NULL;     /* Number of rows in each column        */
  pastix_int_t        totaledgenbr;        /* Total number of edges                */
  pastix_int_t        avredgenbr;          /* Average number of edges              */
  int        baseval_int     = 1; /* Internal baseval, always 1           */
  int        baseval;             /* User baseval                         */
  pastix_int_t        iter, iter2, index;  /* Iterators                            */
  pastix_int_t       *coldist  = NULL;     /* Owner of each column                 */
  int        procnum;             /* Processor number iterator            */
  MPI_Status status;              /* MPI status */
  INTS       currentcol;
  void      *sortptr[2];

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if ( MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_BUILD) ) {
    errorPrint("Do not call ZMURGE_GraphBegin before");
    return MURGE_ERR_ORDER;
  }


  baseval = zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL];


  if ((zmurge_solvers[id]->pastix_data->procnum == root) || (root == -1)) {
    zmurge_solvers[id]->N = N;
  }
  if (root != -1) {
    MPI_Bcast(&zmurge_solvers[id]->N, 1, PASTIX_MPI_INT, root,
              zmurge_solvers[id]->pastix_data->pastix_comm);
    MPI_Bcast(&NNZ, sizeof(INTL), MPI_BYTE, root,
              zmurge_solvers[id]->pastix_data->pastix_comm);
  }
  if ((zmurge_solvers[id]->pastix_data->procnum == root) || (root == -1)) {
    /* Sort Col and Rows with first key, col */
    sortptr[0] = COLS;
    sortptr[1] = ROWS;
#ifdef INTSSIZE64
    qsort2IntAsc(sortptr, NNZ);
#else
    qsort2SmallIntAsc(sortptr, NNZ);
#endif
    /* decide how to distribute the graph */
    MURGE_MEMALLOC(sizecols, zmurge_solvers[id]->N, pastix_int_t);
    memset(sizecols, 0, zmurge_solvers[id]->N*sizeof(pastix_int_t));
    totaledgenbr = 0;
    /* Count how long is each column */
    for (iter = 0; iter < NNZ; iter ++) {
      /* TODO: Dans le cas ou la matrice donnée ne contient que la partie
       triangulaire, il faudra  compter 2 fois les elements
       non diagonaux.
       */
      sizecols[COLS[iter] - 1]++;
      totaledgenbr++;
    }

    avredgenbr = totaledgenbr/zmurge_solvers[id]->pastix_data->procnbr;

    MURGE_MEMALLOC(coldist, zmurge_solvers[id]->N, pastix_int_t);

    for (iter = 0; iter < zmurge_solvers[id]->N; iter++)
      coldist[iter] = zmurge_solvers[id]->pastix_data->procnbr - 1;

    procnum    = 0;
    iter       = 0;

    MURGE_MEMALLOC(localedges, zmurge_solvers[id]->pastix_data->procnbr, pastix_int_t);
    MURGE_MEMALLOC(localn,     zmurge_solvers[id]->pastix_data->procnbr, pastix_int_t);
    memset(localedges, 0, zmurge_solvers[id]->pastix_data->procnbr*sizeof(pastix_int_t));
    memset(localn, 0, zmurge_solvers[id]->pastix_data->procnbr*sizeof(pastix_int_t));

    while (iter < zmurge_solvers[id]->N) {
      localedges[procnum] = 0;
      while (iter < zmurge_solvers[id]->N &&
             (localedges[procnum] < avredgenbr||
              (procnum == zmurge_solvers[id]->pastix_data->procnbr -1))) {
        coldist[iter] = procnum;
        localn[procnum]++;
        localedges[procnum] +=  sizecols[iter];
        iter ++;
      }
      procnum++;
    }

    MURGE_FREE(coldist);

    /* Send data to each processor */

    for (index = 0, procnum = 0;
         procnum < zmurge_solvers[id]->pastix_data->procnbr;
         procnum++) {
      if (procnum != zmurge_solvers[id]->pastix_data->procnum) {
        if (root != -1) {
          sizes[0] = localn[procnum];
          sizes[1] = localedges[procnum];

          /* envoi du nombre de non zeros */
          MPI_Send(sizes,
                   2, PASTIX_MPI_INT, procnum,
                   TAG_SIZE, zmurge_solvers[id]->pastix_data->pastix_comm);
          /* envoi des lignes */
          MPI_Send(&(ROWS[index]),
                   localedges[procnum]*sizeof(INTS), MPI_BYTE, procnum,
                   TAG_ROW, zmurge_solvers[id]->pastix_data->pastix_comm);
          /* envoi des colonnes */
          MPI_Send(&(COLS[index]),
                   localedges[procnum]*sizeof(INTS), MPI_BYTE, procnum,
                   TAG_COL, zmurge_solvers[id]->pastix_data->pastix_comm);
        }
      }
      else
        {
          tmprows = &(ROWS[index]);
          tmpcols = &(COLS[index]);
        }
      index += localedges[procnum];
    }

    zmurge_solvers[id]->n = localn[zmurge_solvers[id]->pastix_data->procnum];
    lnnz = localedges[zmurge_solvers[id]->pastix_data->procnum];
    MURGE_FREE(localn);
    MURGE_FREE(localedges);
  }
  else
    {
      MPI_Recv(sizes,
               2, PASTIX_MPI_INT, root,
               TAG_SIZE, zmurge_solvers[id]->pastix_data->pastix_comm,
               &status);
      zmurge_solvers[id]->n = sizes[0];
      lnnz           = sizes[1];

      MURGE_MEMALLOC(tmprows, lnnz, INTS);
      MPI_Recv(tmprows,
               lnnz*sizeof(INTS), MPI_BYTE, root,
               TAG_ROW, zmurge_solvers[id]->pastix_data->pastix_comm,
               &status);
      MURGE_MEMALLOC(tmpcols, lnnz, INTS);
      MPI_Recv(tmpcols,
               lnnz*sizeof(INTS), MPI_BYTE, root,
               TAG_COL, zmurge_solvers[id]->pastix_data->pastix_comm,
               &status);
    }

  MURGE_MEMALLOC(zmurge_solvers[id]->colptr, zmurge_solvers[id]->n+1, pastix_int_t);
  MURGE_MEMALLOC(zmurge_solvers[id]->l2g,    zmurge_solvers[id]->n  , pastix_int_t);
  /* convert tmpcols/tmprows to CSCd */
  iter2=baseval_int;
  for (iter=0; iter<(zmurge_solvers[id]->n); iter++) {
    zmurge_solvers[id]->colptr[iter] = iter2;
    zmurge_solvers[id]->l2g[iter]    = tmpcols[iter2-baseval_int]
      - baseval + baseval_int;
    currentcol = tmpcols[iter2-baseval_int];
    while (((iter2-baseval) < lnnz) &&
           ((tmpcols[iter2-baseval_int]) == (currentcol))) {
      iter2++;
    }
  }

  if ((zmurge_solvers[id]->pastix_data->procnum != root) && (root != -1))
    MURGE_FREE(tmpcols);

  zmurge_solvers[id]->colptr[zmurge_solvers[id]->n] = iter2;

  if (iter2 != lnnz+baseval) {
    errorPrint("Mauvais nombre d'arrête");
    return MURGE_ERR_PARAMETER;
  }


  MURGE_MEMALLOC(zmurge_solvers[id]->rows, lnnz, pastix_int_t);

  for (iter=0; iter<lnnz; iter++) {
    zmurge_solvers[id]->rows[iter] = tmprows[iter];
  }
  if (!((zmurge_solvers[id]->pastix_data->procnum == root) || (root == -1)))
    MURGE_FREE(tmprows);

  MURGE_STATE_TRUE (zmurge_solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_NODELST_OK);
  return MURGE_SUCCESS;
}

INTS ZMURGE_SetOrdering(INTS id, INTS * permutation) {
  INTS i;
  print_debug(DBG_MURGE, ">> ZMURGE_SetOrdering\n");
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("Build graph before calling ZMURGE_SetOrdering");
    return MURGE_ERR_ORDER;
  }
  if (zmurge_solvers[id]->l2g == NULL) {
    errorPrint("Local to global array is not set");
    return MURGE_ERR_PARAMETER;
  }
  if (permutation == NULL) {
    errorPrint("NULL parameter");
    return MURGE_ERR_PARAMETER;
  }
  MURGE_MEMALLOC(zmurge_solvers[id]->perm, zmurge_solvers[id]->n, pastix_int_t);
  memset(zmurge_solvers[id]->perm, 0, zmurge_solvers[id]->n*sizeof(pastix_int_t));
  for (i = 0; i < zmurge_solvers[id]->n; i++)
    zmurge_solvers[id]->perm[i] = permutation[zmurge_solvers[id]->l2g[i]-1];

  zmurge_solvers[id]->pastix_data->iparm[IPARM_ORDERING] = API_ORDER_PERSONAL;
  zmurge_solvers[id]->pastix_data->iparm[IPARM_LEVEL_OF_FILL] = -1;

  print_debug(DBG_MURGE, "<< ZMURGE_SetOrdering\n");
  return MURGE_SUCCESS;
}

/*******************************************************************************
 * Group: Matrix assembly functions
 */

/* Function: ZMURGE_AssemblySetSequence
 *
 * Create a sequence of entries to build a matrix and store it for being reused.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   coefnbr - Number of entries.
 *   ROWs    - List of rows in the sequence.
 *   COLs    - List of columns in the sequence.
 *   op      - Operation to perform for coefficient which appear
 *             several tim (see <MURGE_ASSEMBLY_OP>).
 *   op2     - Operation to perform when a coefficient is set by
 *             two different processors (see <MURGE_ASSEMBLY_OP>).
 *   mode    - Indicates if user ensure he will respect zmurge_solvers distribution
 *             (see <MURGE_ASSEMBLY_MODE>).
 *   nodes   - 0 entries are entered value by value,
 *             1 entries are entries node by node.
 *   id_seq  - Sequence ID.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 *                         *op*, *mode*, *sym*, or *coefnbr* are not valid.
 */
int ZMURGE_AssemblySetSequence (INTS id, INTL coefnbr, INTS * ROWs, INTS * COLs,
                               INTS op, INTS op2, INTS mode, INTS nodes,
                               INTS * id_seq) {
  z_murge_seq_t * sequence     = NULL;
  INTL         iter;
  ijv_t      **send_ijv      = NULL;
  INTL        *send_ijv_size = NULL;
  INTL        *send_nbr      = NULL;
  INTS         dof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("Graph has to be built before");
    return MURGE_ERR_ORDER;
  }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);


  MURGE_MEMALLOC(sequence, 1, z_murge_seq_t);

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (zmurge_solvers[id]->sequences == NULL) {
    zmurge_solvers[id]->sequences = sequence;
  }
  else
    {
      z_murge_seq_t *last_sequence = zmurge_solvers[id]->sequences;
      while(sequence->next != NULL) {
        last_sequence = sequence->next;
      }
      last_sequence->next = sequence;
    }

  sequence->next         = NULL;
  sequence->indexes      = NULL;
  sequence->recv_indexes = NULL;
  sequence->recv_nbr     = NULL;
  sequence->op_local_entries = op;
  sequence->op_dist_entries  = op;
  sequence->mode  = mode;
  sequence->nodes = (nodes==MURGE_BOOLEAN_FALSE)?API_NO:API_YES;
  sequence->ID = zmurge_solvers[id]->seq_ID++;
  *id_seq = sequence->ID;

  sequence->coefnbr = coefnbr;
  MURGE_MEMALLOC(sequence->indexes, coefnbr, INTL);
  if (sequence->mode == MURGE_ASSEMBLY_FOOL) {
    MURGE_MEMALLOC(send_nbr, zmurge_solvers[id]->pastix_data->procnbr, INTL);
    MURGE_MEMALLOC(sequence->recv_nbr,
                   zmurge_solvers[id]->pastix_data->procnbr,
                   INTL);
    MURGE_MEMALLOC(sequence->recv_indexes,
                   zmurge_solvers[id]->pastix_data->procnbr,
                   INTL*);
    MURGE_MEMALLOC(send_ijv, zmurge_solvers[id]->pastix_data->procnbr, ijv_t*);
    MURGE_MEMALLOC(send_ijv_size, zmurge_solvers[id]->pastix_data->procnbr, INTL);

    for (iter = 0; iter < zmurge_solvers[id]->pastix_data->procnbr; iter++) {
      send_nbr[iter] = 0;
      send_ijv_size[iter] = 1+coefnbr/10;
      MURGE_MEMALLOC(send_ijv[iter], send_ijv_size[iter], ijv_t);
      sequence->recv_indexes[iter] = NULL;
    }
  }

  if (zmurge_solvers[id]->colptr == NULL) {
    /* Need to build a CSC */
    INTS inc = (zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL]==0)?1:0;

    vcsc_init(&zmurge_solvers[id]->vcsc, zmurge_solvers[id]->N, coefnbr/dof, 0, id);
    for (iter = 0; iter < coefnbr; iter++) {
      COEF * VAL = NULL;
      vcsc_add_node(zmurge_solvers[id]->vcsc, COLs[iter]+inc, ROWs[iter]+inc, VAL,
                    MURGE_ASSEMBLY_ADD, id);
    }
    vcsc_to_cscd(zmurge_solvers[id]->vcsc, zmurge_solvers[id]->pastix_data->pastix_comm,
                 &(zmurge_solvers[id]->n), &(zmurge_solvers[id]->colptr),
                 &(zmurge_solvers[id]->rows), NULL,
                 &(zmurge_solvers[id]->l2g), &(zmurge_solvers[id]->g2l),
                 MURGE_ASSEMBLY_ADD, 0, id);
  }

  coefnbr = sequence->coefnbr;

  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_ONLY_PROD)) {
    /* TODO: CHECK That all entries exists in CSC and if not insert it*/
  }

  for (iter = 0; iter < coefnbr; iter++) {
    INTL iter2;
    INTS inc = (zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL]==0)?1:0;
    INTS col = COLs[iter]+inc; /* 1 based */
    INTS row = ROWs[iter]+inc; /* 1 based */
    INTS node_col; /* 0 based */
    INTS node_row; /* 0 based */
    INTS in_node_col; /* 0 based */
    INTS in_node_row; /* 0 based */
    INTS node_col_loc; /* 1 based */

    if (dof > 1 && sequence->nodes == API_NO) {
      node_col     = (col-1 - (col-1)%dof)/dof;
      in_node_col  = (col-1)%dof;
      node_row     = (row-1 - (row-1)%dof)/dof;
      in_node_row  = (row-1)%dof;
    }
    else
      {
        node_col     = col-1;
        in_node_col  = 0;
        node_row     = row-1;
        in_node_row  = 0;
      }

    node_col_loc = zmurge_solvers[id]->g2l[node_col];
    if ( node_col_loc > 0 ) {
      node_col_loc--;
      /* Entry is local */
      for (iter2 = zmurge_solvers[id]->colptr[node_col_loc]-1;
           iter2 < zmurge_solvers[id]->colptr[node_col_loc+1]-1;
           iter2++) {
        if (zmurge_solvers[id]->rows[iter2] == row)
          break;
      }
      if (zmurge_solvers[id]->colptr[node_col_loc+1]-1 == iter2) {
        /* Entry not found in CSC */
        errorPrint("ROW (%ld:%ld) not found in COL (%d:%d) %d",
                   (long)row, node_row, (long)col, node_col,node_col_loc);
        return MURGE_ERR_PARAMETER;
      }
      else
        {
          if (iter2*dof*dof + in_node_col*dof+in_node_row >
              dof*dof*(zmurge_solvers[id]->colptr[zmurge_solvers[id]->n]-1)) {
            return MURGE_ERR_PARAMETER;
          }

          sequence->indexes[iter] = iter2*dof*dof + in_node_col*dof+in_node_row;
        }
    } else {
      /* Entry not local */
      if (sequence->mode == MURGE_ASSEMBLY_RESPECT) {
        errorPrint("COL (%d) is not local (row %d, owner %d)",
                   (long)col+1, (long)row+1, -node_col_loc);
        return MURGE_ERR_PARAMETER;
      } else {
        int owner = -node_col_loc;
        sequence->indexes[iter] = -(owner+1);

        /* => send buffer */
        if (send_nbr[owner] == send_ijv_size[owner]) {
          send_ijv_size[owner] = (INTL)(1.5*send_ijv_size[owner] + 1);
          MURGE_REALLOC(send_ijv[owner], send_ijv_size[owner], ijv_t);
        }
        send_ijv[owner][send_nbr[owner]].i = row;
        send_ijv[owner][send_nbr[owner]].j = col;
        send_nbr[owner]++;
      }
    }
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
    if (node_col == (MURGE_TRACE_COL-1)/dof && node_row == (MURGE_TRACE_ROW-1)/dof)
      fprintf(stdout, "%d, %d index => %d, iter %d node_col_loc %d send_nbr %d\n", node_col, node_row, sequence->indexes[iter], iter, node_col_loc, (node_col_loc>0)?-1:send_nbr[-node_col_loc] );
#endif
  }

  if (sequence->mode != MURGE_ASSEMBLY_RESPECT) {
    MPI_Request * requests;
    ijv_t       * recv_ijv;
    int           size;
    int           lastsize;
    pastix_int_t           iter_coef;

    MPI_Alltoall(send_nbr,           1, MPI_INTL,
                 sequence->recv_nbr, 1, MPI_INTL,
                 zmurge_solvers[id]->pastix_data->pastix_comm);

    MURGE_MEMALLOC(requests, zmurge_solvers[id]->pastix_data->procnbr, MPI_Request);
    for (iter = 0; iter < zmurge_solvers[id]->pastix_data->procnbr; iter++) {
      if (send_nbr[iter] > 0)
        MPI_Isend(send_ijv[iter], send_nbr[iter]*sizeof(ijv_t), MPI_BYTE,
                  iter, TAG_IJV, zmurge_solvers[id]->pastix_data->pastix_comm,
                  &(requests[iter]));
    }

    lastsize = sequence->recv_nbr[0];
    MURGE_MEMALLOC(recv_ijv, lastsize, ijv_t);
    for (iter = 0; iter < zmurge_solvers[id]->pastix_data->procnbr; iter++) {
      MPI_Status status;
      size = sequence->recv_nbr[iter];
      MURGE_MEMALLOC(sequence->recv_indexes[iter], size, INTL);
      if (lastsize < size) {
        MURGE_REALLOC(recv_ijv, size, ijv_t);
        lastsize = size;
      }
      if (size > 0)
        MPI_Recv(recv_ijv, size*sizeof(ijv_t), MPI_BYTE,
                 iter, TAG_IJV, zmurge_solvers[id]->pastix_data->pastix_comm,
                 &status);

      for (iter_coef = 0; iter_coef < size; iter_coef++) {
        INTL iter2;
        INTS col = recv_ijv[iter_coef].j;
        INTS row = recv_ijv[iter_coef].i;
        INTS node_col;
        INTS node_row;
        INTS in_node_col;
        INTS in_node_row;
        INTS node_col_loc;


        if (dof > 1 && sequence->nodes == API_NO) {
          node_col = (col-1-(col-1)%dof)/dof;
          node_row = (row-1-(row-1)%dof)/dof;
          in_node_row = (row-1)%dof;
          in_node_col = (col-1)%col;
        }
        else
          {
            node_col = col-1;
            node_row = row-1;
            in_node_row = 0;
            in_node_col = 0;
          }
        node_col_loc = zmurge_solvers[id]->g2l[node_col];

        if ( node_col_loc > 0 ) {
          /* Entry is local */
          for (iter2 = zmurge_solvers[id]->colptr[node_col_loc-1]-1;
               iter2 < zmurge_solvers[id]->colptr[node_col_loc]-1;
               iter2++) {
            if (zmurge_solvers[id]->rows[iter2] == node_row+1)
              break;
          }
          if (zmurge_solvers[id]->colptr[node_col_loc]-1 == iter2) {
            /* Entry not found in CSC */
            errorPrint("ROW (%ld) not found in COL (%d) %d",
                       (long)row, (long)col, node_col_loc);

            return MURGE_ERR_PARAMETER;
          }
          else
            {
              sequence->recv_indexes[iter][iter_coef] =
                iter2*dof*dof + in_node_col*dof+in_node_row;
            }
        }
        else
          {
            /* Entry not local */
            errorPrint("%s:%d COL (%d) is not local (row %d, owner %d)",
                       __FILE__, __LINE__, (long)col, row, -node_col_loc);
            return MURGE_ERR_PARAMETER;
          }
      }
    }

    MURGE_FREE(recv_ijv);

    for (iter = 0; iter < zmurge_solvers[id]->pastix_data->procnbr; iter++) {
      MPI_Status status;

      if (send_nbr[iter] > 0)
        MPI_Wait(&(requests[iter]), &status);
      MURGE_FREE(send_ijv[iter]);
    }
    MURGE_FREE(send_ijv);
    MURGE_FREE(send_ijv_size);
    MURGE_FREE(requests);
    MURGE_FREE(send_nbr);
  }

  return MURGE_SUCCESS;
}

static inline
int try_complete_reception(INTS    id,
                           z_murge_seq_t * sequence,
                           int     src,
                           int     send_nbr,
                           int     blocksize,
                           pastix_complex64_t * recv_data,
                           pastix_int_t   * recv_cnt) {
  MPI_Status TCR_status;
  int        TCR_flag = 1;

  while (TCR_flag) {
    /* Do we have something to receive ? */
    pastix_int_t        TCR_iter;
    MPI_Iprobe(src, TAG_VAL,
               zmurge_solvers[id]->pastix_data->pastix_comm,
               &TCR_flag, &TCR_status);
    if (TCR_flag) {
      /* Receive and add it */
      int     TCR_src = TCR_status.MPI_SOURCE;
      int     TCR_tag = TCR_status.MPI_TAG;
      pastix_complex64_t * TCR_vals_ptr;
      MPI_Recv(recv_data, send_nbr*blocksize, COMM_FLOAT,
               TCR_src, TCR_tag,
               zmurge_solvers[id]->pastix_data->pastix_comm, &TCR_status);
      TCR_vals_ptr=recv_data;
      for (TCR_iter = 0; TCR_iter < send_nbr; TCR_iter++) {
        pastix_int_t     TCR_index =
          sequence->recv_indexes[TCR_src][recv_cnt[TCR_src]];
        pastix_complex64_t * TCR_node  = &(zmurge_solvers[id]->values[TCR_index]);
        pastix_int_t     TCR_iterdof;
        recv_cnt[TCR_src]++;
        for (TCR_iterdof = 0; TCR_iterdof < blocksize;
             TCR_iterdof++, TCR_node++, TCR_vals_ptr++)
          MURGE_ADD(*TCR_node,
                    *TCR_vals_ptr,
                    sequence->op_local_entries);
      }
    }
  }
  return MURGE_SUCCESS;
}

/*
 * ZMURGE_AssemblyUseSequence
 *
 * Assembly the matrix using a stored sequence.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   id_seq  - Sequence ID.
 *   values  - Values to insert in the CSC.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 *                         *id_seq* or *values* are not valid.
 */
INTS ZMURGE_AssemblyUseSequence(INTS id, INTS id_seq, COEF * values) {
  z_murge_seq_t * sequence;
  INTL          iter;
  INTS          dof;
  pastix_complex64_t       * send_array = NULL;
  pastix_int_t         * send_cnt   = NULL;
  MPI_Request * send_reqs  = NULL;
  MPI_Request * send_reqs2 = NULL;
  pastix_complex64_t       * recv_data  = NULL;
  pastix_int_t         * recv_cnt   = NULL;
  int           send_nbr   = -1;
  int           blocksize;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  sequence = zmurge_solvers[id]->sequences;
  while(sequence != NULL && sequence->ID != id_seq)
    sequence = sequence->next;

  if (sequence == NULL) {
    errorPrint("Sequence %d not found", id_seq);
    sequence = zmurge_solvers[id]->sequences;
    return MURGE_ERR_PARAMETER;
  }

  if (values == NULL) {
    errorPrint("NULL value Pointer");
    return MURGE_ERR_PARAMETER;
  }

  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("Graph has to be built before");
    return MURGE_ERR_ORDER;
  }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (zmurge_solvers[id]->values == NULL) {
    MURGE_MEMALLOC(zmurge_solvers[id]->values,
                   (zmurge_solvers[id]->colptr[zmurge_solvers[id]->n]-1)*dof*dof,
                   pastix_complex64_t);
    memset(zmurge_solvers[id]->values, 0,
           (zmurge_solvers[id]->colptr[zmurge_solvers[id]->n]-1)*dof*dof*sizeof(pastix_complex64_t));
  }

  if (sequence->nodes == API_YES) {
    blocksize = dof*dof;
  }
  else
    {
      blocksize = 1;
    }

  if (sequence->mode != MURGE_ASSEMBLY_RESPECT) {
    pastix_int_t coefnbr, recv_coefnbr;
    coefnbr = sequence->coefnbr;
    MPI_Allreduce(&coefnbr, &recv_coefnbr, 1, PASTIX_MPI_INT, MPI_SUM,
                  zmurge_solvers[id]->pastix_data->pastix_comm);
    send_nbr = recv_coefnbr/100 + 1;

    MURGE_MEMALLOC(send_reqs, zmurge_solvers[id]->pastix_data->procnbr, MPI_Request);
    MURGE_MEMALLOC(send_reqs2, zmurge_solvers[id]->pastix_data->procnbr, MPI_Request);
    MURGE_MEMALLOC(send_array,
                   blocksize*send_nbr*zmurge_solvers[id]->pastix_data->procnbr, pastix_complex64_t);
    MURGE_MEMALLOC(recv_data, blocksize*send_nbr, pastix_complex64_t);
    MURGE_MEMALLOC(recv_cnt,  zmurge_solvers[id]->pastix_data->procnbr, pastix_int_t);
    MURGE_MEMALLOC(send_cnt, zmurge_solvers[id]->pastix_data->procnbr, pastix_int_t);
    for (iter = 0; iter < zmurge_solvers[id]->pastix_data->procnbr; iter++) {
      recv_cnt[iter] = 0;
      send_cnt[iter] = 0;
    }
  }

  for (iter = 0; iter < sequence->coefnbr; iter++) {
    INTL index;
    pastix_complex64_t*node;
    INTS iterdof;

    if (sequence->mode != MURGE_ASSEMBLY_RESPECT) {
      try_complete_reception(id, sequence,
                             MPI_ANY_SOURCE, send_nbr, blocksize,
                             recv_data, recv_cnt);
    }

    index = sequence->indexes[iter];

    if (index>=0) {
      if (index > dof*dof*(zmurge_solvers[id]->colptr[zmurge_solvers[id]->n]-1)) {
        return -1;
      }
      node = &(zmurge_solvers[id]->values[index]);
      for (iterdof = 0; iterdof < blocksize; iterdof++, node++, values++)
        MURGE_ADD(*node, *values, sequence->op_local_entries);
    } else {
      if (sequence->mode == MURGE_ASSEMBLY_RESPECT) {
        errorPrint("Non local entry incompatible with MURGE_ASSEMBLY_RESPECT");
        return MURGE_ERR_PARAMETER;
      }

      if (send_cnt[-index-1] == send_nbr) {
        MPI_Status status;
        int flag = 0;
        while(!flag) {
          MPI_Test(&(send_reqs[-index-1]), &flag, &status);
          if (!flag) {
            try_complete_reception(id, sequence,
                                   MPI_ANY_SOURCE, send_nbr, blocksize,
                                   recv_data, recv_cnt);
          } else {
            send_cnt[-index-1] = 0;
          }
        }
      }

      /* Prepare to send, if send_buff is full send */
      node = &(send_array[blocksize*(send_cnt[-index-1]+send_nbr*(-index-1))]);
      for (iterdof = 0; iterdof < blocksize; iterdof++, node++, values++)
        *node= *values;

      send_cnt[-index-1]++;

      if (send_cnt[-index-1] == send_nbr) {
        MPI_Isend(&(send_array[send_nbr*(-index-1)*blocksize]),
                  send_cnt[-index-1]*blocksize, COMM_FLOAT,
                  -index-1, TAG_VAL,
                  zmurge_solvers[id]->pastix_data->pastix_comm,
                  &(send_reqs[-index-1]));
      }
    }
  }

  if (sequence->mode != MURGE_ASSEMBLY_RESPECT) {
    pastix_int_t * done, done_cnt;
    MURGE_MEMALLOC(done, zmurge_solvers[id]->pastix_data->procnbr, pastix_int_t);
    for (iter = 0; iter < zmurge_solvers[id]->pastix_data->procnbr; iter++)
      done[iter] = API_NO;
    done_cnt = 0;
    for (iter = 0; iter < zmurge_solvers[id]->pastix_data->procnbr; iter++) {
      /* Wait last sends */
      if (send_cnt[iter] == send_nbr) {
        MPI_Status status;
        MPI_Wait(&(send_reqs[iter]), &status);
        send_cnt[iter] = 0;
      }

      /* Send last entries */
      MPI_Isend(&(send_cnt[iter]), 1, PASTIX_MPI_INT,
                iter, TAG_SIZE,
                zmurge_solvers[id]->pastix_data->pastix_comm,
                &(send_reqs[iter]));
      MPI_Isend(&(send_array[blocksize*send_nbr*(iter)]),
                send_cnt[iter]*blocksize, COMM_FLOAT,
                iter, TAG_VAL2,
                zmurge_solvers[id]->pastix_data->pastix_comm,
                &(send_reqs2[iter]));
    }

    while (done_cnt < zmurge_solvers[id]->pastix_data->procnbr) {
      for (iter =0; iter < zmurge_solvers[id]->pastix_data->procnbr; iter++) {
        if (done[iter] == API_NO) {
          try_complete_reception(id, sequence,
                                 iter, send_nbr, blocksize,
                                 recv_data, recv_cnt);

          /* recv last count /entries */
          MPI_Status status;
          pastix_int_t        cnt;
          pastix_int_t        iter2;
          pastix_complex64_t     *myvals_ptr;
          int       flag;
          /* Receive and add it */
          MPI_Iprobe(iter, TAG_SIZE,
                     zmurge_solvers[id]->pastix_data->pastix_comm, &flag, &status);
          if (flag) {
            MPI_Recv(&cnt, 1, PASTIX_MPI_INT,
                     iter, TAG_SIZE,
                     zmurge_solvers[id]->pastix_data->pastix_comm, &status);
            MPI_Recv(recv_data, cnt*blocksize, COMM_FLOAT,
                     iter, TAG_VAL2,
                     zmurge_solvers[id]->pastix_data->pastix_comm, &status);
            myvals_ptr = recv_data;

            for (iter2 = 0; iter2 < cnt; iter2++) {
              INTS iterdof;
              pastix_complex64_t*node;
              pastix_int_t index;

              index = sequence->recv_indexes[iter][recv_cnt[iter]];
              recv_cnt[status.MPI_SOURCE]++;

              node  = &(zmurge_solvers[id]->values[index]);
              for (iterdof = 0; iterdof < dof*dof;
                   iterdof++, node++, myvals_ptr++)
                MURGE_ADD(*node, *myvals_ptr,
                          sequence->op_local_entries);
            }
            done[iter] = API_YES;
            done_cnt++;
          }
        }
      }
    }
    for (iter =0; iter < zmurge_solvers[id]->pastix_data->procnbr; iter++) {
      MPI_Status status;
      MPI_Wait(&(send_reqs[iter]), &(status));
      MPI_Wait(&(send_reqs2[iter]), &(status));
    }

    MURGE_FREE(done);
    MURGE_FREE(send_reqs);
    MURGE_FREE(send_reqs2);
    MURGE_FREE(send_array);
    MURGE_FREE(recv_data);
    MURGE_FREE(recv_cnt);
    MURGE_FREE(send_cnt);
  }
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_VALUES_OK);
  MURGE_DUMP_MATRIX;
  return MURGE_SUCCESS;
}

/*
 * Function: ZMURGE_AssemblyDeleteSequence
 *
 * Destroy an assembly sequence
 *
 *   id      - Solver instance identification number.
 *   id_seq  - Sequence ID.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 *                         *id_seq* is not valid.
 */
INTS ZMURGE_AssemblyDeleteSequence(INTS id, INTS id_seq) {
  z_murge_seq_t * sequence;
  z_murge_seq_t * psequence;
  pastix_int_t iter;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  psequence = NULL;
  sequence  = zmurge_solvers[id]->sequences;
  while(sequence != NULL && sequence->ID != id_seq) {
    psequence = sequence;
    sequence  = sequence->next;
  }

  if (sequence == NULL) {
    errorPrint("Sequence %d not found", id_seq);
    return MURGE_ERR_PARAMETER;
  }

  if (psequence != NULL) {
    psequence->next = sequence->next;
  }
  else
    {
      zmurge_solvers[id]->sequences = sequence->next;
    }
  MURGE_FREE(sequence->indexes);
  MURGE_FREE(sequence->recv_nbr);
  if (NULL != sequence->recv_indexes) {
    for (iter = 0; iter < zmurge_solvers[id]->pastix_data->procnbr; iter++)
      MURGE_FREE(sequence->recv_indexes[iter]);
    MURGE_FREE(sequence->recv_indexes);
  }
  MURGE_FREE(sequence);
  return MURGE_SUCCESS;
}
/*
 * Function: ZMURGE_AssemblyBegin
 *
 * Check that preprocessing has been performed, if not performs it.
 *
 * Allocate ijv structure which will be used to store I,J,v[dof*dof].
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   N       - Global number of nodes in the graph.
 *   coefnbr - Number of coeficients the calling MPI node will add
 *             (in term of coefficients, not node).
 *   op      - Operation to perform for coefficient which appear
 *             several time (see <MURGE_ASSEMBLY_OP>).
 *   op2     - Operation to perform when a coefficient is set by
 *             two different processors (see <MURGE_ASSEMBLY_OP>).
 *   mode    - Indicates if user ensure he will respect zmurge_solvers distribution
 *             (see <MURGE_ASSEMBLY_MODE>).
 *   sym     - Indicates if user will give coefficient in a symmetric way
 *             (ie: only triangullar part) or not.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If graph hasn't been built before.
 *   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 *                         *op*, *mode*, *sym*, or *coefnbr* are not valid.
 */
INTS ZMURGE_AssemblyBegin(INTS id, INTS N, INTL coefnbr, INTS op,
                         INTS op2, INTS mode, INTS sym) {
  int dof;

  print_debug(DBG_MURGE, ">> ZMURGE_AssemblyBegin\n");
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];


  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK)) {
    CHECK_PREPROCESSING(id);

    CHECK_L2G(id);
  }
  else {
    if (mode == MURGE_ASSEMBLY_RESPECT) {
      errorPrint("Can't use MURGE_ASSEMBLY_RESPECT if you didn't build the graph before\n");
    }
  }

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif
  if (coefnbr < 0) {
    zmurge_solvers[id]->dynamic = API_YES;
    coefnbr = -coefnbr;
  }
  else {
#ifdef MURGE_INSERT_DIRECTLY
    if (zmurge_solvers[id]->colptr != NULL) {
      zmurge_solvers[id]->dynamic = API_YES;
      coefnbr = 1+coefnbr/1000;
    }
    else
#endif /* MURGE_INSERT_DIRECTLY */
      {
        zmurge_solvers[id]->dynamic = API_NO;
      }
  }
#ifdef MURGE_INSERT_DIRECTLY
  /* allocate values */
  if (zmurge_solvers[id]->colptr != NULL && zmurge_solvers[id]->values == NULL) {
    MURGE_MEMALLOC(zmurge_solvers[id]->values, dof*dof*(zmurge_solvers[id]->colptr[zmurge_solvers[id]->n]-1), pastix_complex64_t);
    memset(zmurge_solvers[id]->values, 0, (dof*dof*(zmurge_solvers[id]->colptr[zmurge_solvers[id]->n]-1))*sizeof(pastix_complex64_t));
  }
#endif /* MURGE_INSERT_DIRECTLY */
  zmurge_solvers[id]->N           = N;
  zmurge_solvers[id]->coefnbr     = coefnbr;
  zmurge_solvers[id]->nodenbr     = coefnbr/(dof*dof);

  vcsc_init(&zmurge_solvers[id]->vcsc, N, coefnbr/dof*dof, dof, id);

  zmurge_solvers[id]->edgenbr  = coefnbr;
  zmurge_solvers[id]->cnt      = 0;
  zmurge_solvers[id]->cnt_zero = 0;
  zmurge_solvers[id]->cnt_node = 0;
  zmurge_solvers[id]->mode     = mode;
  zmurge_solvers[id]->op       = op;
  zmurge_solvers[id]->op2      = op2;
  zmurge_solvers[id]->sym      = sym;


#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif

  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_VALUES_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_FACTO_OK);

  MURGE_STATE_TRUE(zmurge_solvers[id]->state,  MURGE_MATR_BUILD);
  MURGE_MEMORY_USAGE_PRINT("ZMURGE_AssemblyBegin");
  return MURGE_SUCCESS;
}

/*
 * Function:ZMURGE_AssemblySetValue
 *
 * Check that we are in an assembly section.
 *
 * Check that the number of coefficient entered will not
 * overpass the number of coefficient waited.
 *
 * Adds ROW, COL and value into the zmurge_solvers[id]->tmpijv structure.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   ROW     - Global row number of the coefficient.
 *   COL     - Global column number of the coefficient.
 *   value   - value of the coefficient.
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If we are not in an assembly section.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 */
INTS ZMURGE_AssemblySetValue     (INTS id, INTS ROW, INTS COL, COEF value) {
#ifdef MURGE_CHECK
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_MATR_BUILD))) {
    errorPrint("Need to call ZMURGE_AssemblyBegin first");
    return MURGE_ERR_ORDER;
  }
#endif /* MURGE_CHECK */
  return ZMURGE_AssemblySetValue_(id, ROW, COL, value);
}

static inline
INTS ZMURGE_AssemblySetValue_    (INTS id, INTS ROW, INTS COL, COEF value) {
  int dof;


  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
    COL += 1;
    ROW += 1;
  }

  {
    PASTIX_INT node_col_glob;
    node_col_glob = (COL-1 - (COL-1)%dof)/dof;

    PASTIX_INT node_row_glob;
    node_row_glob = (ROW-1 - (ROW-1)%dof)/dof;
    if ( node_row_glob != node_col_glob ) {
      if (zmurge_solvers[id]->dropmask != NULL) {
        if (zmurge_solvers[id]->dropmask[node_col_glob] &&
            zmurge_solvers[id]->dropmask[node_row_glob]) {
      zmurge_solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
    }
  }

      /* If (node_row_glob, node_col_glob) has to be droped */
      if ( ( zmurge_solvers[id]->droprows != NULL &&
             zmurge_solvers[id]->droprows[node_row_glob] ) ||
           ( zmurge_solvers[id]->dropcols != NULL &&
             zmurge_solvers[id]->dropcols[node_col_glob] ) ) {
        /* If (node_col_glob, node_row_glob) has to be droped */
        if ( ( zmurge_solvers[id]->droprows != NULL &&
               zmurge_solvers[id]->droprows[node_col_glob] ) ||
             ( zmurge_solvers[id]->dropcols != NULL &&
               zmurge_solvers[id]->dropcols[node_row_glob] ) ) {
      zmurge_solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
        } else {
          /* we don't want to break symmetry */
          value = 0;
        }
    }
  }

    if (zmurge_solvers[id]->mode == MURGE_ASSEMBLY_DROPNONLOCAL &&
        NULL != zmurge_solvers[id]->g2l && zmurge_solvers[id]->g2l[node_col_glob] < 1) {
      /* drop non local entries */
      zmurge_solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
    }
  }


  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_ONLY_PROD) && value == 0.0) {
      zmurge_solvers[id]->cnt_zero++;
      return MURGE_SUCCESS;
    }

#ifdef MURGE_CHECK
  if ((COL < 1) || (ROW < 1) ||
      (ROW > (((zmurge_solvers[id]->N)*dof))) ||
      (COL > (((zmurge_solvers[id]->N)*dof)))) {
    errorPrint("COL (%ld) or ROW (%ld) is out of range [1-%ld]",
               (long)COL, (long)ROW, (long)(zmurge_solvers[id]->N*dof));
    return MURGE_ERR_PARAMETER;
  }
#endif /* MURGE_CHECK */

#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
  if (COL == MURGE_TRACE_COL && MURGE_TRACE_ROW == ROW)
    fprintf(stdout, "Setting A(%d,%d) <- %g\n",
            MURGE_TRACE_ROW, MURGE_TRACE_COL, value);
#endif

#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&(zmurge_solvers[id]->mutex_tmpmatrix));
#endif

  if ((zmurge_solvers[id]->dynamic == API_NO) &&
      (zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero +
       (zmurge_solvers[id]->cnt_node )*dof*dof + 1 > zmurge_solvers[id]->coefnbr)) {
    errorPrint("Too many coef added in matrix building session (%ld > %ld)",
               (long)(zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero +
                      (zmurge_solvers[id]->cnt_node )*dof*dof + 1),
               (long)zmurge_solvers[id]->coefnbr);
#ifdef MURGE_THREADSAFE
    pthread_mutex_unlock(&(zmurge_solvers[id]->mutex_tmpmatrix));
#endif
    return MURGE_ERR_ORDER;
  }

  {
    pastix_int_t node_col_glob;
    pastix_int_t node_col_loc;


    node_col_glob = (COL-1 - (COL-1)%dof)/dof;

    if (zmurge_solvers[id]->g2l != NULL) {
      node_col_loc = zmurge_solvers[id]->g2l[node_col_glob];
    }
    else {
      /* If we didn't entered the graph we cannot decide if a column is local or not */
      node_col_loc = -1;
    }
#ifdef MURGE_INSERT_DIRECTLY
    if ( zmurge_solvers[id]->colptr != NULL && node_col_loc > 0 ) {
      pastix_int_t node_row_glob;
      pastix_int_t node_idx;
      int coef_idx;


      node_row_glob = (ROW-1 - (ROW-1)%dof)/dof;

      node_col_loc--;
      /* The column is local we add it into the local CSC */
      for (node_idx = zmurge_solvers[id]->colptr[node_col_loc]-1;
           node_idx < zmurge_solvers[id]->colptr[node_col_loc+1]-1;
           node_idx++)
        if (zmurge_solvers[id]->rows[node_idx]-1 == node_row_glob) break;

      if (node_idx == zmurge_solvers[id]->colptr[node_col_loc+1]-1) {
        if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_ONLY_PROD)) {
          /* we will add it later */
          vcsc_add(zmurge_solvers[id]->vcsc, COL, ROW, value, zmurge_solvers[id]->op, id);
          zmurge_solvers[id]->cnt++;
        }
        else
          {
            errorPrint("ROW (%ld) not found in COL (%d)",
                       (long)ROW, (long)COL);
#  ifdef MURGE_THREADSAFE
            pthread_mutex_unlock(&(zmurge_solvers[id]->mutex_tmpmatrix));
#  endif
            return MURGE_ERR_PARAMETER;
          }
      }
      else
        {
          /* we found the correct node */
          coef_idx = node_idx*dof*dof + ((COL-1)%dof)*dof + (ROW-1)%dof;
          zmurge_solvers[id]->values[coef_idx] =
            func(zmurge_solvers[id]->values[coef_idx], value);
        }
    }
    else
#endif /* MURGE_INSERT_DIRECTLY */
      {
        /* The column has to be sent to the correct CPU */
        if (zmurge_solvers[id]->g2l != NULL) {
          if (node_col_loc <= 0) {
            if (zmurge_solvers[id]->mode == MURGE_ASSEMBLY_RESPECT) {
              errorPrint("Column %d is not local", COL);
              return MURGE_ERR_PARAMETER;
            }
          }
        }

        vcsc_add(zmurge_solvers[id]->vcsc, COL, ROW, value, zmurge_solvers[id]->op, id);

        zmurge_solvers[id]->cnt++;
      }
  }
#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&(zmurge_solvers[id]->mutex_tmpmatrix));
#endif
  return MURGE_SUCCESS;
}


/*
 * Function:ZMURGE_AssemblySetNodeValues
 *
 * Check that we are in an assembly section.
 *
 * Check that the number of coefficient entered will not
 * overpass the number of coefficient waited.
 *
 * Adds ROW, COL and value into the zmurge_solvers[id]->tmpijv structure.
 *
 * Parameters:
 *   id      - Solver instance identification number.
 *   ROW     - Global row number of the coefficient.
 *   COL     - Global column number of the coefficient.
 *   values  - value of the coefficient.(dof^2 matrix)
 *
 * Returns:
 *   MURGE_SUCCESS       - If function runned successfully.
 *   MURGE_ERR_ORDER     - If we are not in an assembly section.
 *   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or
 */
INTS ZMURGE_AssemblySetNodeValues (INTS id, INTS ROW, INTS COL, COEF *values) {
  return ZMURGE_AssemblySetNodeValues_(id, ROW, COL, values);
}

static inline
INTS ZMURGE_AssemblySetNodeValues_ (INTS id, INTS ROW, INTS COL, COEF *values) {
  int dof;

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (dof == 1)
    return ZMURGE_AssemblySetValue_(id, ROW, COL, *values);


  if (zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
    ROW += 1;
    COL += 1;
  }

#ifdef MURGE_CHECK
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_MATR_BUILD))) {
    errorPrint("Need to call ZMURGE_AssemblyBegin first");
    return MURGE_ERR_ORDER;
  }
#ifdef MURGE_THREADSAFE
  pthread_mutex_lock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif
  if (zmurge_solvers[id]->dynamic == API_NO &&
      ( (zmurge_solvers[id]->cnt_node+1)*dof*dof +
        zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero > zmurge_solvers[id]->coefnbr)) {
    errorPrint("Too many coef added in matrix building session (%ld > %ld)",
               (long)((zmurge_solvers[id]->cnt_node+1)*dof*dof + zmurge_solvers[id]->cnt),
               (long)(zmurge_solvers[id]->coefnbr));
#ifdef MURGE_THREADSAFE
    pthread_mutex_unlock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif
    return MURGE_ERR_ORDER;
  }

  if ((COL < 1) || (ROW < 1) ||
      (ROW > ((zmurge_solvers[id]->N))) ||
      (COL > ((zmurge_solvers[id]->N)))) {
    errorPrint("COL (%ld) or ROW (%ld) is out of range [1-%ld]",
               (long)COL, (long)ROW, (long)(zmurge_solvers[id]->N));
    return MURGE_ERR_PARAMETER;
  }
#endif /* MURGE_CHECK */

#if (defined MURGE_TRACE_ROW && defined MURGE_TRACE_COL)
  if ((COL-1)*dof <= MURGE_TRACE_COL-1 &&
      (COL)*dof > MURGE_TRACE_COL-1 &&
      (ROW-1)*dof <= MURGE_TRACE_ROW-1 &&
      ROW*dof > MURGE_TRACE_ROW-1) {

    fprintf(stdout, "Setting A(%d*%d+%d,%d*%d+%d) <- %g\n",
            ROW-1, dof, (MURGE_TRACE_ROW-1-(ROW-1)*dof)+1,
            COL-1, dof, (MURGE_TRACE_COL-1-(COL-1)*dof)+1,
            values[(MURGE_TRACE_COL-1-(COL-1)*dof)*dof +
                   (MURGE_TRACE_ROW-1 - (ROW-1)*dof)]);
  }
#endif

  if (ROW != COL) {
    if (zmurge_solvers[id]->dropmask != NULL)
      if ( zmurge_solvers[id]->dropmask[(COL-1)] &&
           zmurge_solvers[id]->dropmask[(ROW-1)]) {
      zmurge_solvers[id]->cnt_zero+=dof*dof;
      return MURGE_SUCCESS;
    }

    /* If (ROW, COL) has to be droped */
    if ( ( zmurge_solvers[id]->droprows != NULL &&
           zmurge_solvers[id]->droprows[(ROW-1)] ) ||
         ( zmurge_solvers[id]->dropcols != NULL &&
           zmurge_solvers[id]->dropcols[(COL-1)] ) ) {
      /* If (COL, ROW) has to be droped */
      if ( ( zmurge_solvers[id]->droprows != NULL &&
             zmurge_solvers[id]->droprows[(COL-1)] ) ||
           ( zmurge_solvers[id]->dropcols != NULL &&
             zmurge_solvers[id]->dropcols[(ROW-1)] ) ) {
      zmurge_solvers[id]->cnt_zero+=dof*dof;
      return MURGE_SUCCESS;
      } else {
        /* we don't want to break symmetry */
        memset(values, 0, sizeof(COEF)*dof*dof);
    }
    }
  }

  if (zmurge_solvers[id]->mode == MURGE_ASSEMBLY_DROPNONLOCAL &&
      NULL != zmurge_solvers[id]->g2l && zmurge_solvers[id]->g2l[COL-1] < 1) {
    /* drop non local entries */
      zmurge_solvers[id]->cnt_zero+=dof*dof;
      return MURGE_SUCCESS;
    }


  {
    pastix_int_t node_col_loc;
    if (zmurge_solvers[id]->g2l != NULL) {
      node_col_loc = zmurge_solvers[id]->g2l[COL-1];
    }
    else {
      /* If we didn't entered the graph we cannot decide if a column is local or not */
      node_col_loc = -1;
    }

#ifdef MURGE_INSERT_DIRECTLY
    if ( zmurge_solvers[id]->colptr != NULL && node_col_loc > 0 ) {
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
      if ((COL-1)*dof <= MURGE_TRACE_COL-1 &&
          (COL)*dof > MURGE_TRACE_COL-1 &&
          (ROW-1)*dof <= MURGE_TRACE_ROW-1 &&
          ROW*dof > MURGE_TRACE_ROW-1) {
        fprintf(stdout, "The TRACED column is local\n");
      }
#endif
      pastix_complex64_t (*func)(pastix_complex64_t , pastix_complex64_t);
      pastix_int_t node_idx;
      int coef_idx;

      CHOOSE_FUNC(func, zmurge_solvers[id]->op);

      node_col_loc--;
      /* The column is local we add it into the local CSC */
      for (node_idx = zmurge_solvers[id]->colptr[node_col_loc]-1;
           node_idx < zmurge_solvers[id]->colptr[node_col_loc+1]-1;
           node_idx++)
        if (zmurge_solvers[id]->rows[node_idx] == ROW) break;

      if (node_idx == zmurge_solvers[id]->colptr[node_col_loc+1]-1) {
        if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_ONLY_PROD)) {
          vcsc_add_node(zmurge_solvers[id]->vcsc, COL, ROW, values, zmurge_solvers[id]->op, id);
          zmurge_solvers[id]->cnt_node++;
        }
        else
          {
            errorPrint("ROW (%ld) not found in COL (%d)",
                       (long)ROW, (long)COL);
#  ifdef MURGE_THREADSAFE
            pthread_mutex_unlock(&zmurge_solvers[id]->mutex_tmpmatrix);
#  endif
            return MURGE_ERR_PARAMETER;
          }
      }
      /* we found the correct node */
      for ( coef_idx = 0;
            coef_idx < dof*dof;
            coef_idx++)
        zmurge_solvers[id]->values[node_idx*dof*dof + coef_idx] =
          func(zmurge_solvers[id]->values[node_idx*dof*dof+coef_idx],
               values[coef_idx]);
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
      if ((COL-1)*dof <= MURGE_TRACE_COL-1 &&
          (COL)*dof > MURGE_TRACE_COL-1 &&
          (ROW-1)*dof <= MURGE_TRACE_ROW-1 &&
          ROW*dof > MURGE_TRACE_ROW-1) {
        fprintf(stdout, "Setting A(%d*%d+%d,%d*%d+%d) = %g\n",
                ROW-1, dof, (MURGE_TRACE_ROW-1-(ROW-1)*dof)+1,
                COL-1, dof, (MURGE_TRACE_COL-1-(COL-1)*dof)+1,
                zmurge_solvers[id]->values[node_idx*dof*dof +
                                    (MURGE_TRACE_COL-1-(COL-1)*dof)*dof +
                                    (MURGE_TRACE_ROW-1 - (ROW-1)*dof)]);
      }
#endif
    } else
#endif /* MURGE_INSERT_DIRECTLY */
      {
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_COL)
        if ((COL-1)*dof <= MURGE_TRACE_COL-1 &&
            (COL)*dof > MURGE_TRACE_COL-1 &&
            (ROW-1)*dof <= MURGE_TRACE_ROW-1 &&
            ROW*dof > MURGE_TRACE_ROW-1) {
          fprintf(stdout, "The TRACED column is not local (node_col_loc %d)\n", node_col_loc);
        }
#endif

        if (zmurge_solvers[id]->g2l != NULL) {
          if (node_col_loc <= 0) {
            if (zmurge_solvers[id]->mode == MURGE_ASSEMBLY_RESPECT) {
              errorPrint("Column %d is not local", COL);
              return MURGE_ERR_PARAMETER;
            }
          }
        }
        vcsc_add_node(zmurge_solvers[id]->vcsc, COL, ROW, values, zmurge_solvers[id]->op, id);
        zmurge_solvers[id]->cnt_node++;
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
        if ((COL-1)*dof <= MURGE_TRACE_COL-1 &&
            (COL)*dof > MURGE_TRACE_COL-1 &&
            (ROW-1)*dof <= MURGE_TRACE_ROW-1 &&
            ROW*dof > MURGE_TRACE_ROW-1) {
          int index = 0;
          while (index < zmurge_solvers[id]->vcsc.colsizes[COL-1] &&
                 zmurge_solvers[id]->vcsc.rows[COL-1][index] !=  ROW)
            index++;

          fprintf(stdout, "Setting vcsc[%d-1][%d*%d*%d + %d *%d + %d] = %g (vcsc.rows[%d] %d)\n",
                  COL, index, dof, dof, (MURGE_TRACE_COL-1-(COL-1)*dof), dof, 
                  (MURGE_TRACE_ROW-1 - (ROW-1)*dof),
                  zmurge_solvers[id]->vcsc.values[COL-1][index*dof*dof + 
                                                  (MURGE_TRACE_COL-1-(COL-1)*dof)*dof +
                                                  (MURGE_TRACE_ROW-1 - (ROW-1)*dof)],
                  index,
                  zmurge_solvers[id]->vcsc.rows[COL-1][index]);
        }
#endif
      }
  }
#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif
  return MURGE_SUCCESS;
}

INTS ZMURGE_AssemblySetBlockValues(INTS id, INTS nROW, INTS *ROWlist,
                                  INTS nCOL, INTS *COLlist, COEF *values) {
  return ZMURGE_AssemblySetBlockValues_(id, nROW, ROWlist, nCOL, COLlist, values);
}
static inline
INTS ZMURGE_AssemblySetBlockValues_(INTS id, INTS nROW, INTS *ROWlist,
                                   INTS nCOL, INTS *COLlist, COEF *values) {
  pastix_int_t  iter;
  pastix_int_t  iter2;
  INTS ret;
  pastix_int_t  iterv;

#ifdef MURGE_CHECK
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_MATR_BUILD))) {
    errorPrint("Need to call ZMURGE_GraphBegin first");
    return MURGE_ERR_ORDER;
  }
#endif /* MURGE_CHECK */
  iterv = 0;
  for (iter = 0; iter < nCOL; iter ++) {
    for (iter2 = 0; iter2 < nROW; iter2++) {
      if (zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR] == 1) {
        ret = ZMURGE_AssemblySetValue_(id,
                                      ROWlist[iter2],
                                      COLlist[iter],
                                      values[iterv]);
        if (MURGE_SUCCESS != ret)
          return ret;
        iterv++;
      } else {
          if (zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR] < 1)
            return MURGE_ERR_PARAMETER;
        ret = ZMURGE_AssemblySetNodeValues_(id,
                                           ROWlist[iter2],
                                           COLlist[iter] ,
                                           &(values[iterv]));
        if (MURGE_SUCCESS != ret)
          return ret;

        iterv+=zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR]*
          zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
      }
    }
  }
  return MURGE_SUCCESS;
}

INTS ZMURGE_AssemblySetListOfBlockValues(INTS id, INTS nBlocks,
                                        INTS nROW, INTS *ROWlist,
                                        INTS nCOL, INTS *COLlist,
                                        COEF *values) {
  INTS i, ierr;
  int dof;

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  for (i = 0; i < nBlocks; i++) {
    ierr = ZMURGE_AssemblySetBlockValues_(id,
                                         nROW, ROWlist,
                                         nCOL, COLlist,
                                         values);
    if (ierr != MURGE_SUCCESS)
      return ierr;
    ROWlist+=nROW;
    COLlist+=nCOL;
    values+=nROW*nCOL*dof*dof;
  }
  return MURGE_SUCCESS;
}
/*
 Function: ZMURGE_AssemblyEnd

 We have on each proc a part of the matrix in
 two structure, one containing nodes to add
 to the CSCd the other containing simple values.

 We send all data to his owner:
 - We sort our data structures (IJV structures)
 using the "owner" attribute.
 - We send non local data to other processors.

 We merge all data in the node structure.
 - We receive Data and merge node structure with simple
 values one.
 - We look for each coef in node structure, if present we modify the node, if
 not, we search in the CSCd and directly modify it. Else we construct
 a new node and add it.

 We Add this structure to the local CSCd.

 */

INTS ZMURGE_AssemblyEnd(INTS id) {
  int          dof;
#ifdef CENTRALISED
  pastix_int_t         *total_nodelist;
#endif
#ifdef MURGE_TIME
  Clock       clock;
#endif
  CLOCK_INIT;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  print_debug(DBG_MURGE, ">> ZMURGE_AssemblyEnd\n");

  /* Check that we are in a matrix building session */
  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_MATR_BUILD))) {
    errorPrint("Need to call ZMURGE_AssemblyBegin first");
    return MURGE_ERR_ORDER;
  }

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  /* check that we entered the correct number of values */
  if (zmurge_solvers[id]->dynamic == API_NO &&
      zmurge_solvers[id]->cnt_node*dof*dof + zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero != zmurge_solvers[id]->edgenbr) {
    errorPrint("Wrong number of entries  (%ld != %ld) ",
               (long)zmurge_solvers[id]->cnt_node*dof*dof + zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero ,
               (long)zmurge_solvers[id]->edgenbr);
    return MURGE_ERR_ORDER;
  }

  /* some information about skipped zeros entries */
  if (zmurge_solvers[id]->cnt_zero != 0) {
    if (zmurge_solvers[id]->pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NOT) {
      fprintf(stdout,
              "%ld (%.2g %%) zero entries were skipped on proc %ld\n",
              (long)zmurge_solvers[id]->cnt_zero,
              (double)((double)100.0*((double)zmurge_solvers[id]->cnt_zero)/
                       ((double)(zmurge_solvers[id]->cnt_node*dof*dof +
                                 zmurge_solvers[id]->cnt + zmurge_solvers[id]->cnt_zero))),
              (long)zmurge_solvers[id]->pastix_data->procnum);
    }
    else
      {
        if (zmurge_solvers[id]->pastix_data->iparm[IPARM_VERBOSE] > API_VERBOSE_NO) {
          pastix_int_t nz_glob;
          pastix_int_t zeros_glob;
          pastix_int_t nz = (zmurge_solvers[id]->colptr[zmurge_solvers[id]->n]-1)*dof*dof;
          MPI_Reduce( &(zmurge_solvers[id]->cnt_zero), &zeros_glob,
                      1, PASTIX_MPI_INT,
                      MPI_SUM, 0, zmurge_solvers[id]->pastix_data->pastix_comm);
          MPI_Reduce( &nz, &nz_glob,
                      1, PASTIX_MPI_INT,
                      MPI_SUM, 0, zmurge_solvers[id]->pastix_data->pastix_comm);
          if (zmurge_solvers[id]->pastix_data->procnum == 0) {
            fprintf(stdout,
                    "%ld zero entries were skipped"
                    " (from %ld (%.3g%%))\n",
                    (long int)zeros_glob,
                    (long int)nz_glob,
                    100.0*((double)(zeros_glob)/
                           ((double)(nz_glob))));
          }

        }
      }
  }

#ifdef CENTRALISED
  MURGE_MEMALLOC(total_nodelist, zmurge_solvers[id]->N, pastix_int_t);
  for (iter = 0; iter < zmurge_solvers[id]->N; iter++)
    total_nodelist[iter] = iter+1;
#endif

  dof     = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];


  /* If the CSCd is existing it will fill it with new values, else it will create it.
   *
   * NB : VCSC has to fit into CSCd
   */
  vcsc_to_cscd(zmurge_solvers[id]->vcsc, zmurge_solvers[id]->pastix_data->pastix_comm,
               &(zmurge_solvers[id]->n), &(zmurge_solvers[id]->colptr),
               &(zmurge_solvers[id]->rows), &(zmurge_solvers[id]->values),
               &(zmurge_solvers[id]->l2g), &(zmurge_solvers[id]->g2l), zmurge_solvers[id]->op,
               MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_ONLY_PROD), id);
  vcsc_destroy(zmurge_solvers[id]->vcsc, id);

#ifdef MURGE_THREADSAFE
  pthread_mutex_unlock(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif

  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_VALUES_OK);

  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_MATR_BUILD);
  MURGE_DUMP_MATRIX;
  CLOCK_PRINT("ZMURGE_AssemblyEnd");
  return MURGE_SUCCESS;
}

INTS ZMURGE_MatrixReset(INTS id){
  pastix_int_t nbcoef;
  pastix_int_t dof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_ONLY_PROD) &&
      zmurge_solvers[id]->colptr) {
    MURGE_FREE(zmurge_solvers[id]->colptr);
    MURGE_FREE(zmurge_solvers[id]->rows);
    MURGE_FREE(zmurge_solvers[id]->values);
  }
  else {
    if (zmurge_solvers[id]->values != NULL) {
      nbcoef = zmurge_solvers[id]->colptr[zmurge_solvers[id]->n]-1;
      dof    = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

      memset(zmurge_solvers[id]->values, 0, nbcoef*dof*dof*sizeof(pastix_complex64_t));
    }
    MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_VALUES_OK);
  }
  return MURGE_SUCCESS;
}

INTS ZMURGE_MatrixGlobalCSR(INTS id, INTS N, INTL *rowptr, INTS *COLS,
                           COEF *values, INTS root, INTS op, INTS sym) {
  pastix_int_t  dof;
  pastix_int_t  coefnbr;
  pastix_int_t  iter;
  pastix_int_t  iter2;
  INTS ret;
  int baseval;

  CHECK_SOLVER_ID(id);
  if (zmurge_solvers[id]->pastix_data->procnum == 0) {
    errorPrintW("ZMURGE_MatrixGlobalCSR is not optimal with PaStiX, try using ZMURGE_MatrixGlobalCSC instead");
  }
  CHECK_SOLVER_PARAM(id);

  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("Graph has to be built before");
    return MURGE_ERR_ORDER;
  }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);
  baseval = zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  dof     = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  /* Si on a une matrice symetrique autant faire un ZMURGE_MatrixGlobalCSC */
  if (zmurge_solvers[id]->pastix_data->iparm[IPARM_SYM] == API_SYM_YES)
    return  ZMURGE_MatrixGlobalCSC(id, N, rowptr, COLS, values, root, op, sym);

  coefnbr = rowptr[N] - baseval;

  if (zmurge_solvers[id]->pastix_data->procnum == root ||
      (root == -1 && zmurge_solvers[id]->pastix_data->procnum == 0)) {
    ret = ZMURGE_AssemblyBegin(id, N, coefnbr,  op, op,MURGE_ASSEMBLY_FOOL,
                              zmurge_solvers[id]->pastix_data->iparm[IPARM_SYM]);
    if (MURGE_SUCCESS != ret)
      return ret;
    for (iter = 0; iter < N; iter++) {
      for (iter2 = rowptr[iter]; iter2 < rowptr[iter+1]; iter2++) {
        if (dof == 1) {

          ret = ZMURGE_AssemblySetValue_(id,
                                        iter+baseval,
                                        COLS[iter2-baseval],
                                        values[iter2-baseval]);
          if (MURGE_SUCCESS != ret)
            return ret;
        }
        else
          {
            ret = ZMURGE_AssemblySetNodeValues_(id,
                                               iter+baseval,
                                               COLS[iter2-baseval],
                                               &(values[(iter2-baseval)*
                                                        dof*dof]));
            if (MURGE_SUCCESS != ret)
              return ret;
          }
      }
    }
  }
  else
    {
      ret = ZMURGE_AssemblyBegin(id, N, 0, op, op,
                                MURGE_ASSEMBLY_FOOL,
                                zmurge_solvers[id]->pastix_data->iparm[IPARM_SYM]);
      if (MURGE_SUCCESS != ret)
        return ret;
    }

  if (MURGE_SUCCESS != (ret = ZMURGE_AssemblyEnd(id)))
    return ret;

  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_VALUES_OK);
  return MURGE_SUCCESS;
}
/*
 Function: ZMURGE_MatrixGlobalCSC

 Give a CSC on one processor to PaStiX.
 */
INTS ZMURGE_MatrixGlobalCSC(INTS id, INTS N, INTL *COLPTR, INTS *ROWS,
                           COEF *values, INTS root, INTS op, INTS sym) {
  pastix_int_t          *l2g = NULL;
  pastix_int_t           procnum;
  pastix_int_t           localn;
  pastix_int_t          *tmpcolptr;
  pastix_int_t          *tmprows;
  pastix_complex64_t        *tmpvalues;
  pastix_int_t          *tmpcolptr2;
  pastix_int_t          *tmprows2;
  pastix_int_t           tmpn;
  MPI_Status    status;
  pastix_int_t           iter;
  int           dof;
  CSCD_OPERATIONS_t func;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);


  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("Graph has to be built before");
    return MURGE_ERR_ORDER;
  }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);

  CHOOSE_FUNC(func, op);

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  /* Si tout le monde est racine */
  if (root == -1 || root == zmurge_solvers[id]->pastix_data->procnum) {

    if (sizeof(INTS) != sizeof(pastix_int_t)) {
      MURGE_MEMALLOC(tmprows2, (COLPTR[N]-1), pastix_int_t);
      for (iter = 0; iter <  COLPTR[N]-1; iter++) {
        tmprows2[iter] = (pastix_int_t)ROWS[iter];
      }
    }
    else
      {
        tmprows2 = (pastix_int_t*)ROWS;
      }

    if (sizeof(INTL) != sizeof(pastix_int_t)) {
      MURGE_MEMALLOC(tmpcolptr2, N+1, pastix_int_t);
      for (iter = 0; iter <  N+1; iter++) {
        tmpcolptr2[iter] = (pastix_int_t)COLPTR[iter];
      }
    }
    else
      {
        tmpcolptr2 = (pastix_int_t*)COLPTR;
      }


    if (zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
      tmprows2 --;
      values   --;
    }

    MURGE_MEMALLOC(l2g, N, pastix_int_t);

    for (iter = 0; iter < N; iter++)
      l2g[iter] = iter+1;

    /* colptr must be allocated for z_cscd_addlocal_int */
    if (NULL == zmurge_solvers[id]->colptr) {
      MURGE_MEMALLOC(zmurge_solvers[id]->colptr, zmurge_solvers[id]->n+1, pastix_int_t);
      for (iter = 0; iter < zmurge_solvers[id]->n+1; iter++) {
        zmurge_solvers[id]->colptr[iter] = zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
      }
    }

    z_cscd_addlocal_int(zmurge_solvers[id]->n,
                      zmurge_solvers[id]->colptr,
                      zmurge_solvers[id]->rows,
                      zmurge_solvers[id]->values,
                      zmurge_solvers[id]->l2g,
                      N, tmpcolptr2, tmprows2, (pastix_complex64_t*)values, l2g,
                      &tmpn, &tmpcolptr, &tmprows, &tmpvalues, func,
                      dof,  API_YES);
    MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpcolptr), char);
    MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmprows), char);
    MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpvalues), char);
    if (zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
      tmprows2 ++;
      values   ++;
    }

    if (root == -1 && sizeof(INTS) != sizeof(pastix_int_t)) {
      MURGE_FREE(tmprows2);
    }
    if (root == -1 && sizeof(INTL) != sizeof(pastix_int_t)) {
      MURGE_FREE(tmpcolptr2);
    }
    MURGE_FREE(zmurge_solvers[id]->colptr);
    MURGE_FREE(zmurge_solvers[id]->rows);
    MURGE_FREE(zmurge_solvers[id]->values);

    zmurge_solvers[id]->colptr = tmpcolptr;
    zmurge_solvers[id]->rows   = tmprows;
    zmurge_solvers[id]->values = tmpvalues;

  }
  if (root != -1) {
    /* si on est le processeur racine */
    if (root == zmurge_solvers[id]->pastix_data->procnum) {

      /* Pour chaque processeur, on calcule
       la CSCD a ajouter puis on l'envoi.
       */
      for (procnum = 0; procnum < zmurge_solvers[id]->pastix_data->procnbr; procnum++) {
        if (procnum != zmurge_solvers[id]->pastix_data->procnum) {
          MPI_Recv(&localn, 1,  PASTIX_MPI_INT, procnum, TAG_SIZE,
                   zmurge_solvers[id]->pastix_data->pastix_comm, &status);
          MPI_Recv(l2g, localn, PASTIX_MPI_INT, procnum, TAG_L2G,
                   zmurge_solvers[id]->pastix_data->pastix_comm, &status);

          MURGE_MEMALLOC(tmpcolptr, localn + 1, pastix_int_t);

          for (iter = 0; iter < localn+1; iter++) {
            tmpcolptr[iter] = 1;
          }
          if (sizeof(INTS) != sizeof(pastix_int_t)) {
            MURGE_MEMALLOC(tmprows2, (COLPTR[N]-1), pastix_int_t);
            for (iter = 0; iter <  COLPTR[N]-1; iter++) {
              tmprows2[iter] = (pastix_int_t)ROWS[iter];
            }
          }
          else
            {
              tmprows2 = (pastix_int_t*)ROWS;
            }

          if (sizeof(INTL) != sizeof(pastix_int_t)) {
            MURGE_MEMALLOC(tmpcolptr2, N+1, pastix_int_t);
            for (iter = 0; iter <  N+1; iter++) {
              tmpcolptr2[iter] = (pastix_int_t)COLPTR[iter];
            }
          }
          else
            {
              tmpcolptr2 = (pastix_int_t*)ROWS;
            }

          if (zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
            tmprows2 --;
            values   --;
          }
          z_cscd_addlocal_int(localn,
                            tmpcolptr,
                            NULL,
                            NULL,
                            l2g,
                            N, tmpcolptr2, tmprows2, (pastix_complex64_t*)values, l2g,
                            &localn,
                            &tmpcolptr,
                            &tmprows,
                            &tmpvalues, func, dof, API_YES);
	  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpcolptr), char);
	  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmprows), char);
	  MURGE_TRACE_MALLOC(PTR_MEMSIZE(tmpvalues), char);

          if (zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL] == 0) {
            tmprows2 ++;
            values   ++;
          }

          /* On envoi tmpcolptr, tmprows, tmpvalues */
          MPI_Send(tmpcolptr,
                   localn+1, PASTIX_MPI_INT, procnum,
                   TAG_COL, zmurge_solvers[id]->pastix_data->pastix_comm);
          MURGE_FREE(tmpcolptr);

          MPI_Send(tmprows,
                   tmpcolptr[localn - 1], PASTIX_MPI_INT, procnum,
                   TAG_ROW, zmurge_solvers[id]->pastix_data->pastix_comm);
          MURGE_FREE(tmprows);

          MPI_Send(tmpvalues,
                   tmpcolptr[localn - 1], COMM_FLOAT, procnum,
                   TAG_VAL, zmurge_solvers[id]->pastix_data->pastix_comm);
          MURGE_FREE(tmpvalues);

          if (sizeof(INTS) != sizeof(pastix_int_t)) {
            MURGE_FREE(tmprows2);
          }
          if (sizeof(INTL) != sizeof(pastix_int_t)) {
            MURGE_FREE(tmpcolptr2);
          }
        }
        else
          {
            /* La CSCd local a déjà été traitée */
          }

      }
      if (sizeof(INTS) != sizeof(pastix_int_t)) {
        MURGE_FREE(tmprows2);
      }
      if (sizeof(INTL) != sizeof(pastix_int_t)) {
        MURGE_FREE(tmpcolptr2);
      }
    }
    else
      {
        /* Si on est pas la racine, on recoit de la racine la CSCd a ajouter
         et on l'ajoute
         */

        MPI_Send(&zmurge_solvers[id]->n,
                 1, PASTIX_MPI_INT, root,
                 TAG_SIZE, zmurge_solvers[id]->pastix_data->pastix_comm);
        localn = zmurge_solvers[id]->n;
        MPI_Send(zmurge_solvers[id]->l2g,
                 zmurge_solvers[id]->n,
                 PASTIX_MPI_INT, root,
                 TAG_L2G, zmurge_solvers[id]->pastix_data->pastix_comm);

        MPI_Recv(zmurge_solvers[id]->colptr,
                 localn+1, PASTIX_MPI_INT, root,
                 TAG_COL, zmurge_solvers[id]->pastix_data->pastix_comm, &status);

        MPI_Recv(zmurge_solvers[id]->rows,
                 tmpcolptr[localn - 1], PASTIX_MPI_INT, root,
                 TAG_ROW, zmurge_solvers[id]->pastix_data->pastix_comm, &status);

        MPI_Recv(zmurge_solvers[id]->values,
                 tmpcolptr[localn - 1], COMM_FLOAT, root,
                 TAG_VAL, zmurge_solvers[id]->pastix_data->pastix_comm, &status);


      }
  }

  MURGE_DUMP_MATRIX;


  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_VALUES_OK);

  return MURGE_SUCCESS;
}

/*
 Function: ZMURGE_MatrixGlobalIJV

 Add the given global Compress Sparse Column matrix to the matrix.

 Parameters:
 id      - Solver instance identification number.
 N       - Number of edges.
 NNZ     - Number of non zeros.
 ROWS    - Global row number array.
 COLS    - Global column number array.
 values  - values array.
 root    - Root processor for MPI communications.
 op      - Operation to perform if a coefficient appear twice
 (see <MURGE_ASSEMBLY_OP>).
 sym     - Indicates if user will give coefficient in a symmetric way
 (ie: only triangullar part) or not.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range,
 if *root*, *op*, *ROWS* or *COLS* are not valid.

 Fortran interface:
 >
 > SUBROUTINE MURGE_MATRIXGLOBALIJV(ID, N, NNZ, ROWS, COLS, VALUES, &
 >                                & ROOT, OP, SYM, IERROR)
 >   INTS,               INTENT(IN)  :: ID, ROOT, OP, SYM, N
 >   INTL,               INTENT(IN)  :: NNZ
 >   INTS, DIMENSION(0), INTENT(IN)  :: ROWS, COLS
 >   COEF, DIMENSION(0), INTENT(IN)  :: VALUES
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_MATRIXGLOBALIJV
 */
INTS ZMURGE_MatrixGlobalIJV(INTS id, INTS N, INTL NNZ, INTS *ROWS, INTS *COLS,
                           COEF *values, INTS root, INTS op, INTS sym) {
  pastix_int_t        iter;                /* Iterators                            */
  INTS       ret;                 /* Return value                         */
  int        dof;                 /* Number of degree of freedom          */

  CHECK_SOLVER_ID(id);
  if (zmurge_solvers[id]->pastix_data->procnum == 0)
    errorPrintW("ZMURGE_MatrixGlobalIJV is not optimal with PaStiX, try using ZMURGE_MatrixGlobalCSC instead");
  CHECK_SOLVER_PARAM(id);


  if ( MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_BUILD) ) {
    errorPrint("Do not call ZMURGE_GraphBegin before");
    return MURGE_ERR_ORDER;
  }

  if ( MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_MATR_BUILD) ) {
    errorPrint("Do not call ZMURGE_AssemblyBegin before");
    return MURGE_ERR_ORDER;
  }
  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK))) {
    errorPrint("Graph has to be built before");
    return MURGE_ERR_ORDER;
  }

  CHECK_PREPROCESSING(id);

  CHECK_L2G(id);
  dof     = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (zmurge_solvers[id]->pastix_data->procnum == root ||
      (root == -1 && zmurge_solvers[id]->pastix_data->procnum == 0)) {
    ret = ZMURGE_AssemblyBegin(id, N, NNZ, op, op, MURGE_ASSEMBLY_FOOL,
                              zmurge_solvers[id]->pastix_data->iparm[IPARM_SYM]);
    if (MURGE_SUCCESS != ret)
      return ret;

    for (iter = 0; iter < NNZ; iter++) {
      if (dof == 1) {

        ret = ZMURGE_AssemblySetValue_(id,
                                      ROWS[iter],
                                      COLS[iter],
                                      values[iter]);
        if (MURGE_SUCCESS != ret)
          return ret;
      }
      else
        {
          ret = ZMURGE_AssemblySetNodeValues_(id,
                                             ROWS[iter],
                                             COLS[iter],
                                             &(values[iter*dof*dof]));
          if (MURGE_SUCCESS != ret)
            return ret;
        }
    }
  }
  else
    {
      ret = ZMURGE_AssemblyBegin(id, N, 0, op, op,
                                MURGE_ASSEMBLY_FOOL,
                                zmurge_solvers[id]->pastix_data->iparm[IPARM_SYM]);
      if (MURGE_SUCCESS != ret)
        return ret;
    }

  if (MURGE_SUCCESS != (ret = ZMURGE_AssemblyEnd(id)))
    return ret;

  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_VALUES_OK);
  return MURGE_SUCCESS;
}


/*******************************************************************************
 * Group: Filling the right-hand-side member
 */


INTS ZMURGE_SetGlobalRHS(INTS id, COEF *b, INTS root, INTS op) {
  pastix_int_t        iter;
  pastix_int_t        procnum;
  pastix_int_t        localn;
  pastix_int_t        lastn = 0;
  pastix_int_t       *l2g   = NULL;
  pastix_complex64_t     *tmpb  = NULL;
  pastix_complex64_t    (*func)(pastix_complex64_t , pastix_complex64_t);
  MPI_Status status;
  int        dof;
  int        iterdof;
  COEF      *b_recv;
  int        allocated = API_NO;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);

  if (root < -1 || root >= zmurge_solvers[id]->pastix_data->procnbr) {
    errorPrint("Invalid root value");
    return MURGE_ERR_PARAMETER;
  }
  CHECK_L2G(id);

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (root != -1) {
    /* Broadcast and then run algorithm for root = -1 */
    /* TODO : fix mistake and remove BCast            */
    MURGE_MEMALLOC(b_recv, zmurge_solvers[id]->N*dof, pastix_complex64_t);
    if (root ==  zmurge_solvers[id]->pastix_data->procnum) {
      memcpy(b_recv, b, zmurge_solvers[id]->N*dof*sizeof(pastix_complex64_t));
    }
    MPI_Bcast(b_recv, zmurge_solvers[id]->N*dof, COMM_FLOAT, root,
              zmurge_solvers[id]->pastix_data->pastix_comm);
    root = -1;
    allocated = API_YES;
  }
  else {
    b_recv = b;
  }

  if (root == -1) {
    /* Si tous le monde à la racine */
    if (NULL == zmurge_solvers[id]->b) {
      MURGE_MEMALLOC(zmurge_solvers[id]->b, zmurge_solvers[id]->n*dof, pastix_complex64_t);
      memset(zmurge_solvers[id]->b, 0, zmurge_solvers[id]->n*dof*sizeof(pastix_complex64_t));
    }
    for (iter = 0; iter < zmurge_solvers[id]->n; iter++) {
      for (iterdof = 0; iterdof < dof; iterdof++) {

        MURGE_ADD(zmurge_solvers[id]->b[iter*dof +iterdof],
                  b_recv[(zmurge_solvers[id]->l2g[iter]-1)*dof+iterdof],
                  op);
      }
    }
  }
  else
    {
      /* Sinon, on recupère sur la racine tous les loc2globs
       On construit et on envoi tous les right-hand-side locaux
       */
      if (root == zmurge_solvers[id]->pastix_data->procnum) {
        lastn = 0;
        for (procnum = 0;
             procnum < zmurge_solvers[id]->pastix_data->procnbr;
             procnum++) {
          if (procnum != zmurge_solvers[id]->pastix_data->procnum) {

            MPI_Recv(&localn, 1, PASTIX_MPI_INT, procnum, TAG_SIZE,
                     zmurge_solvers[id]->pastix_data->pastix_comm, &status);

            if (lastn < localn) {
              if (NULL != l2g)
                MURGE_FREE(l2g);
              MURGE_MEMALLOC(l2g, localn, pastix_int_t);
              if (tmpb != NULL)
                MURGE_FREE(tmpb);
              MURGE_MEMALLOC(tmpb, localn*dof, pastix_complex64_t);
              lastn = localn;
            }
            MPI_Recv(l2g, localn, PASTIX_MPI_INT, procnum, TAG_L2G,
                     zmurge_solvers[id]->pastix_data->pastix_comm, &status);

            for (iter = 0; iter < localn; iter++) {
              for (iterdof = 0; iterdof < dof; iterdof++) {
                tmpb[iter*dof+iterdof]= b[(l2g[iter]-1)*dof+iterdof];
              }
            }
            MPI_Send(tmpb,
                     localn*dof, COMM_FLOAT, procnum,
                     TAG_VAL, zmurge_solvers[id]->pastix_data->pastix_comm);
          }
        }
        if (NULL != tmpb)
          MURGE_FREE(tmpb);
        if (NULL != l2g)
          MURGE_FREE(l2g);

        /* Le processeur racine construit son bout de RHS */
        if (NULL == zmurge_solvers[id]->b) {
          MURGE_MEMALLOC(zmurge_solvers[id]->b, zmurge_solvers[id]->n*dof, pastix_complex64_t);
          memset(zmurge_solvers[id]->b, 0, zmurge_solvers[id]->n*dof*sizeof(pastix_complex64_t));
        }

        for (iter = 0; iter < zmurge_solvers[id]->n; iter++) {
          for (iterdof = 0; iterdof < dof; iterdof++) {
            MURGE_ADD(zmurge_solvers[id]->b[iter*dof+iterdof],
                      b[(zmurge_solvers[id]->l2g[iter]-1)*dof+iterdof],
                      op);
          }
        }
      }
      else
        {
          /* Sur les procs non racine on recoit simplement le RHS a ajouter */
          MPI_Send(&zmurge_solvers[id]->n,
                   1, PASTIX_MPI_INT, root,
                   TAG_SIZE, zmurge_solvers[id]->pastix_data->pastix_comm);
          MPI_Send(zmurge_solvers[id]->l2g,
                   zmurge_solvers[id]->n,
                   PASTIX_MPI_INT, root,
                   TAG_L2G, zmurge_solvers[id]->pastix_data->pastix_comm);

          MURGE_MEMALLOC(tmpb, zmurge_solvers[id]->n, pastix_complex64_t);
          if (NULL == zmurge_solvers[id]->b) {
            zmurge_solvers[id]->b = tmpb;
          }

          MPI_Recv(tmpb,
                   zmurge_solvers[id]->n*dof, COMM_FLOAT, root,
                   TAG_VAL, zmurge_solvers[id]->pastix_data->pastix_comm, &status);

          if (tmpb != zmurge_solvers[id]->b) {
            for (iter = 0; iter < zmurge_solvers[id]->n*dof; iter++) {
              MURGE_ADD(zmurge_solvers[id]->b[iter],
                        tmpb[iter], op);
            }
            MURGE_FREE(tmpb);
          }
        }
    }
  if (allocated == API_YES) MURGE_FREE(b_recv);
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_RHS_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_SOLVE_DONE);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_REFINE_DONE);
  return MURGE_SUCCESS;
}

INTS ZMURGE_SetLocalRHS (INTS id, COEF *b, INTS op, INTS op2) {
  pastix_int_t        iter;
  pastix_complex64_t    (*func)(pastix_complex64_t , pastix_complex64_t) = NULL;
  int        dof;
#ifdef CENTRALISED
  pastix_int_t        nodenbr;
  pastix_int_t       *intern_nodelist;
  pastix_complex64_t     *tmpb;
#endif

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);
  CHECK_L2G(id);

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

#ifdef CENTRALISED
  nodenbr = z_pastix_getLocalNodeNbr(&(zmurge_solvers[id]->pastix_data));
  MURGE_MEMALLOC(intern_nodelist, nodenbr, pastix_int_t);
  if (NO_ERR != ( z_pastix_getLocalNodeLst(&(zmurge_solvers[id]->pastix_data),
                                         intern_nodelist)))
    return MURGE_ERR_SOLVER;

  if (NULL == zmurge_solvers[id]->b) {
    MURGE_MEMALLOC(zmurge_solvers[id]->b, zmurge_solvers[id]->N*dof, pastix_complex64_t);
    memset(zmurge_solvers[id]->b, 0, zmurge_solvers[id]->N*dof*sizeof(pastix_complex64_t));
  }
  MURGE_MEMALLOC(tmpb, zmurge_solvers[id]->N*dof /* zmurge_solvers[id]->n*dof */, pastix_complex64_t);
  memset(tmpb, 0, zmurge_solvers[id]->N*dof*sizeof(pastix_complex64_t));

  for (iter = 0; iter < nodenbr*dof; iter++) {
    MURGE_ADD(tmpb[(intern_nodelist[(iter-iter%dof)/dof]-1)*dof+iter%dof],
              b[iter], op);
  }
  MPI_Allreduce(tmpb, zmurge_solvers[id]->b, zmurge_solvers[id]->N*dof, COMM_FLOAT, MPI_SUM,
                zmurge_solvers[id]->pastix_data->pastix_comm);

  MURGE_FREE(tmpb);
  MURGE_FREE(intern_nodelist);

#else /* CENTRALISED */
  if (NULL == zmurge_solvers[id]->b) {
    MURGE_MEMALLOC(zmurge_solvers[id]->b, zmurge_solvers[id]->n*dof , pastix_complex64_t);
    memset(zmurge_solvers[id]->b, 0, zmurge_solvers[id]->n*dof*sizeof(pastix_complex64_t));
  }

  for (iter = 0; iter < zmurge_solvers[id]->n*dof; iter++) {
    MURGE_ADD(zmurge_solvers[id]->b[iter],
              b[iter], op);
  }
#endif /* CENTRALISED */

  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_RHS_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_SOLVE_DONE);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_REFINE_DONE);
  return MURGE_SUCCESS;
}

static inline
INTS ZMURGE_SetRHS_local_ (INTS id, INTS n, INTS *coefsidx, COEF *b, INTS op) {
  pastix_int_t             iter;
  pastix_complex64_t         (*func)(pastix_complex64_t , pastix_complex64_t) = NULL;
  pastix_int_t             index;
  int                    baseval;
  int                    dof;
  int                    iterdof;
  z_pastix_data_t *        pastix_data = zmurge_solvers[id]->pastix_data;

  baseval = pastix_data->iparm[IPARM_BASEVAL];
  dof = pastix_data->iparm[IPARM_DOF_NBR];

  for (iter = 0; iter < n; iter++) {
    index = coefsidx[iter]- baseval;
    index = zmurge_solvers[id]->g2l[index] - 1;

    for (iterdof = 0; iterdof < dof; iterdof++) {
      MURGE_ADD(zmurge_solvers[id]->b[index*dof+iterdof],
                b[iter*dof+iterdof], op);
    }
  }
  return MURGE_SUCCESS;
}

INTS ZMURGE_SetRHS      (INTS id, INTS n, INTS *coefsidx, COEF *b, INTS op,
                        INTS op2, INTS mode) {
  pastix_int_t             iter;
  pastix_complex64_t         (*func)(pastix_complex64_t , pastix_complex64_t) = NULL;
  int             baseval;
  int             dof;
  int             iterdof;
  z_pastix_data_t * pastix_data = zmurge_solvers[id]->pastix_data;
  INTS            procnbr     = pastix_data->procnbr;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);
  CHECK_L2G(id);

  baseval = pastix_data->iparm[IPARM_BASEVAL];
  dof = pastix_data->iparm[IPARM_DOF_NBR];

  for (iter = 0; iter < n; iter++)
    if (coefsidx[iter]-baseval >= zmurge_solvers[id]->N) {
      errorPrint("coefsidx[%ld] (%ld) is out of range [1-%ld]",
                 (long)iter, (long)coefsidx[iter]-baseval, (long)(zmurge_solvers[id]->N));
      return MURGE_ERR_PARAMETER;

    }

  if (NULL == zmurge_solvers[id]->b) {
    MURGE_MEMALLOC(zmurge_solvers[id]->b, zmurge_solvers[id]->n*dof , pastix_complex64_t);
    memset(zmurge_solvers[id]->b, 0, zmurge_solvers[id]->n*dof*sizeof(pastix_complex64_t));
  }

  if (mode == MURGE_ASSEMBLY_FOOL) {
    INTS            coefs_rcv_size;
    INTS         *  coefnbr;
    INTS         ** coefs_idx;
    COEF         ** coefs_vals;
    MPI_Request  *  request_cfnbr;
    MPI_Request  *  request_cfidx;
    MPI_Request  *  request_cfvals;

    MURGE_MEMALLOC(coefnbr, procnbr, INTS);
    for (iter = 0; iter <procnbr; iter++)
      coefnbr[iter] = 0;

    /* Count the entries to send to each processor */
    for (iter = 0; iter <n; iter++) {
      pastix_int_t procnum;
      if (zmurge_solvers[id]->g2l[coefsidx[iter]-baseval] > 0)
        procnum = pastix_data->procnum;
      else
        procnum = -zmurge_solvers[id]->g2l[coefsidx[iter]-baseval];

      coefnbr[procnum]++;
    }
    MURGE_MEMALLOC(coefs_idx,  procnbr, INTS*);
    MURGE_MEMALLOC(coefs_vals, procnbr, COEF*);

    for (iter = 0; iter < procnbr; iter++) {
      MURGE_MEMALLOC(coefs_idx[iter],  coefnbr[iter],     INTS);
      MURGE_MEMALLOC(coefs_vals[iter], coefnbr[iter]*dof, COEF);
      coefnbr[iter] = 0;
    }
    /* Prepare the arrays to send to each processors */
    for (iter = 0; iter <n; iter++) {
      pastix_int_t procnum;
      if (zmurge_solvers[id]->g2l[coefsidx[iter]-baseval] > 0)
        procnum = pastix_data->procnum;
      else
        procnum = -zmurge_solvers[id]->g2l[coefsidx[iter]-baseval];

      coefs_idx[procnum][coefnbr[procnum]] = coefsidx[iter];
      for (iterdof = 0; iterdof < dof; iterdof++) {
        coefs_vals[procnum][coefnbr[procnum]*dof+iterdof] =
          b[iter*dof+iterdof];
      }

      coefnbr[procnum]++;
    }

    MURGE_MEMALLOC(request_cfnbr,  procnbr, MPI_Request);
    MURGE_MEMALLOC(request_cfidx,  procnbr, MPI_Request);
    MURGE_MEMALLOC(request_cfvals, procnbr, MPI_Request);

    /* Send the data to the processors */
    for (iter = 0; iter < procnbr; iter++) {
      if (iter == zmurge_solvers[id]->pastix_data->procnum) continue;
      MPI_Isend(&(coefnbr[iter]),    1,                 MPI_INTS,
                iter, TAG_SIZE, pastix_data->pastix_comm,
                &(request_cfnbr[iter]));
      if (coefnbr[iter] > 0) {
        MPI_Isend(coefs_idx[iter],     coefnbr[iter],     MPI_INTS,
                  iter, TAG_ROW, pastix_data->pastix_comm,
                  &(request_cfidx[iter]));
        MPI_Isend(coefs_vals[iter],    coefnbr[iter]*dof, MURGE_MPI_COEF,
                  iter, TAG_VAL, pastix_data->pastix_comm,
                  &(request_cfvals[iter]));
      }
    }

    /* receive the data and run MPI_SetRHS with MURGE_ASSEMBLY_RESPECT */
    {
      coefs_rcv_size = 0;
      INTS       * coefs_idx_rcv  = NULL;
      COEF       * coefs_vals_rcv = NULL;

      for (iter = 0; iter < procnbr; iter++) {
        INTS         coefnbr_rcv;
        MPI_Status   status;
        if (iter == zmurge_solvers[id]->pastix_data->procnum) {
          ZMURGE_SetRHS_local_(id, coefnbr[iter], coefs_idx[iter], coefs_vals[iter],
                              op);
        } else {
          MPI_Recv(&coefnbr_rcv, 1, MPI_INTS, iter, TAG_SIZE,
                   pastix_data->pastix_comm, &status);

          if (coefnbr_rcv > 0) {
            if (coefnbr_rcv > coefs_rcv_size) {
              if (coefs_rcv_size != 0) {
                MURGE_FREE(coefs_idx_rcv);
                MURGE_FREE(coefs_vals_rcv);
              }
              MURGE_MEMALLOC(coefs_idx_rcv,  coefnbr_rcv,         INTS);
              MURGE_MEMALLOC(coefs_vals_rcv, coefnbr_rcv*dof,     COEF);
              coefs_rcv_size = coefnbr_rcv;
            }
            MPI_Recv(coefs_idx_rcv, coefnbr_rcv,     MPI_INTS,
                     iter, TAG_ROW, pastix_data->pastix_comm, &status);
            MPI_Recv(coefs_vals_rcv,coefnbr_rcv*dof, MURGE_MPI_COEF,
                     iter, TAG_VAL, pastix_data->pastix_comm, &status);

            ZMURGE_SetRHS_local_(id, coefnbr_rcv, coefs_idx_rcv, coefs_vals_rcv, op2);
          }
        }
        if (iter == procnbr-1 && coefs_rcv_size != 0) {
          MURGE_FREE(coefs_idx_rcv);
          MURGE_FREE(coefs_vals_rcv);
        }
      }
    }
    /* Now we clean it all */
    for (iter = 0; iter < procnbr; iter++) {
      MPI_Status status;
      if (iter == zmurge_solvers[id]->pastix_data->procnum) continue;
      MPI_Wait(&(request_cfnbr[iter]), &status);
      if (coefnbr[iter] > 0) {
        MPI_Wait(&(request_cfidx[iter]), &status);
        MPI_Wait(&(request_cfvals[iter]), &status);
      }
    }
    MURGE_FREE(request_cfnbr);
    MURGE_FREE(request_cfidx);
    MURGE_FREE(request_cfvals);

    for (iter = 0; iter <procnbr; iter++) {
      MURGE_FREE(coefs_idx[iter]);
      MURGE_FREE(coefs_vals[iter]);
    }
    MURGE_FREE(coefs_idx);
    MURGE_FREE(coefs_vals);
    MURGE_FREE(coefnbr);
  }
  else
    {
      if (!MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_NODELST_OK)) {
        errorPrint("Assembly can be respected only if user asked for node list");
        return MURGE_ERR_ORDER;
      }
      ZMURGE_SetRHS_local_(id, n, coefsidx, b, op);
    }

  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_RHS_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_SOLVE_DONE);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_REFINE_DONE);
  return MURGE_SUCCESS;
}

INTS ZMURGE_RHSReset(INTS id){
  pastix_int_t iter, dof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);
  CHECK_L2G(id);

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (NULL == zmurge_solvers[id]->b)
    MURGE_MEMALLOC(zmurge_solvers[id]->b, zmurge_solvers[id]->n*dof, pastix_complex64_t);

  for (iter = 0; iter < zmurge_solvers[id]->n*dof; iter++)
    zmurge_solvers[id]->b[iter] =0;
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_RHS_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_SOLVE_DONE);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_REFINE_DONE);
  return MURGE_SUCCESS;
}


/*******************************************************************************
 * Group: Getting the solution
 */
static inline
INTS ZMURGE_GetGlobalSolution_(INTS id, COEF *x, INTS root) {
  int    dof;
#ifndef CENTRALISED
  pastix_complex64_t *tmpx = NULL;
  pastix_int_t    iter;
  int    iterdof;
#endif
#ifdef MURGE_TIME
  Clock       clock;
#endif
  CLOCK_INIT;
  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if ((root == zmurge_solvers[id]->pastix_data->procnum ||
       root == -1)
      && NULL == x)
    return MURGE_ERR_PARAMETER;
#ifdef CENTRALISED
  if ((root == zmurge_solvers[id]->pastix_data->procnum ||
       root == -1))
    memcpy(x, zmurge_solvers[id]->b, zmurge_solvers[id]->N*dof*sizeof(pastix_complex64_t));
#else
  MURGE_MEMALLOC(tmpx, zmurge_solvers[id]->N*dof, pastix_complex64_t);

  for (iter = 0; iter < zmurge_solvers[id]->N*dof; iter ++)
    tmpx[iter] = 0;
  for (iter = 0; iter < zmurge_solvers[id]->n; iter ++) {
    for (iterdof = 0; iterdof < dof; iterdof++) {
      tmpx[(zmurge_solvers[id]->l2g[iter]-1)*dof+iterdof] =
        zmurge_solvers[id]->b[iter*dof+iterdof];
    }
  }
  CLOCK_PRINT("MURGE_GetGlobalSolution_ -- copy data");
  if (root == -1) {
  /*   INTS          i; */
  /*   PASTIX_INT   *n_rcv, *disp; */
  /*   PASTIX_INT   *l2g_rcv; */
  /*   PASTIX_FLOAT *rhs_rcv; */
  /*   MURGE_MEMALLOC(n_rcv,   2*zmurge_solvers[id]->pastix_data->procnbr, PASTIX_INT); */
  /*   disp = n_rcv + zmurge_solvers[id]->pastix_data->procnbr; */
  /*   MURGE_MEMALLOC(l2g_rcv, zmurge_solvers[id]->N, PASTIX_INT); */
  /*   MURGE_MEMALLOC(rhs_rcv, zmurge_solvers[id]->N*dof, PASTIX_FLOAT); */
  /*   for (i = 0; i < zmurge_solvers[id]->pastix_data->procnbr; i++) { */
  /*     INTS j,k,d; */
  /*     MPI_Allgather(&(zmurge_solvers[id]->n), 1, COMM_INT, */
  /*                   n_rcv, 1, COMM_INT, */
  /*                   zmurge_solvers[id]->pastix_data->pastix_comm); */
  /*     CLOCK_PRINT("MURGE_GetGlobalSolution_ -- MPI_Allgather"); */
  /*     disp[0] = 0; */
  /*     j = 0; */
  /*     fprintf(stdout, "n_rcv[%d] %d disp[%d] %d\n", j, n_rcv[j], j, disp[j]); */
  /*     for (j = 1; j < zmurge_solvers[id]->pastix_data->procnbr; j++) { */
  /*       disp[j] = disp[j-1] + n_rcv[j-1]; */
  /*       fprintf(stdout, "n_rcv[%d] %d disp[%d] %d\n", j, n_rcv[j], j, disp[j]); */
  /*     } */
  /*     MPI_Allgatherv(zmurge_solvers[id]->l2g, zmurge_solvers[id]->n, COMM_INT, */
  /*                    l2g_rcv, n_rcv, disp, COMM_INT, */
  /*                   zmurge_solvers[id]->pastix_data->pastix_comm); */
  /*     CLOCK_PRINT("MURGE_GetGlobalSolution_ -- MPI_Allgatherv"); */
  /*     for (j = 0; j < zmurge_solvers[id]->pastix_data->procnbr; j++) { */
  /*       n_rcv[j] *=dof; */
  /*       disp[j] *=dof; */
  /*     } */
  /*     MPI_Allgatherv(zmurge_solvers[id]->b, zmurge_solvers[id]->n*dof, COMM_INT, */
  /*                    rhs_rcv, n_rcv, disp, COMM_INT, */
  /*                    zmurge_solvers[id]->pastix_data->pastix_comm); */
  /*     CLOCK_PRINT("MURGE_GetGlobalSolution_ -- MPI_Allgatherv 2"); */
  /*     for (j = 0; j < zmurge_solvers[id]->pastix_data->procnbr; j++) { */
  /*       n_rcv[j] /=dof; */
  /*       disp[j] /=dof; */
  /*     } */
  /*     for (k = 0; k < zmurge_solvers[id]->pastix_data->procnbr; k++) { */
  /*       for (j = 0; j < n_rcv[k]; j++) { */
  /*         for (d=0; d<dof; d++) { */
  /*           x[(l2g_rcv[disp[k]+j]-1)*dof+d] = rhs_rcv[(disp[k]+j)*dof+d]; */
  /*         } */
  /*       } */
  /*     } */
      /* INTS j,d; */
      /* if (i == zmurge_solvers[id]->pastix_data->procnum) { */
      /*   n_rcv   = &(zmurge_solvers[id]->n); */
      /*   l2g_rcv = zmurge_solvers[id]->l2g; */
      /*   rhs_rcv = zmurge_solvers[id]->b; */
      /* } else { */
      /*   n_rcv   = &(n_tmp); */
      /*   l2g_rcv = l2g_tmp; */
      /*   rhs_rcv = rhs_tmp; */
      /* } */
      /* MPI_Bcast(n_rcv, 1, COMM_INT, */
      /*           i, zmurge_solvers[id]->pastix_data->pastix_comm); */
      /* MPI_Bcast(l2g_rcv, *n_rcv, COMM_INT, */
      /*           i, zmurge_solvers[id]->pastix_data->pastix_comm); */
      /* MPI_Bcast(rhs_rcv, (*n_rcv)*dof, COMM_FLOAT, */
      /*           i, zmurge_solvers[id]->pastix_data->pastix_comm); */
      /* for (j = 0; j < *n_rcv; j++) */
      /*   for (d=0; d<dof; d++) */
      /*     x[(l2g_rcv[j]-1)*dof+d] = rhs_rcv[j*dof+d]; */
    /* } */
    MPI_Allreduce(tmpx, x, zmurge_solvers[id]->N*dof, COMM_FLOAT, COMM_SUM,
                  zmurge_solvers[id]->pastix_data->pastix_comm);
  } else {
    /* if (root != zmurge_solvers[id]->pastix_data->procnum) { */
    /*   MPI_Send(&(zmurge_solvers[id]->n), 1, COMM_INT, */
    /*            root, 12345, zmurge_solvers[id]->pastix_data->pastix_comm); */
    /*   MPI_Send(zmurge_solvers[id]->l2g, zmurge_solvers[id]->n, COMM_INT, */
    /*            root, 12346, zmurge_solvers[id]->pastix_data->pastix_comm); */
    /*   MPI_Send(zmurge_solvers[id]->b, zmurge_solvers[id]->n*dof, COMM_FLOAT, */
    /*            root, 12347, zmurge_solvers[id]->pastix_data->pastix_comm); */
    /* } else { */
    /*   PASTIX_INT   *n_rcv, n_tmp, i, j, d; */
    /*   PASTIX_INT   *l2g_rcv, *l2g_tmp; */
    /*   PASTIX_FLOAT *rhs_rcv, *rhs_tmp; */
    /*   MURGE_MEMALLOC(l2g_tmp, zmurge_solvers[id]->N, INTS); */
    /*   MURGE_MEMALLOC(rhs_tmp, zmurge_solvers[id]->N*dof, PASTIX_FLOAT); */
    /*   for (i = 0; i < zmurge_solvers[id]->pastix_data->procnbr; i++) { */
    /*     if (i == root) { */
    /*       n_rcv   = &(zmurge_solvers[id]->n); */
    /*       l2g_rcv = zmurge_solvers[id]->l2g; */
    /*       rhs_rcv = zmurge_solvers[id]->b; */
    /*     } else { */
    /*       MPI_Status status; */
    /*       n_rcv   = &n_tmp; */
    /*       l2g_rcv = l2g_tmp; */
    /*       rhs_rcv = rhs_tmp; */
    /*       MPI_Recv(n_rcv, 1, COMM_INT, */
    /*                i, 12345, zmurge_solvers[id]->pastix_data->pastix_comm, &status); */
    /*       MPI_Recv(l2g_rcv, *n_rcv, COMM_INT, */
    /*                i, 12346, zmurge_solvers[id]->pastix_data->pastix_comm, &status); */
    /*       MPI_Recv(rhs_rcv, (*n_rcv)*dof, COMM_FLOAT, */
    /*                i, 12347, zmurge_solvers[id]->pastix_data->pastix_comm, &status); */
    /*     } */
    /*     for (j = 0; j < *n_rcv; j++) */
    /*       for (d=0; d<dof; d++) */
    /*         x[(l2g_rcv[j]-1)*dof+d] = rhs_rcv[j*dof+d]; */
    /*   } */
    /* } /\* i'm the root *\/ */
      MPI_Reduce(tmpx, x, zmurge_solvers[id]->N*dof, COMM_FLOAT, COMM_SUM, root,
                 zmurge_solvers[id]->pastix_data->pastix_comm);
  } /* root == -1 */

  CLOCK_PRINT("MURGE_GetGlobalSolution_ -- reduction");

  MURGE_FREE(tmpx);
#endif
  CLOCK_PRINT("MURGE_GetGlobalSolution_ -- end");
  return MURGE_SUCCESS;
}

INTS ZMURGE_GetGlobalSolution(INTS id, COEF *x, INTS root) {
#ifdef MURGE_TIME
  Clock       clock;
#endif
  CLOCK_INIT;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CLOCK_PRINT("MURGE_GetGlobalSolution -- CHECK_SOLVER_PARAM");
  CHECK_L2G(id);
  CLOCK_PRINT("MURGE_GetGlobalSolution -- CHECK_L2G");
  CHECK_FACT(id);
  CLOCK_PRINT("MURGE_GetGlobalSolution -- CHECK_FACT");
  ZMURGE_GetGlobalSolution_(id, x, root);
  CLOCK_PRINT("MURGE_GetGlobalSolution -- END");
  return MURGE_SUCCESS;
}

INTS ZMURGE_GetLocalSolution (INTS id, COEF *x) {
  pastix_int_t    iter;
  int    dof;
#ifdef CENTRALISED
  pastix_int_t    nodenbr;
  pastix_int_t   *intern_nodelist;
  int    iterdof;
#endif
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);
  CHECK_FACT(id);

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  if (NULL == x) {
    return MURGE_ERR_PARAMETER;
  }
#ifdef CENTRALISED
  nodenbr = z_pastix_getLocalNodeNbr(&(zmurge_solvers[id]->pastix_data));
  MURGE_MEMALLOC(intern_nodelist, nodenbr, pastix_int_t);
  if (NO_ERR != ( z_pastix_getLocalNodeLst(&(zmurge_solvers[id]->pastix_data),
                                           intern_nodelist)))
    return MURGE_ERR_SOLVER;

  for (iter = 0; iter < nodenbr; iter ++) {
    for (iterdof = 0; iterdof < dof; iterdof++) {
      x[iter*dof+iterdof] = zmurge_solvers[id]->b[(intern_nodelist[iter]-1)*
                                           dof+iterdof];
    }
  }
  MURGE_FREE(intern_nodelist);
#else
  for (iter = 0; iter < zmurge_solvers[id]->n*dof; iter ++) {
    x[iter] = zmurge_solvers[id]->b[iter];
  }
#endif
  return MURGE_SUCCESS;
}

INTS ZMURGE_GetSolution      (INTS id, INTS n, INTS *coefsidx, COEF *x,
                             INTS mode) {
  pastix_int_t    iter;
  COEF  *tmpx = NULL;
  int    dof;
  int    iterdof;
#ifdef MURGE_TIME
  Clock       clock;
#endif
  CLOCK_INIT;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CLOCK_PRINT("MURGE_GetSolution -- CHECK_SOLVER_PARAM");
  CHECK_L2G(id);
  CLOCK_PRINT("MURGE_GetSolution -- CHECK_L2G");
  CHECK_FACT(id);
  CLOCK_PRINT("MURGE_GetSolution -- CHECK_FACT");

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  MURGE_MEMALLOC(tmpx, zmurge_solvers[id]->N*dof, COEF);
  ZMURGE_GetGlobalSolution_(id, tmpx, -1);
  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];
  for (iter = 0; iter < n; iter ++) {
    for (iterdof = 0; iterdof < dof; iterdof++) {
      x[iter*dof+iterdof] = tmpx[(coefsidx[iter]-
                                  zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL])*
                                 dof+iterdof];
    }
  }
  MURGE_FREE(tmpx);
  CLOCK_PRINT("MURGE_GetSolution -- END");


  return MURGE_SUCCESS;
}

/*******************************************************************************
 * Group: Cleaning up this mess
 */

INTS ZMURGE_Clean(INTS id){
  pastix_int_t    * iparm;
  double * dparm;
#ifdef MEMORY_USAGE
  int rank, verb;
  MPI_Comm comm;
  rank =  zmurge_solvers[id]->pastix_data->procnum;
  comm =  zmurge_solvers[id]->pastix_data->pastix_comm;
  if (zmurge_solvers[id]->pastix_data->iparm)
    verb =  zmurge_solvers[id]->pastix_data->iparm[IPARM_VERBOSE];
  else
    verb = API_VERBOSE_NOT;
#endif
  iparm = zmurge_solvers[id]->pastix_data->iparm;
  dparm = zmurge_solvers[id]->pastix_data->dparm;
  if (NULL != zmurge_solvers[id]->colptr)
    MURGE_FREE(zmurge_solvers[id]->colptr);
  if (NULL != zmurge_solvers[id]->rows)
    MURGE_FREE(zmurge_solvers[id]->rows);
  if (NULL != zmurge_solvers[id]->values)
    MURGE_FREE(zmurge_solvers[id]->values);
  if (NULL != zmurge_solvers[id]->b)
    MURGE_FREE(zmurge_solvers[id]->b);
  if (NULL != zmurge_solvers[id]->l2g)
    MURGE_FREE(zmurge_solvers[id]->l2g);
  if (NULL != zmurge_solvers[id]->g2l)
    MURGE_FREE(zmurge_solvers[id]->g2l);
  if (NULL != zmurge_solvers[id]->perm)
    MURGE_FREE(zmurge_solvers[id]->perm);
#ifdef CENTRALISED
  if (NULL != zmurge_solvers[id]->invp)
    MURGE_FREE(zmurge_solvers[id]->invp);
#endif
  if (NULL != zmurge_solvers[id]->dropmask)
    MURGE_FREE(zmurge_solvers[id]->dropmask);
  if (NULL != zmurge_solvers[id]->dropcols)
    MURGE_FREE(zmurge_solvers[id]->dropcols);
  if (NULL != zmurge_solvers[id]->droprows)
    MURGE_FREE(zmurge_solvers[id]->droprows);

  while (NULL != zmurge_solvers[id]->sequences) {
    ZMURGE_AssemblyDeleteSequence(id, zmurge_solvers[id]->sequences->ID);
  }

  if (zmurge_solvers[id]->threadnbr)
    z_stop_threads(id);
  if (NULL != zmurge_solvers[id]->pastix_data)
      z_pastix_task_clean(&(zmurge_solvers[id]->pastix_data),
                          zmurge_solvers[id]->pastix_data->pastix_comm);

#ifndef FORCE_NOSMP
  pthread_mutex_destroy(&(zmurge_solvers[id]->barrier.sync_lock));
  pthread_cond_destroy(&(zmurge_solvers[id]->barrier.sync_cond));
#endif

#ifdef MURGE_THREADSAFE
  pthread_mutex_destroy(&zmurge_solvers[id]->mutex_tmpmatrix);
#endif
  MURGE_FREE_EXT(iparm, IPARM_SIZE*sizeof(pastix_int_t));
  MURGE_FREE_EXT(dparm, DPARM_SIZE*sizeof(double));
  MURGE_MEMORY_USAGE_PRINT2("MURGE_Clean", verb, rank, comm);
  free(zmurge_solvers[id]); zmurge_solvers[id] = NULL;
  /* Hack for clean restart */
  _ZMURGE_InitId(id);
  return MURGE_SUCCESS;
}


INTS ZMURGE_Finalize(){
  pastix_int_t i;

  for (i=0; i< idnbr; i++) {
      ZMURGE_Clean(i);
      if (NULL != zmurge_solvers[i]->pastix_data)
          z_pastix_task_clean(&(zmurge_solvers[i]->pastix_data),
                              zmurge_solvers[i]->pastix_data->pastix_comm);
      free(zmurge_solvers[i]); zmurge_solvers[i] = NULL;
  }

  free(zmurge_solvers); zmurge_solvers = NULL;

  return MURGE_SUCCESS;
}

INTS ZMURGE_GetInfoINT(INTS id,  INTS metric, INTL * value) {
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  pastix_int_t murge_param[1];

  murge_param[MURGE_IINFO_NNZ  - 1024] =  IPARM_NNZEROS;

  if (metric >= 1024) {
    metric = murge_param[metric-1024];
  }

  if (!( metric < IPARM_SIZE )) {
    errorPrint("metric is too big");
    return MURGE_ERR_PARAMETER;
  }

  if (metric < 0) {
    errorPrint("metric is negative");
    return MURGE_ERR_PARAMETER;
  }

  /* TODO : Est-ce qu'on ajoute des tests sur les valeurs rentrées ???? */
  *value = zmurge_solvers[id]->pastix_data->iparm[metric];

  return MURGE_SUCCESS;
}

INTS ZMURGE_GetInfoREAL(INTS id,  INTS metric, REAL * value) {
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  pastix_int_t murge_param[1];

  murge_param[MURGE_RPARAM_EPSILON_ERROR  - 1024] =  DPARM_RELATIVE_ERROR;

  if (metric >= 1024) {
    metric = murge_param[metric-1024];
  }

  if (!( metric < IPARM_SIZE )) {
    errorPrint("metric is too big");
    return MURGE_ERR_PARAMETER;
  }

  if (metric < 0) {
    errorPrint("metric is negative");
    return MURGE_ERR_PARAMETER;
  }

  *value = (REAL)zmurge_solvers[id]->pastix_data->dparm[metric];
  return MURGE_SUCCESS;
}
/*
 Function: MURGE_PrintError

 Print the error message corresponding to ierror
 Parameters:
 error_number  - Error identification number.

 Returns:
 MURGE_ERR_PARAMETER - If ierror does not match an error number
 MURGE_SUCCESS       - If function runned successfully.

 Fortran interface:
 >
 > SUBROUTINE MURGE_PRINTERROR(ERROR_NUMBER, IERROR)
 >   INTS, INTENT(IN)  :: IERROR
 >   INTS, INTENT(OUT) :: ERROR_NUMBER
 > END SUBROUTINE MURGE_PRINTERROR
 */
INTS ZMURGE_PrintError(INTS error_number) {
  return MURGE_SUCCESS;
}

/*
 Function: MURGE_ExitOnError

 Print the error message corresponding to ierror.
 If the ierr is not MURGE_SUCCESS then the program is stopped.

 Parameters:
 ierror         - Error identification number.

 Returns:
 MURGE_SUCCESS   - If function runned successfully,
 stop the program otherwise.

 Fortran interface:
 >
 > SUBROUTINE MURGE_EXITONERROR(ERROR_NUMBER, IERROR)
 >   INTS, INTENT(IN)  :: IERROR
 >   INTS, INTENT(OUT) :: ERROR_NUMBER
 > END SUBROUTINE MURGE_EXITONERROR
 */
INTS ZMURGE_ExitOnError(INTS error_number) {
  if  (error_number == MURGE_SUCCESS)
    return MURGE_SUCCESS;
  else
    exit(1);
}


/*
 Group: Scaling
 */

/*
 Function: ZMURGE_GetGlobalNorm

 Compute the global norm array following a norm rule.

 Must be performed after assembly step.

 Parameters:
 id      - Solver instance identification number.
 norm    - Array of size global column number*dof which will contain
 the norm values
 root    - Indicates which processor will have the norm array
 at the end of the call, -1 for all.
 rule    - Rule to follow to build norm array, see <MURGE_NORM_RULES>

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_GETGLOBALNORM(ID, NORM, ROOT, RULE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, ROOT, RULE
 >   REAL, DIMENSION(0), INTENT(OUT) :: NORM
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_GETGLOBALNORM
 */
INTS ZMURGE_GetGlobalNorm(INTS id, REAL *norm, INTS root, INTS rule) {
  INTS itercol;       /* Each column*/
  pastix_int_t  iterrow;       /* each row entry in each column of the CSCd */
  INTS iterdof_col;   /* each dof on column */
  INTS iterdof_row;   /* each dof on row    */
  INTS scal_idx;    /* Column number of the value */
  pastix_int_t  value_idx;     /* Index of the value */
  REAL*local_norm = NULL;
  INTS dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  print_debug(DBG_MURGE, ">> ZMURGE_GetGlobalNorm\n");
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_NORM_RULE(rule);


  MURGE_MEMALLOC(local_norm, zmurge_solvers[id]->N*dof, REAL);

  for(itercol = 0; itercol <  zmurge_solvers[id]->N*dof; itercol++) {
    local_norm[itercol] = 0;
  }
  for(itercol = 0; itercol <  zmurge_solvers[id]->n; itercol++) {
    for (iterrow = zmurge_solvers[id]->colptr[itercol]-1;
         iterrow < zmurge_solvers[id]->colptr[itercol+1]-1;
         iterrow++) {
      for (iterdof_col = 0;
           iterdof_col < dof;
           iterdof_col++) {
        for (iterdof_row = 0;
             iterdof_row < dof;
             iterdof_row++) {
          value_idx  = iterdof_row + iterdof_col*dof + iterrow*dof*dof;
          switch(rule) {
          case MURGE_NORM_MAX_COL:
            scal_idx = iterdof_col + (zmurge_solvers[id]->l2g[itercol]-1) * dof;
            local_norm[scal_idx] =
              MAX(local_norm[scal_idx],
                  ABS_FLOAT(zmurge_solvers[id]->values[value_idx]));
            break;

          case MURGE_NORM_MAX_ROW:
            scal_idx = iterdof_row + (zmurge_solvers[id]->rows[iterrow]-1) * dof;
            local_norm[scal_idx] =
              MAX(local_norm[scal_idx],
                  ABS_FLOAT(zmurge_solvers[id]->values[value_idx]));
            break;

          case MURGE_NORM_2_COL:
            scal_idx = iterdof_col + (zmurge_solvers[id]->l2g[itercol]-1) * dof;
            local_norm[scal_idx] = local_norm[scal_idx] +
              (REAL)(zmurge_solvers[id]->values[value_idx]*
                     CONJ_FLOAT(zmurge_solvers[id]->values[value_idx]));
            break;

          case MURGE_NORM_2_ROW:
            scal_idx = iterdof_row + (zmurge_solvers[id]->rows[iterrow]-1) * dof;
            local_norm[scal_idx] = local_norm[scal_idx] +
              (REAL)(zmurge_solvers[id]->values[value_idx]*
                     CONJ_FLOAT(zmurge_solvers[id]->values[value_idx]));
            break;

          default:
            errorPrint("Rule not implemented");
            return MURGE_ERR_NOT_IMPLEMENTED;
          }
        }
      }
    }
  }


  if (rule == MURGE_NORM_2_COL ||
      rule == MURGE_NORM_2_ROW) {
    fprintf(stderr, "Reduce on norm\n");
    MPI_Allreduce(local_norm,norm,zmurge_solvers[id]->N*dof,
                  MURGE_MPI_REAL,
                  MPI_SUM,
                  zmurge_solvers[id]->pastix_data->pastix_comm);
    for(itercol = 0; itercol <  zmurge_solvers[id]->N*dof; itercol++) {
      local_norm[itercol] = (REAL)sqrt((double)local_norm[itercol]);
    }
  }
  else
    {
      fprintf(stderr, "Reduce on norm 2\n");
      MPI_Allreduce(local_norm,norm,zmurge_solvers[id]->N*dof,
                    MURGE_MPI_REAL,
                    MPI_MAX,
                    zmurge_solvers[id]->pastix_data->pastix_comm);
    }
  return MURGE_SUCCESS;
}

/*
 Function: ZMURGE_GetLocalNorm

 Compute the local norm array following a norm rule.

 Must be performed after assembly step.

 Parameters:
 id      - Solver instance identification number.
 norm    - Array of size local column number*dof which will contain
 the solution
 rule    - Rule to follow to build norm array, see <MURGE_NORM_RULES>

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_GETLOCALNORM(ID, NORM, RULE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, RULE
 >   REAL, DIMENSION(0), INTENT(OUT) :: NORM
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_GETLOCALNORM
 */
INTS ZMURGE_GetLocalNorm(INTS id, REAL *norm, INTS rule){
  INTS itercol;       /* Each column*/
  pastix_int_t  iterrow;       /* each row entry in each column of the CSCd */
  INTS iterdof_col;   /* each dof on column */
  INTS iterdof_row;   /* each dof on row    */
  INTS column_dof;    /* Column number of the value */
  pastix_int_t  value_idx;     /* Index of the value */
  INTS dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_NORM_RULE(rule);
  CHECK_L2G(id);

  if (rule == MURGE_NORM_MAX_ROW || rule == MURGE_NORM_2_ROW) {
    errorPrint("PaStiX uses column distribution, local norm can't be a norm on rows");
    return MURGE_ERR_PARAMETER;
  }

  for(itercol = 0; itercol <  zmurge_solvers[id]->n*dof; itercol++)
    norm[itercol] = 0;
  for(itercol = 0; itercol <  zmurge_solvers[id]->n; itercol++) {
    for (iterrow = zmurge_solvers[id]->colptr[itercol]-1;
         iterrow < zmurge_solvers[id]->colptr[itercol+1]-1;
         iterrow++) {
      for (iterdof_col = 0;
           iterdof_col < dof;
           iterdof_col++) {
        for (iterdof_row = 0;
             iterdof_row < dof;
             iterdof_row++) {
          column_dof = iterdof_col + itercol * dof;
          value_idx  = iterdof_row + iterdof_col*dof + iterrow*dof*dof;

          if (rule == MURGE_NORM_2_COL) {
            norm[column_dof] = norm[column_dof] +
              (REAL)(zmurge_solvers[id]->values[value_idx]*
                     CONJ_FLOAT(zmurge_solvers[id]->values[value_idx]));
          }
          else
            {
              norm[column_dof] =
                MAX(norm[column_dof],
                    ABS_FLOAT(zmurge_solvers[id]->values[value_idx]));
            }
        }
      }
    }
  }

  for(itercol = 0; itercol <  zmurge_solvers[id]->n*dof; itercol++)
    norm[itercol] = (REAL)sqrt(norm[itercol]);

  return MURGE_SUCCESS;
}


/*
 Function: ZMURGE_GetNorm

 Compute the indicated part of the norm array
 following a norm rule.

 Must be performed after assembly step.


 Parameters:
 id       - Solver instance identification number.
 n        - Number of coefficients user wants to get norm of.
 coefsidx - List of the coefficients user wants to get norm of.
 norm     - Array of size dof*n which will contain
 the solution.
 rule     - Rule to follow to build norm array, see <MURGE_NORM_RULES>
 mode     - Indicates if the user is sure to respect the distribution.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_GETNORM(ID, N, COEFSIDX, NORM, RULE, MODE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, MODE, N, RULE
 >   INTS, DIMENSION(0), INTENT(IN)  :: COEFSIDX
 >   COEF, DIMENSION(0), INTENT(OUT) :: NORM
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_GETNORM
 */
INTS ZMURGE_GetNorm(INTS id,  INTS n, INTS *coefsidx, REAL *norm, INTS rule, INTS mode){
  errorPrint("Not yet implemented");
  return MURGE_ERR_NOT_IMPLEMENTED;

}


/*
 Function: MURGE_ApplyGlobalScaling

 Apply scaling to local unknowns.

 Must be performed after assembly step.

 Parameters:
 id      - Solver instance identification number.
 scal    - Scaling user wants to apply.
 sc_mode - Indicate if the scaling is applied on rows or on columns.
 root    - Indicates which processor that posses the scaling array,
 -1 for all.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_APPLYGLOBALSCALING(ID, SCAL, SC_MODE, ROOT, IERROR)
 >   INTS,               INTENT(IN)  :: ID, ROOT, SC_MODE
 >   REAL, DIMENSION(0), INTENT(OUT) :: SCAL
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_APPLYGLOBALSCALING

 */
INTS MURGE_ApplyGlobalScaling(INTS id, REAL *scal, INTS root, INTS sc_mode){
  INTS itercol;       /* Each column*/
  pastix_int_t  iterrow;       /* each row entry in each column of the CSCd */
  INTS iterdof_col;   /* each dof on column */
  INTS iterdof_row;   /* each dof on row    */
  INTS scal_idx;      /* Scaling array index */
  pastix_int_t  value_idx;     /* Index of the value */
  INTS dof;
  REAL *scaling = NULL;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_SCALING_MODE(sc_mode);
  CHECK_L2G(id);

  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (root == -1) {
    scaling = scal;
  }
  else
    {
      if (root != (zmurge_solvers[id]->pastix_data)->procnum) {
        MURGE_MEMALLOC(scaling, zmurge_solvers[id]->N*dof, REAL);
      }
      else
        {
          scaling = scal;
        }
      MPI_Bcast( scaling, zmurge_solvers[id]->N*dof,
                 MURGE_MPI_REAL,
                 root,
                 zmurge_solvers[id]->pastix_data->pastix_comm);
    }

  for(itercol = 0; itercol <  zmurge_solvers[id]->n; itercol++) {
    for (iterrow = zmurge_solvers[id]->colptr[itercol]-1;
         iterrow < zmurge_solvers[id]->colptr[itercol+1]-1;
         iterrow++) {
      for (iterdof_col = 0;
           iterdof_col < dof;
           iterdof_col++) {
        for (iterdof_row = 0;
             iterdof_row < dof;
             iterdof_row++) {
          if (sc_mode == MURGE_SCAL_COL) {
            scal_idx = iterdof_col + (zmurge_solvers[id]->l2g[itercol]-1) * dof;
          }
          else
            {
              scal_idx = iterdof_row + (zmurge_solvers[id]->rows[iterrow]-1) * dof;
            }
          value_idx  = iterdof_row + iterdof_col*dof + iterrow*dof*dof;
          zmurge_solvers[id]->values[value_idx] =
            zmurge_solvers[id]->values[value_idx] /scaling[scal_idx];
        }
      }
    }
  }
  if (root != -1 && root != (zmurge_solvers[id]->pastix_data)->procnum)
    MURGE_FREE(scaling);

  return MURGE_SUCCESS;
}

/*
 Function: MURGE_ApplyLocalScaling

 Apply the local scaling array on the matrix.

 Must be performed after assembly step.

 Parameters:
 id      - Solver instance identification number.
 scal    - Array of size local column number*dof which will contain
 the solution.
 sc_mode - Indicate if the scaling is applied on rows or on columns.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_APPLYLOCALSCALING(ID, SCAL, SC_MODE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, SC_MODE
 >   REAL, DIMENSION(0), INTENT(OUT) :: SCAL
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_APPLYLOCALSCALING
 */
INTS MURGE_ApplyLocalScaling(INTS id, REAL *scal, INTS sc_mode){
  INTS itercol;       /* Each column*/
  pastix_int_t  iterrow;       /* each row entry in each column of the CSCd */
  INTS iterdof_col;   /* each dof on column */
  INTS iterdof_row;   /* each dof on row    */
  INTS scal_idx;      /* Index in scaling array */
  pastix_int_t  value_idx;     /* Index of the value */
  INTS dof;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_SCALING_MODE(sc_mode);
  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (sc_mode == MURGE_SCAL_ROW) {
    /*
     * Building global to local column number array
     */
    CHECK_L2G(id);
  }
  for(itercol = 0; itercol <  zmurge_solvers[id]->n; itercol++) {
    for (iterrow = zmurge_solvers[id]->colptr[itercol]-1;
         iterrow < zmurge_solvers[id]->colptr[itercol+1]-1;
         iterrow++) {
      for (iterdof_col = 0;
           iterdof_col < dof;
           iterdof_col++) {
        for (iterdof_row = 0;
             iterdof_row < dof;
             iterdof_row++) {

          if (sc_mode == MURGE_SCAL_COL) {
            scal_idx = iterdof_col + itercol * dof;
          }
          else
            {
              scal_idx = iterdof_row + (zmurge_solvers[id]->g2l[zmurge_solvers[id]->rows[iterrow]-1]-1)*dof;
            }
          value_idx  = iterdof_row + iterdof_col*dof + iterrow*dof*dof;
          zmurge_solvers[id]->values[value_idx] =
            zmurge_solvers[id]->values[value_idx] /scal[scal_idx];
        }
      }
    }
  }
  return MURGE_SUCCESS;
}

/*
 Function: MURGE_ApplyScaling

 Apply the scaling array on the indicated part of the matrix

 Must be performed after assembly step.


 Parameters:
 id       - Solver instance identification number.
 n        - Number of coefficients user wants to scale.
 coefsidx - List of the coefficients user wants to scale.
 scal     - Array of size dof*n which will contain
 the solution.
 sc_mode  - Indicate if the scaling is applied on rows or on columns.
 mode     - Indicates if the user is sure to respect the distribution.

 Returns:
 MURGE_SUCCESS       - If function runned successfully.
 MURGE_ERR_PARAMETER - If *id* is not in solver arrays range.
 MURGE_ERR_ORDER     - If the assembly has not been performed.

 Fortran interface:
 >
 > SUBROUTINE MURGE_APPLYSCALING(ID, N, COEFSIDX, SCAL, SC_MODE, MODE, IERROR)
 >   INTS,               INTENT(IN)  :: ID, SC_MODE, MODE, N
 >   INTS, DIMENSION(0), INTENT(IN)  :: COEFSIDX
 >   COEF, DIMENSION(0), INTENT(OUT) :: SCAL
 >   INTS,               INTENT(OUT) :: IERROR
 > END SUBROUTINE MURGE_APPLYSCALING
 */
INTS MURGE_ApplyScaling(INTS id,  INTS n, INTS *coefsidx, REAL *scal,
                        INTS sc_mode, INTS mode){
  INTS itercol;       /* Each column*/
  INTS iterdof_col;   /* each dof on column */
  REAL *scaling = NULL;
  INTS dof;
  INTS baseval;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_SCALING_MODE(sc_mode);
  baseval = zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (mode == MURGE_ASSEMBLY_RESPECT) {
    /*
     * Building global to local column number array
     */
    CHECK_L2G(id);

    MURGE_MEMALLOC(scaling, zmurge_solvers[id]->n*dof, REAL);
    for (itercol = 0; itercol < zmurge_solvers[id]->n*dof; itercol++)
      scaling[itercol] = 1.0;
    for (itercol = 0; itercol < n; itercol++)
      for (iterdof_col =0; iterdof_col < dof; iterdof_col++)
        scaling[(zmurge_solvers[id]->g2l[coefsidx[itercol]-baseval]-1)*dof+iterdof_col] = scal[itercol*dof+iterdof_col];
    MURGE_ApplyLocalScaling(id, scaling, sc_mode);
    MURGE_FREE(scaling);
  }
  else
    {
      REAL * scaling_recv = NULL;
      MURGE_MEMALLOC(scaling, zmurge_solvers[id]->N*dof, REAL);
      for (itercol = 0; itercol < zmurge_solvers[id]->N*dof; itercol++)
        scaling[itercol] = 0.0;
      for (itercol = 0; itercol < n; itercol++)
        for (iterdof_col =0; iterdof_col < dof; iterdof_col++)
          scaling[(coefsidx[itercol]-baseval)*dof+iterdof_col] = scal[itercol*dof+iterdof_col];

      MURGE_MEMALLOC(scaling_recv, zmurge_solvers[id]->N*dof, REAL);
      MPI_Allreduce(scaling, scaling_recv,
                    zmurge_solvers[id]->N*dof,
                    MURGE_MPI_REAL,
                    MPI_SUM,
                    zmurge_solvers[id]->pastix_data->pastix_comm);
      MURGE_FREE(scaling);
      for (itercol = 0; itercol < zmurge_solvers[id]->N*dof; itercol++)
        if (scaling_recv[itercol] == 0.0)
          scaling_recv[itercol] = 1.0;

      for (itercol = 0; itercol < n; itercol++) {
        for (iterdof_col =0; iterdof_col < dof; iterdof_col++) {
          if (scaling[(coefsidx[itercol]-baseval)*dof+iterdof_col] != scal[itercol*dof+iterdof_col]) {
            errorPrint("Multiple entries for the same scaling entry");
            return MURGE_ERR_PARAMETER;
          }
        }
      }


      MURGE_ApplyGlobalScaling(id, scaling_recv, sc_mode, -1);
      MURGE_FREE(scaling_recv);
    }
  return MURGE_SUCCESS;
}





/******************************************************************************
 * Group: Specific PaStiX functions.                                          *
 ******************************************************************************/

/******************************************************************************
 * Function: MURGE_Analyze                                                    *
 *                                                                            *
 * Perform matrix analyze:                                                    *
 *   - Compute a new ordering of the unknows                                  *
 *   - Compute the symbolic factorisation of the matrix                       *
 *   - Distribute column blocks and computation on processors                 *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned succesfuly.                     *
 *   MURGE_ERR_ORDER     - If function the graph is not built.                *
 *   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_Analyze(INTS id){
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_PREPROCESSING(id);

  return MURGE_SUCCESS;
}

/******************************************************************************
 * Function: MURGE_Factorize                                                  *
 *                                                                            *
 * Perform matrix factorization.                                              *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *                                                                            *
 * Returns:                                                                   *
 *   MURGE_SUCCESS       - If function runned succesfuly.                     *
 *   MURGE_ERR_ORDER     - If function the graph is not built.                *
 *   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_Factorize(INTS id){
  z_pastix_data_t   *pastix_data = zmurge_solvers[id]->pastix_data;
  pastix_int_t             *iparm       = pastix_data->iparm;

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);

  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_VALUES_OK))) {
    errorPrint("Need to set values before.");
    return MURGE_ERR_ORDER;
  }
  if (iparm[IPARM_ONLY_RAFF] ==  API_YES) {
    errorPrint("MURGE_Factorize is not compatible with IPARM_ONLY_RAFF == API_YES\n");
    return MURGE_ERR_PARAMETER;
  }

  /* 
   * - fill the internal CSC
   * - delete murge CSC and
   * - perform factorization
   */
  PASTIX_FILLIN_CSC(zmurge_solvers[id]->pastix_data,
		    zmurge_solvers[id]->pastix_data->pastix_comm,
		    zmurge_solvers[id]->n,
		    zmurge_solvers[id]->colptr,
		    zmurge_solvers[id]->rows,
		    zmurge_solvers[id]->values,
		    NULL,
		    0,
		    zmurge_solvers[id]->l2g);
  zmurge_solvers[id]->pastix_data->cscInternFilled = API_YES;

  MURGE_FREE(zmurge_solvers[id]->colptr);
  MURGE_FREE(zmurge_solvers[id]->rows);
  MURGE_FREE(zmurge_solvers[id]->values);

  pastix_data->iparm[IPARM_START_TASK] = API_TASK_NUMFACT;
  iparm[IPARM_END_TASK]                = API_TASK_NUMFACT;

  DPASTIX(&pastix_data,
          pastix_data->pastix_comm,
          zmurge_solvers[id]->n,
          zmurge_solvers[id]->colptr,
          zmurge_solvers[id]->rows,
          zmurge_solvers[id]->values,
          zmurge_solvers[id]->l2g,
          zmurge_solvers[id]->perm,
          NULL,
          NULL,
          zmurge_solvers[id]->nrhs,
          pastix_data->iparm,
          pastix_data->dparm);
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_FACTO_OK);

  return MURGE_SUCCESS;
}

/******************************************************************************
 * Function: MURGE_ForceNoFacto                                               *
 *                                                                            *
 * Prevent Murge from running factorisation even if matrix has changed.       *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 * Returns:                                                                   *
 *   MURGE_SUCCESS                                                            *
 *                                                                            *
 ******************************************************************************/
INTS MURGE_ForceNoFacto(INTS id) {
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_FACTO_OK);
  return MURGE_SUCCESS;
}

/******************************************************************************
 * Function: MURGE_ProductSetLocalNodeNbr                                     *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   n  - Number of local nodes.                                              *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_ProductSetLocalNodeNbr (INTS id, INTS n) {
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  zmurge_solvers[id]->n = n;
  MPI_Allreduce(&zmurge_solvers[id]->n,
                &zmurge_solvers[id]->N, 1, PASTIX_MPI_INT,
                MPI_SUM, zmurge_solvers[id]->pastix_data->pastix_comm);
  return MURGE_SUCCESS;
}


/******************************************************************************
 * Function: MURGE_ProductSetGlobalNodeNbr                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   N  - Number of global nodes.                                             *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_ProductSetGlobalNodeNbr (INTS id, INTS N) {
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  if (MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_NODELST_OK)) {
    errorPrint("%s must be called before MURGE_ProductSetLocalNodeList",
               __FUNCTION__);
    return MURGE_ERR_ORDER;
  }
  zmurge_solvers[id]->N = N;
  return MURGE_SUCCESS;
}
/******************************************************************************
 * Function: MURGE_ProductSetLocalNodeList                                    *
 *                                                                            *
 * Parameters:                                                                *
 *   id  - Solver instance identification number.                             *
 *   l2g - Local to global node numbers.                                      *
 *                                                                            *
 ******************************************************************************/

INTS MURGE_ProductSetLocalNodeList (INTS id, INTS * l2g) {
  INTS i;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  MURGE_MEMALLOC(zmurge_solvers[id]->l2g, zmurge_solvers[id]->n, pastix_int_t);
  for (i = 0; i < zmurge_solvers[id]->n; i++) {
    zmurge_solvers[id]->l2g[i] = l2g[i];
  }

  z_cscd_build_g2l(zmurge_solvers[id]->n,
                 zmurge_solvers[id]->l2g,
                 zmurge_solvers[id]->pastix_data->pastix_comm,
                 &zmurge_solvers[id]->N,
                 &zmurge_solvers[id]->g2l);
  MURGE_TRACE_MALLOC(PTR_MEMSIZE(zmurge_solvers[id]->g2l), char);

  /* No need to work on the graph nor factorize
   when we only perform product */
  MURGE_STATE_TRUE (zmurge_solvers[id]->state, MURGE_ONLY_PROD);
  MURGE_STATE_TRUE (zmurge_solvers[id]->state, MURGE_GRAPH_OK);
  MURGE_STATE_TRUE (zmurge_solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_TRUE (zmurge_solvers[id]->state, MURGE_FACTO_OK);
  MURGE_STATE_TRUE (zmurge_solvers[id]->state, MURGE_NODENBR_OK);
  MURGE_STATE_TRUE (zmurge_solvers[id]->state, MURGE_NODELST_OK);
  MURGE_STATE_FALSE(zmurge_solvers[id]->state, MURGE_VALUES_OK);


  return MURGE_SUCCESS;
}

/******************************************************************************
 * Function: ZMURGE_GetLocalProduct                                            *
 *                                                                            *
 * Perform the product A * X.                                                 *
 *                                                                            *
 * The vector must have been given trough <ZMURGE_SetLocalRHS> or              *
 * <ZMURGE_SetGlobalRHS>.                                                      *
 *                                                                            *
 * Parameters:                                                                *
 *   id - Solver instance identification number.                              *
 *   x  - Array in which the local part of the product will be stored.        *
 * Returns:                                                                   *
 *   MURGE_ERR_ORDER  - If values have not been set.                          *
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
INTS ZMURGE_GetLocalProduct (INTS id, COEF *x) {
  COEF * glob_prod;
  INTS ierr, iter;
  MURGE_MEMALLOC(glob_prod, zmurge_solvers[id]->N, COEF);

  if (MURGE_SUCCESS != (ierr = ZMURGE_GetGlobalProduct(id, glob_prod, -1)))
    return ierr;

  for (iter  = 0; iter < zmurge_solvers[id]->n; iter++) {
    x[iter] = glob_prod[zmurge_solvers[id]->l2g[iter]];
  }
  return MURGE_SUCCESS;
}


/******************************************************************************
 * Function: ZMURGE_GetGlobalProduct                                           *
 *                                                                            *
 * Perform the product A * X.                                                 *
 *                                                                            *
 * The vector must have been given trough <ZMURGE_SetLocalRHS> or              *
 * <ZMURGE_SetGlobalRHS>.                                                      *
 *                                                                            *
 * Parameters:                                                                *
 *   id   - Solver instance identification number.                            *
 *   x    - Array in which the product will be stored.                        *
 *   root - Rank of the process which will own the product at end of call,    *
 *          use -1 for all processes.                                         *
 * Returns:                                                                   *
 *   MURGE_ERR_ORDER  - If values have not been set.                          *
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
INTS ZMURGE_GetGlobalProduct (INTS id, COEF *x, INTS root) {
  INTS dof;

#ifdef MURGE_TIME
  Clock            clock;
#endif

  CLOCK_INIT;
  if (!(MURGE_STATE_ISTRUE(zmurge_solvers[id]->state, MURGE_VALUES_OK))) {
    errorPrint("Need to set values before.");
    return MURGE_ERR_ORDER;
  }

  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);
  CHECK_L2G(id);
  dof = zmurge_solvers[id]->pastix_data->iparm[IPARM_DOF_NBR];

  if (zmurge_solvers[id]->pastix_data->iparm[IPARM_SYM] == API_SYM_YES) {
    errorPrint("Product only available with unsymmetric matrices.");
    return MURGE_ERR_NOT_IMPLEMENTED;
  }

#ifdef MURGE_PRODUCT_CHECK_ZEROS
  {
    INTS itercol, iterrows, row;
    COEF max, sum, norm, critere;
    pastix_int_t cnt, cnt2, cnt_sum, cnt2_sum, nz_glob;
    COEF *values = zmurge_solvers[id]->values;
    INTS baseval, iter, iter2;
    critere = zmurge_solvers[id]->pastix_data->dparm[DPARM_EPSILON_MAGN_CTRL];
    baseval = zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL];

    if (critere < 0.0) {
      critere = -critere;
    }
    else
      {
        max = 0;
        for (itercol = 0; itercol < zmurge_solvers[id]->n; itercol++) {
          for (iter = 0; iter < dof; iter++) {
            sum = 0;
            for (iterrows = zmurge_solvers[id]->colptr[itercol]-baseval;
                 iterrows < zmurge_solvers[id]->colptr[itercol+1]-baseval;
                 iterrows++) {
              row = zmurge_solvers[id]->rows[iterrows]-baseval;
              for (iter2 = 0; iter2 < dof; iter2++) {
                pastix_int_t idx = iterrows*dof*dof+iter*dof+iter2;

                sum = sum + ABS_FLOAT(values[idx]);
              }

            }
            max = MAX(max,sum);
          }
        }

        MPI_Allreduce(&max, &norm, 1, MURGE_MPI_COEF, MPI_MAX,
                      zmurge_solvers[id]->pastix_data->pastix_comm);

        critere = norm*sqrt(critere);
      }
    cnt = 0;
    cnt2 = 0;
    for (itercol = 0; itercol < zmurge_solvers[id]->n; itercol++) {
      for (iter = 0; iter < dof; iter++) {
        sum = 0;
        for (iterrows = zmurge_solvers[id]->colptr[itercol]-baseval;
             iterrows < zmurge_solvers[id]->colptr[itercol+1]-baseval;
             iterrows++) {
          row = zmurge_solvers[id]->rows[iterrows]-baseval;
          for (iter2 = 0; iter2 < dof; iter2++) {
            pastix_int_t idx = iterrows*dof*dof+iter*dof+iter2;
            if (ABS_FLOAT(values[idx]) < critere)
              cnt = cnt + 1;
            if (values[idx] ==  0.0)
              cnt2 = cnt2 + 1;
          }
        }
        max = MAX(max,sum);
      }
    }
    cnt_sum = 0;
    MPI_Reduce(&cnt, &cnt_sum, 1, PASTIX_MPI_INT, MPI_SUM, 0,
               zmurge_solvers[id]->pastix_data->pastix_comm);
    MPI_Reduce(&cnt2, &cnt2_sum, 1, PASTIX_MPI_INT, MPI_SUM, 0,
               zmurge_solvers[id]->pastix_data->pastix_comm);
    cnt = zmurge_solvers[id]->colptr[zmurge_solvers[id]->n]-1;
    MPI_Reduce(&cnt, &nz_glob, 1, PASTIX_MPI_INT, MPI_SUM, 0,
               zmurge_solvers[id]->pastix_data->pastix_comm);
    nz_glob = nz_glob *dof*dof;
    if ((zmurge_solvers[id]->pastix_data)->procnum == 0) {
      fprintf(stdout, "%d zeros in matrix from %d (%.3lg %%) critere :"
              " %.20lg\n",
              cnt_sum, nz_glob,
              (100.0*(double)(cnt_sum)/(double)(nz_glob)), critere);
      fprintf(stdout, "%d real zeros in matrix from %d (%.3lg %%)\n",
              cnt2_sum, nz_glob,
              (100.0*(double)(cnt2_sum)/(double)(nz_glob)));
    }
  }
#endif

  if (zmurge_solvers[id]->threadnbr != 0 &&
      zmurge_solvers[id]->threadnbr != zmurge_solvers[id]->pastix_data->iparm[IPARM_THREAD_NBR])
    z_stop_threads(id);
  if (zmurge_solvers[id]->threadnbr == 0)
    z_start_threads(id);

  /* Wait that the threads are ready then start PRODUCT */
  pthread_mutex_lock(&(zmurge_solvers[id]->mutex_state));
  while (zmurge_solvers[id]->threads_state != MURGE_THREAD_WAIT) {
    pthread_cond_wait(&(zmurge_solvers[id]->cond_state),
                      &(zmurge_solvers[id]->mutex_state));
  }
  zmurge_solvers[id]->threads_state = MURGE_THREAD_PRODUCT;
  pthread_cond_broadcast(&(zmurge_solvers[id]->cond_state));

  /* Wait end of task */
  while (zmurge_solvers[id]->threads_state != MURGE_THREAD_WAIT) {
    pthread_cond_wait(&(zmurge_solvers[id]->cond_state),
                      &(zmurge_solvers[id]->mutex_state));
  }
  pthread_mutex_unlock(&(zmurge_solvers[id]->mutex_state));

  if (zmurge_solvers[id]->threads_data[0].pdata->ret != MURGE_SUCCESS)
    return zmurge_solvers[id]->threads_data[0].pdata->ret;

  if (root == -1) {
    INTS t;
    for (t = 0; t < zmurge_solvers[id]->threadnbr; t++) {
      INTS size = zmurge_solvers[id]->threads_data[t].pdata->last-
        zmurge_solvers[id]->threads_data[t].pdata->first;
      MPI_Allreduce(zmurge_solvers[id]->threads_data[t].pdata->t_prod,
                    x+zmurge_solvers[id]->threads_data[t].pdata->first*dof,
                    size*dof, MURGE_MPI_COEF, MPI_SUM,
                  zmurge_solvers[id]->pastix_data->pastix_comm);
    }
  }
  else {
    INTS t;
    for (t = 0; t < zmurge_solvers[id]->threadnbr; t++) {
      INTS size = zmurge_solvers[id]->threads_data[t].pdata->last-
        zmurge_solvers[id]->threads_data[t].pdata->first;
      MPI_Reduce(zmurge_solvers[id]->threads_data[t].pdata->t_prod,
                 x+zmurge_solvers[id]->threads_data[t].pdata->first*dof,
                 size*dof, MURGE_MPI_COEF, MPI_SUM,
               root, zmurge_solvers[id]->pastix_data->pastix_comm);
  }
  }
  CLOCK_PRINT("MURGE_GetGlobalProduct");

  return MURGE_SUCCESS;
}

/*
 WARNING: NEEDS TO BE CHECKED !
 */
INTS ZMURGE_SetLocalNodeList   (INTS id, INTS nodenbr, INTS *nodelist) {
  pastix_int_t i;
  zmurge_solvers[id]->n = nodenbr;
  /* On detruit colptr et rows, la distribution a changé */
  MURGE_FREE(zmurge_solvers[id]->colptr);
  MURGE_FREE(zmurge_solvers[id]->rows);
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_GRAPH_OK)
    MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_BLEND_OK);
  MURGE_STATE_TRUE(zmurge_solvers[id]->state, MURGE_SYMB_OK);
  MURGE_MEMALLOC(zmurge_solvers[id]->l2g, nodenbr, pastix_int_t);
  for (i = 0; i < nodenbr; i++) {
    zmurge_solvers[id]->l2g[i] = nodelist[i];
  }
  return MURGE_SUCCESS;
}

INTS ZMURGE_GetCommRank(INTS id, int * rank) {
  CHECK_SOLVER_ID(id);
  *rank = (zmurge_solvers[id]->pastix_data)->procnum;
  return MURGE_SUCCESS;

}

INTS ZMURGE_GetCommSize(INTS id, int * size) {
  CHECK_SOLVER_ID(id);
  *size = (zmurge_solvers[id]->pastix_data)->procnbr;
  return MURGE_SUCCESS;
}

INTS ZMURGE_GetOptionINT(INTS id, INTS index, INTS * value) {
  pastix_int_t murge_param[64];
  pastix_int_t * iparm = NULL;
  CHECK_SOLVER_ID(id);
  CHECK_SOLVER_PARAM(id);

  iparm = zmurge_solvers[id]->pastix_data->iparm;

  murge_param[MURGE_IPARAM_BASEVAL       - 1024] =  IPARM_BASEVAL;
  murge_param[MURGE_IPARAM_DOF           - 1024] =  IPARM_DOF_NBR;
  murge_param[MURGE_IPARAM_SYM           - 1024] =  IPARM_SYM;

  if (index == MURGE_IPARAM_SYM) {
    if (iparm[IPARM_SYM] == API_SYM_YES)
      *value = MURGE_BOOLEAN_TRUE;
    else
      *value = MURGE_BOOLEAN_FALSE;
  }
  else
    {
      if (index >= 1024)
        index = murge_param[index-1024];

      *value = iparm[index];
    }
  return MURGE_SUCCESS;
}


INTS ZMURGE_GetComm(INTS id, MPI_Comm * comm) {
  CHECK_SOLVER_ID(id);
  * comm = zmurge_solvers[id]->pastix_data->pastix_comm;
  return MURGE_SUCCESS;
}

INTS ZMURGE_SetDropNodes(INTS id, INTS nodenbr, INTS * dropmask) {
  INTS i;
  CHECK_SOLVER_ID(id);
  MURGE_MEMALLOC(zmurge_solvers[id]->dropmask, nodenbr, char);
  for (i = 0; i < nodenbr; i++) {
    zmurge_solvers[id]->dropmask[i] = (char)dropmask[i];
  }
  return MURGE_SUCCESS;
}


INTS ZMURGE_SetDropRows(INTS id, INTS nodenbr, INTS * droprows) {
  INTS i;
  CHECK_SOLVER_ID(id);
  MURGE_MEMALLOC(zmurge_solvers[id]->droprows, nodenbr, char);
  for (i = 0; i < nodenbr; i++) {
    zmurge_solvers[id]->droprows[i] = (char)droprows[i];
  }
  return MURGE_SUCCESS;
}

INTS ZMURGE_SetDropCols(INTS id, INTS nodenbr, INTS * dropcols) {
  INTS i;
  CHECK_SOLVER_ID(id);
  MURGE_MEMALLOC(zmurge_solvers[id]->dropcols, nodenbr, char);
  for (i = 0; i < nodenbr; i++) {
    zmurge_solvers[id]->dropcols[i] = (char)dropcols[i];
  }
  return MURGE_SUCCESS;
}

INTS MURGE_ColGetNonZerosNbr(INTS id, INTS COL, INTS * nnzNbr) {
  INTS lcol;
  INTS base = zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  CHECK_SOLVER_ID(id);
  CHECK_PREPROCESSING(id);
  COL = COL - base;
  if (COL < 0 || COL > zmurge_solvers[id]->N) {
    errorPrint("invalid column index");
    return MURGE_ERR_PARAMETER;
  }
  lcol = zmurge_solvers[id]->g2l[COL];
  if (lcol > 0)
    *nnzNbr = (INTS)(zmurge_solvers[id]->colptr[lcol] - zmurge_solvers[id]->colptr[lcol-1]);
  else
    *nnzNbr = 0;
  return MURGE_SUCCESS;
}

INTS MURGE_ColGetNonZerosIdx(INTS id, INTS COL, INTS * indexes) {
  INTS lcol;
  INTS base = zmurge_solvers[id]->pastix_data->iparm[IPARM_BASEVAL];
  CHECK_SOLVER_ID(id);
  CHECK_PREPROCESSING(id);
  COL = COL - base;
  if (COL < 1 || COL > zmurge_solvers[id]->N) {
    errorPrint("invalid column index");
    return MURGE_ERR_PARAMETER;
  }
  lcol = zmurge_solvers[id]->g2l[COL];
  if (lcol > 0) {
    INTS i;
    for (i = zmurge_solvers[id]->colptr[lcol-1]-1; i < zmurge_solvers[id]->colptr[lcol]-1; i++)
      indexes[i- zmurge_solvers[id]->colptr[lcol-1]-1] = zmurge_solvers[id]->rows[i];
  }
  return MURGE_SUCCESS;
}
