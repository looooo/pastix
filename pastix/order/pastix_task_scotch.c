/*
 * File: pastix.c
 *
 * PaStiX external functions implementations.
 *
 * Authors:
 *   Mathieu FAVERGE  - faverge@labri.fr
 *   Xavier  LACOSTE  - lacoste@labri.fr
 *   Pierre  RAMET    - ramet@labri.fr
 */

#include "common.h"
#ifdef WITH_SCOTCH
#  ifdef    PASTIX_DISTRIBUTED
#    include <ptscotch.h>
#  else
#    include <scotch.h>
#  endif /* PASTIX_DISTRIBUTED */
#endif /* WITH_SCOTCH */

#include "dof.h"
#include "ftgt.h"
#include "symbol.h"
#include "csc.h"
#include "updown.h"
#include "queue.h"
#include "bulles.h"
#include "solver.h"
#include "assembly.h"
#include "param_blend.h"
#include "order.h"
#include "fax.h"
#include "kass.h"
#include "blend.h"
#include "solverRealloc.h"
#include "sopalin_thread.h"
#include "stack.h"
#include "sopalin3d.h"
#include "sopalin_init.h"
#include "sopalin_option.h"
#include "csc_intern_updown.h"
#include "csc_intern_build.h"
#include "coefinit.h"
#include "out.h"
#include "pastix_internal.h"

#include "csc_utils.h"
#include "cscd_utils.h"
#include "cscd_utils_intern.h"
#include "bordi.h"
#include "sopalin_acces.h"
#include "perf.h"

void global2localperm(pastix_int_t  lN,
                      pastix_int_t *lperm,
                      pastix_int_t *gperm,
                      pastix_int_t *loc2glob);

/* define pour l'affichage */
#define SCOTCH_STRAT_DIRECT                                             \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{rat=0.8,"                                  \
  ""                        "vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{rat=0.8,"                                  \
  ""                        "vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.0},"                       \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{rat=0.8,"                                 \
  ""                         "vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                        "m{rat=0.8,"                                \
  ""                          "vert=100,"                               \
  ""                          "low=h{pass=10},"                         \
  ""                          "asc=f{bal=0.2}};,"                       \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g}}"

#define SCOTCH_STRAT_INCOMP                                             \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.08},"                        \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                       "m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""            "ose=g}}"
#define SCOTCH_STRAT_PERSO                                              \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>%ld)?m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=%ld,cmax=%ld,frat=%f},"                           \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>%ld)?(m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                       "m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=%ld,cmax=%ld,frat=%f},"                           \
  ""      "ose=g}}"

#define PTSCOTCH_STRAT_DIRECT                                           \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{rat=0.8,"                                 \
  ""                        "vert=100,"                                \
  ""                        "low=h{pass=10},"                          \
  ""                        "asc=f{bal=0.2}}|"                         \
  ""                      "m{rat=0.8,"                                 \
  ""                        "vert=100,"                                \
  ""                        "low=h{pass=10},"                          \
  ""                        "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.0},"                         \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{type=h,"                                  \
  ""                         "rat=0.8,"                                 \
  ""                         "vert=100000,"                             \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=08.2}})|"                       \
  ""                       "m{type=h,"                                  \
  ""                         "rat=0.8,"                                 \
  ""                         "vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g}}"
#define PTSCOTCH_STRAT_INCOMP                                           \
  "c{rat=0.7,"                                                          \
  """cpr=n{sep=/(vert>120)?m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}}|"                          \
  ""                      "m{vert=100,"                                 \
  ""                        "low=h{pass=10},"                           \
  ""                        "asc=f{bal=0.2}};,"                         \
  ""      "ole=f{cmin=0,cmax=100000,frat=0.08},"                        \
  ""      "ose=g},"                                                     \
  """unc=n{sep=/(vert>120)?(m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}})|"                        \
  ""                       "m{vert=100,"                                \
  ""                         "low=h{pass=10},"                          \
  ""                         "asc=f{bal=0.2}};,"                        \
  ""      "ole=f{cmin=15,cmax=100000,frat=0.08},"                       \
  ""      "ose=g}}"
#define PTSCOTCH_STRAT_PERSO  "c{rat=0.7,cpr=n{sep=/(vert>%ld)?m{vert=100,low=h{pass=10},asc=f{bal=0.2}}|m{vert=100,low=h{pass=10},asc=f{bal=0.2}};,ole=f{cmin=%ld,cmax=%ld,frat=%f},ose=g},unc=n{sep=/(vert>%ld)?(m{vert=100,low=h{pass=10},asc=f{bal=0.2}})|m{vert=100,low=h{pass=10},asc=f{bal=0.2}};,ole=f{cmin=%ld,cmax=%ld,frat=%f},ose=g}}"


/*******************************************************************************
 *  Section: Macros
 */

/*
  macro: print_onempi

  Print a string using processor 0.
  Uses printf syntax.

  Parameters:
  fmt - Format string (see printf manual).
  ... - Arguments depending on the format string.
*/
#define print_onempi(fmt, ...) if(procnum == 0) fprintf(stdout, fmt, __VA_ARGS__)


/*
  Function: pastix_order_prepare_csc

  Create a copy of user's CSC and prepare it for ordering step.

  Symmetrize the graph and removes diagonal coefficients.

  Parameters:
  pastix_data - PaStiX internal data structure
  n           - Number of column in the CSC.
  colptr      - Start of each column in *rows* array.
  rows        - Row number of each non zeros.
*/
int pastix_order_prepare_csc(pastix_data_t * pastix_data,
                             pastix_int_t             n,
                             pastix_int_t           * colptr,
                             pastix_int_t           * rows)
{
  pastix_int_t * iparm;
  int procnum;

  iparm   = pastix_data->iparm;
  procnum = (int)pastix_data->procnum;
  /* Allocate a copy of col, row */
  pastix_data->bmalcolrow = 1;
  pastix_data->n2   = n;
  MALLOC_INTERN(pastix_data->col2, (n+1),         pastix_int_t);
  MALLOC_INTERN(pastix_data->row2, (colptr[n]-1), pastix_int_t);
  memcpy((void*) pastix_data->col2,(void*)colptr,        (n+1)*sizeof(pastix_int_t));
  memcpy((void*) pastix_data->row2,(void*)rows,  (colptr[n]-1)*sizeof(pastix_int_t));

  /* Symmetrize the graph */
  if (iparm[IPARM_SYM] == API_SYM_YES || iparm[IPARM_SYM] == API_SYM_HER)
    {
      pastix_int_t tmpn;
      pastix_int_t *tmpcol;
      pastix_int_t *tmprow;

      if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
        print_onempi("%s", OUT_SYMGRAPH);
      csc_symgraph_int(n,     pastix_data->col2,   pastix_data->row2,   NULL,
                       &tmpn, &tmpcol, &tmprow, NULL, API_YES);

      memFree_null((pastix_data->col2));
      memFree_null((pastix_data->row2));
      pastix_data->col2 = tmpcol;
      pastix_data->row2 = tmprow;
    }

  /* Remove diagonal coefficient */
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
    print_onempi("%s", OUT_NODIAG);
  csc_noDiag(1, n, pastix_data->col2, pastix_data->row2, NULL);

  return PASTIX_SUCCESS;
}
/*
  Function: pastix_task_scotch

  Execute ordering task, with a centralised graph.

  Free *col2*  and *row2* entries of pastix_data if <pastix_task_scotch>
  has already been called.

  Set *col2*, *row2* and *n2* to a copy of user's CSC.

  Symmetrize this CSC.

  Remove diagonal elements from it.

  Clean last oredering if it exists.
  Depending on *IPARM_ORDERING* :
  - Calls Scotch ordering,
  - Calls Metis ordering,
  - Uses user oredering,
  - Loads oredering stored on disk in a Scotch format.

  Can save computed ordering on disk.

  returns compuited ordering into user arrays.

  Parameters:
  pastix_data - PaStiX data structure.
  pastix_comm - PaStiX MPI communicator.
  n           - size of the matrix/number of vertices.
  colptr      - starting index of each column in row
  row         - row of each element
  perm        - permutation tabular
  invp        - reverse permutation tabular
*/
int pastix_task_scotch(pastix_data_t **pastix_data,
                       MPI_Comm        pastix_comm,
                       pastix_int_t    n,
                       pastix_int_t   *colptr,
                       pastix_int_t   *row,
                       pastix_int_t   *perm,
                       pastix_int_t   *invp)
{
  pastix_int_t              * iparm       = (*pastix_data)->iparm;
  pastix_int_t             * col2;
  pastix_int_t             * row2;
#ifdef WITH_SCOTCH
  SCOTCH_Strat       stratdat;
  char               strat[550];
  pastix_int_t               *colptr_schur  = NULL;
  pastix_int_t               *rows_schur    = NULL;
  pastix_int_t               *perm_schur    = NULL;
  pastix_int_t               *revperm_schur = NULL;
#endif
  Order            * ordemesh;
  double             timer1;
  pastix_int_t                procnum;
  pastix_int_t                iter;
  int                retval     = PASTIX_SUCCESS;
  int                retval_rcv;
  (void)pastix_comm;

  ordemesh  = &((*pastix_data)->ordemesh);
  procnum   =   (*pastix_data)->procnum;

  print_debug(DBG_STEP,"-> pastix_task_scotch\n");
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_onempi("%s", OUT_STEP_ORDER);

  /* Clean ordering/graph if it exists */
  if ((*pastix_data)->malord)
    {
      orderExit(ordemesh);
      (*pastix_data)->malord=0;
    }

  /* Prepare a copy of user's CSC */
  if (!(PASTIX_MASK_ISTRUE(iparm[IPARM_ORDERING], API_ORDER_LOAD)))
      orderPrepareCSC( *pastix_data, n, colptr, row, NULL );

  col2 = (*pastix_data)->col2;
  row2 = (*pastix_data)->row2;

  if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
    print_onempi("%s", OUT_ORDERINIT);

  clockInit(timer1);
  clockStart(timer1);

  switch (iparm[IPARM_ORDERING])
    {
      /*
       * Scotch Ordering
       */
    case API_ORDER_SCOTCH:
#ifndef WITH_SCOTCH
      errorPrint("Scotch ordering needs to compile PaStiX with -DWITH_SCOTCH");
      retval = BADPARAMETER_ERR;
      break;
#else
      orderComputeScotch( *pastix_data );
#endif
      break;


      /*
       *  METIS ordering
       */
    case API_ORDER_METIS:
#ifndef METIS
      errorPrint("Metis ordering needs to compile PaStiX with -DMETIS");
      retval = BADPARAMETER_ERR;
      break;
#else /* METIS */
      {
        pastix_int_t  itervert;
        pastix_int_t  baseval;
        pastix_int_t  opt[8];

        baseval = 1;
        orderInit(ordemesh, n, 0);

        if (sizeof(pastix_int_t) != sizeof(int))
          {
              errorPrint("Inconsistent integer type %lu != %lu\n",
                         sizeof(pastix_int_t), sizeof(SCOTCH_Num));
            retval = INTEGER_TYPE_ERR;
            break;
          }

        if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
          print_onempi("%s", "calling metis...\n");

        /* call METIS and fill ordemesh (provide a partition) */
        opt[OPTION_PTYPE  ] = (iparm[IPARM_DEFAULT_ORDERING]==API_YES)?0:1;

        /* TODO: tester sans cette ligne... 0 if default */
        opt[OPTION_PTYPE  ] = 0;

        opt[OPTION_CTYPE  ] = iparm[IPARM_ORDERING_SWITCH_LEVEL];
        opt[OPTION_ITYPE  ] = iparm[IPARM_ORDERING_CMIN];
        opt[OPTION_RTYPE  ] = iparm[IPARM_ORDERING_CMAX];
        opt[OPTION_DBGLVL ] = iparm[IPARM_ORDERING_FRAT];
        opt[OPTION_OFLAGS ] = iparm[IPARM_STATIC_PIVOTING];
        opt[OPTION_PFACTOR] = iparm[IPARM_METIS_PFACTOR];
        opt[OPTION_NSEPS  ] = iparm[IPARM_NNZEROS];

        /*METIS_NodeND(&n,verttab,edgetab,&baseval,opt,
          ordemesh->permtab,ordemesh->peritab);*/
        METIS_NodeND(&n, col2, row2, &baseval,
                     opt, ordemesh->peritab, ordemesh->permtab);
      }
#endif /* METIS */
      break;

      /*
       * Personal Ordering
       */
    case API_ORDER_PERSONAL:
    {
        orderInit(ordemesh, n, 0);
        memcpy(ordemesh->permtab, perm, n*sizeof(pastix_int_t));
        memcpy(ordemesh->peritab, invp, n*sizeof(pastix_int_t));

        orderLoadFiles( *pastix_data );
    }
    break;

      /*
       * Load ordering with Scotch Format
       */
    case API_ORDER_LOAD:
        orderLoadFiles( *pastix_data );
        break;

    default:
      errorPrint("Ordering not available");
      retval = BADPARAMETER_ERR;
      break;
    }

  fprintf(stderr, "The number of supernodes found is %ld\n", ordemesh->cblknbr );

  MPI_Allreduce(&retval, &retval_rcv, 1, MPI_INT, MPI_MAX, pastix_comm);
  if (retval_rcv != PASTIX_SUCCESS)
    return retval_rcv;

  (*pastix_data)->malord = 1;
  orderBase(ordemesh, 0);

  clockStop(timer1);
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
    print_onempi(TIME_COMPUTE_ORDERING,clockVal(timer1));

  /* Save i/o strategy */
  if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE))
  {
      retval = orderSaveFiles( *pastix_data );
      if (retval != PASTIX_SUCCESS)
        return retval;
  }

  /*
   * Return the ordering to user
   */
  if (iparm[IPARM_ORDERING] != API_ORDER_PERSONAL)
    {
      memcpy(perm, ordemesh->permtab, n*sizeof(pastix_int_t));
      memcpy(invp, ordemesh->peritab, n*sizeof(pastix_int_t));
    }

  iparm[IPARM_START_TASK]++;
  return PASTIX_SUCCESS;
}

#ifdef PASTIX_DISTRIBUTED
/*
  Function: dpastix_order_prepare_cscd

  Create a copy of user's CSCd and prepare it for ordering step.

  Symmetrize the graph and removes diagonal coefficients.

  Symetrize the graph, removes diagonal elements.

  Redistribute the CSCd to be abble to give it to scotch if needed.
  Indeed, PT-Scotch only allows consecutive columns.
  *PTS_permtab* is the permutation tabular from the user's distribution
  to PT-Scotch one, *PTS_peritab* is the reverse permutation.

  Parameters:
  pastix_data  - PaStiX internal data structure
  n            - Number of column in the CSC.
  colptr       - Start of each column in *rows* array.
  rows         - Row number of each non zeros.
  loc2glob     - local to global column number array.
  pastix_comm  - MPI communicator.
*/
int dpastix_order_prepare_cscd(pastix_data_t * pastix_data,
                               pastix_int_t             n,
                               pastix_int_t           * colptr,
                               pastix_int_t           * rows,
                               pastix_int_t           * loc2glob,
                               MPI_Comm        pastix_comm)
{
  pastix_int_t   gN;
  pastix_int_t * iparm;
  pastix_int_t   i;
  int   OK, OKRecv;
  int * alln         = NULL;
  int   nlocal;
  int * displs       = NULL;

  iparm = pastix_data->iparm;
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_CHATTERBOX)
    fprintf(stdout, "> dpastix_order_prepare_cscd\n");

  if (pastix_data->PTS_permtab != NULL)
    {
      memFree_null(pastix_data->PTS_permtab);
      memFree_null(pastix_data->PTS_peritab);
    }
  /* Copy the cscd */
  pastix_data->n2 = n;
  MALLOC_INTERN(pastix_data->col2, n+1, pastix_int_t);
  memcpy(pastix_data->col2, colptr, (n+1)*sizeof(pastix_int_t));
  MALLOC_INTERN(pastix_data->loc2glob2, n, pastix_int_t);
  memcpy((pastix_data->loc2glob2), loc2glob, n*sizeof(pastix_int_t));
  MALLOC_INTERN(pastix_data->row2, colptr[n]-1, pastix_int_t);
  memcpy(pastix_data->row2, rows, (colptr[n]-1)*sizeof(pastix_int_t));

  gN = 0;
  MPI_Allreduce(&n, &gN, 1, PASTIX_MPI_INT, MPI_SUM, pastix_comm);

  /* Symmetrize the graph */
  if (iparm[IPARM_SYM]==API_SYM_YES || iparm[IPARM_SYM] == API_SYM_HER)
    {
      /* Symetric problem */
      /* Build non oriented graph */
      /* build non symmetric csc from symmetric csc */
      /*maillage global*/
      pastix_int_t *tmpia;
      pastix_int_t *tmpja;
      pastix_int_t  tmpn;

      cscd_symgraph_int(pastix_data->n2,   pastix_data->col2,  pastix_data->row2 , NULL,
                        &tmpn, &tmpia, &tmpja, NULL,
                        pastix_data->loc2glob2, pastix_comm, API_YES);

      memFree_null(pastix_data->col2);
      pastix_data->col2 = tmpia;
      memFree_null(pastix_data->row2);
      pastix_data->row2 = tmpja;
      pastix_data->n2   = tmpn;
    }


#  ifdef DEBUG_DPASTIX
  cscd_save(pastix_data->n2,
            pastix_data->col2,
            pastix_data->row2,
            NULL, NULL, pastix_data->loc2glob2,
            iparm[IPARM_DOF_NBR], "cscd_after_sym", pastix_comm);
#  endif


  /* Remove diagonal coefficients */
  cscd_noDiag(pastix_data->n2, pastix_data->col2, pastix_data->row2, NULL, pastix_data->loc2glob2);
#  ifdef DEBUG_DPASTIX
  cscd_save(pastix_data->n2, pastix_data->col2, pastix_data->row2, NULL, NULL,
            pastix_data->loc2glob2, iparm[IPARM_DOF_NBR], "cscd_after_nodiag", pastix_comm);
#  endif
  if (pastix_data->procnbr > 1)
    {
      /* Check if matrix is not allready correct for scotch */
      /* PT-Scotch needs consecutives column blocs */
      OK = 0;
      OKRecv = 0;
      for (i = 0; i < n-1; i++)
        if (loc2glob[i] != loc2glob[i+1] -1)
          OK = 1;

      MPI_Allreduce(&OK, &OKRecv, 1, MPI_INT, MPI_SUM, pastix_comm);
      /* If it is not correct, permut it */
      if (OKRecv != 0)
        {

          /* Correct the cscd for scotch */
          MALLOC_INTERN(alln, pastix_data->procnbr, int);
          nlocal = n;
          MPI_Allgather(&nlocal, 1, MPI_INT,
                        alln, 1, MPI_INT,
                        pastix_comm);

          MALLOC_INTERN(displs, pastix_data->procnbr, int);
          displs[0] = 0;
          for (i = 1; i < pastix_data->procnbr; i++)
            displs[i] = displs[i-1] + alln[i-1];
          for (i = 0; i < n; i++)
            pastix_data->loc2glob2[i] = displs[pastix_data->procnum]+1+i;

          MALLOC_INTERN(pastix_data->PTS_peritab, gN, pastix_int_t);
          MPI_Allgatherv(loc2glob, n, PASTIX_MPI_INT,
                         pastix_data->PTS_peritab, alln, displs, PASTIX_MPI_INT,
                         pastix_comm);

          memFree_null(displs);
          memFree_null(alln);
          MALLOC_INTERN(pastix_data->PTS_permtab, gN, pastix_int_t);
          for (i = 0; i < gN; i++)
            pastix_data->PTS_permtab[pastix_data->PTS_peritab[i]-1] = i+1;

          /* Redistribue la cscd existante */
          for (i = 0; i < (pastix_data->col2)[n] - 1; i++)
            pastix_data->row2[i] = pastix_data->PTS_permtab[(pastix_data->row2)[i]-1];
        }
    }
  if (iparm[IPARM_VERBOSE] > API_VERBOSE_CHATTERBOX)
    fprintf(stdout, "< dpastix_order_prepare_cscd\n");

  return PASTIX_SUCCESS;
}

/*
  Function: dpastix_task_scotch

  Execute ordering task, with a distributed graph.

  In LOAD mode, only load the graph from disk.

  Else, Clean col2, row2, loc2glob2 and grafmesh if ordering task
  has been called before.


  Build the graph and calls PT-Scotch ordering.

  Gather the graph to be abble to perform centralised version of symbolic
  factorisation.

  Save the graph if asked.

  Parameters:
  pastix_data - PaStiX data structure.
  pastix_comm - PaStiX MPI communicator.
  n           - size of the matrix/number of vertices.
  colptr      - starting index of each column in row
  row         - row of each element
  perm        - permutation tabular
  invp        - reverse permutation tabular
  loc2glob    - global index of local columns
*/
int dpastix_task_scotch(pastix_data_t ** pastix_data,
                        MPI_Comm         pastix_comm,
                        pastix_int_t              n,
                        pastix_int_t            * colptr,
                        pastix_int_t            * row,
                        pastix_int_t            * perm,
                        pastix_int_t            * invp,
                        pastix_int_t            * loc2glob)
{
#  ifndef WITH_SCOTCH
  errorPrint("Distributed PaStiX calls only works with -DWITH_SCOTCH");
  return BADPARAMETER_ERR;
#  else
  pastix_int_t              * iparm       = (*pastix_data)->iparm;
  pastix_int_t              * n2;
  pastix_int_t             ** col2;
  pastix_int_t             ** row2;
  pastix_int_t             ** loc2glob2;
  SCOTCH_Graph     * grafmesh;
  SCOTCH_Dordering * ordedat;
  SCOTCH_Ordering    ordering;
  SCOTCH_Dgraph    * dgraph;
  SCOTCH_Strat       stratdat;
  Order            * ordemesh;
  pastix_int_t                gN;
  double             timer1;
  char               strat[550];
  pastix_int_t                i;
  pastix_int_t                procnum;
  int                retval     = PASTIX_SUCCESS;
  int                retval_rcv;

  print_debug(DBG_STEP,"-> pastix_task_scotch\n");

  if (sizeof(pastix_int_t) != sizeof(SCOTCH_Num))
    {
              errorPrint("Inconsistent integer type %lu != %lu\n",
                         sizeof(pastix_int_t), sizeof(SCOTCH_Num));
      return INTEGER_TYPE_ERR;
    }
  ordemesh  = &((*pastix_data)->ordemesh);
  grafmesh  = &(ordemesh->grafmesh);
  ordedat   = &(ordemesh->ordedat);
  dgraph    = &(ordemesh->dgraph);
  procnum   =   (*pastix_data)->procnum;
  n2        = &((*pastix_data)->n2);
  col2      = &((*pastix_data)->col2);
  row2      = &((*pastix_data)->row2);
  loc2glob2 = &((*pastix_data)->loc2glob2);

  if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    print_onempi("%s", OUT_STEP_ORDER);

  /* On nettoie grafmesh et col2/row2 si ils sont déja alloués */
#    ifdef WITH_SCOTCH
  if (ordemesh->malgrf == 1)
    {
      SCOTCH_graphExit(&(ordemesh->grafmesh));
      ordemesh->malgrf = 0;
    }
#    endif
  if ((*pastix_data)->bmalcolrow == 1)
    {
      if ((*col2)      != NULL) memFree_null(*col2);
      if ((*row2)      != NULL) memFree_null(*row2);
      if ((*loc2glob2) != NULL) memFree_null(*loc2glob2);
      (*pastix_data)->bmalcolrow = 0;
    }
  if (iparm[IPARM_ORDERING]  != API_ORDER_LOAD)
    {
      dpastix_order_prepare_cscd(*pastix_data,
                                 n,
                                 colptr,
                                 row,
                                 loc2glob,
                                 pastix_comm);
      /* Allocate a copy of col, row loc2glob */
      (*pastix_data)->bmalcolrow = 1;

    }

  switch (iparm[IPARM_ORDERING])
    {
      /*
       * Scotch Ordering
       */
    case API_ORDER_SCOTCH:
      /* construction du graphe */
      print_debug(DBG_SCOTCH, "> SCOTCH_dgraphInit <\n");
      SCOTCH_dgraphInit(dgraph, pastix_comm);

      print_debug(DBG_SCOTCH, "> SCOTCH_dgraphBuild <\n");
      if ( SCOTCH_dgraphBuild (dgraph,
                               1,              /* baseval */
                               *n2,            /* number of local vertices */
                               *n2,            /* Maximum number of local vertices     */
                               (*col2),
                               NULL,
                               NULL,           /* Local vertex load array (if any)     */
                               NULL,           /* Local vertex label array (if any)    */
                               (*col2)[(*n2)] - 1,
                               (*col2)[(*n2)] - 1,
                               (*row2),        /* Local edge array                     */
                               NULL,           /* Ghost edge array (if any); not const */
                               NULL))
        {
          errorPrint("SCOTCH_dgraphBuild");
          EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }
      ordemesh->malgrf = 1;

      print_debug(DBG_SCOTCH, "> SCOTCH_dgraphCheck <\n");
      if (SCOTCH_dgraphCheck(dgraph))
        {
          errorPrint("pastix : SCOTCH_dgraphCheck");
          EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }

      print_debug(DBG_SCOTCH, "> SCOTCH_stratInit <\n");
      if (SCOTCH_stratInit(&stratdat))
        {
          errorPrint("pastix : SCOTCH_stratInit");
          EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }

      /* TODO : mettre des strategies par défaut */
      if (iparm[IPARM_DEFAULT_ORDERING] == API_YES)
        {
          if (iparm[IPARM_INCOMPLETE] == API_NO)
            sprintf(strat, PTSCOTCH_STRAT_DIRECT);
          else
            sprintf(strat, PTSCOTCH_STRAT_INCOMP);
        }
      else /* personal ordering */
        {
          sprintf(strat, PTSCOTCH_STRAT_PERSO,
                  (long) iparm[IPARM_ORDERING_SWITCH_LEVEL],
                  (long) iparm[IPARM_ORDERING_CMIN],
                  (long) iparm[IPARM_ORDERING_CMAX],
                  ((float)iparm[IPARM_ORDERING_FRAT])/100.,
                  (long) iparm[IPARM_ORDERING_SWITCH_LEVEL],
                  (long) iparm[IPARM_ORDERING_CMIN],
                  (long) iparm[IPARM_ORDERING_CMAX],
                  ((float)iparm[IPARM_ORDERING_FRAT])/100.);

          if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
            print_onempi("Scotch Strategy |%s|\n", strat);
        }

      clockInit(timer1);
      clockStart(timer1);

      /*    print_debug(DBG_SCOTCH, "> SCOTCH_stratDgraphOrder <\n"); */
      /*    if (SCOTCH_stratDgraphOrder(&stratdat, strat)) */
      /*      { */
      /*        errorPrint("pastix : SCOTCH_stratDgraphOrder"); */
      /*        EXIT(MOD_SOPALIN,INTERNAL_ERR); */
      /*      } */
      if (procnum == 0 && iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        errorPrintW("PaStiX works only with PT-Scotch default strategy");

      print_debug(DBG_SCOTCH, "> SCOTCH_dgraphOrderInit <\n");
      if (0 != SCOTCH_dgraphOrderInit(dgraph, ordedat))
        {
          errorPrint("pastix : SCOTCH_dgraphOrderInit");
          EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }

      print_debug(DBG_SCOTCH, "> SCOTCH_dgraphOrderCompute <\n");
      if (0 != SCOTCH_dgraphOrderCompute(dgraph, ordedat, &stratdat))
        {
          errorPrint("pastix : SCOTCH_dgraphOrderCompute");
          EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }

      print_debug(DBG_SCOTCH, "> SCOTCH_stratExit <\n");
      SCOTCH_stratExit(&stratdat);


      /* print_debug(DBG_SCOTCH, "> SCOTCH_dgraphOrderPerm <\n"); */
      /*       if (0 != SCOTCH_dgraphOrderPerm(dgraph, ordedat, perm)) */
      /*  { */
      /*    errorPrint("pastix : SCOTCH_dgraphOrderPerm"); */
      /*    EXIT(MOD_SOPALIN,INTERNAL_ERR); */
      /*  } */

      clockStop(timer1);
      if (iparm[IPARM_VERBOSE] > API_VERBOSE_NOT)
        print_onempi(TIME_COMPUTE_ORDERING,clockVal(timer1));


      /* Clean ordering if it exists */
      if ((*pastix_data)->malord)
        {
          orderExit(ordemesh);
          (*pastix_data)->malord=0;
        }

      if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
        print_onempi("%s", OUT_ORDERINIT);
      gN = 0;
      MPI_Allreduce(&n, &gN, 1, PASTIX_MPI_INT, MPI_SUM, pastix_comm);

      orderInit(ordemesh, gN, gN);
      (*pastix_data)->malord=1;

      for (i=0;i<gN+1;i++)
        ordemesh->rangtab[i] = 0;

      SCOTCH_dgraphCorderInit (dgraph,
                               &ordering,
                               (SCOTCH_Num *)ordemesh->permtab,
                               (SCOTCH_Num *)ordemesh->peritab,
                               &ordemesh->cblknbr,
                               ordemesh->rangtab,
                               NULL);

      if (procnum == 0)
        {
          SCOTCH_dgraphOrderGather (dgraph, ordedat, &ordering);
        }
      else
        {
          SCOTCH_dgraphOrderGather (dgraph, ordedat, NULL);
        }
      MPI_Bcast(&ordemesh->cblknbr, 1                   , PASTIX_MPI_INT, 0, pastix_comm);
      MPI_Bcast(ordemesh->rangtab, (ordemesh->cblknbr+1), PASTIX_MPI_INT, 0, pastix_comm);
      MPI_Bcast(ordemesh->permtab, gN                   , PASTIX_MPI_INT, 0, pastix_comm);
      MPI_Bcast(ordemesh->peritab, gN                   , PASTIX_MPI_INT, 0, pastix_comm);

      global2localperm(n, perm, ((*pastix_data)->ordemesh).permtab, loc2glob);
      /* Gathering graph */
      print_debug(DBG_SCOTCH, "> SCOTCH_dgraphGather <\n");
      SCOTCH_dgraphGather (dgraph, &(ordemesh->grafmesh));
      SCOTCH_dgraphCorderExit(dgraph, &ordering);
      SCOTCH_dgraphOrderExit(dgraph, ordedat);
      SCOTCH_dgraphExit(dgraph);

      SCOTCH_graphBase(&(ordemesh->grafmesh), 0);
      orderBase(ordemesh, 0);

      if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE))
        {
          pastix_int_t nsave;
          pastix_int_t * colptrsave = NULL;
          pastix_int_t * rowsave    = NULL;
          if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE_CSC))
            {
              cscd2csc_int((*pastix_data)->n2, (*pastix_data)->col2, (*pastix_data)->row2, NULL,
                           NULL, NULL, NULL,
                           &nsave, &colptrsave, &rowsave, NULL,
                           NULL, NULL, NULL,
                           (*pastix_data)->loc2glob2, pastix_comm, iparm[IPARM_DOF_NBR], API_YES);
            }

          retval = orderSaveFiles( *pastix_data );
          if (PASTIX_MASK_ISTRUE(iparm[IPARM_IO_STRATEGY], API_IO_SAVE_CSC))
          {
              memFree_null(colptrsave);
              memFree_null(rowsave);
          }
          if (retval != PASTIX_SUCCESS)
              break;
        }

#    ifdef FORGET_PARTITION
      {
        pastix_int_t iter;
        memFree_null(ordemesh->rangtab);
        MALLOC_INTERN(ordemesh->rangtab, n+1, pastix_int_t);
        for (iter=0;iter<n+1;iter++)
          ordemesh->rangtab[iter]=iter;
        ordemesh->cblknbr=n;
      }
#    endif

      break;
    case API_ORDER_PERSONAL:
      {
        pastix_int_t iter;
        pastix_int_t * tmpperm = NULL;
        gN = 0;
        MPI_Allreduce(&n, &gN, 1, PASTIX_MPI_INT, MPI_SUM, pastix_comm);
        MALLOC_INTERN(tmpperm, gN, pastix_int_t);
        /* Clean ordering if it exists */
        if ((*pastix_data)->malord)
          {
            orderExit(ordemesh);
            (*pastix_data)->malord=0;
          }

        if (iparm[IPARM_VERBOSE] > API_VERBOSE_YES)
          print_onempi("%s", OUT_ORDERINIT);

        gN = 0;
        MPI_Allreduce(&n, &gN, 1, PASTIX_MPI_INT, MPI_SUM, pastix_comm);

        orderInit(ordemesh, gN, gN);
        (*pastix_data)->malord=1;

        memset(ordemesh->permtab, 0, gN   *sizeof(pastix_int_t));

        memset(tmpperm, 0, gN*sizeof(pastix_int_t));
        for (iter = 0; iter < n; iter++)
          tmpperm[loc2glob[iter]-1] = perm[iter]-1;

        MPI_Allreduce(tmpperm, ordemesh->permtab, gN, PASTIX_MPI_INT, MPI_SUM, pastix_comm);
        for (iter = 0; iter < gN; iter++)
          ordemesh->peritab[ordemesh->permtab[iter]] = iter;

        break;
      }
    case API_ORDER_LOAD:
    {
        orderLoadFiles( *pastix_data );
        break;
    }

    default:
      errorPrint("Ordering not available with distributed interface");
      retval = BADPARAMETER_ERR;
      break;
    }

  MPI_Allreduce(&retval, &retval_rcv, 1, MPI_INT, MPI_MAX, pastix_comm);
  if (retval_rcv != PASTIX_SUCCESS)
    return retval_rcv;

  (*pastix_data)->malord=1;

  iparm[IPARM_START_TASK]++;
#  endif /* WITH_SCOTCH */
  return PASTIX_SUCCESS;
}
#endif /* PASTIX_DISTRIBUTED */
