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
    Order            * ordemesh;
    double             timer1;
    pastix_int_t                procnum;
    int                retval     = PASTIX_SUCCESS;
    int                retval_rcv;
    (void)pastix_comm;

    procnum  = (*pastix_data)->procnum;

    print_debug(DBG_STEP,"-> pastix_task_scotch\n");
    if (iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
        print_onempi("%s", OUT_STEP_ORDER);

    /* Clean ordering if it exists */
    if ((*pastix_data)->ordemesh != NULL) {
        orderExit((*pastix_data)->ordemesh);
    } else {
        MALLOC_INTERN( (*pastix_data)->ordemesh, 1, Order );
    }
    ordemesh = (*pastix_data)->ordemesh;
    orderInit( ordemesh, 0, 0 );

    /* Prepare a copy of user's CSC */
    if (!(PASTIX_MASK_ISTRUE(iparm[IPARM_ORDERING], API_ORDER_LOAD)))
        orderPrepareCSC( *pastix_data, n, colptr, row, NULL );

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
        {
            orderComputeScotch( *pastix_data, (*pastix_data)->csc );
        }
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

        assert( 0 );
        orderLoadFiles( *pastix_data );
    }
    break;

    /*
     * Load ordering with Scotch Format
     */
    case API_ORDER_LOAD:
        assert( 0 );
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
