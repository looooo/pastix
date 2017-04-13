/**
 *
 * @copyright 2004-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author David Goudin
 * @author Pascal Henon
 * @author Francois Pellegrini
 *
 **/
#ifndef _PASTIX_H_
#define _PASTIX_H_

#include "pastix/config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include "pastix/api.h"
#include "pastix/datatypes.h"

#include <math.h>
#if defined(HAVE_MPI)
#include <mpi.h>

#else
#include "pastix/nompi.h"
#endif
#include "pastix/old_api.h"

#define PASTIX_SUCESS  0

/**
 * Main structure of the pastix solver associated to a given problem
 */
struct pastix_data_s;
typedef struct pastix_data_s pastix_data_t;

struct pastix_graph_s;
typedef struct pastix_graph_s pastix_graph_t;

struct pastix_order_s;
typedef struct pastix_order_s Order;

struct pastix_spm_s;
typedef struct pastix_spm_s pastix_spm_t;
typedef struct pastix_spm_s pastix_csc_t;

struct SolverMatrix_;
typedef struct SolverMatrix_ SolverMatrix;

enum pastix_driver_e;
typedef enum pastix_driver_e pastix_driver_t;

/**
 * Some define for old pastix compatibility
 */
#define pastix_float_t void
#define API_SYM_YES PastixSymmetric
#define API_SYM_HER PastixHermitian
#define API_SYM_NO PastixGeneral

/*
 * Group: Main PaStiX functions
 */
/*
 * Function: pastix
 *
 * Computes steps of the resolution of Ax=b linear system,
 * using direct methods.
 *
 * The matrix is given in CSC format.
 *
 * Parameters:
 *   pastix_data - Data used for a step by step execution.
 *   pastix_comm - MPI communicator which compute the resolution.
 *   n           - Size of the system.
 *   colptr      - Tabular containing the start of each column in row
 *                 and avals tabulars.
 *   row         - Tabular containing the row number for each element
 *                 sorted by column.
 *   avals       - Tabular containing the values of each element
 *                 sorted by column.
 *   perm        - Permutation tabular for the renumerotation of the unknowns.
 *   invp        - Reverse permutation tabular for the renumerotation
 *                 of the unknowns.
 *   b           - Right hand side vector(s).
 *   rhs         - Number of right hand side vector(s).
 *   iparm       - Integer parameters given to pastix.
 *   dparm       - Double parameters given to pâstix.
 *
 * About: Example
 *
 *   from file <simple.c> :
 *
 *   > /\*******************************************\/
 *   > /\*    Check Matrix format                  *\/
 *   > /\*******************************************\/
 *   > /\*
 *   >  * Matrix needs :
 *   >  *    - to be in fortran numbering
 *   >  *    - to have only the lower triangular part in symmetric case
 *   >  *    - to have a graph with a symmetric structure in unsymmetric case
 *   >  *\/
 *   > mat_type = API_SYM_NO;
 *   > if (MTX_ISSYM(type)) mat_type = API_SYM_YES;
 *   > if (MTX_ISHER(type)) mat_type = API_SYM_HER;
 *   > pastix_checkMatrix( MPI_COMM_WORLD, verbosemode,
 *   >                     mat_sym,
 *   >                     API_YES,
 *   >                     ncol, &colptr, &rows, &values, NULL);
 *   >
 *   > /\*******************************************\/
 *   > /\* Initialize parameters to default values *\/
 *   > /\*******************************************\/
 *   > iparm[IPARM_MODIFY_PARAMETER] = API_NO;
 *   > pastix(&pastix_data, MPI_COMM_WORLD,
 *   >        ncol, colptr, rows, values,
 *   >        perm, invp, rhs, 1, iparm, dparm);
 *   >
 *   > /\*******************************************\/
 *   > /\*       Customize some parameters         *\/
 *   > /\*******************************************\/
 *   > iparm[IPARM_THREAD_NBR] = nbthread;
 *   > iparm[IPARM_SYM] = mat_type;
 *   > switch (mat_type)
 *   >   {
 *   >     case API_SYM_YES:
 *   >       iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
 *   >       break;
 *   >     case API_SYM_HER:
 *   >       iparm[IPARM_FACTORIZATION] = API_FACT_LDLH;
 *   >       break;
 *   >     default:
 *   >       iparm[IPARM_FACTORIZATION] = API_FACT_LU;
 *   >   }
 *   > iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
 *   > iparm[IPARM_END_TASK]            = API_TASK_CLEAN;
 *   >
 *   > /\*******************************************\/
 *   > /\*           Save the rhs                  *\/
 *   > /\*    (it will be replaced by solution)    *\/
 *   > /\*******************************************\/
 *   > rhssaved = malloc(ncol*sizeof(pastix_complex64_t));
 *   > memcpy(rhssaved, rhs, ncol*sizeof(pastix_complex64_t));
 *   >
 *   > /\*******************************************\/
 *   > /\*           Call pastix                   *\/
 *   > /\*******************************************\/
 *   > perm = malloc(ncol*sizeof(pastix_int_t));
 *   > invp = malloc(ncol*sizeof(pastix_int_t));
 *   >
 *   > pastix(&pastix_data, MPI_COMM_WORLD,
 *   >  ncol, colptr, rows, values,
 *   >  perm, invp, rhs, 1, iparm, dparm);
 */
void pastix(pastix_data_t **pastix_data,  MPI_Comm pastix_comm, pastix_int_t n,
            pastix_int_t *colptr, pastix_int_t *row, void *avals, pastix_int_t *perm,
            pastix_int_t *invp, void *b, pastix_int_t nrhs, pastix_int_t *iparm,
            double *dparm);

void pastixInitParam( pastix_int_t *iparm,
                      double       *dparm );
void pastixInit( pastix_data_t **pastix_data,
                 MPI_Comm        pastix_comm,
                 pastix_int_t   *iparm,
                 double         *dparm );

void pastixFinalize( pastix_data_t **pastix_data,
                     MPI_Comm        pastix_comm,
                     pastix_int_t   *iparm,
                     double         *dparm );

/**
 * Main steps of the solver
 */
int pastix_task_analyze( pastix_data_t      *pastix_data,
                         pastix_spm_t       *spm );
int pastix_task_numfact( pastix_data_t      *pastix_data,
                         const pastix_spm_t *spm );
int pastix_task_solve  ( pastix_data_t      *pastix_data,
                         const pastix_spm_t *spm,
                         int                 nrhs,
                         void               *b,
                         int                 ldb );
int pastix_task_refine ( pastix_data_t *pastix_data,
                         void          *x,
                         pastix_int_t   rhsnbr,
                         void          *b);

/**
 * Analyze subtasks
 */
int pastix_subtask_order     ( pastix_data_t      *pastix_data,
                               const pastix_spm_t *spm,
                               pastix_int_t       *perm,
                               pastix_int_t       *invp );
int pastix_subtask_symbfact  ( pastix_data_t      *pastix_data,
                               pastix_int_t       *perm,
                               pastix_int_t       *invp );
int pastix_subtask_reordering( pastix_data_t      *pastix_data );
int pastix_subtask_blend     ( pastix_data_t      *pastix_data );

/**
 * Numerical factorization subtasks
 */
int pastix_subtask_spm2bcsc  ( pastix_data_t      *pastix_data,
                               const pastix_spm_t *spm );
int pastix_subtask_bcsc2ctab ( pastix_data_t      *pastix_data,
                               const pastix_spm_t *spm );
int pastix_subtask_sopalin   ( pastix_data_t      *pastix_data,
                               const pastix_spm_t *spm );

void
pastix_setSchurUnknownList(pastix_data_t *pastix_data,
                           pastix_int_t   n,
                           pastix_int_t  *list);

#endif /* _PASTIX_H_ */
