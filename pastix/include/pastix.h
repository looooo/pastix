/**
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.2.0
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
#include "pastix/datatypes.h"
#include "pastix/api.h"
#include <math.h>
#if defined(HAVE_MPI)
#include <mpi.h>
#else
#include "pastix/nompi.h"
#endif

/** ****************************************************************************
 *
 *  PaStiX constants - Compatible with CBLAS & LAPACK
 *  The naming and numbering is consistent with:
 *
 *    1) CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz),
 *    2) C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/).
 *
 **/
#define PastixNoTrans       111
#define PastixTrans         112
#define PastixConjTrans     113

#define PastixGeneral       111
#define PastixSymmetric     112
#define PastixHermitian     113

/** ****************************************************************************
 * Sparse matrix format
 **/
typedef enum pastix_fmttype_e {
    PastixCSC = 0,
    PastixCSR = 1,
    PastixIJV = 2
} pastix_fmttype_t;


#define PASTIX_SUCESS  0

struct pastix_data_s;
typedef struct pastix_data_s pastix_data_t;

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
void pastix(pastix_data_t **pastix_data, MPI_Comm pastix_comm,
            pastix_int_t n, pastix_int_t *colptr, pastix_int_t *row,
            pastix_complex64_t *avals, pastix_int_t *perm, pastix_int_t *invp, pastix_complex64_t *b, pastix_int_t rhs,
            pastix_int_t *iparm, double *dparm);

/*
 * Function: dpastix
 *
 *   Computes steps of the resolution of Ax=b linear system,
 *   using direct methods.
 *   Here the matrix is given distributed.
 *
 *   The matrix is given in CSCD format.
 *
 *   Parameters:
 *      pastix_data - Data used for a step by step execution.
 *      pastix_comm - MPI communicator which compute the resolution.
 *      n           - Size of the system.
 *      colptr      - Tabular containing the start of each column in
 *                    *row* and *avals* tabulars.
 *      row         - Tabular containing the row number for each element
 *                    sorted by column.
 *      avals       - Tabular containing the values of each element
 *                    sorted by column.
 *      loc2glob    - Global column number of the local columns.
 *      perm        - Permutation tabular for the renumerotation
 *                    of the unknowns.
 *      invp        - Reverse permutation tabular for the renumerotation
 *                    of the unknowns.
 *      b           - Right hand side vector(s).
 *      rhs         - Number of right hand side vector(s).
 *      iparm       - Integer parameters given to pastix.
 *      dparm       - Double parameters given to pâstix.
 *
 * About: Example
 *
 *   from file <simple_dist.c> :
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
 *   > pastix_checkMatrix( MPI_COMM_WORLD, verbosemode,
 *   >                     (MTX_ISSYM(type) ? API_SYM_YES : API_SYM_NO),
 *   >                     API_YES,
 *   >                     ncol, &colptr, &rows, &values, &loc2glob, 1);
 *   >
 *   > /\*******************************************\/
 *   > /\* Initialize parameters to default values *\/
 *   > /\*******************************************\/
 *   > iparm[IPARM_MODIFY_PARAMETER] = API_NO;
 *   > dpastix(&pastix_data, MPI_COMM_WORLD,
 *   >         ncol, colptr, rows, values, loc2glob,
 *   >         perm, invp, rhs, 1, iparm, dparm);
 *   >
 *   > /\*******************************************\/
 *   > /\*       Customize some parameters         *\/
 *   > /\*******************************************\/
 *   > iparm[IPARM_THREAD_NBR] = nbthread;
 *   > if (MTX_ISSYM(type))
 *   >   {
 *   >     iparm[IPARM_SYM]           = API_SYM_YES;
 *   >     iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
 *   >   }
 *   > else
 *   >   {
 *   >     iparm[IPARM_SYM]           = API_SYM_NO;
 *   >     iparm[IPARM_FACTORIZATION] = API_FACT_LU;
 *   >   }
 *   > iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
 *   >
 *   > iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
 *   > iparm[IPARM_END_TASK]            = API_TASK_BLEND;
 *   >
 *   > /\*******************************************\/
 *   > /\*           Call pastix                   *\/
 *   > /\*******************************************\/
 *   > perm = malloc(ncol*sizeof(pastix_int_t));
 *   > /\* No need to allocate invp in dpastix *\/
 *   >
 *   > dpastix(&pastix_data, MPI_COMM_WORLD,
 *   >         ncol, colptr, rows, NULL, loc2glob,
 *   >         perm, NULL, NULL, 1, iparm, dparm);
 *   >
 *   > /\* Redistributing the matrix *\/
 *   >
 *   > ncol2 = pastix_getLocalNodeNbr(&pastix_data);
 *   >
 *   > if (NULL == (loc2glob2 = malloc(ncol2 * sizeof(pastix_int_t))))
 *   >   {
 *   >     fprintf(stderr, "Malloc error\n");
 *   >     return EXIT_FAILURE;
 *   >   }
 *   >
 *   > pastix_getLocalNodeLst(&pastix_data, loc2glob2);
 *   >
 *   > if (EXIT_SUCCESS != cscd_redispatch(ncol,   colptr,   rows,
 *   >                                     values,    rhs,  loc2glob,
 *   >                                     ncol2, &colptr2, &rows2,
 *   >                                     &values2, &rhs2, loc2glob2,
 *   >                                     MPI_COMM_WORLD))
 *   >   return EXIT_FAILURE;
 *   >
 *   > free(colptr);
 *   > free(rows);
 *   > free(values);
 *   > free(rhs);
 *   > free(loc2glob);
 *   > free(perm);
 *   >
 *   > iparm[IPARM_START_TASK]          = API_TASK_NUMFACT;
 *   > iparm[IPARM_END_TASK]            = API_TASK_CLEAN;
 *   >
 *   >
 *   > dpastix(&pastix_data, MPI_COMM_WORLD,
 *   >         ncol2, colptr2, rows2, values2, loc2glob2,
 *   >         perm, invp, rhs2, 1, iparm, dparm);
 *   >
 */
void dpastix(pastix_data_t **pastix_data, MPI_Comm pastix_comm,
             pastix_int_t n, pastix_int_t *colptr, pastix_int_t *row,
             pastix_complex64_t *avals, pastix_int_t * loc2glob, pastix_int_t *perm, pastix_int_t *invp,
             pastix_complex64_t *b, pastix_int_t rhs, pastix_int_t *iparm, double *dparm);

/*
 * Group: Thread functions
 */

/*
  Function: pastix_bindThreads

  Set bindtab in pastix_data, it gives for each thread the CPU to bind in to.
  bindtab follows this organisation :

  bindtab[threadnum] = cpu to set thread threadnum.

  Parameters:
    pastix_data - Structure de donnée pour l'utilisation step by step
    thrdnbr     - Nombre de threads / Taille du tableau
    bindtab     - Tableau de correspondance entre chaque thread et coeur de la machine
*/
void pastix_bindThreads(pastix_data_t *pastix_data, pastix_int_t thrdnbr, pastix_int_t *bindtab);

/*
 * Group: Checking the matrix.
 */

/*
 * Function: pastix_checkMatrix
 *
 * Check the matrix :
 * - Renumbers in Fortran numerotation (base 1) if needed (base 0)
 * - Check that the matrix contains no doubles,  with flagcor == API_YES,
 *   correct it.
 * - Can scale the matrix if compiled with -DMC64 -DSCALING (untested)
 * - Checks the symetry of the graph in non symmetric mode.
 *   With non distributed matrices, with flagcor == API_YES,
 *   correct the matrix.
 * - sort the CSC.
 *
 * Parameters:
 *   pastix_comm - PaStiX MPI communicator
 *   verb        - Level of prints (API_VERBOSE_[NOT|NO|YES])
 *   flagsym     - Indicate if the given matrix is symetric
 *                 (API_SYM_YES or API_SYM_NO)
 *   flagcor     - Indicate if we permit the function to reallocate the matrix.
 *   n           - Number of local columns.
 *   colptr      - First element of each row in *row* and *avals*.
 *   row         - Row of each element of the matrix.
 *   avals       - Value of each element of the matrix.
 *   loc2glob    - Global column number of local columns
 *                 (NULL if not distributed).
 *   dof         - Number of degrees of freedom.
 */
pastix_int_t pastix_checkMatrix(MPI_Comm pastix_comm, pastix_int_t verb, pastix_int_t flagsym, pastix_int_t flagcor,
                                pastix_int_t n, pastix_int_t **colptr, pastix_int_t **row, pastix_complex64_t **avals,
                                pastix_int_t **loc2glob, pastix_int_t dof);

/*
 * Group: Getting solver distribution.
 */

/*
  Function: pastix_getLocalNodeNbr

  Return the node number in the new distribution computed by analyze step.
  Needs analyze step to be runned with pastix_data before.

  Parameters:
    pastix_data - Data used for a step by step execution.

  Returns:
    Number of local nodes/columns in new distribution.
 */
pastix_int_t pastix_getLocalNodeNbr(pastix_data_t ** pastix_data);

/*
  Function: pastix_getLocalNodeLst

  Fill in nodelst with the list of local nodes/columns.
  Needs nodelst to be allocated with nodenbr*sizeof(pastix_int_t),
  where nodenbr has been computed by <pastix_getLocalNodeNbr>.

  Parameters:
    pastix_data - Data used for a step by step execution.
    nodelst     - An array where to write the list of local nodes/columns.
 */
pastix_int_t pastix_getLocalNodeLst(pastix_data_t ** pastix_data, pastix_int_t * nodelst);

/*
  Function: pastix_getLocalUnknownNbr

  Return the unknown number in the new distribution computed by analyze step.
  Needs analyze step to be runned with pastix_data before.

  Parameters:
    pastix_data - Data used for a step by step execution.

  Returns:
    Number of local unknowns/columns in new distribution.
 */
pastix_int_t pastix_getLocalUnknownNbr(pastix_data_t ** pastix_data);

/*
  Function: pastix_getLocalUnknownLst

  Fill in unknownlst with the list of local unknowns/clumns.
  Needs unknownlst to be allocated with unknownnbr*sizeof(pastix_int_t),
  where unknownnbr has been computed by <pastix_getLocalUnknownNbr>.

  Parameters:
    pastix_data - Data used for a step by step execution.
    unknownlst     - An array where to write the list of local unknowns/columns.
 */
pastix_int_t pastix_getLocalUnknownLst(pastix_data_t ** pastix_data, pastix_int_t * unknownlst);

/*
 * Group: About the Schur complement.
 */


/*
  Function: pastix_setSchurUnknownList

  Set the list of unknowns to isolate at the end
  of the matrix via permutations.

  Has to be called if using IPARM_SCHUR = API_YES.

  Parameters:
    pastix_data - Data used for a step by step execution.
    n           - Number of unknowns.
    list        - List of unknowns.
*/
pastix_int_t pastix_setSchurUnknownList(pastix_data_t * pastix_data,
             pastix_int_t  n,
             pastix_int_t *list);

/*
  Function: pastix_getSchurLocalNodeNbr

  Compute the number of nodes in the local part of the Schur.

  Parameters:
    pastix_data - Common data structure for PaStiX calls.
    nodeNbr     - (out) Number of nodes in Schur (local).

  Returns:
    NO_ERR      - For the moment

  TODO: Error management.
*/
pastix_int_t pastix_getSchurLocalNodeNbr(pastix_data_t * pastix_data, pastix_int_t * nodeNbr);

/*
  Function: pastix_getSchurLocalUnkownNbr

  Compute the number of unknowns in the local part of the Schur.

  Parameters:
    pastix_data - Common data structure for PaStiX calls.
    unknownNbr  - (out) Number of unknowns in Schur (local).

  Returns:
    NO_ERR      - For the moment

  TODO: Error management.
*/
pastix_int_t pastix_getSchurLocalUnkownNbr(pastix_data_t * pastix_data,
                                  pastix_int_t * unknownNbr);

/*
  Function: pastix_getSchurLocalNodeList

  Compute the list of nodes in the local part of the Schur.

  Parameters:
    pastix_data - Common data structure for PaStiX calls.
    nodes     - (out) Nodes in Schur (local).

  Returns:
    NO_ERR      - For the moment

  TODO: Error management.
*/
pastix_int_t pastix_getSchurLocalNodeList(pastix_data_t * pastix_data, pastix_int_t * nodes);

/*
  Function: pastix_getSchurLocalUnkownList

  Compute the list of unknowns in the local part of the Schur.

  Parameters:
    pastix_data - Common data structure for PaStiX calls.
    unknowns    - (out) Unknowns in Schur (local).

  Returns:
    NO_ERR      - For the moment

  TODO: Error management.
*/
pastix_int_t pastix_getSchurLocalUnknownList(pastix_data_t * pastix_data,
                                    pastix_int_t * unknowns);

/*
  Function: pastix_setSchurArray

  Give user memory area to store Schur in PaStiX.

  Parameters:
    pastix_data - Common data structure for PaStiX calls.
    array       - Memory area to store the Schur.

  Returns:
    NO_ERR      - For the moment

  TODO: Error management.
*/
pastix_int_t pastix_setSchurArray(pastix_data_t * pastix_data, pastix_complex64_t * array);

/*
  Function: pastix_getSchur

  Get the Schur complement from PaStiX.

  Schur complement is a dense block in a
  column scheme.

  Parameters:
    pastix_data - Data used for a step by step execution.
    schur - Array to fill-in with Schur complement.

*/
pastix_int_t pastix_getSchur(pastix_data_t * pastix_data,
                             pastix_complex64_t * schur);

/*
 * Group: About parameters.
 */

/*
  Function: pastix_initParam

  Sets default parameters for iparm and dparm

  Parameters:
    iparm - tabular of IPARM_SIZE integer parameters.
    dparm - tabular of DPARM_SIZE double parameters.
*/
void pastix_initParam(pastix_int_t    *iparm,
                      double *dparm);

/*
 * Group: About Scaling.
 *
 * Working on it...
 */

/*
 * Function: pastix_unscale
 *
 * Unscale the factorized SolverMatrix.
 *
 * Unscale the factorized SolverMatrix in pastix_data,
 * according to pastix_data->scalerowtab,
 * pastix_data->iscalerowtab, pastix_data->scalecoltab and
 * pastix_data->iscalecoltab
 *
 * (Elements in pastix_data->iscaleXXXtab are supposed to be the inverse of
 * elements in pastix_data->scaleXXXtab).
 *
 * Parameters:
 *   pastix_data - Common data structure for PaStiX calls.
 *   sym         - API_YES if the matrix is symmetric
 *                   ( only pastix_data->scalerowtab and
`*                     pastix_data->iscalerowtab
 *                     are used in this case),
 *                 API_NO otherwise
 */
void pastix_unscale ( pastix_data_t *pastix_data, pastix_int_t sym);

unsigned long pastix_getMemoryUsage(void);
unsigned long pastix_getMaxMemoryUsage(void);

struct pastix_graph_s;
typedef struct pastix_graph_s pastix_graph_t;

struct Order_;
typedef struct Order_ Order;

struct pastix_spm_s;
typedef struct pastix_spm_s pastix_csc_t;

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

int pastix_task_order(pastix_data_t *pastix_data,
                      pastix_int_t   n,
                      const pastix_int_t *colptr,
                      const pastix_int_t *row,
                      const pastix_int_t *loc2glob,
                      pastix_int_t  *perm,
                      pastix_int_t  *invp);

int pastix_task_symbfact(pastix_data_t *pastix_data,
                         pastix_int_t  *perm,
                         pastix_int_t  *invp);

int pastix_task_blend( pastix_data_t *pastix_data );

int pastix_task_sopalin( pastix_data_t *pastix_data,
                          pastix_csc_t  *csc );

#endif /* _PASTIX_H_ */
