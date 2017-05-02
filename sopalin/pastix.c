/**
 *
 * @file pastix.c
 *
 * PaStiX main interface for compatibility with former releases
 *
 * @copyright 2011-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathias Hastaran
 * @date 2011-11-11
 *
 **/
#include "common.h"
#include "spm.h"

/**
 *******************************************************************************
 *
 * @brief PaStiX main solver function.
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix data structure of the solver to store the state of the
 *          solver at every call.
 *
 * @param[in] pastix_comm
 *          The MPI communicator to use for the distributed solver.
 *
 * @param[in] n
 *          The size of the sparse matrix system to solve.
 *
 * @param[inout] colptr
 *          The pointer to the column index of the sparse matrix in the CSC
 *          format.
 *          On exit, the base value of the array might have changed, and/or the
 *          pointer might have been freed if IPARM_FREE_CSCUSER is set, and the
 *          factorization step is performed.
 *
 * @param[inout] row
 *          The pointer to the row array of the sparse matrix in the CSC
 *          format.
 *          On exit, the base value of the array might have changed, and/or the
 *          pointer might have been freed if IPARM_FREE_CSCUSER is set, and the
 *          factorization step is performed.
 *
 * @param[inout] avals
 *          The pointer to the values array of the sparse matrix in the CSC
 *          format.
 *          On exit, the pointer might have been freed if IPARM_FREE_CSCUSER is
 *          set, and the factorization step is performed.
 *
 * @param[inout] perm
 *          The pointer to the permutation array.
 *          On entry: the pointer might be allocated to store the generated
 *          permutation on exit, or to provide the user permutation.
 *          On exit, the permutation used by the solver is returned if perm is
 *          not NULL.
 *
 * @param[inout] invp
 *          The pointer to the inverse permutation array.
 *          On entry: the pointer might be allocated to store the generated
 *          inverse permutation on exit, or to provide the user permutation.
 *          On exit, the inverse permutation used by the solver is returned if
 *          invp is not NULL.
 *
 * @param[inout] b
 *          Array of size n -by- nrhs
 *          On entry, contains the nrhs vectors of the problem.
 *          On exit, contains the nrhs solution vectors of the problem.
 *
 * @param[in] nrhs
 *          The number of right hand side in the problem.
 *
 * @param[inout] iparm
 *          Array of size IPARM_SIZE
 *          On entry, contains all the integer parameters of the solver.
 *          On exit, the aray is updated with integer outputs of the solver.
 *
 * @param[inout] dparm
 *          Array of size DPARM_SIZE
 *          On entry, contains all the double parameters of the solver.
 *          On exit, the aray is updated with double outputs of the solver.
 *
 *******************************************************************************
 *
 * Example:
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
 *
 *******************************************************************************
 */
void
pastix( pastix_data_t **pastix_data_ptr,
        MPI_Comm        pastix_comm,
        pastix_int_t    n,
        pastix_int_t   *colptr,
        pastix_int_t   *row,
        void           *avals,
        pastix_int_t   *perm,
        pastix_int_t   *invp,
        void           *b,
        pastix_int_t    nrhs,
        pastix_int_t   *iparm,
        double         *dparm )
{
    pastix_data_t *pastix_data;
    pastix_spm_t *spm = NULL;
    int ret;
    size_t size;

    iparm[IPARM_ERROR_NUMBER] = PASTIX_SUCCESS;

    /*
     * Initialize iparm/dparm to default values and exit when called with
     * IPARM_MODIFY_PARAMETER set to API_NO
     */
    if (iparm[IPARM_MODIFY_PARAMETER] == API_NO)
    {
        pastixInitParam(iparm, dparm);
        iparm[IPARM_MODIFY_PARAMETER] = API_YES;
        return;
    }

    /*
     * Initialization step
     * Create the pastix_data structure and initialize the runtimes
     */
    if (iparm[IPARM_END_TASK] < API_TASK_INIT) {
        return;
    }

    if (iparm[IPARM_START_TASK] == API_TASK_INIT) {
        if (*pastix_data_ptr != NULL)
        {
            /*
             * Let's consider the user want to restart pastix with different
             * parameters
             */
            if (iparm[IPARM_VERBOSE] > 1) {
                pastix_print( 0, 0, "WARNING: PaStiX schedulers restarted\n" );
            }
            pastixFinalize( pastix_data_ptr, pastix_comm, iparm, dparm );
        }

        /*
         * Initialize pastix_data structure, and start scheduler(s)
         */
        pastixInit( pastix_data_ptr, pastix_comm, iparm, dparm );

        iparm[IPARM_START_TASK]++;
    }

    /*
     * Return now if only initialization is required
     */
    if (iparm[IPARM_END_TASK] < API_TASK_ORDERING) {
        return;
    }

    pastix_data = *pastix_data_ptr;

    /*
     * Initialize the internal spm structure.
     * We perform if only if starting step is lower than numerical
     * factorization, because further steps are using the internal bcsc for
     * computations with A.
     */
    if ( iparm[IPARM_START_TASK] <= API_TASK_NUMFACT) {
        if ( (pastix_data->csc != NULL) &&
             ((pastix_data->csc->n      != n)                       ||
              (pastix_data->csc->nnz    != (colptr[n] - colptr[0])) ||
              (pastix_data->csc->colptr != colptr)                  ||
              (pastix_data->csc->rowptr != row)) )
        {
            /*
             * This is a new csc, we need to delete the old one stored in pastix_data, and create a new one
             * We do not use spmExit, because the user allocated the fields.
             */
            memFree_null( pastix_data->csc );
        }

        if ( pastix_data->csc == NULL )
        {
            spm = malloc(sizeof( pastix_spm_t ));
            spmInit( spm );

            /*
             * Check and set the matrix type
             */
            if (iparm[IPARM_MTX_TYPE] == -1 ) {
                printf("Pastix old interface: you have to use --iparm iparm_mtx_type\n");
                iparm[IPARM_ERROR_NUMBER] = PASTIX_ERR_BADPARAMETER;
                spmExit( spm );
                free(spm);
                return  ;
            } else {
                spm->mtxtype = iparm[IPARM_MTX_TYPE];
            }

            spm->flttype = iparm[IPARM_FLOAT];
            spm->fmttype = PastixCSC;

            if ( iparm[IPARM_DOF_NBR] < 1 ) {
                fprintf(stderr,
                        "pastix: Variadic dofs are not supported in old pastix interface."
                        "        Please switch to the new interface to use this feature\n");
                memFree_null( spm );
                iparm[IPARM_ERROR_NUMBER] = PASTIX_ERR_NOTIMPLEMENTED;
                spmExit( spm );
                free(spm);
                return;
            }
            spm->dof  = iparm[IPARM_DOF_NBR];

            spm->n    = n;
            spm->nnz  = colptr[n] - colptr[0];

            spm->colptr = colptr;
            spm->rowptr = row;
            spm->values = avals;

            spmUpdateComputedFields( spm );

            pastix_data->csc = spm;
        }
        else {
            spm = (pastix_spm_t*)(pastix_data->csc);
        }

        /*
         * Update value field if given only at numerical steps
         */
        if ( spm->values == NULL ) {
            spm->values = avals;
        }
    }
    else {
        spm = (pastix_spm_t*)(pastix_data->csc);
    }

    /*
     * Ordering
     */
    if (iparm[IPARM_START_TASK] == API_TASK_ORDERING)
    {
        ret = pastix_subtask_order( pastix_data, spm, perm, invp );
        if (PASTIX_SUCCESS != ret)
        {
            iparm[IPARM_ERROR_NUMBER] = ret;
            return;
        }
        iparm[IPARM_START_TASK]++;
    }

    /*
     * Symbolic factorization
     */
    if (iparm[IPARM_END_TASK] < API_TASK_SYMBFACT) {
        return;
    }

    if (iparm[IPARM_START_TASK] == API_TASK_SYMBFACT)
    {
        ret = pastix_subtask_symbfact( pastix_data, perm, invp );
        if (PASTIX_SUCCESS != ret)
        {
            iparm[IPARM_ERROR_NUMBER] = ret;
            return;
        }
        iparm[IPARM_START_TASK]++;
    }

    /*
     * Analyze step
     */
    if (iparm[IPARM_END_TASK] < API_TASK_ANALYSE) {
        return;
    }

    if (iparm[IPARM_START_TASK] == API_TASK_ANALYSE)
    {
        ret = pastix_subtask_blend( pastix_data );
        if (PASTIX_SUCCESS != ret)
        {
            iparm[IPARM_ERROR_NUMBER] = ret;
            return;
        }
        iparm[IPARM_START_TASK]++;
    }

    /*
     * Numerical factorisation
     */
    if (iparm[IPARM_END_TASK] < API_TASK_NUMFACT) {
        return;
    }

    if (iparm[IPARM_START_TASK] == API_TASK_NUMFACT)
    {
        ret = pastix_task_numfact( pastix_data, spm );
        if (PASTIX_SUCCESS != ret) {
            iparm[IPARM_ERROR_NUMBER] = ret;
            return;
        }
        iparm[IPARM_START_TASK]++;
    }

    /*
     * Solve
     */
    if (iparm[IPARM_END_TASK] < API_TASK_SOLVE) {
        return;
    }

    if (iparm[IPARM_START_TASK] == API_TASK_SOLVE) {
        void *tmp;
        size = pastix_size_of( spm->flttype ) * spm->n;
        tmp = malloc(size);

        /*
         * Backup the initial b if we need to perform an iterative
         * refinement after the solve step
         */
        if (iparm[IPARM_END_TASK] > API_TASK_SOLVE) {
            memcpy(tmp, b, size);
            pastix_data->b = tmp;
        }
        pastix_task_solve( pastix_data, spm, nrhs, b, spm->n );
        iparm[IPARM_START_TASK]++;

        /*
         * Backup the first solution x0 if the user wants to come back later for
         * iterative refinement
         */
        if (iparm[IPARM_END_TASK] == API_TASK_SOLVE) {
            memcpy(tmp, b, size);
            pastix_data->x0 = tmp;
        }
    }

    /*
     * Refinement
     */
    if (iparm[IPARM_END_TASK] < API_TASK_REFINE) {
        return;
    }

    if (iparm[IPARM_START_TASK] == API_TASK_REFINE) {
        void *refineB  = pastix_data->b;
        void *refineX0 = pastix_data->x0;
        size = pastix_size_of( spm->flttype ) * spm->n;
        if ( !refineB ) {
            if ( !refineX0 ) {
                /*
                 * Neither b or x0 have been saved.
                 * Then, we need to start with x0 as a null vector. For that, we
                 * backup the original b, and we use the given b as x in the
                 * refinement step to store the solution in place.
                 */
                /* refineB = malloc(size); */
                /* memcpy(refineB, b, size); */

                /* refineX0 = b; */
                /* memset(refineX0, 0, size); */
                /* exit(0);  */
                /*
                 * Neither b and x0 have been saved, this should never happen.
                 */
                fprintf(stderr, "Neither b and x0 have been saved, this should never happen\n");
                iparm[IPARM_ERROR_NUMBER] = PASTIX_ERR_UNKNOWN;
                return;
            }
            else {
                /*
                 * x0 is saved, but not b. It means that we exit the pastix
                 * function call between the solve and refinemnet
                 * step. Therefor, b holds the original b.
                 */
                refineB = b;
            }
        }
        else {
            if ( !refineX0 ) {
                /*
                 * b is saved, but not x0. It means that we did not exit the
                 * pastix function call between solve and refinement steps.
                 * Therefor, b holds the initial solution x0 from the solve step.
                 */
                refineX0 = b;
            }
            else {
                /*
                 * Both x0 and b are saved. This should never happen.
                 */
                fprintf(stderr, "Both b and x0 are defined, this should never happen\n");
                iparm[IPARM_ERROR_NUMBER] = PASTIX_ERR_UNKNOWN;
                return;
            }
        }
        pastix_task_refine( pastix_data, refineX0, nrhs, refineB );
        iparm[IPARM_START_TASK]++;

        /*
         * Let's return the solution to the user
         */
        if ( b != refineX0 ) {
            memcpy(b, refineB, size);
        }

        if ( pastix_data->x0 ) {
            free( pastix_data->x0 );
            pastix_data->x0 = NULL;
        }
        if ( pastix_data->b ) {
            free( pastix_data->b );
            pastix_data->b = NULL;
        }
    }

    /*
     * Cleaning
     */
    if (iparm[IPARM_END_TASK] < API_TASK_CLEAN) {
        return;
    }

    if (iparm[IPARM_START_TASK] == API_TASK_CLEAN) {
        pastixFinalize( pastix_data_ptr, pastix_comm, iparm, dparm );
        iparm[IPARM_START_TASK]++;
    }
}
