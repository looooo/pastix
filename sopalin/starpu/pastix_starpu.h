/**
 *
 * @file pastix_starpu.h
 *
 * StarPU support for the numerical factorization and solve of PaStiX.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2013-06-24
 *
 * @addtogroup pastix_starpu
 * @{
 *   This module describes the functionnality provided by the runtime system
 *   StarPU for the numerical factorization and solve.
 *
 **/
#ifndef _PASTIX_STARPU_H_
#define _PASTIX_STARPU_H_

#include <starpu.h>

#if defined(PASTIX_WITH_MPI)
#define starpu_insert_task starpu_mpi_insert_task
#define pastix_codelet(_codelet_) MPI_COMM_WORLD, _codelet_
#else
#define pastix_codelet(_codelet_) _codelet_
#endif

/**
 * @brief StarPU descriptor stucture for the sparse matrix.
 */
typedef struct starpu_sparse_matrix_desc_s {
    int             typesze;    /**< Arithmetic size                                                                 */
    int             mtxtype;    /**< Matrix structure: PastixGeneral, PastixSymmetric or PastixHermitian.            */
    SolverMatrix   *solvmtx;    /**< Solver matrix structure that describes the problem and stores the original data */
    void          **d_blocktab; /**< Pointer to GPU arrays that contains frownum,lrownum of each block for Fermi     */
} starpu_sparse_matrix_desc_t;

/**
 * @brief StarPU descriptor for the vectors linked to a given sparse matrix.
 */
typedef struct starpu_sparse_vector_desc_t {
    int             typesze; /**< Arithmetic size                                                                 */
    SolverMatrix   *solvmtx; /**< Solver matrix structure that describes the problem and stores the original data */
} starpu_sparse_vector_desc_t;

void starpu_sparse_matrix_init( starpu_sparse_matrix_desc_t *desc,
                                SolverMatrix *solvmtx,
                                int typesize, int mtxtype,
                                int nodes, int myrank );
void starpu_sparse_matrix_destroy( starpu_sparse_matrix_desc_t *desc );

void starpu_sparse_vector_init( starpu_sparse_vector_desc_t *desc,
                                int typesze, int nodes, int myrank );
void starpu_sparse_vector_destroy( starpu_sparse_vector_desc_t *desc );

void pastix_starpu_init( pastix_data_t *pastix,
                         int *argc, char **argv[] );
void pastix_starpu_finalize( pastix_data_t *pastix );

#endif /* _PASTIX_STARPU_H_ */

/**
 *@}
 */
