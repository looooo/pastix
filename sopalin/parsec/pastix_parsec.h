/**
 *
 * @file pastix_parsec.h
 *
 * PaRSEC support for the numerical factorization and solve of PaStiX.
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2013-06-24
 *
 * @addtogroup pastix_parsec
 * @{
 *   This module describes the functionnality provided by the runtime system
 *   PaRSEC for the numerical factorization and solve.
 *
 **/
#ifndef _PASTIX_PARSEC_H_
#define _PASTIX_PARSEC_H_

#include <parsec/data_distribution.h>

/**
 * @brief PaRSEC descriptor stucture for the sparse matrix.
 */
typedef struct parsec_sparse_matrix_desc_s {
    parsec_ddesc_t  super;      /**< Every PaRSEC descriptors must inherit from parsec_desc_t                        */
    int             typesze;    /**< Arithmetic size                                                                 */
    int             mtxtype;    /**< Matrix structure: PastixGeneral, PastixSymmetric or PastixHermitian.            */
    SolverMatrix   *solvmtx;    /**< Solver matrix structure that describes the problem and stores the original data */
    void          **d_blocktab; /**< Pointer to GPU arrays that contains frownum,lrownum of each block for Fermi     */
} parsec_sparse_matrix_desc_t;

/**
 * @brief PaRSEC descriptor for the vectors linked to a given sparse matrix.
 */
typedef struct parsec_sparse_vector_desc_t {
    parsec_ddesc_t   super;  /**< Every PaRSEC descriptors must inherit from parsec_desc_t                        */
    int             typesze; /**< Arithmetic size                                                                 */
    SolverMatrix   *solvmtx; /**< Solver matrix structure that describes the problem and stores the original data */
} parsec_sparse_vector_desc_t;

void parsec_sparse_matrix_init( SolverMatrix *solvmtx,
                                int typesize, int mtxtype,
                                int nodes, int myrank );
void parsec_sparse_matrix_destroy( parsec_sparse_matrix_desc_t *desc );

void parsec_sparse_vector_init( parsec_sparse_vector_desc_t *desc,
                                int typesze, int nodes, int myrank );
void parsec_sparse_vector_destroy( parsec_sparse_vector_desc_t *desc );

void pastix_parsec_init( pastix_data_t *pastix,
                         int *argc, char **argv[] );
void pastix_parsec_finalize( pastix_data_t *pastix );

#endif /* _PASTIX_PARSEC_H_ */

/**
 *@}
 */
