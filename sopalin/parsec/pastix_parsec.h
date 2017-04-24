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
typedef struct sparse_matrix_desc_s {
    parsec_ddesc_t  super;      /**< Every PaRSEC descriptors must inherit from parsec_desc_t                        */
    int             typesze;    /**< Arithmetic size                                                                 */
    int             mtxtype;    /**< Matrix structure: PastixGeneral, PastixSymmetric or PastixHermitian.            */
    SolverMatrix   *solvmtx;    /**< Solver matrix structure that describes the problem and stores the original data */
    void          **d_blocktab; /**< Pointer to GPU arrays that contains frownum,lrownum of each block for Fermi     */
} sparse_matrix_desc_t;

/**
 * @brief PaRSEC descriptor for the vectors linked to a given sparse matrix.
 */
typedef struct sparse_vector_desc_t {
    parsec_ddesc_t   super;  /**< Every PaRSEC descriptors must inherit from parsec_desc_t                        */
    int             typesze; /**< Arithmetic size                                                                 */
    SolverMatrix   *solvmtx; /**< Solver matrix structure that describes the problem and stores the original data */
} sparse_vector_desc_t;

void sparse_matrix_init( sparse_matrix_desc_t *desc,
                         SolverMatrix *solvmtx,
                         int typesize, int mtxtype,
                         int nodes, int myrank );
void sparse_matrix_destroy( sparse_matrix_desc_t *desc );

void sparse_vector_init( sparse_vector_desc_t *desc,
                         int typesze, int nodes, int myrank);
void sparse_vector_destroy( sparse_vector_desc_t *desc );

void pastix_parsec_init( pastix_data_t *pastix,
                         int *argc, char **argv[] );
void pastix_parsec_finalize( pastix_data_t *pastix );

#endif /* _PASTIX_PARSEC_H_ */

/**
 *@}
 */
