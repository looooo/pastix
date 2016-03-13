/**
 *
 * @file sparse-matrix.h
 *
 * @author Mathieu Faverge
 * @date 2011-03-01
 * @precisions normal z -> c d s
 *
 **/

#ifndef _SPARSE_MATRIX_H_
#define _SPARSE_MATRIX_H_

#include <dague/data_distribution.h>

typedef struct sparse_matrix_desc_s {
    dague_ddesc_t         super;
    dague_data_t        **data_map;
    int                   typesze; /* Type size                           */
    SolverMatrix         *solvmtx;
} sparse_matrix_desc_t;

typedef struct sparse_vector_desc_t {
    dague_ddesc_t         super;
    int                   typesze; /* Type size                           */
    SolverMatrix         *solvmtx;
} sparse_vector_desc_t;

void sparse_matrix_init( sparse_matrix_desc_t *desc,
                         SolverMatrix *solvmtx,
                         int typesize, int nodes, int myrank );
void sparse_matrix_destroy( sparse_matrix_desc_t *desc );

void sparse_vector_init( sparse_vector_desc_t *desc,
                         int typesze, int nodes, int myrank);
void sparse_vector_destroy( sparse_vector_desc_t *desc );

void pastix_parsec_init( pastix_data_t *pastix,
                         int *argc, char **argv[] );
void pastix_parsec_finalize( pastix_data_t *pastix );

#endif /* _SPARSE_MATRIX_H_ */
