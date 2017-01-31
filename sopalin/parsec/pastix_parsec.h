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

#include <parsec/data_distribution.h>

typedef struct sparse_matrix_desc_s {
    parsec_ddesc_t   super;
    parsec_data_t  **datamap_cblk;
    parsec_data_t  **datamap_blok;
    int             typesze;   /*< Type size                                                            */
    int             mtxtype;   /*< Matrix structure: PastixGeneral, PastixSymmetric or PastixHermitian. */
    SolverMatrix   *solvmtx;
    void          **d_blocktab;
} sparse_matrix_desc_t;

typedef struct sparse_vector_desc_t {
    parsec_ddesc_t         super;
    int                   typesze; /* Type size                           */
    SolverMatrix         *solvmtx;
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

#endif /* _SPARSE_MATRIX_H_ */
