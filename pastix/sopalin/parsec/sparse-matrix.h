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

#include <stdarg.h>
#include <stdint.h>
#if defined(HAVE_CUDA)
#include <cuda.h>
#endif

/* #if !defined(PASTIX_STR_H) && !defined(_PASTIX_H_) */
/* struct pastix_data_t; */
/* typedef struct pastix_data_t pastix_data_t; */
/* #endif */

typedef struct sparse_matrix_desc_s {
    dague_ddesc_t         super;
    dague_data_t        **data_map;
    int                   typesze; /* Type size                           */
#if defined(HAVE_CUDA)
    CUdeviceptr          *d_blocktab;
#endif
    SolverMatrix         *solvmtx;
    pastix_int_t          gpu_limit;
} sparse_matrix_desc_t;

typedef struct sparse_vector_desc_t {
    dague_ddesc_t         super;
    int                   typesze; /* Type size                           */
    SolverMatrix         *solvmtx;
} sparse_vector_desc_t;

typedef struct sparse_context_s {
    int          format;     /* Matrix file format                         */
    int          factotype;
    int          coresnbr;   /* Number of cores to use for Pastix          */
    int          verbose;    /* Level of verbose                           */
    char        *matrixname; /* Filename to get the matrix                 */
    char        *rhsname;    /* Filename to get the matrix                 */
    char        *ordername;  /* Filename where the ordering is stored      */
    char        *symbname;   /* Filename where the symbol matrix is stored */
    char        *type;       /* Type of the matrix                         */
    char        *rhstype;    /* Type of the RHS                            */
    pastix_int_t  n;          /* Number of unknowns/columns/rows            */
    pastix_int_t  nnz;        /* Number of non-zero values in the input matrix */
    pastix_int_t *colptr;     /* Vector of size N+1 storing the starting point of each column in the array rows */
    pastix_int_t *rows;       /* Indices of the rows present in each column */
    void        *values;     /* Values of the matrix                       */
    void        *rhs;        /* Right Hand Side                            */
    pastix_int_t *permtab;    /* vector of permutation                      */
    pastix_int_t *peritab;    /* vector of inverse permutation              */
    pastix_int_t  iparm[IPARM_SIZE];
    double       dparm[DPARM_SIZE];
    sparse_matrix_desc_t *desc; /* Pointer to symbol matrix structure */
    sparse_vector_desc_t *rhsdesc; /* Pointer to symbol matrix structure */
} sparse_context_t;

void sparse_matrix_init( sparse_matrix_desc_t *desc,
                         SolverMatrix *solvmtx, int typesize, int nodes, int myrank);
void sparse_matrix_destroy( sparse_matrix_desc_t *desc );

void sparse_vector_init( sparse_vector_desc_t *desc,
                         int typesze, int nodes, int myrank);
void sparse_vector_destroy( sparse_vector_desc_t *desc );

double sparse_matrix_zrdmtx( sparse_context_t *dspctxt );
double sparse_matrix_crdmtx( sparse_context_t *dspctxt );
double sparse_matrix_drdmtx( sparse_context_t *dspctxt );
double sparse_matrix_srdmtx( sparse_context_t *dspctxt );

void sparse_matrix_zcheck( sparse_context_t *dspctxt );
void sparse_matrix_ccheck( sparse_context_t *dspctxt );
void sparse_matrix_dcheck( sparse_context_t *dspctxt );
void sparse_matrix_scheck( sparse_context_t *dspctxt );

void sparse_matrix_zclean( sparse_context_t *dspctxt );
void sparse_matrix_cclean( sparse_context_t *dspctxt );
void sparse_matrix_dclean( sparse_context_t *dspctxt );
void sparse_matrix_sclean( sparse_context_t *dspctxt );

void sparse_vector_zinit( sparse_context_t *dspctxt );
void sparse_vector_cinit( sparse_context_t *dspctxt );
void sparse_vector_dinit( sparse_context_t *dspctxt );
void sparse_vector_sinit( sparse_context_t *dspctxt );
void sparse_vector_zfinalize( sparse_context_t *dspctxt );
void sparse_vector_cfinalize( sparse_context_t *dspctxt );
void sparse_vector_dfinalize( sparse_context_t *dspctxt );
void sparse_vector_sfinalize( sparse_context_t *dspctxt );

#endif /* _SPARSE_MATRIX_H_ */
