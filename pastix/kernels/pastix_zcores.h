/*
 * Copyright (c) 2010      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 */
#ifndef _PASTIX_ZCORES_H_
#define _PASTIX_ZCORES_H_

#define pastix_cblk_lock( cblk_ )
#define pastix_cblk_unlock( cblk_ )

#if defined(NAPA_SOPALIN)
static inline int is_block_inside_fblock( SolverBlok *blok, SolverBlok *fblok ) {
    return (((blok->frownum >= fblok->frownum) &&
             (blok->lrownum <= fblok->lrownum)) ||
            ((blok->frownum <= fblok->frownum) &&
             (blok->lrownum >= fblok->lrownum)) ||
            ((blok->frownum <= fblok->frownum) &&
             (blok->lrownum >= fblok->frownum)) ||
            ((blok->frownum <= fblok->lrownum) &&
             (blok->lrownum >= fblok->lrownum)));
}
#else
static inline int is_block_inside_fblock( SolverBlok *blok, SolverBlok *fblok ) {
    return ((blok->frownum >= fblok->frownum) &&
            (blok->lrownum <= fblok->lrownum));
}
#endif /* defined(NAPA_SOPALIN) */

void core_zgetro(int m, int n,
                 pastix_complex64_t *A, int lda,
                 pastix_complex64_t *B, int ldb);

int core_zgeadd(int trans, int M, int N, pastix_complex64_t alpha,
                const pastix_complex64_t *A, int LDA,
                      pastix_complex64_t *B, int LDB);

int core_zgemdm(int transA, int transB,
                int M, int N, int K,
                pastix_complex64_t alpha, pastix_complex64_t *A, int LDA,
                pastix_complex64_t *B, int LDB,
                pastix_complex64_t beta, pastix_complex64_t *C, int LDC,
                pastix_complex64_t *D, int incD,
                pastix_complex64_t *WORK, int LWORK);


int core_zgetrfsp1d( SolverMatrix *datacode,
                     pastix_int_t c,
                     pastix_complex64_t *L,
                     pastix_complex64_t *U,
                     double criteria );

void core_zgetrfsp1d_gemm( SolverMatrix *datacode,
                           pastix_int_t cblknum,
                           pastix_int_t bloknum,
                           pastix_int_t fcblknum,
                           pastix_complex64_t *L,
                           pastix_complex64_t *U,
                           pastix_complex64_t *Cl,
                           pastix_complex64_t *Cu,
                           pastix_complex64_t *work );

int core_zhetrfsp1d( SolverMatrix *datacode,
                     pastix_int_t c,
                     pastix_complex64_t *L,
                     double criteria,
                     pastix_complex64_t *work );

void core_zhetrfsp1d_gemm( SolverMatrix *datacode,
                           pastix_int_t cblknum,
                           pastix_int_t bloknum,
                           pastix_int_t fcblknum,
                           pastix_complex64_t *L,
                           pastix_complex64_t *C,
                           pastix_complex64_t *work1,
                           pastix_complex64_t *work2 );

int core_zpotrfsp1d( SolverMatrix *datacode,
                     pastix_int_t c,
                     pastix_complex64_t *L,
                     double criteria );

void core_zpotrfsp1d_gemm( SolverMatrix *datacode,
                           pastix_int_t cblknum,
                           pastix_int_t bloknum,
                           pastix_int_t fcblknum,
                           pastix_complex64_t *L,
                           pastix_complex64_t *C,
                           pastix_complex64_t *work);

int core_zsytrfsp1d( SolverMatrix *datacode,
                     pastix_int_t c,
                     pastix_complex64_t *L,
                     double criteria,
                     pastix_complex64_t *work );

void core_zsytrfsp1d_gemm( SolverMatrix *datacode,
                           pastix_int_t cblknum,
                           pastix_int_t bloknum,
                           pastix_int_t fcblknum,
                           pastix_complex64_t *L,
                           pastix_complex64_t *C,
                           pastix_complex64_t *work1,
                           pastix_complex64_t *work2 );

#endif /* _CORE_Z_H_ */
