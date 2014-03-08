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

int core_zgeadd(int M, int N, pastix_complex64_t alpha,
                const pastix_complex64_t *A, int LDA,
                      pastix_complex64_t *B, int LDB);

void core_zaxpyt(int m, int n, pastix_complex64_t alpha,
                 pastix_complex64_t *A, int lda,
                 pastix_complex64_t *B, int ldb);

void core_zgetrfsp1d(pastix_complex64_t *L,
                     pastix_complex64_t *U,
                     SolverMatrix *datacode,
                     pastix_int_t c,
                     double criteria);

void core_zgetrfsp1d_gemm(pastix_int_t cblknum,
                          pastix_int_t bloknum,
                          pastix_int_t fcblknum,
                          pastix_complex64_t *L,
                          pastix_complex64_t *U,
                          pastix_complex64_t *Cl,
                          pastix_complex64_t *Cu,
                          pastix_complex64_t *work,
                          SolverMatrix *datacode);


void core_zhetrfsp1d(pastix_complex64_t *L,
                     pastix_complex64_t *work,
                     SolverMatrix *datacode,
                     pastix_int_t c,
                     double criteria);

void core_zhetrfsp1d_gemm(pastix_int_t cblknum,
                          pastix_int_t bloknum,
                          pastix_int_t fcblknum,
                          pastix_complex64_t *L,
                          pastix_complex64_t *C,
                          pastix_complex64_t *work1,
                          pastix_complex64_t *work2,
                          SolverMatrix *datacode);

int core_zpotrfsp1d(SolverMatrix *datacode,
                    pastix_int_t c,
                    pastix_complex64_t *L,
                    double criteria);

void core_zpotrfsp1d_gemm(SolverMatrix *datacode,
                          pastix_int_t cblknum,
                          pastix_int_t bloknum,
                          pastix_int_t fcblknum,
                          pastix_complex64_t *L,
                          pastix_complex64_t *C,
                          pastix_complex64_t *work);

void core_zsytrfsp1d(pastix_complex64_t *L,
                     pastix_complex64_t *work,
                     SolverMatrix *datacode,
                     pastix_int_t c,
                     double criteria);

void core_zsytrfsp1d_gemm(pastix_int_t cblknum,
                          pastix_int_t bloknum,
                          pastix_int_t fcblknum,
                          pastix_complex64_t *L,
                          pastix_complex64_t *C,
                          pastix_complex64_t *work1,
                          pastix_complex64_t *work2,
                          SolverMatrix *datacode);

#endif /* _CORE_Z_H_ */
