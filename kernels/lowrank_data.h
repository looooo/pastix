/**
 * @file lowrank_data.h
 *
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Gregoire Pichon
 * @date 2016-11
 *
 **/
#ifndef _LOWRANK_DATA_H_
#define _LOWRANK_DATA_H_

void core_sge2lr_SVD_interface( double tol, pastix_int_t m, pastix_int_t n,
                                void *A, pastix_int_t lda,
                                void *Alr );
void core_dge2lr_SVD_interface( double tol, pastix_int_t m, pastix_int_t n,
                                void *A, pastix_int_t lda,
                                void *Alr );
void core_cge2lr_SVD_interface( double tol, pastix_int_t m, pastix_int_t n,
                                void *A, pastix_int_t lda,
                                void *Alr );
void core_zge2lr_SVD_interface( double tol, pastix_int_t m, pastix_int_t n,
                                void *A, pastix_int_t lda,
                                void *Alr );
void core_sge2lr_RRQR_interface( double tol, pastix_int_t m, pastix_int_t n,
                                 void *A, pastix_int_t lda,
                                 void *Alr );
void core_dge2lr_RRQR_interface( double tol, pastix_int_t m, pastix_int_t n,
                                 void *A, pastix_int_t lda,
                                 void *Alr );
void core_cge2lr_RRQR_interface( double tol, pastix_int_t m, pastix_int_t n,
                                 void *A, pastix_int_t lda,
                                 void *Alr );
void core_zge2lr_RRQR_interface( double tol, pastix_int_t m, pastix_int_t n,
                                 void *A, pastix_int_t lda,
                                 void *Alr );

int core_srradd_SVD_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_drradd_SVD_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_crradd_SVD_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_zrradd_SVD_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_srradd_RRQR_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_drradd_RRQR_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_crradd_RRQR_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);
int core_zrradd_RRQR_interface( double tol, int transA1, void *alpha,
                               pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                               pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                               pastix_int_t offx, pastix_int_t offy);

#endif /* _LOWRANK_DATA_H */

