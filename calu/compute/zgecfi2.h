/**
 *
 * @file zgecfi2.h
 *
 *  PLASMA InPlaceTransformation module
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 *  This work is the implementation of an inplace transformation 
 *  based on the GKK algorithm by Gustavson, Karlsson, Kagstrom 
 *  and its fortran implementation.
 *
 * @version 2.5.0
 * @author Mathieu Faverge
 * @date 2010-11-15
 *
 * @precisions normal z -> c d s
 *
 **/

#ifndef ZGECFI2_H
#define ZGECFI2_H

#define ipt_call( name, m1, n1, mb, nb ) \
  ipt_z##name(plasma, (m1),     (n1),     (A+A11), (mb),     (nb),     sequence, request); \
  ipt_z##name(plasma, (m1),     (n-(n1)), (A+A12), (mb),     (n-(n1)), sequence, request); \
  ipt_z##name(plasma, (m-(m1)), (n1),     (A+A21), (m-(m1)), (nb),     sequence, request); \
  ipt_z##name(plasma, (m-(m1)), (n-(n1)), (A+A22), (m-(m1)), (n-(n1)), sequence, request);

#define ipt_cal2( name, m1, n1, mb, nb ) \
  ipt_z##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m1),     (n1),     (A+A11), (mb),     (nb),     sequence, request); \
  ipt_z##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m1),     (n-(n1)), (A+A12), (mb),     (n-(n1)), sequence, request); \
  ipt_z##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m-(m1)), (n1),     (A+A21), (m-(m1)), (nb),     sequence, request); \
  ipt_z##name(plasma, PlasmaIPT_NoDep, PlasmaIPT_NoDep, (m-(m1)), (n-(n1)), (A+A22), (m-(m1)), (n-(n1)), sequence, request);

/* one transformation */
#define ipt_zrm2rrrb(  plasma, m, n, A, mb, nb, seq, req) ipt_zcm2ccrb(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
#define ipt_zrrrb2rm(  plasma, m, n, A, mb, nb, seq, req) ipt_zccrb2cm(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
int ipt_zcm2ccrb  (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zccrb2cm  (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

int ipt_zccrb2crrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zcrrb2ccrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zrcrb2rrrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zrrrb2rcrb(plasma_context_t *plasma, PLASMA_enum idep, PLASMA_enum odep, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

#define ipt_zcrrb2rrrb(plasma, m, n, A, mb, nb, seq, req) ipt_zccrb2rcrb((plasma), (m), (n), (A), (mb), (nb), (seq), (req));
#define ipt_zrcrb2ccrb(plasma, m, n, A, mb, nb, seq, req) ipt_zccrb2rcrb((plasma), (n), (m), (A), (nb), (mb), (seq), (req));
#define ipt_zrrrb2crrb(plasma, m, n, A, mb, nb, seq, req) ipt_zccrb2rcrb((plasma), (n), (m), (A), (nb), (mb), (seq), (req));
int ipt_zccrb2rcrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

/* 2 transformations */
#define ipt_zrm2crrb(  plasma, m, n, A, mb, nb, seq, req) ipt_zcm2rcrb(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
#define ipt_zcrrb2rm(  plasma, m, n, A, mb, nb, seq, req) ipt_zrcrb2cm(  (plasma), (n), (m), (A), (nb), (mb), (seq), (req));
int ipt_zcm2rcrb  (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zrcrb2cm  (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

int ipt_zccrb2rrrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zrrrb2ccrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zcrrb2rcrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zrcrb2crrb(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

int ipt_zcm2crrb  (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zcrrb2cm  (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zrcrb2rm  (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zrm2rcrb  (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

/* 3 transformations */
int ipt_zcm2rrrb  (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zrrrb2cm  (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zccrb2rm  (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zrm2ccrb  (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);

/* 4 transformations */
int ipt_zcm2rm    (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);
int ipt_zrm2cm    (plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req);


int ipt_zpanel2all(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req); 
int ipt_zall2panel(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req); 
int ipt_zpanel2tile(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req); 
int ipt_ztile2panel(plasma_context_t *plasma, int m, int n, PLASMA_Complex64_t *A, int mb, int nb, PLASMA_sequence *seq, PLASMA_request *req); 
#endif /* ZGECFI2_H */
