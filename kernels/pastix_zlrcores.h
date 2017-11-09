/**
 * @file pastix_zcores.h
 *
 * PaStiX kernel header.
 *
 * @copyright 2011-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2017-10-05
 * @precisions normal z -> c d s
 *
 */
#ifndef _pastix_zlrcores_h_
#define _pastix_zlrcores_h_

#include "pastix_lowrank.h"

/**
 *
 * @addtogroup kernel_lr
 * @{
 *    This module contains all the low-rank kernels working on pastix_lr_t
 *    matrix representations.
 *
 *    @name PastixComplex64 low-rank kernels
 *    @{
 */
void core_zlralloc( pastix_int_t M, pastix_int_t N, pastix_int_t rkmax, pastix_lrblock_t *A );
void core_zlrfree ( pastix_lrblock_t *A );
int  core_zlrsze  ( int copy, pastix_int_t M, pastix_int_t N, pastix_lrblock_t *A, int newrk, int newrkmax, pastix_int_t rklimit );
int  core_zlr2ge  ( pastix_trans_t trans, pastix_int_t M, pastix_int_t N, const pastix_lrblock_t *Alr, pastix_complex64_t *A, pastix_int_t lda );

pastix_int_t core_zlr_rklimit( pastix_int_t M, pastix_int_t N );
void core_zlrcpy  ( const pastix_lr_t *lowrank,
                    pastix_trans_t transA, pastix_complex64_t alpha,
                    pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                    pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                    pastix_int_t offx, pastix_int_t offy );
void core_zlrconcatenate_u( pastix_complex64_t alpha,
                            pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                            pastix_int_t M2,                        pastix_lrblock_t *B,
                            pastix_int_t offx,
                            pastix_complex64_t *u1u2 );
void core_zlrconcatenate_v( pastix_trans_t transA1, pastix_complex64_t alpha,
                            pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                                             pastix_int_t N2,       pastix_lrblock_t *B,
                            pastix_int_t offy,
                            pastix_complex64_t *v1v2 );
pastix_fixdbl_t core_zlrorthu( pastix_trans_t transV, pastix_int_t M,  pastix_int_t N, pastix_int_t K,
                               pastix_complex64_t *U, pastix_int_t ldu,
                               pastix_complex64_t *V, pastix_int_t ldv );


typedef struct core_zlrmm_s {
    const pastix_lr_t      *lowrank;
    pastix_trans_t          transA;
    pastix_trans_t          transB;
    pastix_int_t            M, N, K;
    pastix_int_t            Cm, Cn;
    pastix_int_t            offx, offy;
    pastix_complex64_t      alpha;
    const pastix_lrblock_t *A;
    const pastix_lrblock_t *B;
    pastix_complex64_t      beta;
    pastix_lrblock_t       *C;
    pastix_complex64_t     *work;
    pastix_int_t            lwork;
    pastix_atomic_lock_t   *lock;
} core_zlrmm_t;

#define PASTE_CORE_ZLRMM_PARAMS(_a_)                   \
    const pastix_lr_t      *lowrank = (_a_)->lowrank;  \
    pastix_trans_t          transA  = (_a_)->transA;   \
    pastix_trans_t          transB  = (_a_)->transB;   \
    pastix_int_t            M       = (_a_)->M;        \
    pastix_int_t            N       = (_a_)->N;        \
    pastix_int_t            K       = (_a_)->K;        \
    pastix_int_t            Cm      = (_a_)->Cm;       \
    pastix_int_t            Cn      = (_a_)->Cn;       \
    pastix_int_t            offx    = (_a_)->offx;     \
    pastix_int_t            offy    = (_a_)->offy;     \
    pastix_complex64_t      alpha   = (_a_)->alpha;    \
    const pastix_lrblock_t *A       = (_a_)->A;        \
    const pastix_lrblock_t *B       = (_a_)->B;        \
    pastix_complex64_t      beta    = (_a_)->beta;     \
    pastix_lrblock_t       *C       = (_a_)->C;        \
    pastix_complex64_t     *work    = (_a_)->work;     \
    pastix_int_t            lwork   = (_a_)->lwork;    \
    pastix_atomic_lock_t   *lock    = (_a_)->lock;

#define PASTE_CORE_ZLRMM_VOID                   \
    (void)lowrank;                              \
    (void)transA;                               \
    (void)transB;                               \
    (void)M;                                    \
    (void)N;                                    \
    (void)K;                                    \
    (void)Cm;                                   \
    (void)Cn;                                   \
    (void)offx;                                 \
    (void)offy;                                 \
    (void)alpha;                                \
    (void)A;                                    \
    (void)B;                                    \
    (void)beta;                                 \
    (void)C;                                    \
    (void)work;                                 \
    (void)lwork;                                \
    (void)lock

static inline pastix_complex64_t *
core_zlrmm_resize_ws( core_zlrmm_t *params,
                      ssize_t newsize )
{
    if ( params->lwork < newsize ) {
        if ( params->lwork == -1 ) {
            params->work  = malloc( newsize * sizeof(pastix_complex64_t) );
        }
        else {
            params->work  = realloc( params->work, newsize * sizeof(pastix_complex64_t) );
            params->lwork = newsize;
        }
    }
    return params->work;
}

//void core_zlrmm( core_zlrmm_t *params );
pastix_fixdbl_t
core_zfrfr2lr( core_zlrmm_t     *params,
               pastix_lrblock_t *AB,
               int              *infomask,
               pastix_int_t      Kmax );

pastix_fixdbl_t
core_zfrlr2lr( core_zlrmm_t     *params,
               pastix_lrblock_t *AB,
               int              *infomask,
               pastix_int_t      Brkmin );

pastix_fixdbl_t
core_zlrfr2lr( core_zlrmm_t     *params,
               pastix_lrblock_t *AB,
               int              *infomask,
               pastix_int_t      Arkmin );

pastix_fixdbl_t
core_zlrlr2lr( core_zlrmm_t     *params,
               pastix_lrblock_t *AB,
               int              *infomask );

void core_zlrmm_Cfr( core_zlrmm_t *params );
void core_zlrmm_Clr( core_zlrmm_t *params );
void core_zlrmm_Cnull( core_zlrmm_t *params );

void core_zlrmm( const pastix_lr_t *lowrank,
                 pastix_trans_t transA, pastix_trans_t transB,
                 pastix_int_t M, pastix_int_t N, pastix_int_t K,
                 pastix_int_t Cm, pastix_int_t Cn,
                 pastix_int_t offx, pastix_int_t offy,
                 pastix_complex64_t alpha, const pastix_lrblock_t *A,
                 const pastix_lrblock_t *B,
                 pastix_complex64_t beta,        pastix_lrblock_t *C,
                 pastix_complex64_t *work, pastix_int_t ldwork,
                 pastix_atomic_lock_t *lock );

/**
 *     @}
 * @}
 *
 * @addtogroup kernel_lr_svd
 * @{
 *    This is the SVD implementation of the low-rank kernels based on the LAPACK
 *    GESVD function.
 *
 *    @name PastixComplex64 SVD low-rank kernels
 *    @{
 */
pastix_fixdbl_t core_zge2lr_SVD( double tol, pastix_int_t rklimit, pastix_int_t m, pastix_int_t n,
                                 const pastix_complex64_t *A, pastix_int_t lda,
                                 pastix_lrblock_t *Alr );
int  core_zrradd_SVD( const pastix_lr_t *lowrank, pastix_trans_t transA1, pastix_complex64_t alpha,
                      pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                      pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                      pastix_int_t offx, pastix_int_t offy );

pastix_fixdbl_t core_zge2lr_SVD_interface( pastix_fixdbl_t tol, pastix_int_t rklimit,
                                           pastix_int_t m, pastix_int_t n,
                                           const void *Aptr, pastix_int_t lda, void *Alr );
int  core_zrradd_SVD_interface( const pastix_lr_t *lowrank, pastix_trans_t transA1, const void *alphaptr,
                                pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                                pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                                pastix_int_t offx, pastix_int_t offy );

/**
 *     @}
 * @}
 *
 * @addtogroup kernel_lr_rrqr
 * @{
 *    This is the rank-revealing QR implementation of the low-rank kernels based
 *    on the modified LAPACK GEQP3 function.
 *
 *    @name PastixComplex64 RRQR low-rank kernels
 *    @{
 */
pastix_fixdbl_t core_zge2lr_RRQR( double tol, pastix_int_t rklimit, pastix_int_t m, pastix_int_t n,
                                  const pastix_complex64_t *A, pastix_int_t lda,
                                  pastix_lrblock_t *Alr );
int  core_zrradd_RRQR( const pastix_lr_t *lowrank, pastix_trans_t transA1, pastix_complex64_t alpha,
                       pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                       pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                       pastix_int_t offx, pastix_int_t offy );


pastix_fixdbl_t core_zge2lr_RRQR_interface( pastix_fixdbl_t tol, pastix_int_t rklimit,
                                            pastix_int_t m, pastix_int_t n,
                                            const void *Aptr, pastix_int_t lda, void *Alr );

int  core_zrradd_RRQR_interface( const pastix_lr_t *lowrank, pastix_trans_t transA1, const void *alphaptr,
                                 pastix_int_t M1, pastix_int_t N1, const pastix_lrblock_t *A,
                                 pastix_int_t M2, pastix_int_t N2,       pastix_lrblock_t *B,
                                 pastix_int_t offx, pastix_int_t offy );
/**
 *     @}
 * @}
 *
 */
#endif /* _pastix_zlrcores_h_ */
