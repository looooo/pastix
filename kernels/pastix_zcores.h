/**
 * @file pastix_zcores.h
 *
 * PaStiX kernel header.
 *
 * @copyright 2011-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Esragul Korkmaz
 * @author Gregoire Pichon
 * @author Tony Delarue
 * @author Alycia Lisito
 * @date 2023-03-30
 * @precisions normal z -> c d s
 *
 */
#ifndef _pastix_zcores_h_
#define _pastix_zcores_h_

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#define pastix_cblk_lock( cblk_ )    pastix_atomic_lock( &((cblk_)->lock) )
#define pastix_cblk_unlock( cblk_ )  pastix_atomic_unlock( &((cblk_)->lock) )
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/**
 * @addtogroup kernel_blas_lapack
 * @{
 *    This module contains all the BLAS and LAPACK-like kernels that are working
 *    on lapack layout matrices.
 *
 *    @name PastixComplex64 BLAS kernels
 *    @{
 */
void core_zplrnt( int                    m,
                  int                    n,
                  pastix_complex64_t    *A,
                  int                    lda,
                  int                    gM,
                  int                    m0,
                  int                    n0,
                  unsigned long long int seed );
void core_zgetmo( int                       m,
                  int                       n,
                  const pastix_complex64_t *A,
                  int                       lda,
                  pastix_complex64_t       *B,
                  int                       ldb );
int core_zgeadd( pastix_trans_t            trans,
                 pastix_int_t              M,
                 pastix_int_t              N,
                 pastix_complex64_t        alpha,
                 const pastix_complex64_t *A,
                 pastix_int_t              LDA,
                 pastix_complex64_t        beta,
                 pastix_complex64_t       *B,
                 pastix_int_t              LDB );
int core_zgemdm( pastix_trans_t            transA,
                 pastix_trans_t            transB,
                 int                       M,
                 int                       N,
                 int                       K,
                 pastix_complex64_t        alpha,
                 const pastix_complex64_t *A,
                 int                       LDA,
                 const pastix_complex64_t *B,
                 int                       LDB,
                 pastix_complex64_t        beta,
                 pastix_complex64_t       *C,
                 int                       LDC,
                 const pastix_complex64_t *D,
                 int                       incD,
                 pastix_complex64_t       *WORK,
                 int                       LWORK );
int core_zpqrcp( double              tol,
                 pastix_int_t        maxrank,
                 int                 full_update,
                 pastix_int_t        nb,
                 pastix_int_t        m,
                 pastix_int_t        n,
                 pastix_complex64_t *A,
                 pastix_int_t        lda,
                 pastix_int_t       *jpvt,
                 pastix_complex64_t *tau,
                 pastix_complex64_t *work,
                 pastix_int_t        lwork,
                 double             *rwork );
int core_zrqrcp( double              tol,
                 pastix_int_t        maxrank,
                 int                 refine,
                 pastix_int_t        nb,
                 pastix_int_t        m,
                 pastix_int_t        n,
                 pastix_complex64_t *A,
                 pastix_int_t        lda,
                 pastix_int_t       *jpvt,
                 pastix_complex64_t *tau,
                 pastix_complex64_t *work,
                 pastix_int_t        lwork,
                 double             *rwork );
int core_zrqrrt( double              tol,
                 pastix_int_t        maxrank,
                 pastix_int_t        nb,
                 pastix_int_t        m,
                 pastix_int_t        n,
                 pastix_complex64_t *A,
                 pastix_int_t        lda,
                 pastix_complex64_t *tau,
                 pastix_complex64_t *B,
                 pastix_int_t        ldb,
                 pastix_complex64_t *tau_b,
                 pastix_complex64_t *work,
                 pastix_int_t        lwork,
                 double              normA );
int core_ztqrcp( double              tol,
                 pastix_int_t        maxrank,
                 int                 unused,
                 pastix_int_t        nb,
                 pastix_int_t        m,
                 pastix_int_t        n,
                 pastix_complex64_t *A,
                 pastix_int_t        lda,
                 pastix_int_t       *jpvt,
                 pastix_complex64_t *tau,
                 pastix_complex64_t *work,
                 pastix_int_t        lwork,
                 double             *rwork );
int core_ztradd( pastix_uplo_t             uplo,
                 pastix_trans_t            trans,
                 pastix_int_t              M,
                 pastix_int_t              N,
                 pastix_complex64_t        alpha,
                 const pastix_complex64_t *A,
                 pastix_int_t              LDA,
                 pastix_complex64_t        beta,
                 pastix_complex64_t       *B,
                 pastix_int_t              LDB);
int core_zscalo( pastix_trans_t            trans,
                 pastix_int_t              M,
                 pastix_int_t              N,
                 const pastix_complex64_t *A,
                 pastix_int_t              lda,
                 const pastix_complex64_t *D,
                 pastix_int_t              ldd,
                 pastix_complex64_t       *B,
                 pastix_int_t              ldb );

/**
 *    @}
 *    @name PastixComplex64 Othogonalization kernels for low-rank updates
 *    @{
 */
pastix_fixdbl_t core_zlrorthu_fullqr( pastix_int_t        M,
                                      pastix_int_t        N,
                                      pastix_int_t        rank,
                                      pastix_complex64_t *U,
                                      pastix_int_t        ldu,
                                      pastix_complex64_t *V,
                                      pastix_int_t        ldv );
pastix_fixdbl_t core_zlrorthu_partialqr( pastix_int_t        M,
                                         pastix_int_t        N,
                                         pastix_int_t        r1,
                                         pastix_int_t       *r2ptr,
                                         pastix_int_t        offx,
                                         pastix_int_t        offy,
                                         pastix_complex64_t *U,
                                         pastix_int_t        ldu,
                                         pastix_complex64_t *V,
                                         pastix_int_t        ldv );
pastix_fixdbl_t core_zlrorthu_cgs( pastix_int_t        M1,
                                   pastix_int_t        N1,
                                   pastix_int_t        M2,
                                   pastix_int_t        N2,
                                   pastix_int_t        r1,
                                   pastix_int_t       *r2ptr,
                                   pastix_int_t        offx,
                                   pastix_int_t        offy,
                                   pastix_complex64_t *U,
                                   pastix_int_t        ldu,
                                   pastix_complex64_t *V,
                                   pastix_int_t        ldv );

/**
 *    @}
 *    @name PastixComplex64 LAPACK kernels
 *    @{
 */
void core_zpotrfsp( pastix_int_t        n,
                    pastix_complex64_t *A,
                    pastix_int_t        lda,
                    pastix_int_t       *nbpivots,
                    double              criterion );
void core_zpxtrfsp( pastix_int_t        n,
                    pastix_complex64_t *A,
                    pastix_int_t        lda,
                    pastix_int_t       *nbpivots,
                    double              criterion );
void core_zgetrfsp( pastix_int_t        n,
                    pastix_complex64_t *A,
                    pastix_int_t        lda,
                    pastix_int_t       *nbpivots,
                    double              criterion );
#if defined(PRECISION_z) || defined(PRECISION_c)
void core_zhetrfsp( pastix_int_t        n,
                    pastix_complex64_t *A,
                    pastix_int_t        lda,
                    pastix_int_t       *nbpivots,
                    double              criterion );
#endif
void core_zsytrfsp( pastix_int_t        n,
                    pastix_complex64_t *A,
                    pastix_int_t        lda,
                    pastix_int_t       *nbpivots,
                    double              criterion );

/**
 *    @}
 * @}
 *
 * @addtogroup kernel_fact
 * @{
 *    This module contains all the kernel working at the solver matrix structure
 *    level for the numerical factorization step.
 *
 *    @name PastixComplex64 cblk-BLAS CPU kernels
 *    @{
 */

int cpucblk_zgeaddsp1d( const SolverCblk         *cblk1,
                        SolverCblk               *cblk2,
                        const pastix_complex64_t *L1,
                        pastix_complex64_t       *L2,
                        const pastix_complex64_t *U1,
                        pastix_complex64_t       *U2 );

pastix_fixdbl_t cpucblk_zgemmsp( pastix_coefside_t   sideA,
                                 pastix_trans_t      trans,
                                 const SolverCblk   *cblk,
                                 const SolverBlok   *blok,
                                 SolverCblk         *fcblk,
                                 const void         *A,
                                 const void         *B,
                                 void               *C,
                                 pastix_complex64_t *work,
                                 pastix_int_t        lwork,
                                 const pastix_lr_t  *lowrank );
void cpucblk_ztrsmsp( pastix_side_t      side,
                      pastix_uplo_t      uplo,
                      pastix_trans_t     trans,
                      pastix_diag_t      diag,
                      const SolverCblk  *cblk,
                      const void        *A,
                      void              *C,
                      const pastix_lr_t *lowrank );
void cpucblk_zscalo ( pastix_trans_t     trans,
                      SolverCblk        *cblk,
                      void              *dataL,
                      void              *dataLD );

pastix_fixdbl_t cpublok_zgemmsp( pastix_trans_t     trans,
                                 const SolverCblk  *cblk,
                                 SolverCblk        *fcblk,
                                 pastix_int_t       blok_mk,
                                 pastix_int_t       blok_nk,
                                 pastix_int_t       blok_mn,
                                 const void        *A,
                                 const void        *B,
                                 void              *C,
                                 const pastix_lr_t *lowrank );
pastix_fixdbl_t cpublok_ztrsmsp( pastix_side_t      side,
                                 pastix_uplo_t      uplo,
                                 pastix_trans_t     trans,
                                 pastix_diag_t      diag,
                                 const SolverCblk  *cblk,
                                 pastix_int_t       blok_m,
                                 const void        *A,
                                 void              *C,
                                 const pastix_lr_t *lowrank );
void cpublok_zscalo ( pastix_trans_t  trans,
                      SolverCblk     *cblk,
                      pastix_int_t    blok_m,
                      const void     *A,
                      const void     *dataD,
                      void           *dataB );

/**
 *    @}
 *    @name PastixComplex64 cblk LU kernels
 *    @{
 */
int cpucblk_zgetrfsp1d_getrf( SolverMatrix *solvmtx,
                              SolverCblk   *cblk,
                              void         *L,
                              void         *U );
int cpucblk_zgetrfsp1d_panel( SolverMatrix *solvmtx,
                              SolverCblk   *cblk,
                              void         *L,
                              void         *U );
int cpucblk_zgetrfsp1d      ( SolverMatrix       *solvmtx,
                              SolverCblk         *cblk,
                              pastix_complex64_t *work,
                              pastix_int_t        lwork );

/**
 *    @}
 *    @name PastixComplex64 cblk Cholesky kernels
 *    @{
 */
int cpucblk_zpotrfsp1d_potrf( SolverMatrix *solvmtx,
                              SolverCblk   *cblk,
                              void         *dataL );
int cpucblk_zpotrfsp1d_panel( SolverMatrix *solvmtx,
                              SolverCblk   *cblk,
                              void         *L );
int cpucblk_zpotrfsp1d( SolverMatrix       *solvmtx,
                        SolverCblk         *cblk,
                        pastix_complex64_t *work,
                        pastix_int_t        lwork );

/**
 *    @}
 */

#if defined(PRECISION_z) || defined(PRECISION_c)
 /**
 *    @name PastixComplex64 cblk LDL^h kernels
 *    @{
 */
int cpucblk_zhetrfsp1d_hetrf( SolverMatrix *solvmtx,
                              SolverCblk   *cblk,
                              void         *L );
int cpucblk_zhetrfsp1d_panel( SolverMatrix *solvmtx,
                              SolverCblk   *cblk,
                              void         *L,
                              void         *DLh );
int cpucblk_zhetrfsp1d( SolverMatrix       *solvmtx,
                        SolverCblk         *cblk,
                        pastix_complex64_t *work1,
                        pastix_complex64_t *work2,
                        pastix_int_t        lwork );

/**
 *    @}
 *    @name PastixComplex64 cblk LL^t kernels
 *    @{
 */
int cpucblk_zpxtrfsp1d_pxtrf( SolverMatrix *solvmtx,
                              SolverCblk   *cblk,
                              void         *dataL );
int cpucblk_zpxtrfsp1d_panel( SolverMatrix *solvmtx,
                              SolverCblk   *cblk,
                              void         *dataL );
int cpucblk_zpxtrfsp1d( SolverMatrix       *solvmtx,
                        SolverCblk         *cblk,
                        pastix_complex64_t *work,
                        pastix_int_t        lwork );

/**
 *    @}
 */
#endif

 /**
 *    @name PastixComplex64 cblk LDL^t kernels
 *    @{
 */
int cpucblk_zsytrfsp1d_sytrf( SolverMatrix *solvmtx,
                              SolverCblk   *cblk,
                              void         *dataL );
int cpucblk_zsytrfsp1d_panel( SolverMatrix *solvmtx,
                              SolverCblk   *cblk,
                              void         *L,
                              void         *DLt );
int cpucblk_zsytrfsp1d( SolverMatrix       *solvmtx,
                        SolverCblk         *cblk,
                        pastix_complex64_t *Dlt,
                        pastix_complex64_t *work,
                        pastix_int_t        lwork );

/**
 *    @}
 *    @name PastixComplex64 initialization and additionnal routines
 *    @{
 */
void cpucblk_zalloc_lrws( const SolverCblk   *cblk,
                          pastix_lrblock_t   *lrblok,
                          pastix_complex64_t *ws );
void cpucblk_zalloc_lr( pastix_coefside_t  side,
                        SolverCblk        *cblk,
                        int                rkmax );
void cpucblk_zalloc_fr( pastix_coefside_t  side,
                        SolverCblk        *cblk );
void cpucblk_zalloc( pastix_coefside_t  side,
                     SolverCblk        *cblk );
void cpucblk_zfree( pastix_coefside_t  side,
                    SolverCblk        *cblk );
void cpucblk_zfillin( pastix_coefside_t    side,
                      const SolverMatrix  *solvmtx,
                      const pastix_bcsc_t *bcsc,
                      pastix_int_t         itercblk );
void cpucblk_zinit( pastix_coefside_t    side,
                    const SolverMatrix  *solvmtx,
                    const pastix_bcsc_t *bcsc,
                    pastix_int_t         itercblk,
                    const char          *directory );
void cpucblk_zgetschur( const SolverCblk   *cblk,
                        int                 upper_part,
                        pastix_complex64_t *S,
                        pastix_int_t        lds );
void cpucblk_zdump( pastix_coefside_t  side,
                    const SolverCblk  *cblk,
                    FILE              *stream );
int cpucblk_zdiff( pastix_coefside_t  side,
                   const SolverCblk  *cblkA,
                   SolverCblk        *cblkB );
void cpucblk_zadd( pastix_coefside_t  side,
                   double             alpha,
                   const SolverCblk  *cblkA,
                   SolverCblk        *cblkB,
                   const pastix_lr_t *lowrank );

/**
 *    @}
 *    @name PastixComplex64 MPI routines
 *    @{
 */
int cpucblk_zincoming_deps( int                mt_flag,
                            pastix_coefside_t  side,
                            SolverMatrix      *solvmtx,
                            SolverCblk        *cblk );
void cpucblk_zrelease_deps( pastix_coefside_t  side,
                            SolverMatrix      *solvmtx,
                            const SolverCblk  *cblk,
                            SolverCblk        *fcbk );
void cpucblk_zrequest_cleanup( pastix_coefside_t  side,
                               pastix_int_t       sched,
                               SolverMatrix      *solvmtx );
void cpucblk_zupdate_reqtab( SolverMatrix *solvmtx );
#if defined( PASTIX_WITH_MPI )
void cpucblk_zmpi_progress( pastix_coefside_t  side,
                            SolverMatrix      *solvmtx,
                            int                threadid );
void cpucblk_zisend_rhs_bwd( SolverMatrix *solvmtx,
                             pastix_rhs_t  rhsb,
                             SolverCblk   *cblk );
#endif
void cpucblk_zmpi_rhs_fwd_progress( const args_solve_t *enums,
                                    SolverMatrix       *solvmtx,
                                    pastix_rhs_t        rhsb,
                                    int                 threadid );
void cpucblk_zrelease_rhs_fwd_deps( const args_solve_t *enums,
                                    SolverMatrix       *solvmtx,
                                    pastix_rhs_t        rhsb,
                                    const SolverCblk   *cblk,
                                    SolverCblk         *fcbk );
int cpucblk_zincoming_rhs_fwd_deps( int                 rank,
                                    const args_solve_t *enums,
                                    SolverMatrix       *solvmtx,
                                    SolverCblk         *cblk,
                                    pastix_rhs_t        rhsb );
void cpucblk_zrequest_rhs_fwd_cleanup( const args_solve_t *enums,
                                       pastix_int_t        sched,
                                       SolverMatrix       *solvmtx,
                                       pastix_rhs_t        rhsb );

void cpucblk_zmpi_rhs_bwd_progress( const args_solve_t *enums,
                                    SolverMatrix       *solvmtx,
                                    pastix_rhs_t        rhsb,
                                    int                 threadid );
void cpucblk_zrelease_rhs_bwd_deps( const args_solve_t *enums,
                                    SolverMatrix       *solvmtx,
                                    pastix_rhs_t        rhsb,
                                    const SolverCblk   *cblk,
                                    SolverCblk         *fcbk );
int cpucblk_zincoming_rhs_bwd_deps( int                 rank,
                                    const args_solve_t *enums,
                                    SolverMatrix       *solvmtx,
                                    SolverCblk         *cblk,
                                    pastix_rhs_t        rhsb );
void cpucblk_zrequest_rhs_bwd_cleanup( const args_solve_t *enums,
                                       pastix_int_t        sched,
                                       SolverMatrix       *solvmtx,
                                       pastix_rhs_t        rhsb );
void cpucblk_zsend_rhs_forward( const SolverMatrix *solvmtx,
                                SolverCblk         *cblk,
                                pastix_rhs_t        b );
void cpucblk_zrecv_rhs_forward( const SolverMatrix *solvmtx,
                                SolverCblk         *cblk,
                                pastix_complex64_t *work,
                                pastix_rhs_t        b );
void cpucblk_zsend_rhs_backward( const SolverMatrix *solvmtx,
                                 SolverCblk         *cblk,
                                 pastix_rhs_t        b );
void cpucblk_zrecv_rhs_backward( const SolverMatrix *solvmtx,
                                 SolverCblk         *cblk,
                                 pastix_rhs_t        b );

/**
 *    @}
 *    @name PastixComplex64 compression/uncompression routines
 *    @{
 */
pastix_fixdbl_t cpublok_zcompress( const pastix_lr_t *lowrank,
                                   pastix_int_t        M,
                                   pastix_int_t        N,
                                   pastix_lrblock_t   *blok );
pastix_int_t cpucblk_zcompress( const SolverMatrix *solvmtx,
                                pastix_coefside_t   side,
                                int                 max_ilulvl,
                                SolverCblk         *cblk );
void cpucblk_zuncompress( pastix_coefside_t  side,
                          SolverCblk        *cblk );
void cpucblk_zmemory( pastix_coefside_t   side,
                      const SolverMatrix *solvmtx,
                      SolverCblk         *cblk,
                      pastix_int_t       *orig,
                      pastix_int_t       *gain );

/**
 *    @}
 * @}
 *
 * @addtogroup kernel_solve
 * @{
 *    This module contains all the kernel working on the solver matrix structure
 *    for the solve step.
 *
 */

void solve_blok_ztrsm( pastix_side_t       side,
                       pastix_uplo_t       uplo,
                       pastix_trans_t      trans,
                       pastix_diag_t       diag,
                       const SolverCblk   *cblk,
                       int                 nrhs,
                       const void         *dataA,
                       pastix_complex64_t *b,
                       int                 ldb );
void solve_blok_zgemm( pastix_side_t             side,
                       pastix_trans_t            trans,
                       pastix_int_t              nrhs,
                       const SolverCblk         *cblk,
                       const SolverBlok         *blok,
                       SolverCblk               *fcbk,
                       const void               *dataA,
                       const pastix_complex64_t *B,
                       pastix_int_t              ldb,
                       pastix_complex64_t       *C,
                       pastix_int_t              ldc );

void solve_cblk_ztrsmsp_forward( const args_solve_t *enums,
                                 SolverMatrix       *datacode,
                                 const SolverCblk   *cblk,
                                 pastix_rhs_t        b );
void solve_cblk_ztrsmsp_backward( const args_solve_t *enums,
                                  SolverMatrix       *datacode,
                                  SolverCblk         *cblk,
                                  pastix_rhs_t        b );

void solve_cblk_zdiag( const SolverCblk   *cblk,
                       int                 nrhs,
                       pastix_complex64_t *b,
                       int                 ldb,
                       pastix_complex64_t *work );
/**
 * @}
 *
 * @addtogroup kernel_fact_null
 * @{
 *    This module contains the three terms update functions for the LDL^t and
 *    LDL^h factorizations.
 *
 */
#if defined(PRECISION_z) || defined(PRECISION_c)
void core_zhetrfsp1d_gemm( const SolverCblk         *cblk,
                           const SolverBlok         *blok,
                           SolverCblk               *fcblk,
                           const pastix_complex64_t *L,
                           pastix_complex64_t       *C,
                           pastix_complex64_t       *work );
#endif
void core_zsytrfsp1d_gemm( const SolverCblk         *cblk,
                           const SolverBlok         *blok,
                           SolverCblk               *fcblk,
                           const pastix_complex64_t *L,
                           pastix_complex64_t       *C,
                           pastix_complex64_t       *work );

int
cpucblk_zpotrfsp1dplus( SolverMatrix *solvmtx,
                        SolverCblk   *cblk );
void
cpucblk_zpotrfsp1dplus_update( SolverMatrix       *solvmtx,
                               SolverBlok         *blok,
                               pastix_complex64_t *work,
                               pastix_int_t        lwork );
int
cpucblk_zgetrfsp1dplus( SolverMatrix *solvmtx,
                        SolverCblk   *cblk );
void
cpucblk_zgetrfsp1dplus_update( SolverMatrix       *solvmtx,
                               SolverBlok         *blok,
                               pastix_complex64_t *work,
                               pastix_int_t        lwork );

/**
 * @}
 */

#endif /* _pastix_zcores_h_ */
