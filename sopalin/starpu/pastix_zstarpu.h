/**
 *
 * @file pastix_zstarpu.h
 *
 * Pastix StarPU codelets header
 *
 * @copyright 2016-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Vincent Bridonneau
 * @date 2022-10-17
 *
 * @precisions normal z -> z c d s
 *
 **/
#ifndef _pastix_zstarpu_h_
#define _pastix_zstarpu_h_

void starpu_task_cblk_zgetrfsp( sopalin_data_t *sopalin_data,
                                SolverCblk     *cblk,
                                int             prio );
void starpu_task_cblk_zhetrfsp( sopalin_data_t *sopalin_data,
                                SolverCblk     *cblk,
                                int             prio );
void starpu_task_cblk_zpotrfsp( sopalin_data_t *sopalin_data,
                                SolverCblk     *cblk,
                                int             prio );
void starpu_task_cblk_zpxtrfsp( sopalin_data_t *sopalin_data,
                                SolverCblk     *cblk,
                                int             prio );
void starpu_task_cblk_zsytrfsp( sopalin_data_t *sopalin_data,
                                SolverCblk     *cblk,
                                int             prio );

void starpu_task_blok_zgetrf( sopalin_data_t *sopalin_data,
                              SolverCblk     *cblk,
                              int             prio );
void starpu_task_blok_zhetrf( sopalin_data_t *sopalin_data,
                              SolverCblk     *cblk,
                              int             prio );
void starpu_task_blok_zpotrf( sopalin_data_t *sopalin_data,
                              SolverCblk     *cblk,
                              int             prio );
void starpu_task_blok_zpxtrf( sopalin_data_t *sopalin_data,
                              SolverCblk     *cblk,
                              int             prio );
void starpu_task_blok_zsytrf( sopalin_data_t *sopalin_data,
                              SolverCblk     *cblk,
                              int             prio );
void starpu_stask_cblk_zdiag( sopalin_data_t *sopalin_data,
                              SolverCblk     *cblk,
                              int             prio );
void starpu_stask_blok_zgemm( sopalin_data_t   *sopalin_data,
                              pastix_coefside_t coef,
                              pastix_side_t     side,
                              pastix_trans_t    trans,
                              const SolverCblk *cblk,
                              const SolverBlok *blok,
                              SolverCblk       *fcbk,
                              pastix_int_t      prio );
void starpu_stask_blok_ztrsm( sopalin_data_t   *sopalin_data,
                              pastix_coefside_t coef,
                              pastix_side_t     side,
                              pastix_uplo_t     uplo,
                              pastix_trans_t    trans,
                              pastix_diag_t     diag,
                              const SolverCblk *cblk,
                              pastix_int_t      prio );

void starpu_task_cblk_zgemmsp( sopalin_data_t   *sopalin_data,
                               pastix_coefside_t sideA,
                               pastix_coefside_t sideB,
                               pastix_trans_t    trans,
                               const SolverCblk *cblk,
                               const SolverBlok *blok,
                               SolverCblk       *fcblk,
                               int               prio );
void starpu_task_blok_zgemmsp( sopalin_data_t   *sopalin_data,
                               pastix_coefside_t sideA,
                               pastix_coefside_t sideB,
                               pastix_trans_t    trans,
                               const SolverCblk *cblk,
                               SolverCblk       *fcblk,
                               const SolverBlok *blokA,
                               const SolverBlok *blokB,
                               int               prio );

void starpu_task_blok_ztrsmsp( sopalin_data_t   *sopalin_data,
                               pastix_coefside_t coef,
                               pastix_side_t     side,
                               pastix_uplo_t     uplo,
                               pastix_trans_t    trans,
                               pastix_diag_t     diag,
                               const SolverCblk *cblk,
                               SolverBlok       *blok,
                               int               prio );

void starpu_task_blok_zscalo( sopalin_data_t   *sopalin_data,
                              pastix_trans_t  trans,
                              SolverCblk     *cblk,
                              SolverBlok     *blok,
                              int             prio );

void starpu_zdiag ( pastix_data_t      *pastix_data,
                    sopalin_data_t     *sopalin_data,
                    int                 nrhs,
                    pastix_complex64_t *b,
                    int                 ldb );
void starpu_zpotrf( pastix_data_t  *pastix_data,
                    sopalin_data_t *sopalin_data );
void starpu_zpxtrf( pastix_data_t  *pastix_data,
                    sopalin_data_t *sopalin_data );
void starpu_zgetrf( pastix_data_t  *pastix_data,
                    sopalin_data_t *sopalin_data );
void starpu_zhetrf( pastix_data_t  *pastix_data,
                    sopalin_data_t *sopalin_data );
void starpu_zsytrf( pastix_data_t  *pastix_data,
                    sopalin_data_t *sopalin_data );
void starpu_ztrsm ( pastix_data_t  *pastix_data,
                    enums_trsm_t   *enums,
                    sopalin_data_t *sopalin_data,
                    pastix_rhs_t    b );

#endif /* _pastix_zstarpu_h_ */
