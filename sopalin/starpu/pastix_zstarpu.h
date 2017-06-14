/**
 *
 * @file pastix_zstarpu.c
 *
 * Pastix StarPU codelets header
 *
 * @copyright 2016-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @date 2013-06-24
 *
 * @precisions normal z -> z c d s
 *
 **/

void
starpu_task_cblk_zpotrfsp1d_panel( sopalin_data_t *sopalin_data,
                                   SolverCblk     *cblk );

void
starpu_task_cblk_zgetrfsp1d_panel( sopalin_data_t *sopalin_data,
                                   SolverCblk     *cblk );

void
starpu_task_cblk_zgemmsp( pastix_coefside_t sideA,
                          pastix_coefside_t sideB,
                          pastix_trans_t    trans,
                          const SolverCblk *cblk,
                          const SolverBlok *blok,
                          SolverCblk       *fcblk,
                          sopalin_data_t   *sopalin_data );

void
starpu_task_blok_zpotrf( sopalin_data_t *sopalin_data,
                         SolverCblk     *cblk );

void
starpu_task_blok_zgetrf( sopalin_data_t *sopalin_data,
                         SolverCblk     *cblk );

void
starpu_task_blok_ztrsmsp( pastix_coefside_t coef,
                          pastix_side_t     side,
                          pastix_uplo_t     uplo,
                          pastix_trans_t    trans,
                          pastix_diag_t     diag,
                          const SolverCblk *cblk,
                          SolverBlok       *blok,
                          sopalin_data_t   *sopalin_data );

void
starpu_task_blok_zgemmsp( pastix_coefside_t sideA,
                          pastix_coefside_t sideB,
                          pastix_trans_t    trans,
                          const SolverCblk *cblk,
                          SolverCblk       *fcblk,
                          const SolverBlok *blokA,
                          const SolverBlok *blokB,
                          sopalin_data_t   *sopalin_data );

void
starpu_zpotrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data );

void
starpu_zgetrf( pastix_data_t  *pastix_data,
               sopalin_data_t *sopalin_data );


