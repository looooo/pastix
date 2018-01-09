/**
 *
 * @file laplacian.h
 *
 * @copyright 2011-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2011-11-11
 *
 **/
#ifndef _laplacian_h_
#define _laplacian_h_

void z_spmLaplacian_7points( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3, pastix_fixdbl_t alpha, pastix_fixdbl_t beta );
void c_spmLaplacian_7points( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3, pastix_fixdbl_t alpha, pastix_fixdbl_t beta );
void d_spmLaplacian_7points( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3, pastix_fixdbl_t alpha, pastix_fixdbl_t beta );
void s_spmLaplacian_7points( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3, pastix_fixdbl_t alpha, pastix_fixdbl_t beta );
void p_spmLaplacian_7points( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3, pastix_fixdbl_t alpha, pastix_fixdbl_t beta );

void z_spmExtendedLaplacian2D( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2 );
void c_spmExtendedLaplacian2D( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2 );
void d_spmExtendedLaplacian2D( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2 );
void s_spmExtendedLaplacian2D( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2 );
void p_spmExtendedLaplacian2D( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2 );

void z_spmExtendedLaplacian3D( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );
void c_spmExtendedLaplacian3D( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );
void d_spmExtendedLaplacian3D( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );
void s_spmExtendedLaplacian3D( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );
void p_spmExtendedLaplacian3D( pastix_spm_t *spm, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );

int laplacian_parse_info( const char   *filename,
                          pastix_spm_t *spm,
                          pastix_int_t *dim1,
                          pastix_int_t *dim2,
                          pastix_int_t *dim3,
                          double       *alpha,
                          double       *beta );

#endif /* _laplacian_h_ */
