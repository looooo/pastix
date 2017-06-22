/**
 * @file laplacian.h
 *
 *  $COPYRIGHTS$
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2011-11-11
 *
 **/
#ifndef _LAPLACIAN_H_
#define _LAPLACIAN_H_

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

#endif /* _LAPLACIAN_H_ */
