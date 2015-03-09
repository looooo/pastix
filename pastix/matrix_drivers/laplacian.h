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

void z_spmLaplacian1D( pastix_csc_t *csc, pastix_int_t dim1 );
void c_spmLaplacian1D( pastix_csc_t *csc, pastix_int_t dim1 );
void d_spmLaplacian1D( pastix_csc_t *csc, pastix_int_t dim1 );
void s_spmLaplacian1D( pastix_csc_t *csc, pastix_int_t dim1 );
void p_spmLaplacian1D( pastix_csc_t *csc, pastix_int_t dim1 );

void z_spmLaplacian2D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2 );
void c_spmLaplacian2D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2 );
void d_spmLaplacian2D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2 );
void s_spmLaplacian2D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2 );
void p_spmLaplacian2D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2 );

void z_spmLaplacian3D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );
void c_spmLaplacian3D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );
void d_spmLaplacian3D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );
void s_spmLaplacian3D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );
void p_spmLaplacian3D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );

void z_spmExtendedLaplacian2D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2 );
void c_spmExtendedLaplacian2D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2 );
void d_spmExtendedLaplacian2D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2 );
void s_spmExtendedLaplacian2D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2 );
void p_spmExtendedLaplacian2D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2 );

void z_spmExtendedLaplacian3D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );
void c_spmExtendedLaplacian3D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );
void d_spmExtendedLaplacian3D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );
void s_spmExtendedLaplacian3D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );
void p_spmExtendedLaplacian3D( pastix_csc_t *csc, pastix_int_t dim1, pastix_int_t dim2, pastix_int_t dim3 );

#endif /* _LAPLACIAN_H_ */
