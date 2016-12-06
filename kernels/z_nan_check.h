/**
 * @file z_nan_check.h
 *
 * Copyright (c) 2016      Inria. All rights reserved.
 *
 * @precisions normal z -> s d c
 */
#ifndef _Z_NAN_CHECK_H_
#define _Z_NAN_CHECK_H_

#define PASTIX_DEBUG_LR

//#define PASTIX_LR_CHECKNAN
#if defined(PASTIX_LR_CHECKNAN)
#define LAPACKE_zlacpy_work LAPACKE_zlacpy
#define LAPACKE_zlaset_work LAPACKE_zlaset

#define LAPACKE_zunmlq_work( _layout_, _side_, _trans_, _m_, _n_, _k_, _a_, _lda_, _tau_, _c_, _ldc_, _w_, _ldw_ ) \
    LAPACKE_zunmlq( _layout_, _side_, _trans_, _m_, _n_, _k_, _a_, _lda_, _tau_, _c_, _ldc_ )
#define LAPACKE_zunmqr_work( _layout_, _side_, _trans_, _m_, _n_, _k_, _a_, _lda_, _tau_, _c_, _ldc_, _w_, _ldw_ ) \
    LAPACKE_zunmqr( _layout_, _side_, _trans_, _m_, _n_, _k_, _a_, _lda_, _tau_, _c_, _ldc_ )

#define LAPACKE_zgeqrf_work( _layout_, _m_, _n_, _a_, _lda_, _tau_, _w_, _ldw_ ) \
    LAPACKE_zgeqrf( _layout_, _m_, _n_, _a_, _lda_, _tau_ )
#define LAPACKE_zgelqf_work( _layout_, _m_, _n_, _a_, _lda_, _tau_, _w_, _ldw_ ) \
    LAPACKE_zgelqf( _layout_, _m_, _n_, _a_, _lda_, _tau_ )

#if defined(PRECISION_z) || defined(PRECISION_c)
#define MYLAPACKE_zgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ ) \
    LAPACKE_zgesvd( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, (double*)(_w_) )
#else
#define MYLAPACKE_zgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ ) \
    LAPACKE_zgesvd( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, (double*)(_w_) )
#endif

#else

#if defined(PRECISION_z) || defined(PRECISION_c)
#define MYLAPACKE_zgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ ) \
    LAPACKE_zgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ )
#else
#define MYLAPACKE_zgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_, _rw_ ) \
    LAPACKE_zgesvd_work( _layout_, _jobu_, jobv_, _m_, _n_, _a_, _lda_, _s_, _u_, _ldu_, _v_, _ldv_, _w_, _ldw_ )
#endif

#endif /* defined(PASTIX_LR_CHECKNAN) */

#endif /* _Z_NAN_CHECK_H_ */
