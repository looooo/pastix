/**
 *
 * @file z_refine_functions.h
 *
 * PaStiX refinement functions implementations.
 *
 * @copyright 2015-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/

#ifndef _z_refine_functions_h_
#define _z_refine_functions_h_

void *
z_Pastix_Malloc( size_t size );

void
z_Pastix_Free( void *x );


void
z_Pastix_Verbose( double t0, double tf, double err, pastix_int_t nb_iters );

void
z_Pastix_End( pastix_data_t *pastix_data, pastix_complex64_t err,
              pastix_int_t nb_iters, double tf,
              void *x, pastix_complex64_t *gmresx );

void
z_Pastix_X( pastix_data_t *pastix_data, void *x, pastix_complex64_t *gmresx );

pastix_int_t
z_Pastix_n( pastix_data_t *pastix_data );

void
z_Pastix_B( void *b, pastix_complex64_t *refineb, pastix_int_t n );

pastix_complex64_t
z_Pastix_Eps( pastix_data_t *pastix_data );

pastix_int_t
z_Pastix_Itermax( pastix_data_t *pastix_data );

pastix_int_t
z_Pastix_Krylov_Space( pastix_data_t *pastix_data );

pastix_complex64_t
z_Pastix_Norm2( pastix_complex64_t *x, pastix_int_t n );

void
z_Pastix_Precond( pastix_data_t *pastix_data, pastix_complex64_t *s, pastix_complex64_t *d );

void
z_Pastix_Scal( pastix_int_t n, pastix_complex64_t alpha, pastix_complex64_t *x );

#if defined(PRECISION_z) || defined(PRECISION_c)
void
z_Pastix_Dotc( pastix_int_t n, pastix_complex64_t *x,
               pastix_complex64_t *y, pastix_complex64_t *r );
#endif

void
z_Pastix_Dotu( pastix_int_t n, pastix_complex64_t *x,
               pastix_complex64_t *y, pastix_complex64_t *r );

void
z_Pastix_Ax( pastix_bcsc_t *bcsc, pastix_complex64_t *x, pastix_complex64_t *r );

void
z_Pastix_bMAx( pastix_bcsc_t *bcsc, pastix_complex64_t *b,
               pastix_complex64_t *x, pastix_complex64_t *r );

void
z_Pastix_BYPX( pastix_int_t n, pastix_complex64_t *beta,
               pastix_complex64_t *y, pastix_complex64_t *x );

void
z_Pastix_AXPY( pastix_int_t n, double coeff, pastix_complex64_t *alpha,
               pastix_complex64_t *x, pastix_complex64_t *y );

pastix_int_t
z_Pastix_me( void *arg );

struct z_solver
{
    pastix_complex64_t* (* Synchro)(void *, void *, int);
    void* (* Malloc)(size_t);
    void (* Free)(void*);

    void (* Verbose)(double, double, double, pastix_int_t);
    void (* End)(pastix_data_t *, pastix_complex64_t, pastix_int_t,
                 double, void*, pastix_complex64_t*);
    void (* X)(pastix_data_t *, void *, pastix_complex64_t *);
    pastix_int_t (* N)(pastix_data_t *);
    void (* B)(void *, pastix_complex64_t *, pastix_int_t);
    pastix_complex64_t (* Eps)(pastix_data_t *);
    pastix_int_t (* Itermax)(pastix_data_t *);
    pastix_int_t (* Krylov_Space)(pastix_data_t *);
    pastix_int_t (* me)(void *);

    pastix_complex64_t (* Norm)(pastix_complex64_t *, pastix_int_t);
    void (* Precond)(pastix_data_t *, pastix_complex64_t *, pastix_complex64_t *);

    void (* Scal)(pastix_int_t, pastix_complex64_t, pastix_complex64_t *);
    void (* Dotc)(pastix_int_t, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *);
    void (* Ax)(pastix_bcsc_t *, pastix_complex64_t *, pastix_complex64_t *);

    void (* bMAx)(pastix_bcsc_t *, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *);
    void (* BYPX)(pastix_int_t, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *);
    void (* AXPY)(pastix_int_t, double, pastix_complex64_t *,
                  pastix_complex64_t *, pastix_complex64_t *);
};

void z_Pastix_Solveur(struct z_solver *);


void z_gmres_smp   ( pastix_data_t *pastix_data, void *x, void *b );
void z_grad_smp    ( pastix_data_t *pastix_data, void *x, void *b );
void z_pivot_smp   ( pastix_data_t *pastix_data, void *x, void *b );
void z_bicgstab_smp( pastix_data_t *pastix_data, void *x, void *b );

#endif /* _z_refine_functions_h_ */
