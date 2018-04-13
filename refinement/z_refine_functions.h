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

pastix_complex64_t
z_Pastix_Dotc( pastix_int_t n,
               const pastix_complex64_t *x,
               const pastix_complex64_t *y );

struct z_solver
{
    pastix_complex64_t* (* Synchro)(void *, void *, int);
    void* (* Malloc)(size_t);
    void  (* Free)(void*);

    void  (* Verbose)(double, double, double, pastix_int_t);
    void  (* End)(pastix_data_t *, pastix_complex64_t, pastix_int_t,
                 double, void*, pastix_complex64_t*);

    void (* X)(pastix_data_t *, void *, pastix_complex64_t *);
    pastix_int_t (* N)(pastix_data_t *);
    void (* B)(const pastix_complex64_t *, pastix_complex64_t *, pastix_int_t);
    pastix_complex64_t (* Eps)(pastix_data_t *);
    pastix_int_t (* Itermax)(pastix_data_t *);
    pastix_int_t (* Krylov_Space)(pastix_data_t *);
    pastix_int_t (* me)(void *);

    double (* Norm)(pastix_int_t, const pastix_complex64_t *);
    void (* Precond)(pastix_data_t *, pastix_complex64_t *);

    void (* Scal)(pastix_int_t, pastix_complex64_t, pastix_complex64_t *);
    void (* Ax)(pastix_bcsc_t *, pastix_complex64_t *, pastix_complex64_t *);

    void (* bMAx)(pastix_bcsc_t *, const pastix_complex64_t *, const pastix_complex64_t *, pastix_complex64_t *);
    void (* BYPX)(pastix_int_t, pastix_complex64_t *, pastix_complex64_t *, pastix_complex64_t *);
    void (* AXPY)(pastix_int_t, pastix_complex64_t,  const pastix_complex64_t *, pastix_complex64_t *);

    void   (*output_oneiter)(double, double, double, pastix_int_t);
    void   (*scal)( pastix_int_t, pastix_complex64_t, pastix_complex64_t * );
    pastix_complex64_t (*dot) ( pastix_int_t, const pastix_complex64_t *, const pastix_complex64_t * );
    void   (*copy)( pastix_int_t, const pastix_complex64_t *, pastix_complex64_t * );
    void   (*axpy)( pastix_int_t, pastix_complex64_t, const pastix_complex64_t *, pastix_complex64_t *);
    void   (*spmv)( pastix_data_t *, pastix_complex64_t, const pastix_complex64_t *, pastix_complex64_t, pastix_complex64_t * );
    void   (*spsv)( pastix_data_t *, pastix_complex64_t * );
    double (*norm)( pastix_int_t, const pastix_complex64_t * );
    void   (*trsv)( pastix_data_t *, pastix_complex64_t * );
    void   (*gemv)( pastix_data_t *, pastix_complex64_t * );
};

void z_Pastix_Solveur(struct z_solver *);


void z_gmres_smp   ( pastix_data_t *pastix_data, void *x, void *b );
void z_grad_smp    ( pastix_data_t *pastix_data, void *x, void *b );
void z_pivot_smp   ( pastix_data_t *pastix_data, void *x, void *b );
void z_bicgstab_smp( pastix_data_t *pastix_data, void *x, void *b );

#endif /* _z_refine_functions_h_ */
