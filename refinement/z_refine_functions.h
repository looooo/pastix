/**
 *
 * @file z_refine_functions.h
 *
 * PaStiX refinement functions implementations.
 *
 * @copyright 2015-2023 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.3.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Gregoire Pichon
 * @author Theophile Terraz
 * @author Vincent Bridonneau
 * @date 2022-09-24
 * @precisions normal z -> c d s
 *
 * @addtogroup pastix_dev_refine
 * @{
 *   @brief TODO.
 *
 **/

#ifndef _z_refine_functions_h_
#define _z_refine_functions_h_

#include "common.h"
#include "bcsc/bcsc.h"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * @brief TODO
 */
struct z_solver
{
    pastix_int_t    (* getN   )   (pastix_data_t *);
    pastix_fixdbl_t (* getEps )   (pastix_data_t *);
    pastix_int_t    (* getImax)   (pastix_data_t *);
    pastix_int_t    (* getRestart)(pastix_data_t *);

    void* (*malloc)(size_t);
    void  (*free)(void*);

    void   (*output_oneiter)( pastix_fixdbl_t, pastix_fixdbl_t, double, pastix_int_t);
    void   (*output_final)( pastix_data_t *, pastix_complex64_t, pastix_int_t,
                            pastix_fixdbl_t, void*, pastix_complex64_t*);

    void   (*scal)( pastix_data_t *, pastix_int_t, pastix_complex64_t, pastix_complex64_t * );
    pastix_complex64_t (*dot) ( pastix_data_t *, pastix_int_t, const pastix_complex64_t *, const pastix_complex64_t * );
    void   (*copy)( pastix_data_t *, pastix_int_t, const pastix_complex64_t *, pastix_complex64_t * );
    void   (*axpy)( pastix_data_t *, pastix_int_t, pastix_complex64_t, const pastix_complex64_t *, pastix_complex64_t *);
    void   (*spmv)( const pastix_data_t *, pastix_trans_t, pastix_complex64_t, const pastix_complex64_t *, pastix_complex64_t, pastix_complex64_t * );
    void   (*spsv)( pastix_data_t *, pastix_complex64_t *, pastix_complex32_t * );
    double (*norm)( pastix_data_t *, pastix_int_t, const pastix_complex64_t * );
    void   (*gemv)( pastix_data_t *, pastix_int_t, pastix_int_t,
                    pastix_complex64_t, const pastix_complex64_t *, pastix_int_t,
                    const pastix_complex64_t *, pastix_complex64_t, pastix_complex64_t *);
};
#endif /* DOXYGEN_SHOULD_SKIP_THIS */

void z_refine_init(struct z_solver *, pastix_data_t*);

pastix_int_t z_gmres_smp   ( pastix_data_t *pastix_data, pastix_rhs_t xp, pastix_rhs_t bp );
pastix_int_t z_grad_smp    ( pastix_data_t *pastix_data, pastix_rhs_t xp, pastix_rhs_t bp );
pastix_int_t z_pivot_smp   ( pastix_data_t *pastix_data, pastix_rhs_t xp, pastix_rhs_t bp );
pastix_int_t z_bicgstab_smp( pastix_data_t *pastix_data, pastix_rhs_t xp, pastix_rhs_t bp );

/**
 * @}
 */
#endif /* _z_refine_functions_h_ */
