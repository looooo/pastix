/**
 *
 * @file starpu_zgetrfsp_cpu.c
 *
 *  PaStiX kernel routines for StarPU
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2011-11-11
 * @precisions normal z -> c d s
 *
 **/
#include <assert.h>
#include "starpu_zdefines.h"
#include "starpu_zsubmit.h"

#include "common.h"
#include "z_sopalin3d.h"
#include "z_solver.h"
#include "pastix_zcores.h"
#include "sopalin_acces.h"

static pastix_complex64_t zone  =  1.;
static pastix_complex64_t mzone = -1.;
static pastix_complex64_t zzero =  0.;


/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zgetrfsp1d_getrf_cpu - Computes the LU factorization of one panel.
 *
 *******************************************************************************
 *
 * @param[in, out] buffers
 *          Array of 2 buffers:
 *            -# L : The pointer to the lower matrix storing the coefficients
 *                   of the panel. Must be of size cblk.stride -by- cblk.width
 *            -# U : The pointer to the upper matrix storing the coefficients
 *                   of the panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in, out] _args
 *          Package of arguments to unpack with starpu_codelet_unpack_args():
 *            -# sopalin_data : common sopalin structure.
 *            -# cblk : Pointer to the structure representing the panel to
 *                      factorize in the cblktab array. Next column blok must
 *                      be accessible through cblk[1].
 *
 *
 *******************************************************************************/
void starpu_zgetrfsp1d_getrf_cpu(void * buffers[], void * _args)
{
    z_Sopalin_Data_t     *sopalin_data;
    z_SolverCblk         *cblk;
    pastix_complex64_t *L      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *U      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
    pastix_int_t        stride = STARPU_MATRIX_GET_LD(buffers[0]);
    int                 me     = starpu_worker_get_id();

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk);

    sopalin_data->thread_data[me]->nbpivot +=
        core_zgetrfsp1d_getrf(cblk, L, U,
                              sopalin_data->critere);
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zgetrfsp1d_trsm_cpu - Apply all the trsm updates on one panel.
 *
 *******************************************************************************
 *
 * @param[in, out] buffers
 *          Array of 2 buffers:
 *            -# L : The pointer to the lower matrix storing the coefficients
 *                   of the panel. Must be of size cblk.stride -by- cblk.width
 *            -# U : The pointer to the upper matrix storing the coefficients
 *                   of the panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in, out] _args
 *          Package of arguments to unpack with starpu_codelet_unpack_args():
 *            -# sopalin_data : common sopalin structure.
 *            -# cblk : Pointer to the structure representing the panel to
 *                      factorize in the cblktab array. Next column blok must
 *                      be accessible through cblk[1].
 *
 *
 *******************************************************************************/
void starpu_zgetrfsp1d_trsm_cpu(void * buffers[], void * _args)
{
    z_Sopalin_Data_t     *sopalin_data;
    z_SolverCblk         *cblk;
    pastix_complex64_t *L      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *U      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
    pastix_int_t        stride = STARPU_MATRIX_GET_LD(buffers[0]);

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk);

    core_zgetrfsp1d_trsm(cblk,
                         L,
                         U);

    SUBMIT_GEMMS_IF_NEEDED;
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zgetrfsp1d_cpu - Computes the LU factorization of one panel and apply
 * all the trsm updates to this panel.
 *
 *******************************************************************************
 *
 * @param[in, out] buffers
 *          Array of 2 buffers:
 *            -# L : The pointer to the lower matrix storing the coefficients
 *                   of the panel. Must be of size cblk.stride -by- cblk.width
 *            -# U : The pointer to the upper matrix storing the coefficients
 *                   of the panel. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in, out] _args
 *          Package of arguments to unpack with starpu_codelet_unpack_args():
 *            -# sopalin_data : common sopalin structure.
 *            -# cblk : Pointer to the structure representing the panel to
 *                      factorize in the cblktab array. Next column blok must
 *                      be accessible through cblk[1].
 *
 *
 *******************************************************************************/
void starpu_zgetrfsp1d_cpu(void * buffers[], void * _args)
{
    z_Sopalin_Data_t     *sopalin_data;
    z_SolverCblk         *cblk;
    pastix_complex64_t *L      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *U      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
    pastix_int_t        stride = STARPU_MATRIX_GET_LD(buffers[0]);
    int                 me     = starpu_worker_get_id();

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk);

    sopalin_data->thread_data[me]->nbpivot +=
        core_zgetrfsp1d(cblk,
                        L,
                        U,
                        sopalin_data->critere);

    SUBMIT_GEMMS_IF_NEEDED;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zgetrfsp1d_gemm_cpu - Computes the Cholesky factorization of one panel
 * and apply all the trsm updates to this panel.
 *
 *******************************************************************************
 *
 * @param[in, out] buffers
 *          Array of 5 buffers:
 *            -# L : The pointer to the lower matrix storing the coefficients
 *                   of the panel. Must be of size cblk.stride -by- cblk.width
 *            -# Cl : The pointer to the lower matrix storing the coefficients
 *                    of the updated panel. Must be of size cblk.stride -by-
 *                    cblk.width
 *            -# U : The pointer to the upper matrix storing the coefficients
 *                   of the panel. Must be of size cblk.stride -by- cblk.width
 *            -# Cu : The pointer to the upper matrix storing the coefficients
 *                    of the updated panel. Must be of size cblk.stride -by-
 *                    cblk.width
 *            -# work : Scratch buffer of size SOLV_COEFMAX.
 *
 * @param[in, out] _args
 *          Package of arguments to unpack with starpu_codelet_unpack_args():
 *            -# sopalin_data : common sopalin structure.
 *            -# cblk : Pointer to the structure representing the panel to
 *                      factorize in the cblktab array. Next column blok must
 *                      be accessible through cblk[1].
 *            -# blok : The pointer to the data structure that describes the
 *                      blok from which we compute the contributions.
 *            -# fcblk : The pointer to the data structure that describes the
 *                       panel on which we compute the contributions. Next
 *                       column blok must be accessible through fcblk[1].
 *
 *
 *******************************************************************************/

void starpu_zgetrfsp1d_gemm_cpu(void * buffers[], void * _args)
{
    z_Sopalin_Data_t     *sopalin_data;
    z_SolverCblk         *cblk;
    z_SolverBlok         *blok;
    z_SolverCblk         *fcblk;
    pastix_int_t        stride   = STARPU_MATRIX_GET_LD(buffers[0]);
    pastix_complex64_t *L    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *Cl   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
    pastix_complex64_t *U    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[2]);
    pastix_complex64_t *Cu   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[3]);
    pastix_complex64_t *work = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[4]);

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk, &blok, &fcblk);

    core_zgetrfsp1d_gemm(cblk,
                         blok,
                         fcblk,
                         L, U, Cl, Cu,
                         work);
    SUBMIT_TRF_IF_NEEDED;
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zgetrfsp1d_geadd_cpu - Computes the addition of a fanin column block
 * into the destination column block.
 *
 *******************************************************************************
 *
 * @param[in, out] buffers
 *          Array of 4 buffers:
 *            -# L : The pointer to the lower matrix storing the coefficients
 *                   of the panel. Must be of size cblk1.stride -by- cblk1.width
 *            -# Cl : The pointer to the lower matrix storing the coefficients
 *                    of the updated panel. Must be of size cblk2.stride -by-
 *                    cblk2.width
 *            -# U : The pointer to the upper matrix storing the coefficients
 *                   of the panel. Must be of size cblk1.stride -by- cblk1.width
 *            -# Cu : The pointer to the upper matrix storing the coefficients
 *                    of the updated panel. Must be of size cblk2.stride -by-
 *                    cblk2.width
 *
 * @param[in, out] _args
 *          Package of arguments to unpack with starpu_codelet_unpack_args():
 *            -# sopalin_data : common sopalin structure.
 *            -# cblk1 : Pointer to the structure representing the panel to
 *                       factorize in the cblktab array. Next column blok must
 *                       be accessible through cblk1[1].
 *            -# cblk2 : Pointer to the structure representing the panel to
 *                       factorize in the cblktab array. Next column blok must
 *                       be accessible through cblk2[1].
 *
 *******************************************************************************/
void
starpu_zgetrfsp1d_geadd_cpu( void * buffers[],
                             void * _args )
{
    z_Sopalin_Data_t     *sopalin_data;
    z_SolverCblk *cblk1, *cblk2;
    pastix_complex64_t *L    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *Cl   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
    pastix_complex64_t *U    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[2]);
    pastix_complex64_t *Cu   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[3]);
    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk1, &cblk2);

    core_zgeaddsp1d(cblk1,
                    cblk2,
                    L, Cl,
                    U, Cu);
}
