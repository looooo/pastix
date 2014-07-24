/**
 *
 * @file starpu_zpotrfsp_cpu.c
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
#include <starpu.h>

#include "common.h"
#include "sopalin3d.h"
#include "z_solver.h"
#include "pastix_zcores.h"
#include "sopalin_acces.h"
#include "starpu_zdefines.h"
#include "starpu_zsubmit.h"

static pastix_complex64_t zone  =  1.;
static pastix_complex64_t mzone = -1.;
static pastix_complex64_t zzero =  0.;


/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zpotrfsp1d_potrf_cpu - Computes the LLt factorization of one panel.
 *
 *******************************************************************************
 *
 * @param[in, out] buffers
 *          Array of 2 buffers:
 *            -# L : The pointer to the lower matrix storing the coefficients
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
void starpu_zpotrfsp1d_potrf_cpu(void * buffers[], void * _args)
{
    Sopalin_Data_t     *sopalin_data;
    d_SolverCblk         *cblk;
    pastix_complex64_t *L      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_int_t        stride = STARPU_MATRIX_GET_LD(buffers[0]);
    int                 me     = starpu_worker_get_id();

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk);

    sopalin_data->thread_data[me]->nbpivot +=
        core_zpotrfsp1d_potrf(cblk, L,
                              sopalin_data->critere);
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zpotrfsp1d_trsm_cpu - Apply all the trsm updates on one panel.
 *
 *******************************************************************************
 *
 * @param[in, out] buffers
 *          Array of 2 buffers:
 *            -# L : The pointer to the lower matrix storing the coefficients
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
void starpu_zpotrfsp1d_trsm_cpu(void * buffers[], void * _args)
{
    Sopalin_Data_t     *sopalin_data;
    d_SolverCblk         *cblk;
    pastix_complex64_t *L      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_int_t        stride = STARPU_MATRIX_GET_LD(buffers[0]);

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk);

    core_zpotrfsp1d_trsm(cblk, L);

    SUBMIT_GEMMS_IF_NEEDED;
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zpotrfsp1d_cpu - Computes the LU factorization of one panel and apply
 * all the trsm updates to this panel.
 *
 *******************************************************************************
 *
 * @param[in, out] buffers
 *          Array of 2 buffers:
 *            -# L : The pointer to the lower matrix storing the coefficients
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
void starpu_zpotrfsp1d_cpu(void * buffers[], void * _args)
{
    Sopalin_Data_t     *sopalin_data;
    d_SolverCblk         *cblk;
    pastix_complex64_t *L      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_int_t        stride = STARPU_MATRIX_GET_LD(buffers[0]);
    int                 me     = starpu_worker_get_id();

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk);

    sopalin_data->thread_data[me]->nbpivot +=
        core_zpotrfsp1d(cblk,
                        L,
                        sopalin_data->critere);

    SUBMIT_GEMMS_IF_NEEDED;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zpotrfsp1d_gemm_cpu - Computes the Cholesky factorization of one panel
 * and apply all the trsm updates to this panel.
 *
 *******************************************************************************
 *
 * @param[in, out] buffers
 *          Array of 2 buffers:
 *            -# L : The pointer to the lower matrix storing the coefficients
 *                   of the panel. Must be of size cblk.stride -by- cblk.width
 *            -# Cl : The pointer to the lower matrix storing the coefficients
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
void starpu_zpotrfsp1d_gemm_cpu(void * buffers[], void * _args)
{
    Sopalin_Data_t     *sopalin_data;
    d_SolverCblk         *cblk;
    d_SolverBlok         *blok;
    d_SolverCblk         *fcblk;
    pastix_int_t        stride   = STARPU_MATRIX_GET_LD(buffers[0]);
    pastix_complex64_t *L    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *Cl   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
    pastix_complex64_t *work = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[2]);

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk, &blok, &fcblk);

    core_zpotrfsp1d_gemm(cblk,
                         blok,
                         fcblk,
                         L, Cl,
                         work);
    SUBMIT_TRF_IF_NEEDED;
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zpotrfsp1d_geadd_cpu - Computes the addition of a fanin column block
 * into the destination column block.
 *
 *******************************************************************************
 *
 * @param[in, out] buffers
 *          Array of 2 buffers:
 *            -# L : The pointer to the lower matrix storing the coefficients
 *                   of the panel. Must be of size cblk1.stride -by- cblk1.width
 *            -# Cl : The pointer to the lower matrix storing the coefficients
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
starpu_zpotrfsp1d_geadd_cpu(void * buffers[], void * _args) {
    Sopalin_Data_t * sopalin_data;
    d_SolverCblk *cblk1, *cblk2;
    pastix_complex64_t *L    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *Cl   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk1, &cblk2);

    core_zgeaddsp1d(cblk1,
                    cblk2,
                    L, Cl,
                    NULL, NULL);
}
