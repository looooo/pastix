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
 * starpu_zsytrfsp1d_sytrf_cpu - Computes the LDLt factorization of one panel.
 *
 *******************************************************************************
 *
 * @param[in, out] buffers
 *          Array of 2 buffers:
 *            -# L : The pointer to the lower matrix storing the coefficients
 *                   of the panel. Must be of size cblk.stride -by- cblk.width
 *            -# work : Scratch buffers of size SOLV_COEFMAX.
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
void starpu_zsytrfsp1d_sytrf_cpu(void * buffers[], void * _args)
{
    Sopalin_Data_t     *sopalin_data;
    z_SolverCblk         *cblk;
    pastix_complex64_t *L      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *work   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_int_t        stride = STARPU_MATRIX_GET_LD(buffers[0]);
    int                 me     = starpu_worker_get_id();

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk);

    sopalin_data->thread_data[me]->nbpivot +=
        core_zsytrfsp1d_sytrf(cblk, L,
                              sopalin_data->critere,
                              work);
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zsytrfsp1d_trsm_cpu - Apply all the trsm updates on one panel.
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
void starpu_zsytrfsp1d_trsm_cpu(void * buffers[], void * _args)
{
    Sopalin_Data_t     *sopalin_data;
    z_SolverCblk         *cblk;
    pastix_complex64_t *L      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_int_t        stride = STARPU_MATRIX_GET_LD(buffers[0]);

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk);

    core_zsytrfsp1d_trsm(cblk, L);

    SUBMIT_GEMMS_IF_NEEDED;
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zsytrfsp1d_cpu - Computes the LU factorization of one panel and apply
 * all the trsm updates to this panel.
 *
 *******************************************************************************
 *
 * @param[in, out] buffers
 *          Array of 2 buffers:
 *            -# L : The pointer to the lower matrix storing the coefficients
 *                   of the panel. Must be of size cblk.stride -by- cblk.width
 *            -# work : Scratch buffers of size SOLV_COEFMAX.
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
void starpu_zsytrfsp1d_cpu(void * buffers[], void * _args)
{
    Sopalin_Data_t     *sopalin_data;
    z_SolverCblk         *cblk;
    pastix_complex64_t *L      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *work   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_int_t        stride = STARPU_MATRIX_GET_LD(buffers[0]);
    int                 me     = starpu_worker_get_id();

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk);

    sopalin_data->thread_data[me]->nbpivot +=
        core_zsytrfsp1d(cblk,
                        L,
                        sopalin_data->critere,
                        work);

    SUBMIT_GEMMS_IF_NEEDED;
#ifdef PASTIX_DUMP_CBLK
    {
        char name[256];
        sprintf(name, "cblk_%d_after_trf_trsm", cblk->gcblknum);
        z_cblk_save(cblk, name, L);
    }
#endif
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zsytrfsp1d_gemm_cpu - Computes the Cholesky factorization of one panel
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
 *            -# work + work2 : Scratch buffers of size SOLV_COEFMAX.
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

void starpu_zsytrfsp1d_gemm_cpu(void * buffers[], void * _args)
{
    Sopalin_Data_t     *sopalin_data;
    z_SolverCblk         *cblk;
    z_SolverBlok         *blok;
    z_SolverCblk         *fcblk;
    pastix_complex64_t *L    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *Cl   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
    pastix_complex64_t *work = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[2]);
    pastix_int_t        ldw  = (pastix_int_t)STARPU_VECTOR_GET_NX(buffers[2])/2;
    pastix_complex64_t *work2= work + ldw;

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk, &blok, &fcblk);
    assert(z_is_block_inside_fblock(blok, fcblk->fblokptr));
    assert(cblk->stride == STARPU_MATRIX_GET_LD(buffers[0]));
    assert(z_cblk_colnbr(cblk)  == STARPU_MATRIX_GET_NY(buffers[0]));
    core_zsytrfsp1d_gemm(cblk,
                         blok,
                         fcblk,
                         L, Cl,
                         work,
                         work2);
    SUBMIT_TRF_IF_NEEDED;
#ifdef PASTIX_DUMP_CBLK
    {
        char name[256];
        sprintf(name, "cblk_%d_after_gemm_%d_%d_%d_on_%d", fcblk->gcblknum,
                cblk->gcblknum, blok - cblk->fblokptr, fcblk->gcblknum,
                sopalin_data->datacode->clustnum);
        z_cblk_save(fcblk, name, Cl);
    }
#endif
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zsytrfsp1d_geadd_cpu - Computes the addition of a fanin column block
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
starpu_zsytrfsp1d_syadd_cpu(void * buffers[], void * _args) {
    Sopalin_Data_t *sopalin_data;
    z_SolverCblk *cblk1, *cblk2;
    pastix_complex64_t *L    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *Cl   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk1, &cblk2);
#ifdef PASTIX_DUMP_CBLK
    {
        char name[256];
        sprintf(name, "cblk_dst_%d_before_add_from_%d", cblk2->gcblknum,
                fcblk_getorigin(sopalin_data->datacode, cblk1));
        z_cblk_save(cblk2, name, Cl);
    }
    {
        char name[256];
        sprintf(name, "cblk_src_%d_before_add_from_%d", cblk1->gcblknum,
                fcblk_getorigin(sopalin_data->datacode, cblk1));
        z_cblk_save(cblk1, name, L);
    }
#endif
    core_zgeaddsp1d(cblk1,
                    cblk2,
                    L, Cl,
                    NULL, NULL);
#ifdef PASTIX_DUMP_CBLK
    {
        char name[256];
        sprintf(name, "cblk_%d_after_add_from_%d", cblk2->gcblknum,
                fcblk_getorigin(sopalin_data->datacode, cblk1));
        z_cblk_save(cblk2, name, Cl);
    }
#endif
}
