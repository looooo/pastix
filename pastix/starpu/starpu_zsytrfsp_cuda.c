/**
 *
 * @file starpu_zpotrfsp_cuda.c
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
#include "z_sopalin3d.h"
#include "z_solver.h"
#include "pastix_zcores.h"
#include "sopalin_acces.h"
#include "starpu_defines.h"
#include "starpu_zsubmit.h"
#include "pastix_cuda_helper.h"
#include "sparse_zgemm_fermi.h"
static pastix_complex64_t zone  =  1.;
static pastix_complex64_t mzone = -1.;
static pastix_complex64_t zzero =  0.;


/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zsytrfsp1d_sytrf_cuda - Computes the LDLt factorization of one panel.
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
void starpu_zsytrfsp1d_sytrf_cuda(void * buffers[], void * _args)
{
    z_Sopalin_Data_t     *sopalin_data;
    z_SolverCblk         *cblk;
    pastix_complex64_t *L      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *work   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_int_t        stride = STARPU_MATRIX_GET_LD(buffers[0]);
    int                 me     = starpu_worker_get_id();

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk);

    fprintf(stderr, "NOT YET IMPLEMENTED");
    /* sopalin_data->thread_data[me]->nbpivot += */
    /*     core_zsytrfsp1d_sytrf(cblk, L, */
    /*                           sopalin_data->critere, */
    /*                           work); */
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zsytrfsp1d_trsm_cuda - Apply all the trsm updates on one panel.
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
void starpu_zsytrfsp1d_trsm_cuda(void * buffers[], void * _args)
{
    z_Sopalin_Data_t     *sopalin_data;
    z_SolverCblk         *cblk;
    pastix_complex64_t *L      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_int_t        stride = STARPU_MATRIX_GET_LD(buffers[0]);

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk);

    fprintf(stderr, "NOT YET IMPLEMENTED");
    /* core_zsytrfsp1d_trsm(cblk, L); */

    SUBMIT_GEMMS_IF_NEEDED;
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zsytrfsp1d_cuda - Computes the LU factorization of one panel and apply
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
void starpu_zsytrfsp1d_cuda(void * buffers[], void * _args)
{
    z_Sopalin_Data_t     *sopalin_data;
    z_SolverCblk         *cblk;
    pastix_complex64_t *L      = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *work   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_int_t        stride = STARPU_MATRIX_GET_LD(buffers[0]);
    int                 me     = starpu_worker_get_id();

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk);

    fprintf(stderr, "NOT YET IMPLEMENTED");
    /* sopalin_data->thread_data[me]->nbpivot += */
    /*     core_zsytrfsp1d(cblk, */
    /*                     L, */
    /*                     sopalin_data->critere, */
    /*                     work); */
    SUBMIT_GEMMS_IF_NEEDED;
#ifdef PASTIX_DUMP_CBLK
    {
        char name[256];
        sprintf(name, "cblk_%d_after_trf_trsm", cblk->gcblknum);
        cblk_save(cblk, name, L);
    }
#endif
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zsytrfsp1d_gemm_cuda - Computes the Cholesky factorization of one panel
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

void starpu_zsytrfsp1d_gemm_cuda(void * buffers[], void * _args)
{
    z_Sopalin_Data_t     *sopalin_data;
    z_SolverMatrix       *datacode;
    z_SolverCblk         *cblk;
    z_SolverBlok         *blok;
    z_SolverCblk         *fcblk;
    pastix_int_t dimi, dimj, dima, dimb;
    pastix_complex64_t *Aik, *Aij;
    pastix_complex64_t *L    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *Cl   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
    pastix_complex64_t *work = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[2]);
    pastix_int_t        ldw  = (pastix_int_t)STARPU_VECTOR_GET_NX(buffers[2])/2;
    int                *all_blocktab = (int*)STARPU_VECTOR_GET_PTR(buffers[3]);
    int                *blocktab, *fblocktab;
    pastix_complex64_t *work2= work + ldw;

    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk, &blok, &fcblk);
    datacode = sopalin_data->datacode;
    assert(z_is_block_inside_fblock(blok, fcblk->fblokptr));
    assert(cblk->stride == STARPU_MATRIX_GET_LD(buffers[0]));
    assert(cblk_colnbr(cblk)  == STARPU_MATRIX_GET_NY(buffers[0]));
    fblocktab = &(all_blocktab[2*blok_getnum(datacode,
                                             fcblk->fblokptr)]);
    if (cblk_ishalo(datacode, fcblk))
        blocktab  = &(all_blocktab[2*(SYMB_BLOKNBR + hblok_getnum(datacode, blok))]);
    else
        if (cblk_islocal(datacode, fcblk))
            blocktab  = &(all_blocktab[2*(blok_getnum(datacode, blok))]);

    /* Aik starts at the First blok */
    Aik = L + blok->coefind;
    dimj = blok_rownbr(blok);
    dimi = cblk->stride - blok->coefind;

    {
        CU_FLOAT cu_alpha = CU_FLOAT_INIT(1.0,  0.0);
        CU_FLOAT cu_beta  = CU_FLOAT_INIT(0.0,  0.0);
        GENERATE_SM_VERSION_NAME(gemdm)('n', 't',
                                        (int)dimi, (int)dimj, (int)dima,
                                        cu_alpha,
                                        (CU_FLOAT*)Aik, (int)cblk->stride,
                                        (CU_FLOAT*)L,   (int)cblk->stride+1,
                                        (CU_FLOAT*)Aik, (int)cblk->stride,
                                        cu_beta,
                                        (CU_FLOAT*)Cl,  (int)fcblk->stride,
                                        cblk_bloknbr(cblk),
                                        blocktab,
                                        cblk_bloknbr(fcblk),
                                        fblocktab,
                                        starpu_cuda_get_local_stream());
        cudaStreamSynchronize(starpu_cuda_get_local_stream());
    }
    SUBMIT_TRF_IF_NEEDED;
#ifdef PASTIX_DUMP_CBLK
    {
        char name[256];
        sprintf(name, "cblk_%d_after_gemm_%d_%d_%d_on_%d", fcblk->gcblknum,
                cblk->gcblknum, blok - cblk->fblokptr, fcblk->gcblknum,
                sopalin_data->datacode->clustnum);
        cblk_save(fcblk, name, Cl);
    }
#endif
}


/**
 *******************************************************************************
 *
 * @ingroup pastix_starpu_kernel
 *
 * starpu_zsytrfsp1d_geadd_cuda - Computes the addition of a fanin column block
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
starpu_zsytrfsp1d_syadd_cuda(void * buffers[], void * _args) {
    z_Sopalin_Data_t *sopalin_data;
    z_SolverCblk *cblk1, *cblk2;
    pastix_complex64_t *L    = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[0]);
    pastix_complex64_t *Cl   = (pastix_complex64_t*)STARPU_MATRIX_GET_PTR(buffers[1]);
    starpu_codelet_unpack_args(_args, &sopalin_data, &cblk1, &cblk2);
#ifdef PASTIX_DUMP_CBLK
    {
        char name[256];
        sprintf(name, "cblk_dst_%d_before_add_from_%d", cblk2->gcblknum,
                fcblk_getorigin(sopalin_data->datacode, cblk1));
        cblk_save(cblk2, name, Cl);
    }
    {
        char name[256];
        sprintf(name, "cblk_src_%d_before_add_from_%d", cblk1->gcblknum,
                fcblk_getorigin(sopalin_data->datacode, cblk1));
        cblk_save(cblk1, name, L);
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
        cblk_save(cblk2, name, Cl);
    }
#endif
}
