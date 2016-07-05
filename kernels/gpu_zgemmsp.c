/**
 *
 * @file gpu_zgemmsp.c
 *
 *  PaStiX kernel routines
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
#include "common.h"
#include <cblas.h>
#include "blend/solver.h"
#include "kernels/pastix_zcores.h"
#include "kernels/pastix_cuda.h"
#include <cublas.h>

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * gpu_zgemmsp_1d1d - Computes the updates associated to one off-diagonal block
 * between two cblk stored as 1D block columns.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          If uplo == PastixLower, the contribution of:
 *          (block .. (cblk[1].fblokptr-1)) -by- block is computed and added to
 *          C, otherwise the contribution:
 *          (block+1 .. (cblk[1].fblokptr-1)) -by- block is computed and added
 *          to C.
 *          The pointer to the data structure that describes the panel from
 *          which we compute the contributions. Next column blok must be
 *          accessible through cblk[1].
 *
 * @param[in] trans
 *          Specify the transposition used for the B matrix. It has to be either
 *          PastixTrans or PastixConjTrans.
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok
 *          The block from which we compute the contributions.
 *
 * @param[in] fcblk
 *          The pointer to the data structure that describes the panel on which
 *          we compute the contributions. The C pointer must be one of the
 *          oceftab from this fcblk. Next column blok must be accessible through
 *          fcblk[1].
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in] B The pointer to the coeftab of the cblk.lcoeftab matrix storing
 *          the coefficients of the panel, if Symmetric/Hermitian cases or if
 *          upper part is computed; cblk.ucoeftab otherwise. Must be of size
 *          cblk.stride -by- cblk.width
 *
 * @param[in,out] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 * @param[in] work
 *          Temporary memory buffer.
 *
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal
 *          block.
 *
 *******************************************************************************/
void
gpu_zgemmsp( int uplo, int trans,
             const SolverCblk      *cblk,
             const SolverBlok      *blok,
                   SolverCblk      *fcblk,
             const cuDoubleComplex *A,
             const cuDoubleComplex *B,
                   cuDoubleComplex *C,
                   cudaStream_t stream )
{
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex mzone = make_cuDoubleComplex(-1., 0.);
    cuDoubleComplex zone  = make_cuDoubleComplex( 1., 0.);
#else
    double mzone = -1.;
    double zone  =  1.;
#endif
    gemm_params_t params;
    const SolverBlok *iterblok;
    const SolverBlok *fblok;
    const SolverBlok *lblok;

    pastix_int_t stride, stridef, indblok;
    pastix_int_t N, K, max_m = 0;
    int i, shift, count, ldb;

    shift = (uplo == PastixUpper) ? 1 : 0;

    stride  = cblk->stride;
    stridef = fcblk->stride;
    K = cblk_colnbr( cblk );

    /* First blok */
    indblok = blok->coefind;

    N = blok_rownbr( blok );

    /* Move B to the right pointer */
    B = B + indblok;
    ldb = (cblk->cblktype & CBLK_SPLIT) ? N : stride;

    /* Get the first block of the distant panel */
    fblok = fcblk->fblokptr;

    /* Get the last block to stop the iteration */
    lblok = cblk[1].fblokptr;
    count = (lblok - blok) - shift;

    for (iterblok=blok+shift, i=0; iterblok<lblok; iterblok++, i++) {
        /* Find facing blok */
        while (!is_block_inside_fblock( iterblok, fblok ))
        {
            fblok++;
            assert( fblok < fcblk[1].fblokptr );
        }

        stridef = (fcblk->cblktype  & CBLK_SPLIT) ? blok_rownbr( fblok ) : stridef;
        params.p[i].M    = blok_rownbr( iterblok );
        params.p[i].Aptr = A + iterblok->coefind;
        params.p[i].lda  = (cblk->cblktype  & CBLK_SPLIT) ? params.p[i].M : stride;
        params.p[i].Cptr = C +
            fblok->coefind + iterblok->frownum - fblok->frownum +
            (blok->frownum - fcblk->fcolnum) * stridef;
        params.p[i].ldc  = stridef;

        max_m = pastix_imax( max_m, params.p[i].M);

        if (i+1 == MAX_BATCH_COUNT) {
            pastix_zgemm_vbatched_nt(
                trans, N, K,
                /* alpha  */  mzone,
                /* B      */  B, ldb,
                /* beta   */  zone,
                max_m, MAX_BATCH_COUNT,
                stream, params );

            /* Restart the loop */
            i = -1;
            count -= MAX_BATCH_COUNT;
            max_m = 0;
        }
    }

    if (count > 0) {
        pastix_zgemm_vbatched_nt(
            trans, N, K,
            /* alpha  */  mzone,
            /* B      */  B, ldb,
            /* beta   */  zone,
            max_m, count,
            stream, params );
    }
}

static char transstr[3] = { 'N', 'T', 'C' };

void
gpu_zgemmsp_2d2dsub( int trans,
                     pastix_int_t blok_mk,
                     pastix_int_t blok_kn,
                     pastix_int_t blok_mn,
                     const SolverCblk      *cblk,
                           SolverCblk      *fcblk,
                     const cuDoubleComplex *A,
                     const cuDoubleComplex *B,
                           cuDoubleComplex *C,
                           cudaStream_t stream )
{
#if defined(PRECISION_z) || defined(PRECISION_c)
    cuDoubleComplex mzone = make_cuDoubleComplex(-1., 0.);
    cuDoubleComplex zone  = make_cuDoubleComplex( 1., 0.);
#else
    double mzone = -1.;
    double zone  =  1.;
#endif
    const SolverBlok *blokA, *blokB, *blokC;
    const SolverBlok *bA, *bB, *bC;
    const SolverBlok *fblokK, *lblokK;
    const SolverBlok *fblokN, *lblokN;

    const cuDoubleComplex *Aptr, *Bptr;
    cuDoubleComplex *Cptr;
    pastix_int_t M, N, K, lda, ldb, ldc, cblk_n, cblk_m;
    size_t offsetA, offsetB, offsetC;

    /* cblk is stored in 1D and fcblk in 2D */
    assert( cblk->cblktype  & CBLK_SPLIT );
    assert( fcblk->cblktype & CBLK_SPLIT );

    /**
     * Blocs on column K
     */
    fblokK = cblk[0].fblokptr;
    lblokK = cblk[1].fblokptr;

    blokB = fblokK + blok_kn;
    offsetB = blokB->coefind;
    cblk_n = blokB->fcblknm;

    blokA = fblokK + blok_mk;
    offsetA = blokA->coefind;
    cblk_m = blokA->fcblknm;

    /**
     * Blocs on column N
     */
    fblokN = fcblk[0].fblokptr;
    lblokN = fcblk[1].fblokptr;

    blokC = fblokN + blok_mn;
    offsetC = blokC->coefind;
    assert( blokC->fcblknm == cblk_m );

    K = cblk_colnbr( cblk );

    cublasSetKernelStream( stream );
    bC = blokC;
    for (bA = blokA; (bA < lblokK) && (bA->fcblknm == cblk_m); bA++) {
        M = blok_rownbr(bA);
        Aptr = A + bA->coefind - offsetA;
        lda = M;

        /* Find facing C blok */
        while (!is_block_inside_fblock( bA, bC )) {
            bC++;
            assert( bC < lblokN );
        }

        Cptr = C + bC->coefind - offsetC;
        ldc = blok_rownbr(bC);

        for (bB = blokB; (bB <= bA) && (bB->fcblknm == cblk_n); bB++) {
            N = blok_rownbr( bB );
            Bptr = B + bB->coefind - offsetB;
            ldb = N;

            cublasZgemm( 'N', transstr[trans - PastixNoTrans],
                         M, N, K,
                         mzone, Aptr, lda,
                                Bptr, ldb,
                          zone, Cptr + (bA->frownum - bC->frownum)
                                     + (bB->frownum - fcblk->fcolnum) * ldc, ldc );
        }
    }
}

