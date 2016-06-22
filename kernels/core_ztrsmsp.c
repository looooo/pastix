/**
 *
 * @file core_ztrsmsp.c
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
#include "pastix_zcores.h"

static pastix_complex64_t zone  =  1.;

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_ztrsmsp_1d1d - Apply all the trsm updates on a panel stored in 1D
 * layout.
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
static inline int
core_ztrsmsp_1d( int side, int uplo, int trans, int diag,
                       SolverCblk         *cblk,
                 const pastix_complex64_t *A,
                       pastix_complex64_t *C )
{
    SolverBlok *fblok;
    pastix_int_t M, N, lda;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    lda   = cblk->stride;
    fblok = cblk->fblokptr;  /* The diagonal block */

    /* vertical dimension */
    M = lda - N;

    /* if there is an extra-diagonal bloc in column block */
    assert( fblok + 1 < cblk[1].fblokptr );
    assert( blok_rownbr( fblok) == N );

    /* first extra-diagonal bloc in column block address */
    C = C + fblok[1].coefind;

    cblas_ztrsm(CblasColMajor,
                side, uplo, trans, diag,
                M, N,
                CBLAS_SADDR(zone), A, lda,
                                   C, lda);

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_ztrsmsp_2d - Computes the updates associated to one off-diagonal block
 * between two cblk stored in 2D.
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
 * @param[in] tol
 *          Tolerance for low-rank compression kernels
 *
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal
 *          block.
 *
 *******************************************************************************/
static inline int
core_ztrsmsp_2d( int side, int uplo, int trans, int diag,
                       SolverCblk         *cblk,
                 const pastix_complex64_t *A,
                       pastix_complex64_t *C )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda, ldc;
    pastix_complex64_t *blokC;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */
    lda   = blok_rownbr( fblok );

    assert( blok_rownbr(fblok) == N );

    for (blok=fblok+1; blok<lblok; blok++) {

        blokC = C + blok->coefind;
        M   = blok_rownbr(blok);
        ldc = M;

        cblas_ztrsm(CblasColMajor,
                    side, uplo, trans, diag,
                    M, N,
                    CBLAS_SADDR(zone), A, lda,
                                       blokC, ldc);
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * core_ztrsmsp - Computes the updates associated to one off-diagonal block.
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
 * @param[in] tol
 *          Tolerance for low-rank compression kernels
 *
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal
 *          block.
 *
 *******************************************************************************/
void core_ztrsmsp( int side, int uplo, int trans, int diag,
                         SolverCblk         *cblk,
                   const pastix_complex64_t *A,
                         pastix_complex64_t *C )
{
    if (  cblk[0].fblokptr + 1 < cblk[1].fblokptr ) {
        if ( cblk->cblktype & CBLK_SPLIT ) {
            core_ztrsmsp_2d( side, uplo, trans, diag,
                             cblk, A, C );
        }
        else {
            core_ztrsmsp_1d( side, uplo, trans, diag,
                             cblk, A, C );
        }
    }
}

