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
#include "cblas.h"
#include "blend/solver.h"
#include "pastix_zcores.h"

static pastix_complex64_t  zone =  1.;
static pastix_complex64_t mzone = -1.;

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
 * @param[in] side
 *
 * @param[in] uplo
 *
 * @param[in] trans
 *
 * @param[in] diag
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and C pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
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
    assert(!(cblk->cblktype & CBLK_SPLIT));

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
 * @param[in] side
 *
 * @param[in] uplo
 *
 * @param[in] trans
 *
 * @param[in] diag
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and C pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
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
    assert( cblk->cblktype & CBLK_SPLIT );

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
 * core_ztrsmsp_2d - Computes the updates associated to one off-diagonal block
 * between two cblk stored in 2D.
 *
 *******************************************************************************
 *
 * @param[in] side
 *
 * @param[in] uplo
 *
 * @param[in] trans
 *
 * @param[in] diag
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and C pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] blok_m
 *          Index of the first off-diagonal block in cblk that is solved. The
 *          TRSM is also applied to all the folowing blocks which are facing the
 *          same diagonal block
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
 *
 *******************************************************************************
 *
 * @return
 *          The number of static pivoting during factorization of the diagonal
 *          block.
 *
 *******************************************************************************/
int
core_ztrsmsp_2dsub( int side, int uplo, int trans, int diag,
                          SolverCblk         *cblk,
                          pastix_int_t        blok_m,
                    const pastix_complex64_t *A,
                          pastix_complex64_t *C )
{
    const SolverBlok *fblok, *lblok, *blok;
    pastix_int_t M, N, lda, ldc, offset, cblk_m;
    pastix_complex64_t *Cptr;

    N     = cblk->lcolnum - cblk->fcolnum + 1;
    fblok = cblk[0].fblokptr;  /* The diagonal block */
    lblok = cblk[1].fblokptr;  /* The diagonal block of the next cblk */
    lda   = blok_rownbr( fblok );

    assert( blok_rownbr(fblok) == N );
    assert( cblk->cblktype & CBLK_SPLIT );

    blok   = fblok + blok_m;
    offset = blok->coefind;
    cblk_m = blok->fcblknm;

    for (; (blok < lblok) && (blok->fcblknm == cblk_m); blok++) {

        Cptr = C + blok->coefind - offset;
        M   = blok_rownbr(blok);
        ldc = M;

        cblas_ztrsm(CblasColMajor,
                    side, uplo, trans, diag,
                    M, N,
                    CBLAS_SADDR(zone), A, lda,
                                       Cptr, ldc);
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
 * @param[in] side
 *
 * @param[in] uplo
 *
 * @param[in] trans
 *
 * @param[in] diag
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] A
 *          The pointer to the coeftab of the cblk.lcoeftab matrix storing the
 *          coefficients of the panel when the Lower part is computed,
 *          cblk.ucoeftab otherwise. Must be of size cblk.stride -by- cblk.width
 *
 * @param[in,out] C
 *          The pointer to the fcblk.lcoeftab if the lower part is computed,
 *          fcblk.ucoeftab otherwise.
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



/**
 *******************************************************************************
 *
 * @ingroup pastix_kernel
 *
 * solve_ztrsmsp - Apply the solve related to one cblk to all the right hand side.
 *
 *******************************************************************************
 *
 * @param[in] side
 *
 * @param[in] uplo
 *
 * @param[in] trans
 *
 * @param[in] diag
 *
 * @param[in] datacode
 *
 * @param[in] cblk
 *          The cblk structure to which block belongs to. The A and B pointers
 *          must be the coeftab of this column block.
 *          Next column blok must be accessible through cblk[1].
 *
 * @param[in] nrhs
 *          The number of right hand side.
 *
 * @param[in,out] b
 *          The pointer to vectors of the right hand side
 *
 * @param[in] ldb
 *          The leading dimension of b
 *
 *******************************************************************************
 *
 * @return
 *
 *******************************************************************************/
void solve_ztrsmsp( int side, int uplo, int trans, int diag,
                    SolverMatrix *datacode, SolverCblk *cblk,
                    int nrhs, pastix_complex64_t *b, int ldb )
{
    SolverCblk *fcbk;
    SolverBlok *blok;
    pastix_complex64_t *A;
    pastix_int_t j, tempm, tempn, lda;

    tempn = cblk->lcolnum - cblk->fcolnum + 1;
    lda = (cblk->cblktype & CBLK_SPLIT) ? blok_rownbr( cblk->fblokptr ) : cblk->stride;

    /*
     *  Left / Upper / NoTrans
     */
    if (side == PastixLeft) {
        if (uplo == PastixUpper) {
            A = (pastix_complex64_t*)(cblk->ucoeftab);

            /*  We store U^t, so we swap uplo and trans */
            if (trans == PastixNoTrans) {

                /* Solve the diagonal block */
                cblas_ztrsm(
                    CblasColMajor, CblasLeft, CblasLower,
                    CblasTrans, (enum CBLAS_DIAG)diag,
                    tempn, nrhs,
                    CBLAS_SADDR(zone), A, lda,
                                       b + cblk->lcolidx, ldb );

                /* Apply the update */
                for (j = cblk[1].brownum-1; j>=cblk[0].brownum; j-- ) {
                    blok = datacode->bloktab + datacode->browtab[j];
                    fcbk = datacode->cblktab + blok->lcblknm;
                    tempm = fcbk->lcolnum - fcbk->fcolnum + 1;
                    tempn = blok->lrownum - blok->frownum + 1;
                    A   = (pastix_complex64_t*)(fcbk->ucoeftab);
                    lda = (fcbk->cblktype & CBLK_SPLIT) ? tempn : fcbk->stride;

                    cblas_zgemm(
                        CblasColMajor, CblasTrans, CblasNoTrans,
                        tempm, nrhs, tempn,
                        CBLAS_SADDR(mzone), A + blok->coefind, lda,
                                            b + cblk->lcolidx + blok->frownum - cblk->fcolnum, ldb,
                        CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
                }
            }
            /*
             *  Left / Upper / [Conj]Trans
             */
            else {
                assert(0 /* Not implemented */);
            }
        }
        else {
            A = (pastix_complex64_t*)(cblk->lcoeftab);

            /*
             *  Left / Lower / NoTrans
             */
            if (trans == PastixNoTrans) {

                /* In sequential */
                assert( cblk->fcolnum == cblk->lcolidx );

                /* Solve the diagonal block */
                cblas_ztrsm(
                    CblasColMajor, CblasLeft, CblasLower,
                    CblasNoTrans, (enum CBLAS_DIAG)diag,
                    tempn, nrhs,
                    CBLAS_SADDR(zone), A, lda,
                                       b + cblk->lcolidx, ldb );

                /* Apply the update */
                for (blok = cblk[0].fblokptr+1; blok < cblk[1].fblokptr; blok++ ) {
                    fcbk  = datacode->cblktab + blok->fcblknm;
                    tempm = blok->lrownum - blok->frownum + 1;

                    assert( blok->frownum >= fcbk->fcolnum );
                    assert( tempm <= (fcbk->lcolnum - fcbk->fcolnum + 1));

                    lda = (cblk->cblktype & CBLK_SPLIT) ? tempm : cblk->stride;

                    cblas_zgemm(
                        CblasColMajor, CblasNoTrans, CblasNoTrans,
                        tempm, nrhs, tempn,
                        CBLAS_SADDR(mzone), A + blok->coefind, lda,
                                            b + cblk->lcolidx, ldb,
                        CBLAS_SADDR(zone),  b + fcbk->lcolidx + blok->frownum - fcbk->fcolnum, ldb );
                }
            }
            /*
             *  Left / Lower / [Conj]Trans
             */
            else {
                /* Solve the diagonal block */
                cblas_ztrsm(
                    CblasColMajor, CblasLeft, CblasLower,
                    (enum CBLAS_TRANSPOSE)trans, (enum CBLAS_DIAG)diag,
                    tempn, nrhs,
                    CBLAS_SADDR(zone), A, lda,
                                       b + cblk->lcolidx, ldb );

                /* Apply the update */
                for (j = cblk[1].brownum-1; j>=cblk[0].brownum; j-- ) {
                    blok = datacode->bloktab + datacode->browtab[j];
                    fcbk = datacode->cblktab + blok->lcblknm;
                    tempm = fcbk->lcolnum - fcbk->fcolnum + 1;
                    tempn = blok->lrownum - blok->frownum + 1;
                    A   = (pastix_complex64_t*)(fcbk->lcoeftab);
                    lda = (fcbk->cblktype & CBLK_SPLIT) ? tempn : fcbk->stride;

                    cblas_zgemm(
                        CblasColMajor, (enum CBLAS_TRANSPOSE)trans, CblasNoTrans,
                        tempm, nrhs, tempn,
                        CBLAS_SADDR(mzone), A + blok->coefind, lda,
                                            b + cblk->lcolidx + blok->frownum - cblk->fcolnum, ldb,
                        CBLAS_SADDR(zone),  b + fcbk->lcolidx, ldb );
                }
            }
        }
    }
    /**
     * Right
     */
    else {
        assert(0 /* Not implemented */);
    }
}


